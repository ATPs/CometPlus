from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import sys
from typing import Dict, List, Optional, Tuple

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from .peptide_index_digestion import (
    iter_fasta,
    iter_peptides,
    iter_protein_sequences,
    normalize_sequence,
    peptide_flanks,
)
from .peptide_index_variants import (
    SequenceTask,
    VariantContext,
    VariantTaskRow,
    apply_set_masses,
    build_aa_masses,
    enumerate_sequence_variants,
    enumerate_variant_chunk,
    get_static_mods,
    init_variant_worker,
    resolve_thread_count,
)
from .pyLoadParameters import CometParams


def build_tables(
    params: CometParams,
    fasta_paths: Optional[List[str]] = None,
    protein_sequences: Optional[List[str]] = None,
    max_proteins: Optional[int] = None,
    progress: bool = False,
    progress_every: int = 1000,
    use_protein_name: bool = False,
    threads: int = 1,
    unimod_variable_map: Optional[Dict[int, str]] = None,
    unimod_fixed_map: Optional[Dict[str, str]] = None,
) -> Dict[str, List[Tuple]]:
    _ = progress_every  # Backward compatibility: retained in signature, no periodic logging.
    thread_count = resolve_thread_count(threads)
    residue_mods, term_mods = get_static_mods(params)
    aa_masses, h2o_mass = build_aa_masses(bool(params.parsed.get("mass_type_parent", 1)))
    apply_set_masses(aa_masses, params)
    var_mods = sorted(params.variable_mods, key=lambda mod: mod.index)
    stderr_isatty = getattr(sys.stderr, "isatty", lambda: False)()
    use_tqdm = bool(progress and tqdm is not None and stderr_isatty)

    enzyme = params.enzymes["search"]
    enzyme2 = params.enzymes["search2"] if params.parsed.get("search_enzyme2_number", 0) else None
    num_termini = int(params.parsed.get("num_enzyme_termini", 2))
    max_missed = int(params.parsed.get("allowed_missed_cleavage", 2))
    min_len, max_len = params.parsed.get("peptide_length_range", (5, 50))
    min_mass, max_mass = params.parsed.get("digest_mass_range", (600.0, 5000.0))
    clip_nterm_methionine = int(params.parsed.get("clip_nterm_methionine", 0))
    max_var_mods = int(params.parsed.get("max_variable_mods_in_peptide", 5))
    require_var_mod = int(params.parsed.get("require_variable_mod", 0))

    protein_rows: List[Tuple] = []
    seq_info: Dict[str, Dict] = {}
    peptide_records = 0

    def progress_log(message: str) -> None:
        if progress and not use_tqdm:
            print(message, file=sys.stderr, flush=True)

    def record_peptide(
        pep_seq: str,
        abs_start: int,
        abs_end: int,
        pr_start: int,
        pr_end: int,
        prev_aa: str,
        next_aa: str,
        protein_len: int,
        nterm_offset: int,
        protein_id,
        offset: int,
    ) -> None:
        nonlocal peptide_records
        peptide_records += 1
        entry = seq_info.get(pep_seq)
        if entry is None:
            entry = {
                "protein_ids": set(),
                "primary": (
                    offset,
                    protein_id,
                    abs_start,
                    abs_end,
                    prev_aa,
                    next_aa,
                    protein_len,
                    nterm_offset,
                ),
                "locations": set(),
            }
            seq_info[pep_seq] = entry
        entry["protein_ids"].add(protein_id)
        entry["locations"].add((protein_id, pr_start, pr_end))
        primary = entry["primary"]
        if offset < primary[0]:
            entry["primary"] = (
                offset,
                protein_id,
                abs_start,
                abs_end,
                prev_aa,
                next_aa,
                protein_len,
                nterm_offset,
            )

    protein_count = 0
    max_proteins_label = f" (max {max_proteins})" if max_proteins else ""
    digestion_bar = None
    if progress:
        if use_tqdm:
            digestion_total = None
            if protein_sequences:
                digestion_total = len(protein_sequences)
                if max_proteins is not None:
                    digestion_total = min(digestion_total, max_proteins)
            digestion_bar = tqdm(
                total=digestion_total,
                desc="Index: digest proteins",
                unit="protein",
                file=sys.stderr,
                leave=False,
            )
        else:
            progress_log(f"Index: starting protein digestion{max_proteins_label}")
    if protein_sequences:
        protein_iter = iter_protein_sequences(protein_sequences)
    else:
        fasta_paths = fasta_paths or []
        protein_iter = (
            (header, seq, offset)
            for fasta_path in fasta_paths
            for header, seq, offset in iter_fasta(fasta_path)
        )
    for header, seq_raw, offset in protein_iter:
        if max_proteins is not None and protein_count >= max_proteins:
            break
        protein_count += 1
        if use_protein_name:
            protein_id = header.split()[0] if header else str(protein_count)
        else:
            protein_id = protein_count
        seq = normalize_sequence(seq_raw)
        pr_seq = seq.replace("*", "")
        protein_rows.append((protein_id, pr_seq, offset, header))
        if not seq:
            if digestion_bar is not None:
                digestion_bar.update(1)
            continue
        protein_starts_with_m = seq.startswith("M")
        raw_pos = 0
        pr_pos = 0
        for segment in seq.split("*"):
            segment_len = len(segment)
            if segment_len == 0:
                raw_pos += 1
                continue
            segment_start = raw_pos
            pr_segment_start = pr_pos
            for start, end in iter_peptides(
                segment,
                enzyme,
                enzyme2,
                num_termini,
                max_missed,
                min_len,
                max_len,
            ):
                pep_seq = segment[start : end + 1]
                abs_start = segment_start + start
                abs_end = segment_start + end
                pr_start = pr_segment_start + start + 1
                pr_end = pr_segment_start + end + 1
                prev_aa, next_aa = peptide_flanks(seq, abs_start, abs_end)
                record_peptide(
                    pep_seq,
                    abs_start,
                    abs_end,
                    pr_start,
                    pr_end,
                    prev_aa,
                    next_aa,
                    len(seq),
                    0,
                    protein_id,
                    offset,
                )
            raw_pos += segment_len + 1
            pr_pos += segment_len

        if clip_nterm_methionine and protein_starts_with_m and len(seq) > 1:
            clipped_seq = seq[1:]
            clipped_segment = clipped_seq.split("*", 1)[0]
            if clipped_segment:
                for start, end in iter_peptides(
                    clipped_segment,
                    enzyme,
                    enzyme2,
                    num_termini,
                    max_missed,
                    min_len,
                    max_len,
                ):
                    if start != 0:
                        continue
                    pep_seq = clipped_segment[start : end + 1]
                    abs_start = start + 1
                    abs_end = end + 1
                    pr_start = start + 2
                    pr_end = end + 2
                    _, next_aa = peptide_flanks(seq, abs_start, abs_end)
                    record_peptide(
                        pep_seq,
                        abs_start,
                        abs_end,
                        pr_start,
                        pr_end,
                        "-",
                        next_aa,
                        len(seq),
                        1,
                        protein_id,
                        offset,
                    )

        if digestion_bar is not None:
            digestion_bar.update(1)

    if digestion_bar is not None:
        digestion_bar.close()
    progress_log(
        f"Index: finished proteins {protein_count}; "
        f"peptides seen {peptide_records}; unique sequences {len(seq_info)}"
    )

    # Assign sequence IDs deterministically
    sequences = sorted(seq_info.keys())
    seq_id_map: Dict[str, int] = {seq: idx + 1 for idx, seq in enumerate(sequences)}

    peptide_sequence_rows: List[Tuple] = []
    peptide_sequence_protein_rows: List[Tuple] = []
    peptide_protein_location_rows: List[Tuple] = []

    for seq in sequences:
        entry = seq_info[seq]
        _, primary_protein_id, *_ = entry["primary"]
        seq_id = seq_id_map[seq]
        peptide_sequence_rows.append((seq_id, seq, len(seq), primary_protein_id))
        for pid in sorted(entry["protein_ids"]):
            peptide_sequence_protein_rows.append((seq_id, pid))
        for protein_id, pep_start, pep_end in entry["locations"]:
            peptide_protein_location_rows.append((protein_id, seq_id, pep_start, pep_end, seq))

    def _protein_id_sort_key(protein_id) -> Tuple[int, object]:
        if isinstance(protein_id, int):
            return (0, protein_id)
        return (1, str(protein_id))

    peptide_protein_location_rows.sort(
        key=lambda row: (
            _protein_id_sort_key(row[0]),
            row[2],
            row[3],
            row[1],
            row[4],
        )
    )

    # Variable mod table rows
    variable_mod_rows: List[Tuple] = []
    for mod in var_mods:
        variable_mod_rows.append(
            (
                mod.index,
                mod.residues,
                mod.mass,
                mod.binary_mod,
                mod.min_per_pep,
                mod.max_per_pep,
                mod.term_distance,
                mod.which_term,
                mod.require_this_mod,
                mod.neutral_loss1,
                mod.neutral_loss2,
                unimod_variable_map.get(mod.index) if unimod_variable_map is not None else None,
            )
        )

    static_mod_rows: List[Tuple] = []
    for residue, delta in sorted(residue_mods.items()):
        static_mod_rows.append(
            (
                residue,
                delta,
                "residue",
                unimod_fixed_map.get(f"add_{residue}") if unimod_fixed_map is not None else None,
            )
        )
    if abs(term_mods["Nterm_peptide"]) > 1e-12:
        static_mod_rows.append(
            (
                "-",
                term_mods["Nterm_peptide"],
                "N-term",
                unimod_fixed_map.get("add_Nterm_peptide") if unimod_fixed_map is not None else None,
            )
        )
    if abs(term_mods["Cterm_peptide"]) > 1e-12:
        static_mod_rows.append(
            (
                "-",
                term_mods["Cterm_peptide"],
                "C-term",
                unimod_fixed_map.get("add_Cterm_peptide") if unimod_fixed_map is not None else None,
            )
        )
    if abs(term_mods["Nterm_protein"]) > 1e-12:
        static_mod_rows.append(
            (
                "-",
                term_mods["Nterm_protein"],
                "protein N-term",
                unimod_fixed_map.get("add_Nterm_protein") if unimod_fixed_map is not None else None,
            )
        )
    if abs(term_mods["Cterm_protein"]) > 1e-12:
        static_mod_rows.append(
            (
                "-",
                term_mods["Cterm_protein"],
                "protein C-term",
                unimod_fixed_map.get("add_Cterm_protein") if unimod_fixed_map is not None else None,
            )
        )

    sequence_tasks: List[SequenceTask] = []
    for seq in sequences:
        entry = seq_info[seq]
        (
            _offset,
            _primary_protein_id,
            start,
            end,
            prev_aa,
            next_aa,
            protein_len,
            nterm_offset,
        ) = entry["primary"]
        seq_id = seq_id_map[seq]
        sequence_tasks.append(
            (
                seq_id,
                seq,
                start,
                end,
                prev_aa,
                next_aa,
                protein_len,
                nterm_offset,
            )
        )

    peptide_variant_rows: List[Tuple] = []
    peptide_variant_mod_rows: List[Tuple] = []
    variant_id = 0

    total_sequences = len(sequence_tasks)
    variant_bar = None
    if total_sequences:
        if use_tqdm:
            variant_bar = tqdm(
                total=total_sequences,
                desc="Index: enumerate variants",
                unit="sequence",
                file=sys.stderr,
                leave=False,
            )
        else:
            progress_log(f"Index: starting variant enumeration for {total_sequences} sequences")
    variant_context = VariantContext(
        aa_masses=aa_masses,
        h2o_mass=h2o_mass,
        residue_mods=residue_mods,
        term_mods=term_mods,
        var_mods=var_mods,
        max_var_mods=max_var_mods,
        require_var_mod=require_var_mod,
        min_mass=min_mass,
        max_mass=max_mass,
        unimod_variable_map=unimod_variable_map,
        unimod_fixed_map=unimod_fixed_map,
    )

    variant_task_rows: List[VariantTaskRow] = []
    if thread_count > 1 and total_sequences > 1:
        chunk_size = max(1, min(500, total_sequences // (thread_count * 8) or 1))
        futures = []
        with ProcessPoolExecutor(
            max_workers=thread_count,
            initializer=init_variant_worker,
            initargs=(variant_context,),
        ) as executor:
            for chunk_start in range(0, total_sequences, chunk_size):
                chunk = sequence_tasks[chunk_start : chunk_start + chunk_size]
                futures.append((len(chunk), executor.submit(enumerate_variant_chunk, chunk)))

            for chunk_len, future in futures:
                chunk_rows = future.result()
                variant_task_rows.extend(chunk_rows)
                if variant_bar is not None:
                    variant_bar.update(chunk_len)
    else:
        for task in sequence_tasks:
            variant_task_rows.extend(enumerate_sequence_variants(task, variant_context))
            if variant_bar is not None:
                variant_bar.update(1)

    if variant_bar is not None:
        variant_bar.close()
    if total_sequences:
        progress_log(
            f"Index: finished variant enumeration for {total_sequences} sequences; "
            f"variants {len(variant_task_rows)}"
        )

    for (
        seq_id,
        mh_plus,
        prev_aa,
        next_aa,
        sites_text,
        sites_unimod_text,
        total_mods,
        fixed_sites_text,
        fixed_sites_unimod_text,
        fixed_mod_count,
        mass_bin10,
        sites,
    ) in variant_task_rows:
        variant_id += 1
        pep_seq = sequences[seq_id - 1]
        peptide_variant_rows.append(
            (
                variant_id,
                seq_id,
                pep_seq,
                mh_plus,
                prev_aa,
                next_aa,
                sites_text,
                sites_unimod_text,
                total_mods,
                fixed_sites_text,
                fixed_sites_unimod_text,
                fixed_mod_count,
                mass_bin10,
            )
        )
        for pos, mod_idx in enumerate(sites):
            if mod_idx:
                peptide_variant_mod_rows.append(
                    (
                        variant_id,
                        pos,
                        mod_idx,
                        unimod_variable_map.get(mod_idx)
                        if unimod_variable_map is not None
                        else None,
                    )
                )

    # Sort variants by mass then sequence for reproducibility
    peptide_variant_rows.sort(key=lambda row: (row[3], row[2], row[6]))

    return {
        "protein": protein_rows,
        "peptide_sequence": peptide_sequence_rows,
        "peptide_sequence_protein": peptide_sequence_protein_rows,
        "peptide_protein_location": peptide_protein_location_rows,
        "peptide_variant": peptide_variant_rows,
        "peptide_variant_mod": peptide_variant_mod_rows,
        "static_mod": static_mod_rows,
        "variable_mod": variable_mod_rows,
    }


__all__ = ["build_tables"]
