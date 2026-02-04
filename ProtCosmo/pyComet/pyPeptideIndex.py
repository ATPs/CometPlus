#!/usr/bin/env python3
"""Generate Comet-like peptide index tables from FASTA + comet.params.

This script parses comet.params, digests FASTA sequences, enumerates variable
mods, and writes TSV tables that mirror the suggested relational layout in
study/01-db-creation-suggestion.md.
"""

from __future__ import annotations

import argparse
import datetime
import os
import sys
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Tuple, Sequence

sys.path.append(os.path.dirname(__file__))

from pyLoadParameters import CometParams, EnzymeDefinition, VarMod, parse_comet_params


PROTON_MASS = 1.00727646688
HYDROGEN_MONO = 1.00782503223
OXYGEN_MONO = 15.99491461957


class CometHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def build_aa_masses(mono: bool) -> Tuple[Dict[str, float], float]:
    if mono:
        h = HYDROGEN_MONO
        o = OXYGEN_MONO
        c = 12.0
        n = 14.00307400443
        s = 31.9720711744
        se = 79.9165218
    else:
        h = 1.00794
        o = 15.9994
        c = 12.0107
        n = 14.0067
        s = 32.065
        se = 78.96

    h2o = h + h + o
    masses = {
        "G": c * 2 + h * 3 + n + o,
        "A": c * 3 + h * 5 + n + o,
        "S": c * 3 + h * 5 + n + o * 2,
        "P": c * 5 + h * 7 + n + o,
        "V": c * 5 + h * 9 + n + o,
        "T": c * 4 + h * 7 + n + o * 2,
        "C": c * 3 + h * 5 + n + o + s,
        "U": c * 3 + h * 5 + n + o + se,
        "L": c * 6 + h * 11 + n + o,
        "I": c * 6 + h * 11 + n + o,
        "N": c * 4 + h * 6 + n * 2 + o * 2,
        "D": c * 4 + h * 5 + n + o * 3,
        "Q": c * 5 + h * 8 + n * 2 + o * 2,
        "K": c * 6 + h * 12 + n * 2 + o,
        "E": c * 5 + h * 7 + n + o * 3,
        "M": c * 5 + h * 9 + n + o + s,
        "H": c * 6 + h * 7 + n * 3 + o,
        "F": c * 9 + h * 9 + n + o,
        "R": c * 6 + h * 12 + n * 4 + o,
        "Y": c * 9 + h * 9 + n + o * 2,
        "W": c * 11 + h * 10 + n * 2 + o,
        "O": c * 12 + h * 19 + n * 3 + o * 2,
        "B": 0.0,
        "J": 0.0,
        "X": 0.0,
        "Z": 0.0,
    }
    return masses, h2o


def apply_set_masses(masses: Dict[str, float], params: CometParams) -> None:
    for key, value in params.parsed.items():
        if not key.startswith("set_"):
            continue
        if not isinstance(value, (int, float)):
            continue
        if abs(value) < 1e-12:
            continue
        residue = key.split("_", 2)[1].upper()[:1]
        if residue in masses:
            masses[residue] = float(value)


def get_static_mods(params: CometParams) -> Tuple[Dict[str, float], Dict[str, float]]:
    residue_mods: Dict[str, float] = {}
    term_mods = {
        "Nterm_peptide": float(params.parsed.get("add_Nterm_peptide", 0.0)),
        "Cterm_peptide": float(params.parsed.get("add_Cterm_peptide", 0.0)),
        "Nterm_protein": float(params.parsed.get("add_Nterm_protein", 0.0)),
        "Cterm_protein": float(params.parsed.get("add_Cterm_protein", 0.0)),
    }
    for key, value in params.parsed.items():
        if not key.startswith("add_"):
            continue
        if key in (
            "add_Nterm_peptide",
            "add_Cterm_peptide",
            "add_Nterm_protein",
            "add_Cterm_protein",
        ):
            continue
        if not isinstance(value, (int, float)):
            continue
        if abs(value) < 1e-12:
            continue
        residue = key.split("_", 2)[1].upper()[:1]
        residue_mods[residue] = float(value)
    return residue_mods, term_mods


def iter_fasta(path: str) -> Iterable[Tuple[str, str, int]]:
    header = None
    seq_chunks: List[str] = []
    offset = 0
    with open(path, "rb") as handle:
        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                if header is not None:
                    yield header, "".join(seq_chunks), offset
                header = line[1:].decode("utf-8", errors="replace").strip()
                header = header.replace("\t", " ").replace("\r", " ").replace("\n", " ")
                seq_chunks = []
                offset = pos
            else:
                seq_chunks.append(line.decode("utf-8", errors="replace").strip())
        if header is not None:
            yield header, "".join(seq_chunks), offset


def iter_protein_sequences(sequences: Sequence[str]) -> Iterable[Tuple[str, str, int]]:
    offset = 0
    for index, seq in enumerate(sequences, start=1):
        header = f"protein_{index}"
        yield header, seq, offset
        offset += len(seq) + 1


def normalize_sequence(seq: str) -> str:
    seq = seq.upper()
    seq = "".join(ch for ch in seq if ch.isalpha() or ch == "*")
    return seq


def is_cleavage_site(seq: str, index: int, enzyme: EnzymeDefinition) -> bool:
    if enzyme.is_no_enzyme():
        return True
    if index < 0 or index >= len(seq):
        return False
    c_current = seq[index]
    if enzyme.offset == 0:
        flank = seq[index - 1] if index - 1 >= 0 else None
    else:
        flank = seq[index + 1] if index + 1 < len(seq) else None
    if c_current not in enzyme.break_aa:
        return False
    if flank is None:
        return True
    if flank in enzyme.no_break_aa:
        return False
    return True


def check_enzyme_termini(
    seq: str,
    start: int,
    end: int,
    enzyme: EnzymeDefinition,
    enzyme2: Optional[EnzymeDefinition],
    num_termini: int,
) -> bool:
    if enzyme.is_no_enzyme() and (enzyme2 is None or enzyme2.is_no_enzyme()):
        return True

    def _start_cleavage(e: EnzymeDefinition) -> bool:
        if start == 0 or seq[start - 1] == "*":
            return True
        idx = start - e.offset
        flank = start - 1 + e.offset
        if idx < 0 or flank < 0 or flank >= len(seq):
            return False
        return seq[idx] in e.break_aa and seq[flank] not in e.no_break_aa

    def _end_cleavage(e: EnzymeDefinition) -> bool:
        if end == len(seq) - 1 or seq[end + 1] == "*":
            return True
        idx = end + 1 - e.offset
        flank = end + e.offset
        if idx < 0 or flank < 0 or flank >= len(seq):
            return False
        return seq[idx] in e.break_aa and seq[flank] not in e.no_break_aa

    begin_ok = _start_cleavage(enzyme)
    end_ok = _end_cleavage(enzyme)

    if enzyme2 and not enzyme2.is_no_enzyme():
        if not begin_ok:
            begin_ok = _start_cleavage(enzyme2)
        if not end_ok:
            end_ok = _end_cleavage(enzyme2)

    if num_termini == 2:
        return begin_ok and end_ok
    if num_termini == 1:
        return begin_ok or end_ok
    if num_termini == 8:
        return begin_ok
    if num_termini == 9:
        return end_ok
    return True


def count_missed_cleavages(
    seq: str,
    start: int,
    end: int,
    enzyme: EnzymeDefinition,
    enzyme2: Optional[EnzymeDefinition],
) -> int:
    if enzyme.is_no_enzyme() and (enzyme2 is None or enzyme2.is_no_enzyme()):
        return 0

    if enzyme.offset == 0:
        begin_ref = start + 1
        end_ref = end
    else:
        begin_ref = start
        end_ref = end - 1

    missed = 0
    for i in range(begin_ref, end_ref + 1):
        if i < 0 or i >= len(seq):
            continue
        if is_cleavage_site(seq, i, enzyme):
            if (enzyme.offset == 1 and i != end) or (enzyme.offset == 0 and i != start):
                missed += 1
                continue
        if enzyme2 and not enzyme2.is_no_enzyme() and is_cleavage_site(seq, i, enzyme2):
            if (enzyme2.offset == 1 and i != end) or (enzyme2.offset == 0 and i != start):
                missed += 1
    return missed


def cleavage_positions(seq: str, enzyme: EnzymeDefinition) -> List[int]:
    if enzyme.is_no_enzyme():
        return list(range(0, len(seq) + 1))
    cuts = [0]
    if enzyme.offset == 1:
        for i in range(len(seq) - 1):
            if is_cleavage_site(seq, i, enzyme):
                cuts.append(i + 1)
    else:
        for i in range(1, len(seq)):
            if is_cleavage_site(seq, i, enzyme):
                cuts.append(i)
    cuts.append(len(seq))
    return sorted(set(cuts))


def combined_cleavage_positions(
    seq: str, enzyme: EnzymeDefinition, enzyme2: Optional[EnzymeDefinition]
) -> List[int]:
    cuts = set(cleavage_positions(seq, enzyme))
    if enzyme2 and not enzyme2.is_no_enzyme():
        cuts.update(cleavage_positions(seq, enzyme2))
    return sorted(cuts)


def iter_peptides(
    seq: str,
    enzyme: EnzymeDefinition,
    enzyme2: Optional[EnzymeDefinition],
    num_termini: int,
    max_missed: int,
    min_len: int,
    max_len: int,
) -> Iterable[Tuple[int, int]]:
    if enzyme.is_no_enzyme() and (enzyme2 is None or enzyme2.is_no_enzyme()):
        for start in range(len(seq)):
            for end in range(start, min(len(seq), start + max_len)):
                length = end - start + 1
                if length < min_len:
                    continue
                yield start, end
        return

    if num_termini == 2:
        cuts = combined_cleavage_positions(seq, enzyme, enzyme2)
        for i, start in enumerate(cuts[:-1]):
            for j in range(i + 1, min(len(cuts), i + max_missed + 2)):
                end = cuts[j] - 1
                if end < start:
                    continue
                length = end - start + 1
                if length < min_len or length > max_len:
                    continue
                yield start, end
        return

    for start in range(len(seq)):
        for end in range(start, min(len(seq), start + max_len)):
            length = end - start + 1
            if length < min_len:
                continue
            if not check_enzyme_termini(seq, start, end, enzyme, enzyme2, num_termini):
                continue
            missed = count_missed_cleavages(seq, start, end, enzyme, enzyme2)
            if missed > max_missed:
                continue
            yield start, end


def peptide_flanks(seq: str, start: int, end: int) -> Tuple[str, str]:
    prev_aa = "-" if start == 0 else seq[start - 1]
    next_aa = "-" if end == len(seq) - 1 else seq[end + 1]
    if prev_aa == "*":
        prev_aa = "-"
    if next_aa == "*":
        next_aa = "-"
    return prev_aa, next_aa


def build_var_mod_candidates(
    seq: str,
    var_mods: List[VarMod],
    protein_start: int,
    protein_length: int,
    nterm_offset: int,
) -> List[List[int]]:
    candidates: List[List[int]] = []
    pep_len = len(seq)
    protein_end = protein_start + pep_len - 1
    protein_nterm = protein_start == nterm_offset
    protein_cterm = protein_end == protein_length - 1

    for mod in var_mods:
        positions: List[int] = []
        residues = mod.residues
        residue_chars = residues.replace("n", "").replace("c", "")
        for i, residue in enumerate(seq):
            if residue not in residue_chars:
                continue
            if mod.term_distance >= 0:
                if mod.which_term == 0 and (protein_start + i - nterm_offset) > mod.term_distance:
                    continue
                if mod.which_term == 1 and (protein_length - 1 - (protein_start + i)) > mod.term_distance:
                    continue
                if mod.which_term == 2 and i > mod.term_distance:
                    continue
                if mod.which_term == 3 and (pep_len - 1 - i) > mod.term_distance:
                    continue
            positions.append(i)

        if "n" in residues:
            ok = True
            if mod.which_term == 0 and not protein_nterm:
                ok = False
            if mod.which_term == 1 and not protein_cterm:
                ok = False
            if mod.which_term == 2 and not True:
                ok = False
            if mod.which_term == 3 and not True:
                ok = False
            if ok:
                positions.append(pep_len)  # N-term index

        if "c" in residues:
            ok = True
            if mod.which_term == 0 and not protein_nterm:
                ok = False
            if mod.which_term == 1 and not protein_cterm:
                ok = False
            if mod.which_term == 2 and not True:
                ok = False
            if mod.which_term == 3 and not True:
                ok = False
            if ok:
                positions.append(pep_len + 1)  # C-term index

        candidates.append(sorted(set(positions)))

    return candidates


def enumerate_var_mods(
    seq: str,
    var_mods: List[VarMod],
    candidates: List[List[int]],
    max_mods: int,
) -> List[Tuple[List[int], int, float]]:
    pep_len = len(seq)
    sites = [0] * (pep_len + 2)
    results: List[Tuple[List[int], int, float]] = []
    counts = [0] * len(var_mods)

    def backtrack(idx: int, total_mods: int, mass_delta: float) -> None:
        if idx == len(var_mods):
            for i, mod in enumerate(var_mods):
                if counts[i] < mod.min_per_pep:
                    return
                if mod.require_this_mod > 0 and counts[i] == 0:
                    return
            results.append((sites.copy(), total_mods, mass_delta))
            return

        mod = var_mods[idx]
        positions = candidates[idx]
        if mod.mass == 0.0 or mod.residues == "-":
            backtrack(idx + 1, total_mods, mass_delta)
            return

        max_for_mod = min(mod.max_per_pep, len(positions))
        if mod.binary_mod:
            max_for_mod = min(max_for_mod, 1)
        min_for_mod = min(mod.min_per_pep, max_for_mod)

        for k in range(min_for_mod, max_for_mod + 1):
            if total_mods + k > max_mods:
                break
            if k == 0:
                counts[idx] = 0
                backtrack(idx + 1, total_mods, mass_delta)
                continue
            for combo in combinations(positions, k):
                if any(sites[pos] != 0 for pos in combo):
                    continue
                for pos in combo:
                    sites[pos] = mod.index
                counts[idx] = k
                backtrack(idx + 1, total_mods + k, mass_delta + mod.mass * k)
                for pos in combo:
                    sites[pos] = 0
        counts[idx] = 0

    backtrack(0, 0, 0.0)
    return results


def format_var_mod_sites(sites: List[int]) -> str:
    pairs = []
    for pos, mod_idx in enumerate(sites):
        if mod_idx:
            pairs.append(f"{pos}:{mod_idx}")
    return ";".join(pairs)


def compute_peptide_mass(
    seq: str,
    aa_masses: Dict[str, float],
    h2o_mass: float,
    residue_mods: Dict[str, float],
    term_mods: Dict[str, float],
    is_protein_nterm: bool,
    is_protein_cterm: bool,
) -> Optional[float]:
    mass = h2o_mass + PROTON_MASS
    for residue in seq:
        base = aa_masses.get(residue)
        if base is None or base == 0.0:
            return None
        mass += base + residue_mods.get(residue, 0.0)
    mass += term_mods["Nterm_peptide"] + term_mods["Cterm_peptide"]
    if is_protein_nterm:
        mass += term_mods["Nterm_protein"]
    if is_protein_cterm:
        mass += term_mods["Cterm_protein"]
    return mass


def build_tables(
    params: CometParams,
    fasta_paths: Optional[List[str]] = None,
    protein_sequences: Optional[List[str]] = None,
    max_proteins: Optional[int] = None,
    progress: bool = False,
    progress_every: int = 1000,
    use_protein_name: bool = False,
) -> Dict[str, List[Tuple]]:
    residue_mods, term_mods = get_static_mods(params)
    aa_masses, h2o_mass = build_aa_masses(bool(params.parsed.get("mass_type_parent", 1)))
    apply_set_masses(aa_masses, params)
    var_mods = sorted(params.variable_mods, key=lambda mod: mod.index)

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
        if progress:
            print(message, file=sys.stderr, flush=True)

    def record_peptide(
        pep_seq: str,
        abs_start: int,
        abs_end: int,
        prev_aa: str,
        next_aa: str,
        protein_len: int,
        nterm_offset: int,
        protein_id: int,
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
            }
            seq_info[pep_seq] = entry
        entry["protein_ids"].add(protein_id)
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
        protein_rows.append((protein_id, offset, header))
        if not seq:
            continue
        protein_starts_with_m = seq.startswith("M")
        pos = 0
        for segment in seq.split("*"):
            segment_len = len(segment)
            if segment_len == 0:
                pos += 1
                continue
            segment_start = pos
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
                prev_aa, next_aa = peptide_flanks(seq, abs_start, abs_end)
                record_peptide(
                    pep_seq,
                    abs_start,
                    abs_end,
                    prev_aa,
                    next_aa,
                    len(seq),
                    0,
                    protein_id,
                    offset,
                )
            pos += segment_len + 1

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
                    _, next_aa = peptide_flanks(seq, abs_start, abs_end)
                    record_peptide(
                        pep_seq,
                        abs_start,
                        abs_end,
                        "-",
                        next_aa,
                        len(seq),
                        1,
                        protein_id,
                        offset,
                    )

        if progress and protein_count % progress_every == 0:
            progress_log(
                f"Index: processed {protein_count} proteins; "
                f"peptides seen {peptide_records}; unique sequences {len(seq_info)}"
            )

    progress_log(
        f"Index: finished proteins {protein_count}; "
        f"peptides seen {peptide_records}; unique sequences {len(seq_info)}"
    )

    # Assign sequence IDs deterministically
    sequences = sorted(seq_info.keys())
    seq_id_map: Dict[str, int] = {seq: idx + 1 for idx, seq in enumerate(sequences)}

    peptide_sequence_rows: List[Tuple] = []
    peptide_sequence_protein_rows: List[Tuple] = []

    for seq in sequences:
        entry = seq_info[seq]
        _, primary_protein_id, *_ = entry["primary"]
        seq_id = seq_id_map[seq]
        peptide_sequence_rows.append((seq_id, seq, len(seq), primary_protein_id))
        for pid in sorted(entry["protein_ids"]):
            peptide_sequence_protein_rows.append((seq_id, pid))

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
            )
        )

    static_mod_rows: List[Tuple] = []
    for residue, delta in sorted(residue_mods.items()):
        static_mod_rows.append((residue, delta, "residue"))
    if abs(term_mods["Nterm_peptide"]) > 1e-12:
        static_mod_rows.append(("-", term_mods["Nterm_peptide"], "N-term"))
    if abs(term_mods["Cterm_peptide"]) > 1e-12:
        static_mod_rows.append(("-", term_mods["Cterm_peptide"], "C-term"))
    if abs(term_mods["Nterm_protein"]) > 1e-12:
        static_mod_rows.append(("-", term_mods["Nterm_protein"], "protein N-term"))
    if abs(term_mods["Cterm_protein"]) > 1e-12:
        static_mod_rows.append(("-", term_mods["Cterm_protein"], "protein C-term"))

    peptide_variant_rows: List[Tuple] = []
    peptide_variant_mod_rows: List[Tuple] = []
    variant_id = 0

    total_sequences = len(sequences)
    if total_sequences:
        progress_log(f"Index: enumerating variants for {total_sequences} sequences")

    for seq_idx, seq in enumerate(sequences, start=1):
        entry = seq_info[seq]
        (
            offset,
            primary_protein_id,
            start,
            end,
            prev_aa,
            next_aa,
            protein_len,
            nterm_offset,
        ) = entry["primary"]
        is_protein_nterm = start == nterm_offset
        is_protein_cterm = end == protein_len - 1

        base_mass = compute_peptide_mass(
            seq,
            aa_masses,
            h2o_mass,
            residue_mods,
            term_mods,
            is_protein_nterm,
            is_protein_cterm,
        )
        if base_mass is None:
            continue

        candidates = build_var_mod_candidates(
            seq,
            var_mods,
            start,
            protein_len,
            nterm_offset,
        )
        mod_combos = enumerate_var_mods(seq, var_mods, candidates, max_var_mods)
        for sites, total_mods, delta_mass in mod_combos:
            if require_var_mod and total_mods == 0:
                continue
            mh_plus = base_mass + delta_mass
            if mh_plus < min_mass or mh_plus > max_mass:
                continue
            variant_id += 1
            mass_bin10 = int(mh_plus * 10.0)
            sites_text = format_var_mod_sites(sites)
            seq_id = seq_id_map[seq]
            peptide_variant_rows.append(
                (
                    variant_id,
                    seq_id,
                    mh_plus,
                    prev_aa,
                    next_aa,
                    sites_text,
                    total_mods,
                    mass_bin10,
                )
            )
            for pos, mod_idx in enumerate(sites):
                if mod_idx:
                    peptide_variant_mod_rows.append((variant_id, pos, mod_idx))

        if progress and seq_idx % progress_every == 0:
            progress_log(
                f"Index: variants for {seq_idx}/{total_sequences} sequences; "
                f"variants {variant_id}"
            )

    # Sort variants by mass then sequence for reproducibility
    peptide_variant_rows.sort(key=lambda row: (row[2], sequences[row[1] - 1], row[5]))

    return {
        "protein": protein_rows,
        "peptide_sequence": peptide_sequence_rows,
        "peptide_sequence_protein": peptide_sequence_protein_rows,
        "peptide_variant": peptide_variant_rows,
        "peptide_variant_mod": peptide_variant_mod_rows,
        "static_mod": static_mod_rows,
        "variable_mod": variable_mod_rows,
    }


def generate_peptide_index_tables(
    params_path: str,
    fasta_paths: Optional[List[str]] = None,
    protein_sequences: Optional[List[str]] = None,
    max_proteins: Optional[int] = None,
    progress: bool = False,
    use_protein_name: bool = False,
) -> Tuple[Dict[str, List[Tuple]], Tuple]:
    params = parse_comet_params(params_path)
    tables = build_tables(
        params,
        fasta_paths=fasta_paths,
        protein_sequences=protein_sequences,
        max_proteins=max_proteins,
        progress=progress,
        use_protein_name=use_protein_name,
    )
    run_id = 1
    created_at = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    database_label = " ".join(fasta_paths or []) if fasta_paths else "<inline_sequences>"
    index_run_row = (
        run_id,
        params.version,
        params.params_path,
        database_label,
        float(params.parsed.get("digest_mass_range", (0.0, 0.0))[0]),
        float(params.parsed.get("digest_mass_range", (0.0, 0.0))[1]),
        int(params.parsed.get("peptide_length_range", (0, 0))[0]),
        int(params.parsed.get("peptide_length_range", (0, 0))[1]),
        created_at,
        params.to_json(),
    )
    return tables, index_run_row


def format_field(value: Optional[object]) -> str:
    if value is None:
        return "\\N"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_tsv(path: str, rows: List[Tuple]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        for row in rows:
            handle.write("\t".join(format_field(value) for value in row) + "\n")


def main() -> None:
    description = (
        "Generate TSV tables that mirror Comet's peptide index layout.\n"
        "The parser follows Comet.cpp::LoadParameters rules (comments at '#',\n"
        "variable_modNN requires 8 fields, and [COMET_ENZYME_INFO] is parsed at EOF).\n\n"
        "Key comet.params settings used here:\n"
        "  database_name             FASTA path(s) unless overridden by -D/--database\n"
        "  search_enzyme_number      selects enzyme row in [COMET_ENZYME_INFO]\n"
        "  search_enzyme2_number     optional second enzyme (0 disables)\n"
        "  num_enzyme_termini        1=semi, 2=fully (default), 8 N-term only, 9 C-term only\n"
        "  allowed_missed_cleavage   max internal cleavage sites (default 2)\n"
        "  peptide_length_range      peptide length filter (default 5..50)\n"
        "  digest_mass_range         MH+ mass filter (default 600..5000)\n"
        "  max_variable_mods_in_peptide  total var mods per peptide (default 5)\n"
        "  require_variable_mod      1=keep only modified variants\n"
        "  clip_nterm_methionine     1=allow peptides starting at position 2 if protein starts with M\n"
        "  mass_type_parent          1=mono (default), 0=average\n"
        "  add_* / add_Nterm_*        static mods\n"
        "  variable_modNN            <mass> <residues> <binary> <min|max> <term_dist> <which_term> <required> <neutral_loss>\n"
        "                           - binary=1 enforces at most one site per peptide\n"
        "                           - term_dist=-1 disables distance filtering\n"
        "                           - which_term: 0=protein N, 1=protein C, 2=peptide N, 3=peptide C\n\n"
        "Example:\n"
        "  python pyPeptideIndex.py -P comet.params -D db.fasta -N out/comet\n"
        "  python pyPeptideIndex.py -P comet.params --protein ACDEFGHIK -N out/comet\n\n"
        "Outputs (TSV, no headers):\n"
        "  <prefix>.index_run.tsv\n"
        "  <prefix>.static_mod.tsv\n"
        "  <prefix>.variable_mod.tsv\n"
        "  <prefix>.protein.tsv\n"
        "  <prefix>.peptide_sequence.tsv\n"
        "  <prefix>.peptide_sequence_protein.tsv\n"
        "  <prefix>.peptide_variant.tsv\n"
        "  <prefix>.peptide_variant_mod.tsv\n"
    )
    epilog = (
        "Examples:\n"
        "  python pyPeptideIndex.py --params comet.params\n"
        "  python pyPeptideIndex.py -P comet.params -D db.fasta -N out/comet\n"
        "  python pyPeptideIndex.py -D db1.fasta -D db2.fasta --prefix results/comet\n"
    )
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=CometHelpFormatter,
    )
    parser.add_argument(
        "-P",
        "--params",
        default="comet.params",
        help=(
            "Comet params file. Used for database_name, enzyme rules, peptide_length_range,\n"
            "digest_mass_range, static/variable mods, mass_type_parent, and related settings."
        ),
    )
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-D",
        "--database",
        action="append",
        nargs="+",
        help=(
            "FASTA file(s). Overrides database_name from the params file.\n"
            "Repeat -D or pass multiple paths in one -D to include multiple FASTA files."
        ),
    )
    input_group.add_argument(
        "--protein",
        action="append",
        nargs="+",
        help=(
            "Protein sequence(s) supplied directly on the command line.\n"
            "Repeat --protein or pass multiple sequences in one --protein.\n"
            "Synthetic headers are assigned as protein_1, protein_2, ... and file offsets are synthetic."
        ),
    )
    parser.add_argument(
        "-N",
        "--prefix",
        default="comet",
        help="Output basename (can include a directory).",
    )
    parser.add_argument(
        "--max-record",
        type=int,
        default=None,
        help="Maximum number of proteins to process (default: all).",
    )
    parser.add_argument(
        "--use-protein-name",
        action="store_true",
        help="Use the protein name (first token in FASTA header) as protein_id.",
    )
    args = parser.parse_args()

    params = parse_comet_params(args.params)
    fasta_paths: List[str] = []
    protein_sequences: List[str] = []
    if args.database:
        for group in args.database:
            fasta_paths.extend(group)
    elif args.protein:
        for group in args.protein:
            protein_sequences.extend(group)
    else:
        db_name = params.parsed.get("database_name", "")
        if db_name:
            fasta_paths = db_name.split()

    if not fasta_paths and not protein_sequences:
        raise SystemExit(
            "No FASTA database or protein sequences specified; use --database, --protein, or set database_name in params."
        )

    output_dir = os.path.dirname(args.prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    tables = build_tables(
        params,
        fasta_paths=fasta_paths,
        protein_sequences=protein_sequences,
        max_proteins=args.max_record,
        progress=True,
        use_protein_name=args.use_protein_name,
    )

    run_id = 1
    created_at = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    database_label = " ".join(fasta_paths) if fasta_paths else "<inline_sequences>"
    index_run_row = (
        run_id,
        params.version,
        params.params_path,
        database_label,
        float(params.parsed.get("digest_mass_range", (0.0, 0.0))[0]),
        float(params.parsed.get("digest_mass_range", (0.0, 0.0))[1]),
        int(params.parsed.get("peptide_length_range", (0, 0))[0]),
        int(params.parsed.get("peptide_length_range", (0, 0))[1]),
        created_at,
        params.to_json(),
    )

    write_tsv(f"{args.prefix}.index_run.tsv", [index_run_row])
    write_tsv(
        f"{args.prefix}.static_mod.tsv",
        [(run_id, *row) for row in tables["static_mod"]],
    )
    write_tsv(
        f"{args.prefix}.variable_mod.tsv",
        [(run_id, *row) for row in tables["variable_mod"]],
    )
    write_tsv(
        f"{args.prefix}.protein.tsv",
        [(run_id, *row) for row in tables["protein"]],
    )
    write_tsv(
        f"{args.prefix}.peptide_sequence.tsv",
        [(run_id, *row) for row in tables["peptide_sequence"]],
    )
    write_tsv(
        f"{args.prefix}.peptide_sequence_protein.tsv",
        [(run_id, *row) for row in tables["peptide_sequence_protein"]],
    )
    write_tsv(
        f"{args.prefix}.peptide_variant.tsv",
        [(run_id, *row) for row in tables["peptide_variant"]],
    )
    write_tsv(
        f"{args.prefix}.peptide_variant_mod.tsv",
        [(run_id, *row) for row in tables["peptide_variant_mod"]],
    )


if __name__ == "__main__":
    main()


__all__ = [
    "build_tables",
    "generate_peptide_index_tables",
    "write_tsv",
    "parse_comet_params",
]
