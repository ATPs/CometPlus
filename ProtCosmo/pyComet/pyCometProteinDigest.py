#!/usr/bin/env python3
"""Digest proteins and/or peptide inputs into internal novel peptide TSV rows."""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

sys.path.append(os.path.dirname(__file__))

from utils.peptide_index_digestion import iter_fasta, iter_peptides, normalize_sequence
from utils.pyLoadParameters import CometParams, parse_comet_params


def normalize_peptide_token(text: str) -> str:
    return "".join(ch.upper() for ch in text if ch.isalpha())


def normalize_header(text: str, fallback: str) -> str:
    value = text.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()
    if not value:
        return fallback
    return " ".join(value.split())


def first_header_token(text: str, fallback: str) -> str:
    normalized = normalize_header(text, fallback)
    parts = normalized.split()
    if not parts:
        return fallback
    return parts[0]


def flatten_groups(groups: Optional[Sequence[Sequence[str]]]) -> List[str]:
    values: List[str] = []
    for group in groups or []:
        values.extend(group)
    return values


def iter_tokens_from_text(text: str) -> Iterator[str]:
    chars: List[str] = []
    for ch in text:
        if ch == "," or ch.isspace():
            chars.append(" ")
        else:
            chars.append(ch)
    for token in "".join(chars).split():
        peptide = normalize_peptide_token(token)
        if peptide:
            yield peptide


def is_fasta_like_file(path: str) -> bool:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                return True
    return False


def iter_peptides_from_fasta_file(path: str) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    seq_chunks: List[str] = []

    def flush_current() -> Iterator[Tuple[str, str]]:
        if header is None:
            return
        sequence = normalize_peptide_token("".join(seq_chunks))
        if sequence:
            yield sequence, normalize_header(header, sequence)

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                for row in flush_current():
                    yield row
                header = line[1:]
                seq_chunks = []
            else:
                if header is None:
                    continue
                seq_chunks.append(line)

    for row in flush_current():
        yield row


def iter_peptides_from_text_file(path: str) -> Iterator[Tuple[str, str]]:
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            for peptide in iter_tokens_from_text(line):
                yield peptide, peptide


def iter_peptides_from_input(value: str) -> Iterator[Tuple[str, str]]:
    if os.path.isfile(value):
        if is_fasta_like_file(value):
            yield from iter_peptides_from_fasta_file(value)
        else:
            yield from iter_peptides_from_text_file(value)
        return

    for peptide in iter_tokens_from_text(value):
        yield peptide, peptide


def iter_database_peptides(
    fasta_paths: Sequence[str],
    params: CometParams,
    max_record: Optional[int],
) -> Iterator[Tuple[str, str]]:
    enzyme = params.enzymes["search"]
    enzyme2 = params.enzymes["search2"] if params.parsed.get("search_enzyme2_number", 0) else None
    num_termini = int(params.parsed.get("num_enzyme_termini", 2))
    max_missed = int(params.parsed.get("allowed_missed_cleavage", 2))
    min_len, max_len = params.parsed.get("peptide_length_range", (5, 50))
    clip_nterm_methionine = int(params.parsed.get("clip_nterm_methionine", 0))

    protein_count = 0
    for fasta_path in fasta_paths:
        for header, seq_raw, _offset in iter_fasta(fasta_path):
            if max_record is not None and protein_count >= max_record:
                return
            protein_count += 1

            protein_id = first_header_token(header, str(protein_count))
            sequence = normalize_sequence(seq_raw)
            if not sequence:
                continue

            protein_starts_with_m = sequence.startswith("M")
            for segment in sequence.split("*"):
                if not segment:
                    continue
                for start, end in iter_peptides(
                    segment,
                    enzyme,
                    enzyme2,
                    num_termini,
                    max_missed,
                    min_len,
                    max_len,
                ):
                    peptide = segment[start : end + 1]
                    if peptide:
                        yield peptide, protein_id

            if clip_nterm_methionine and protein_starts_with_m and len(sequence) > 1:
                clipped_seq = sequence[1:]
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
                        peptide = clipped_segment[start : end + 1]
                        if peptide:
                            yield peptide, protein_id


def looks_like_sqlite(path: str) -> bool:
    try:
        with open(path, "rb") as handle:
            magic = handle.read(16)
        return magic == b"SQLite format 3\x00"
    except OSError:
        return False


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def load_peptide_exclude_tsv(path: str) -> Set[str]:
    excluded: Set[str] = set()
    pep_seq_index: Optional[int] = None

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("#"):
                header = [field.strip() for field in line[1:].split("\t")]
                lowered = [field.lower() for field in header]
                if "pep_seq" in lowered:
                    pep_seq_index = lowered.index("pep_seq")
                continue

            fields = line.split("\t")
            if pep_seq_index is None:
                lowered = [field.strip().lower() for field in fields]
                if "pep_seq" in lowered:
                    pep_seq_index = lowered.index("pep_seq")
                    continue
                pep_seq_index = 0

            if pep_seq_index >= len(fields):
                continue

            peptide = normalize_peptide_token(fields[pep_seq_index].strip())
            if peptide:
                excluded.add(peptide)

    return excluded


def load_peptide_exclude_sqlite(path: str) -> Set[str]:
    excluded: Set[str] = set()
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    cur = conn.cursor()

    try:
        cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name;"
        )
        table_names = [row[0] for row in cur.fetchall()]
        if len(table_names) != 1:
            raise ValueError(
                f"peptide-exclude sqlite must contain exactly one table; found {len(table_names)}"
            )

        table_name = table_names[0]
        quoted_table = quote_identifier(table_name)

        cur.execute(f"PRAGMA table_info({quoted_table});")
        columns = [str(row[1]).lower() for row in cur.fetchall()]
        if "pep_seq" not in columns:
            raise ValueError(
                f"peptide-exclude sqlite table '{table_name}' must include column 'pep_seq'"
            )

        cur.execute(f"SELECT pep_seq FROM {quoted_table};")
        for row in cur.fetchall():
            peptide = normalize_peptide_token(str(row[0]) if row and row[0] is not None else "")
            if peptide:
                excluded.add(peptide)
    finally:
        conn.close()

    return excluded


def load_peptide_exclude(path: Optional[str]) -> Set[str]:
    if not path:
        return set()
    if not os.path.isfile(path):
        raise FileNotFoundError(f"peptide-exclude file does not exist: {path}")

    if path.lower().endswith(".sqlite") or looks_like_sqlite(path):
        return load_peptide_exclude_sqlite(path)
    return load_peptide_exclude_tsv(path)


def add_record(
    peptide_map: Dict[str, List[str]],
    peptide_protein_sets: Dict[str, Set[str]],
    peptide: str,
    protein_id: str,
    excluded: Set[str],
) -> bool:
    norm_peptide = normalize_peptide_token(peptide)
    if not norm_peptide:
        return False
    if norm_peptide in excluded:
        return False

    norm_protein_id = normalize_header(protein_id, norm_peptide)

    proteins = peptide_map.get(norm_peptide)
    if proteins is None:
        peptide_map[norm_peptide] = [norm_protein_id]
        peptide_protein_sets[norm_peptide] = {norm_protein_id}
        return True

    protein_set = peptide_protein_sets[norm_peptide]
    if norm_protein_id not in protein_set:
        protein_set.add(norm_protein_id)
        proteins.append(norm_protein_id)
    return True


def write_output(path: str, peptide_map: Dict[str, List[str]]) -> None:
    output_dir = os.path.dirname(path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("peptide\tpeptide_id\tprotein_id\n")
        next_id = 1
        for peptide, protein_ids in peptide_map.items():
            peptide_id = f"COMETPLUS_NOVEL_{next_id}"
            protein_field = ";".join(protein_ids)
            handle.write(f"{peptide}\t{peptide_id}\t{protein_field}\n")
            next_id += 1


def build_output_rows(args: argparse.Namespace) -> Dict[str, List[str]]:
    peptide_map: Dict[str, List[str]] = {}
    peptide_protein_sets: Dict[str, Set[str]] = {}

    excluded = load_peptide_exclude(args.peptide_exclude)

    database_paths = flatten_groups(args.database)
    if database_paths:
        params = parse_comet_params(args.params)
        for peptide, protein_id in iter_database_peptides(database_paths, params, args.max_record):
            add_record(peptide_map, peptide_protein_sets, peptide, protein_id, excluded)

    peptide_inputs = flatten_groups(args.peptide)
    for peptide_input in peptide_inputs:
        for peptide, protein_id in iter_peptides_from_input(peptide_input):
            add_record(peptide_map, peptide_protein_sets, peptide, protein_id, excluded)

    return peptide_map


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Digest proteins from FASTA using comet.params rules and/or import peptide inputs, "
            "then write a 3-column internal novel peptide TSV."
        )
    )
    parser.add_argument(
        "-D",
        "--database",
        action="append",
        nargs="+",
        default=None,
        help="FASTA database file(s). Repeat -D or pass multiple files in one -D.",
    )
    parser.add_argument(
        "--max-record",
        type=int,
        default=None,
        help="Maximum number of proteins to read from --database inputs.",
    )
    parser.add_argument(
        "-P",
        "--params",
        default="comet.params",
        help="Comet params file used for digestion settings.",
    )
    parser.add_argument(
        "--peptide",
        action="append",
        nargs="+",
        default=None,
        help=(
            "Peptide input item(s): single peptide, comma/whitespace token list, text file, or FASTA file. "
            "Repeat --peptide as needed."
        ),
    )
    parser.add_argument(
        "--peptide-exclude",
        default=None,
        help=(
            "Optional peptide exclusion file: TSV (column pep_seq) or SQLite (single table with column pep_seq)."
        ),
    )
    parser.add_argument(
        "-O",
        "--output",
        required=True,
        help="Output TSV path.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    has_database = bool(flatten_groups(args.database))
    has_peptide = bool(flatten_groups(args.peptide))
    if not has_database and not has_peptide:
        parser.error("At least one of --database or --peptide must be provided.")

    if args.max_record is not None and args.max_record <= 0:
        parser.error("--max-record must be > 0")

    peptide_map = build_output_rows(args)
    write_output(args.output, peptide_map)

    print(f"output={args.output}")
    print(f"unique_peptides={len(peptide_map)}")


if __name__ == "__main__":
    main()
