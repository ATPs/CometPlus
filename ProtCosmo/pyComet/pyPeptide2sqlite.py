#!/usr/bin/env python3
"""Convert pyPeptideIndex peptide_sequence TSV to SQLite."""

from __future__ import annotations

import argparse
import os
import sqlite3
from typing import Iterator, Optional


def _resolve_output_path(input_path: str, output_path: Optional[str]) -> str:
    if output_path:
        return output_path
    return input_path + ".sqlite"


def _iter_pep_seq(tsv_path: str) -> Iterator[str]:
    pep_seq_index: Optional[int] = None

    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as handle:
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
                pep_seq_index = 2 if len(fields) >= 3 else 0

            if pep_seq_index >= len(fields):
                continue

            pep_seq = fields[pep_seq_index].strip().upper()
            if pep_seq and pep_seq != "\\N":
                yield pep_seq


def tsv_to_sqlite(input_path: str, output_path: str, batch_size: int = 10000) -> None:
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"input file does not exist: {input_path}")

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if os.path.exists(output_path):
        os.remove(output_path)

    conn = sqlite3.connect(output_path)
    cur = conn.cursor()

    try:
        cur.execute("PRAGMA journal_mode=OFF;")
        cur.execute("PRAGMA synchronous=OFF;")
        cur.execute("PRAGMA temp_store=MEMORY;")

        cur.execute(
            "CREATE TABLE IF NOT EXISTS peptides (pep_seq TEXT PRIMARY KEY) WITHOUT ROWID;"
        )

        cur.execute("BEGIN;")

        batch = []
        total_rows = 0
        for pep_seq in _iter_pep_seq(input_path):
            batch.append((pep_seq,))
            total_rows += 1
            if len(batch) >= batch_size:
                cur.executemany(
                    "INSERT OR IGNORE INTO peptides(pep_seq) VALUES (?);",
                    batch,
                )
                batch.clear()

        if batch:
            cur.executemany(
                "INSERT OR IGNORE INTO peptides(pep_seq) VALUES (?);",
                batch,
            )

        conn.commit()

        cur.execute("SELECT COUNT(*) FROM peptides;")
        unique_rows = cur.fetchone()[0]
        print(f"input={input_path}")
        print(f"output={output_path}")
        print(f"rows_read={total_rows}")
        print(f"rows_inserted={unique_rows}")
    finally:
        conn.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert a pyPeptideIndex peptide_sequence TSV file to SQLite. "
            "The output table is peptides(pep_seq TEXT PRIMARY KEY)."
        )
    )
    parser.add_argument("--input", required=True, help="Input peptide_sequence TSV file.")
    parser.add_argument(
        "--output",
        default=None,
        help="Output SQLite file path (default: <input>.sqlite).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Number of rows per executemany batch.",
    )
    args = parser.parse_args()

    if args.batch_size <= 0:
        parser.error("--batch-size must be > 0")

    output_path = _resolve_output_path(args.input, args.output)
    tsv_to_sqlite(args.input, output_path, batch_size=args.batch_size)


if __name__ == "__main__":
    main()
