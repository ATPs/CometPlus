from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple

from .pyLoadParameters import EnzymeDefinition


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


__all__ = [
    "iter_fasta",
    "iter_protein_sequences",
    "normalize_sequence",
    "is_cleavage_site",
    "check_enzyme_termini",
    "count_missed_cleavages",
    "cleavage_positions",
    "combined_cleavage_positions",
    "iter_peptides",
    "peptide_flanks",
]

