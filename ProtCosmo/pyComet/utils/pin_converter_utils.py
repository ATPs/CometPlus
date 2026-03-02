from __future__ import annotations

import csv
import math
import re
import sys
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from .unimod_matcher import normalize_unimod_id


_FLOAT_TOKEN_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


@dataclass(frozen=True)
class FlexUnimodEntry:
    mono_mass_delta: float
    residues: str
    unimod_id: str
    unimod_name: str
    row_index: int


class WarningSink:
    """De-duplicate warning messages and print once to stderr."""

    def __init__(self) -> None:
        self._seen: set[str] = set()

    def warn(self, message: str) -> None:
        if message in self._seen:
            return
        self._seen.add(message)
        print(f"warning: {message}", file=sys.stderr)


def _detect_delimiter(unimod_path: str) -> str:
    with open(unimod_path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" in line and line.count("\t") >= line.count(","):
                return "\t"
            return ","
    return ","


def load_unimod_entries_flexible(unimod_path: str) -> List[FlexUnimodEntry]:
    delimiter = _detect_delimiter(unimod_path)
    entries: List[FlexUnimodEntry] = []
    with open(unimod_path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if not row:
                continue
            first = row[0].strip()
            if not first or first.startswith("#"):
                continue
            if len(row) < 3:
                continue
            try:
                mass_delta = float(first)
            except ValueError:
                continue
            residues = row[1].strip()
            unimod_id = normalize_unimod_id(row[2])
            if not residues or not unimod_id:
                continue
            unimod_name = " ".join(part.strip() for part in row[3:] if part.strip())
            entries.append(
                FlexUnimodEntry(
                    mono_mass_delta=mass_delta,
                    residues=residues,
                    unimod_id=unimod_id,
                    unimod_name=unimod_name,
                    row_index=len(entries),
                )
            )
    return entries


def _residue_tokens(value: str) -> set[str]:
    tokens: set[str] = set()
    for ch in value.strip():
        if ch == "n" or ch == "c":
            tokens.add(ch)
        elif ch.isalpha():
            tokens.add(ch.upper())
    return tokens


def _ppm_error(observed_mass: float, reference_mass: float) -> float:
    if abs(reference_mass) <= 1e-12:
        return math.inf
    return abs(observed_mass - reference_mass) / abs(reference_mass) * 1e6


def _select_best_unimod(
    entries: Sequence[FlexUnimodEntry],
    mass_delta: float,
    residue_code: str,
    ppm_tolerance: float,
) -> Optional[FlexUnimodEntry]:
    observed_tokens = _residue_tokens(residue_code)
    if not observed_tokens:
        return None

    candidates: List[Tuple[float, int, FlexUnimodEntry]] = []
    for entry in entries:
        entry_tokens = _residue_tokens(entry.residues)
        if not observed_tokens.issubset(entry_tokens):
            continue
        ppm = _ppm_error(mass_delta, entry.mono_mass_delta)
        candidates.append((ppm, entry.row_index, entry))

    if not candidates:
        return None

    candidates.sort(key=lambda item: (item[0], item[1]))
    best_ppm, _best_order, best_entry = candidates[0]
    if best_ppm <= ppm_tolerance:
        return best_entry
    return None


def _dot_positions_outside_brackets(text: str) -> List[int]:
    positions: List[int] = []
    depth = 0
    for idx, ch in enumerate(text):
        if ch == "[":
            depth += 1
        elif ch == "]":
            if depth > 0:
                depth -= 1
        elif ch == "." and depth == 0:
            positions.append(idx)
    return positions


def _strip_pin_flanks(peptide: str) -> str:
    dots = _dot_positions_outside_brackets(peptide)
    if len(dots) < 2:
        return peptide
    first = dots[0]
    second = dots[1]
    if second <= first:
        return peptide
    return peptide[first + 1 : second]


def _find_closing_bracket(text: str, open_pos: int) -> int:
    if open_pos < 0 or open_pos >= len(text) or text[open_pos] != "[":
        return -1
    return text.find("]", open_pos + 1)


def _normalize_mod_tag(
    mod_tag: str,
    residue_code: str,
    entries: Sequence[FlexUnimodEntry],
    ppm_tolerance: float,
    warning_sink: WarningSink,
    context: str,
) -> str:
    if not (mod_tag.startswith("[") and mod_tag.endswith("]")):
        warning_sink.warn(f"invalid mod token '{mod_tag}' in peptide '{context}'")
        return mod_tag

    mass_text = mod_tag[1:-1].strip()
    if not _FLOAT_TOKEN_RE.match(mass_text):
        warning_sink.warn(
            f"non-numeric mod mass '{mass_text}' for residue '{residue_code}' in peptide '{context}'"
        )
        return mod_tag

    mass_delta = float(mass_text)
    matched = _select_best_unimod(entries, mass_delta, residue_code, ppm_tolerance)
    if matched is None:
        warning_sink.warn(
            f"cannot map mass '{mass_text}' at residue '{residue_code}' in peptide '{context}'"
        )
        return mod_tag
    return f"[U:{matched.unimod_id}]"


def _extract_nterm_mod(core: str) -> Tuple[str, Optional[str], Optional[str]]:
    if not core.startswith("n["):
        return core, None, None
    close = _find_closing_bracket(core, 1)
    if close < 0:
        return core, None, None
    return core[close + 1 :], core[: close + 1], core[1 : close + 1]


def _extract_cterm_mod(core: str) -> Tuple[str, Optional[str], Optional[str]]:
    if not core.endswith("]"):
        return core, None, None
    open_pos = core.rfind("[")
    if open_pos <= 0:
        return core, None, None
    if core[open_pos - 1] != "c":
        return core, None, None
    return core[: open_pos - 1], core[open_pos - 1 :], core[open_pos:]


def _normalize_core_peptide(
    core: str,
    entries: Sequence[FlexUnimodEntry],
    ppm_tolerance: float,
    warning_sink: WarningSink,
    original_peptide: str,
) -> str:
    body = core
    nterm_original: Optional[str]
    nterm_bracket: Optional[str]
    body, nterm_original, nterm_bracket = _extract_nterm_mod(body)

    cterm_original: Optional[str]
    cterm_bracket: Optional[str]
    body, cterm_original, cterm_bracket = _extract_cterm_mod(body)

    out_parts: List[str] = []
    idx = 0
    while idx < len(body):
        ch = body[idx]
        out_parts.append(ch)

        if ch.isalpha() and idx + 1 < len(body) and body[idx + 1] == "[":
            close = _find_closing_bracket(body, idx + 1)
            if close > idx + 1:
                mod_tag = body[idx + 1 : close + 1]
                replaced = _normalize_mod_tag(
                    mod_tag,
                    ch.upper(),
                    entries,
                    ppm_tolerance,
                    warning_sink,
                    original_peptide,
                )
                out_parts.append(replaced)
                idx = close + 1
                continue
        idx += 1

    normalized_body = "".join(out_parts)

    prefix = ""
    suffix = ""
    if nterm_original is not None and nterm_bracket is not None:
        mapped = _normalize_mod_tag(
            nterm_bracket,
            "n",
            entries,
            ppm_tolerance,
            warning_sink,
            original_peptide,
        )
        if mapped.startswith("[U:"):
            prefix = f"{mapped}-"
        else:
            prefix = f"{nterm_original}-"

    if cterm_original is not None and cterm_bracket is not None:
        mapped = _normalize_mod_tag(
            cterm_bracket,
            "c",
            entries,
            ppm_tolerance,
            warning_sink,
            original_peptide,
        )
        if mapped.startswith("[U:"):
            suffix = f"-{mapped}"
        else:
            suffix = f"-{cterm_original}"

    return prefix + normalized_body + suffix


def normalize_pin_peptide_to_unimod(
    peptide_text: str,
    entries: Sequence[FlexUnimodEntry],
    ppm_tolerance: float,
    warning_sink: WarningSink,
) -> str:
    peptide = peptide_text.strip()
    if not peptide:
        return peptide
    core = _strip_pin_flanks(peptide)
    return _normalize_core_peptide(core, entries, ppm_tolerance, warning_sink, peptide)


__all__ = [
    "FlexUnimodEntry",
    "WarningSink",
    "load_unimod_entries_flexible",
    "normalize_pin_peptide_to_unimod",
]
