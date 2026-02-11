from __future__ import annotations

from dataclasses import dataclass
import csv
import math
import re
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .pyLoadParameters import CometParams


_ZERO_EPS = 1e-12
_UNIMOD_ID_RE = re.compile(r"(\d+)$")


@dataclass(frozen=True)
class UnimodEntry:
    mono_mass_delta: float
    residues: str
    unimod_id: str
    unimod_name: str
    row_index: int


@dataclass(frozen=True)
class MissingUnimodMatch:
    mod_type: str
    key: str
    mass_delta: float
    residues: str
    reason: str
    nearest_unimod_id: Optional[str] = None
    nearest_unimod_name: Optional[str] = None
    nearest_ppm_error: Optional[float] = None

    def to_message(self) -> str:
        prefix = (
            f"{self.mod_type} {self.key}: mass={self.mass_delta:.6f}, "
            f"residues={self.residues!r}, reason={self.reason}"
        )
        if self.nearest_unimod_id is None:
            return prefix
        return (
            f"{prefix}; nearest=UNIMOD:{self.nearest_unimod_id} "
            f"({self.nearest_unimod_name}) at {self.nearest_ppm_error:.3f} ppm"
        )


@dataclass(frozen=True)
class UnimodResolution:
    variable_mod_to_unimod: Dict[int, str]
    fixed_mod_to_unimod: Dict[str, str]
    missing: List[MissingUnimodMatch]


@dataclass(frozen=True)
class _ObservedMod:
    key: str
    mass_delta: float
    residues: str


def normalize_unimod_id(raw_value: str) -> str:
    value = raw_value.strip()
    if not value:
        return ""
    if ":" in value:
        value = value.split(":", 1)[1].strip()
    match = _UNIMOD_ID_RE.search(value)
    if not match:
        return ""
    return str(int(match.group(1)))


def _residue_tokens(value: str) -> Set[str]:
    tokens: Set[str] = set()
    for ch in value.strip():
        if ch == "n" or ch == "c":
            tokens.add(ch)
        elif ch.isalpha():
            tokens.add(ch.upper())
    return tokens


def _ppm_error(observed_mass: float, reference_mass: float) -> float:
    if abs(reference_mass) <= _ZERO_EPS:
        return math.inf
    return abs(observed_mass - reference_mass) / abs(reference_mass) * 1e6


def load_unimod_entries(unimod_path: str) -> List[UnimodEntry]:
    entries: List[UnimodEntry] = []
    with open(unimod_path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
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
            unimod_name = ",".join(row[3:]).strip() if len(row) > 3 else ""
            entries.append(
                UnimodEntry(
                    mono_mass_delta=mass_delta,
                    residues=residues,
                    unimod_id=unimod_id,
                    unimod_name=unimod_name,
                    row_index=len(entries),
                )
            )
    return entries


def _iter_active_variable_mods(params: CometParams) -> List[_ObservedMod]:
    rows: List[_ObservedMod] = []
    for mod in sorted(params.variable_mods, key=lambda item: item.index):
        if abs(float(mod.mass)) <= _ZERO_EPS:
            continue
        residues = mod.residues.strip()
        if residues == "-":
            continue
        rows.append(
            _ObservedMod(
                key=f"variable_mod{mod.index:02d}",
                mass_delta=float(mod.mass),
                residues=residues,
            )
        )
    return rows


def _iter_active_fixed_mods(params: CometParams) -> List[_ObservedMod]:
    rows: Dict[str, _ObservedMod] = {}
    terminal_keys = {
        "add_Nterm_peptide": "n",
        "add_Nterm_protein": "n",
        "add_Cterm_peptide": "c",
        "add_Cterm_protein": "c",
    }

    for key, value in params.parsed.items():
        if not key.startswith("add_"):
            continue
        if not isinstance(value, (int, float)):
            continue
        mass_delta = float(value)
        if abs(mass_delta) <= _ZERO_EPS:
            continue
        if key in terminal_keys:
            rows[key] = _ObservedMod(
                key=key,
                mass_delta=mass_delta,
                residues=terminal_keys[key],
            )
            continue
        residue = key.split("_", 2)[1].upper()[:1]
        if not residue:
            continue
        mapped_key = f"add_{residue}"
        rows[mapped_key] = _ObservedMod(
            key=mapped_key,
            mass_delta=mass_delta,
            residues=residue,
        )
    return [rows[key] for key in sorted(rows.keys())]


def _select_best_unimod(
    entries: Sequence[UnimodEntry],
    mass_delta: float,
    residues: str,
    ppm_tolerance: float,
) -> Tuple[Optional[UnimodEntry], Optional[float], Optional[UnimodEntry], Optional[float], str]:
    observed_tokens = _residue_tokens(residues)
    if not observed_tokens:
        return None, None, None, None, "no residue token in mod definition"

    candidates: List[Tuple[float, int, UnimodEntry]] = []
    for entry in entries:
        entry_tokens = _residue_tokens(entry.residues)
        if not observed_tokens.issubset(entry_tokens):
            continue
        ppm_error = _ppm_error(mass_delta, entry.mono_mass_delta)
        candidates.append((ppm_error, entry.row_index, entry))

    if not candidates:
        return None, None, None, None, "no residue-compatible UniMod entry"

    candidates.sort(key=lambda item: (item[0], item[1]))
    nearest_ppm, _nearest_order, nearest_entry = candidates[0]
    if nearest_ppm <= ppm_tolerance:
        return nearest_entry, nearest_ppm, nearest_entry, nearest_ppm, ""
    return (
        None,
        None,
        nearest_entry,
        nearest_ppm,
        "nearest residue-compatible UniMod exceeds ppm tolerance",
    )


def resolve_params_unimod_ids(
    params: CometParams,
    unimod_path: str,
    ppm_tolerance: float = 10.0,
) -> UnimodResolution:
    entries = load_unimod_entries(unimod_path)

    variable_mod_to_unimod: Dict[int, str] = {}
    fixed_mod_to_unimod: Dict[str, str] = {}
    missing: List[MissingUnimodMatch] = []

    for observed in _iter_active_variable_mods(params):
        matched, _matched_ppm, nearest, nearest_ppm, reason = _select_best_unimod(
            entries,
            observed.mass_delta,
            observed.residues,
            ppm_tolerance,
        )
        if matched is not None:
            mod_index = int(observed.key[-2:])
            variable_mod_to_unimod[mod_index] = matched.unimod_id
            continue
        missing.append(
            MissingUnimodMatch(
                mod_type="variable_mod",
                key=observed.key,
                mass_delta=observed.mass_delta,
                residues=observed.residues,
                reason=reason,
                nearest_unimod_id=nearest.unimod_id if nearest else None,
                nearest_unimod_name=nearest.unimod_name if nearest else None,
                nearest_ppm_error=nearest_ppm,
            )
        )

    for observed in _iter_active_fixed_mods(params):
        matched, _matched_ppm, nearest, nearest_ppm, reason = _select_best_unimod(
            entries,
            observed.mass_delta,
            observed.residues,
            ppm_tolerance,
        )
        if matched is not None:
            fixed_mod_to_unimod[observed.key] = matched.unimod_id
            continue
        missing.append(
            MissingUnimodMatch(
                mod_type="fixed_mod",
                key=observed.key,
                mass_delta=observed.mass_delta,
                residues=observed.residues,
                reason=reason,
                nearest_unimod_id=nearest.unimod_id if nearest else None,
                nearest_unimod_name=nearest.unimod_name if nearest else None,
                nearest_ppm_error=nearest_ppm,
            )
        )

    return UnimodResolution(
        variable_mod_to_unimod=variable_mod_to_unimod,
        fixed_mod_to_unimod=fixed_mod_to_unimod,
        missing=missing,
    )


__all__ = [
    "UnimodEntry",
    "MissingUnimodMatch",
    "UnimodResolution",
    "normalize_unimod_id",
    "load_unimod_entries",
    "resolve_params_unimod_ids",
]
