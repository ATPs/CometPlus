#!/usr/bin/env python3
"""Comet params parser (Python).

This module mirrors the key parsing behavior in Comet.cpp::LoadParameters():
- strip comments at '#'
- only parse lines with '='
- variable_modNN requires 8 fields
- [COMET_ENZYME_INFO] block is parsed at EOF

The resulting CometParams object is intended to be reused by other modules.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple, Any
import json
import re


DEFAULTS = {
    "database_name": "",
    "search_enzyme_number": 1,
    "search_enzyme2_number": 0,
    "sample_enzyme_number": 1,
    "num_enzyme_termini": 2,
    "allowed_missed_cleavage": 2,
    "clip_nterm_methionine": 0,
    "digest_mass_range": (600.0, 5000.0),
    "peptide_length_range": (5, 50),
    "max_variable_mods_in_peptide": 5,
    "require_variable_mod": 0,
    "mass_type_parent": 1,
}


@dataclass
class VarMod:
    index: int
    mass: float
    residues: str
    binary_mod: int
    min_per_pep: int
    max_per_pep: int
    term_distance: int
    which_term: int
    require_this_mod: int
    neutral_loss1: float
    neutral_loss2: float


@dataclass
class EnzymeDefinition:
    number: int
    name: str
    offset: int
    break_aa: str
    no_break_aa: str
    raw_line: str

    def is_no_enzyme(self) -> bool:
        name = self.name.lower()
        if name in ("cut_everywhere", "no_enzyme", "no_enzyme_selected", "no_enzyme_selected"):
            return True
        if self.break_aa == "-" and self.no_break_aa == "-":
            return True
        return False


@dataclass
class CometParams:
    params_path: str
    version: str
    raw_params: Dict[str, str]
    parsed: Dict[str, Any]
    variable_mods: List[VarMod]
    enzymes: Dict[str, EnzymeDefinition]
    enzyme_table: List[EnzymeDefinition]

    def to_json(self) -> str:
        payload = {
            "params_path": self.params_path,
            "version": self.version,
            "raw_params": dict(self.raw_params),
            "parsed": dict(self.parsed),
            "variable_mods": [asdict(x) for x in self.variable_mods],
            "enzymes": {k: asdict(v) for k, v in self.enzymes.items()},
        }
        return json.dumps(payload, sort_keys=True)


_VAR_MOD_RE = re.compile(r"^variable_mod(\d\d)$")


def _trim_whitespace(value: str) -> str:
    return value.strip()


def _parse_int(value: str) -> int:
    try:
        return int(value.strip().split()[0])
    except ValueError:
        return 0


def _parse_float(value: str) -> float:
    try:
        return float(value.strip().split()[0])
    except ValueError:
        return 0.0


def _parse_int_range(value: str) -> Tuple[int, int]:
    parts = value.strip().split()
    if len(parts) >= 2:
        return int(parts[0]), int(parts[1])
    if len(parts) == 1:
        return int(parts[0]), int(parts[0])
    return 0, 0


def _parse_float_range(value: str) -> Tuple[float, float]:
    parts = value.strip().split()
    if len(parts) >= 2:
        return float(parts[0]), float(parts[1])
    if len(parts) == 1:
        return float(parts[0]), float(parts[0])
    return 0.0, 0.0


def _parse_mass_list(value: str) -> List[float]:
    values: List[float] = []
    for token in value.strip().split():
        try:
            val = float(token)
        except ValueError:
            continue
        if val >= 0.0:
            values.append(val)
    values.sort()
    return values


def _parse_variable_mod(index: int, value: str) -> VarMod:
    fields = value.strip().split()
    if len(fields) != 8:
        raise ValueError(
            f"variable_mod{index:02d} expects 8 fields; got {len(fields)}: {value!r}"
        )
    mass = float(fields[0])
    residues = fields[1]
    binary_mod = int(fields[2])
    min_max = fields[3]
    term_distance = int(fields[4])
    which_term = int(fields[5])
    require_this_mod = int(fields[6])
    neutral_loss = fields[7]

    if "," in min_max:
        min_per, max_per = min_max.split(",", 1)
        try:
            min_per_pep = int(min_per)
        except ValueError:
            min_per_pep = 0
        try:
            max_per_pep = int(max_per)
        except ValueError:
            max_per_pep = 0
    else:
        min_per_pep = 0
        try:
            max_per_pep = int(min_max)
        except ValueError:
            max_per_pep = 0

    if "," in neutral_loss:
        nl1, nl2 = neutral_loss.split(",", 1)
        try:
            neutral_loss1 = float(nl1)
        except ValueError:
            neutral_loss1 = 0.0
        try:
            neutral_loss2 = float(nl2)
        except ValueError:
            neutral_loss2 = 0.0
    else:
        try:
            neutral_loss1 = float(neutral_loss)
        except ValueError:
            neutral_loss1 = 0.0
        neutral_loss2 = 0.0

    return VarMod(
        index=index,
        mass=mass,
        residues=residues,
        binary_mod=binary_mod,
        min_per_pep=min_per_pep,
        max_per_pep=max_per_pep,
        term_distance=term_distance,
        which_term=which_term,
        require_this_mod=require_this_mod,
        neutral_loss1=neutral_loss1,
        neutral_loss2=neutral_loss2,
    )


def _parse_enzyme_line(line: str) -> Optional[EnzymeDefinition]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    match = re.match(r"^(\d+)\.\s+(\S+)\s+(-?\d+)\s+(\S+)\s+(\S+)", line)
    if not match:
        return None
    number = int(match.group(1))
    name = match.group(2)
    offset = int(match.group(3))
    break_aa = match.group(4)
    no_break_aa = match.group(5)
    return EnzymeDefinition(
        number=number,
        name=name,
        offset=offset,
        break_aa=break_aa,
        no_break_aa=no_break_aa,
        raw_line=line,
    )


def parse_comet_params(params_path: str) -> CometParams:
    version = ""
    raw_params: Dict[str, str] = {}
    parsed: Dict[str, Any] = {}
    variable_mods: List[VarMod] = []
    enzyme_table: List[EnzymeDefinition] = []

    with open(params_path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    # Version check pass
    for line in lines:
        if line.startswith("# comet_version "):
            version = line.strip()[len("# comet_version "):].strip()
            break
    if not version:
        raise ValueError("Missing '# comet_version' entry in params file")

    in_enzyme_block = False
    for line in lines:
        if line.startswith("[COMET_ENZYME_INFO]"):
            in_enzyme_block = True
            continue
        if in_enzyme_block:
            enzyme = _parse_enzyme_line(line)
            if enzyme:
                enzyme_table.append(enzyme)
            continue

        # Strip comments
        if "#" in line:
            line = line.split("#", 1)[0]
        if "=" not in line:
            continue

        left, right = line.split("=", 1)
        param_name = left.strip().split()[0]
        param_val = right.strip()
        raw_params[param_name] = param_val

        var_match = _VAR_MOD_RE.match(param_name)
        if var_match:
            index = int(var_match.group(1))
            variable_mods.append(_parse_variable_mod(index, param_val))
            continue

        if param_name in ("database_name", "peff_obo", "spectral_library_name"):
            parsed[param_name] = _trim_whitespace(param_val)
        elif param_name in (
            "allowed_missed_cleavage",
            "clip_nterm_methionine",
            "num_enzyme_termini",
            "search_enzyme_number",
            "search_enzyme2_number",
            "sample_enzyme_number",
            "max_variable_mods_in_peptide",
            "require_variable_mod",
            "mass_type_parent",
            "mass_type_fragment",
            "equal_I_and_L",
        ):
            parsed[param_name] = _parse_int(param_val)
        elif param_name in ("peptide_length_range",):
            parsed[param_name] = _parse_int_range(param_val)
        elif param_name in ("digest_mass_range", "ms1_mass_range", "clear_mz_range"):
            parsed[param_name] = _parse_float_range(param_val)
        elif param_name in ("mass_offsets", "precursor_NL_ions"):
            parsed[param_name] = _parse_mass_list(param_val)
        elif param_name.startswith("add_") or param_name.startswith("set_"):
            parsed[param_name] = _parse_float(param_val)
        else:
            # store parsed as raw string for anything else
            parsed[param_name] = param_val

    # Apply defaults
    for key, value in DEFAULTS.items():
        if key not in parsed:
            parsed[key] = value

    search_enzyme_number = int(parsed.get("search_enzyme_number", DEFAULTS["search_enzyme_number"]))
    search_enzyme2_number = int(parsed.get("search_enzyme2_number", DEFAULTS["search_enzyme2_number"]))
    sample_enzyme_number = int(parsed.get("sample_enzyme_number", DEFAULTS["sample_enzyme_number"]))
    allowed_missed = int(parsed.get("allowed_missed_cleavage", DEFAULTS["allowed_missed_cleavage"]))

    enzyme_map: Dict[int, EnzymeDefinition] = {e.number: e for e in enzyme_table}

    def _get_enzyme(number: int) -> EnzymeDefinition:
        enzyme = enzyme_map.get(number)
        if not enzyme:
            raise ValueError(f"Enzyme number {number} missing from [COMET_ENZYME_INFO] block")
        return enzyme

    if search_enzyme2_number != 0 and search_enzyme2_number not in enzyme_map:
        raise ValueError(
            f"Enzyme number {search_enzyme2_number} missing from [COMET_ENZYME_INFO] block"
        )

    enzymes = {
        "search": _get_enzyme(search_enzyme_number),
        "search2": _get_enzyme(search_enzyme2_number)
        if search_enzyme2_number in enzyme_map
        else EnzymeDefinition(
            number=0,
            name="Cut_everywhere",
            offset=0,
            break_aa="-",
            no_break_aa="-",
            raw_line="0.  Cut_everywhere 0 - -",
        ),
        "sample": _get_enzyme(sample_enzyme_number),
    }
    parsed["allowed_missed_cleavage"] = allowed_missed

    return CometParams(
        params_path=params_path,
        version=version,
        raw_params=raw_params,
        parsed=parsed,
        variable_mods=variable_mods,
        enzymes=enzymes,
        enzyme_table=enzyme_table,
    )


__all__ = [
    "VarMod",
    "EnzymeDefinition",
    "CometParams",
    "parse_comet_params",
]
