#!/usr/bin/env python3
"""Match active Comet params modifications to UniMod IDs."""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Tuple

sys.path.append(os.path.dirname(__file__))

from utils.pyLoadParameters import CometParams, parse_comet_params
from utils.unimod_matcher import resolve_params_unimod_ids


_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_UNIMOD_PATH = os.path.join(_SCRIPT_DIR, "resource", "common_unimod.csv")
_ZERO_EPS = 1e-12


def _normalize_unimod_arg(unimod_path: str) -> str:
    value = unimod_path.strip()
    if not value:
        return ""
    if value.lower() == "none":
        return ""
    return os.path.expanduser(value)


def _iter_active_variable_mods(params: CometParams) -> List[Tuple[int, float, str]]:
    rows: List[Tuple[int, float, str]] = []
    for mod in sorted(params.variable_mods, key=lambda item: item.index):
        if abs(float(mod.mass)) <= _ZERO_EPS:
            continue
        residues = mod.residues.strip()
        if residues == "-":
            continue
        rows.append((mod.index, float(mod.mass), residues))
    return rows


def _iter_active_fixed_mods(params: CometParams) -> List[Tuple[str, float, str]]:
    rows: Dict[str, Tuple[str, float, str]] = {}
    term_to_residue = {
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
        if key in term_to_residue:
            rows[key] = (key, mass_delta, term_to_residue[key])
            continue
        residue = key.split("_", 2)[1].upper()[:1]
        if not residue:
            continue
        mapped_key = f"add_{residue}"
        rows[mapped_key] = (mapped_key, mass_delta, residue)
    return [rows[key] for key in sorted(rows.keys())]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Resolve active Comet parameter modifications to UniMod IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-P",
        "--params",
        default="comet.params",
        help="Comet params file path.",
    )
    parser.add_argument(
        "--unimod",
        default=_DEFAULT_UNIMOD_PATH,
        help=(
            "UniMod CSV mapping file. Use 'None' to disable matching.\n"
            "If missing/None, prints a message and exits successfully."
        ),
    )
    parser.add_argument(
        "--ppm",
        type=float,
        default=10.0,
        help="PPM tolerance for mapping modifications to UniMod entries.",
    )
    args = parser.parse_args()

    if args.ppm < 0.0:
        parser.error("--ppm must be >= 0")

    params = parse_comet_params(args.params)
    unimod_path = _normalize_unimod_arg(args.unimod)
    if not unimod_path:
        print("cannot match unimod: --unimod is None; continue normally.")
        raise SystemExit(0)
    if not os.path.isfile(unimod_path):
        print(f"cannot match unimod: file does not exist: {unimod_path}; continue normally.")
        raise SystemExit(0)

    resolution = resolve_params_unimod_ids(params, unimod_path, args.ppm)
    if resolution.missing:
        print(
            f"cannot match unimod: {len(resolution.missing)} active modification(s) were not matched:"
        )
        for missing in resolution.missing:
            print(f"  - {missing.to_message()}")
        raise SystemExit(1)

    print("variable modifications:")
    variable_count = 0
    for index, mass_delta, residues in _iter_active_variable_mods(params):
        unimod_id = resolution.variable_mod_to_unimod.get(index)
        print(f"  variable_mod{index:02d}\tmass={mass_delta:.6f}\tresidues={residues}\tunimod={unimod_id}")
        variable_count += 1

    print("fixed modifications:")
    fixed_count = 0
    for key, mass_delta, residues in _iter_active_fixed_mods(params):
        unimod_id = resolution.fixed_mod_to_unimod.get(key)
        print(f"  {key}\tmass={mass_delta:.6f}\tresidues={residues}\tunimod={unimod_id}")
        fixed_count += 1

    print(f"matched variable mods: {variable_count}")
    print(f"matched fixed mods: {fixed_count}")


if __name__ == "__main__":
    main()
