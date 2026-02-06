from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import os
from typing import Dict, List, Optional, Sequence, Tuple

from .pyLoadParameters import CometParams, VarMod


PROTON_MASS = 1.00727646688
HYDROGEN_MONO = 1.00782503223
OXYGEN_MONO = 15.99491461957

SequenceTask = Tuple[int, str, int, int, str, str, int, int]
VariantTaskRow = Tuple[int, float, str, str, str, int, int, List[int]]


@dataclass(frozen=True)
class VariantContext:
    aa_masses: Dict[str, float]
    h2o_mass: float
    residue_mods: Dict[str, float]
    term_mods: Dict[str, float]
    var_mods: List[VarMod]
    max_var_mods: int
    require_var_mod: int
    min_mass: float
    max_mass: float


_VARIANT_WORKER_CONTEXT: Optional[VariantContext] = None


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


def resolve_thread_count(thread_count: int) -> int:
    if thread_count < 0:
        raise ValueError("thread count must be >= 0")
    if thread_count == 0:
        return max(1, os.cpu_count() or 1)
    return thread_count


def enumerate_sequence_variants(
    task: SequenceTask,
    context: VariantContext,
) -> List[VariantTaskRow]:
    (
        seq_id,
        seq,
        start,
        end,
        prev_aa,
        next_aa,
        protein_len,
        nterm_offset,
    ) = task

    is_protein_nterm = start == nterm_offset
    is_protein_cterm = end == protein_len - 1

    base_mass = compute_peptide_mass(
        seq,
        context.aa_masses,
        context.h2o_mass,
        context.residue_mods,
        context.term_mods,
        is_protein_nterm,
        is_protein_cterm,
    )
    if base_mass is None:
        return []

    candidates = build_var_mod_candidates(
        seq,
        context.var_mods,
        start,
        protein_len,
        nterm_offset,
    )
    mod_combos = enumerate_var_mods(
        seq,
        context.var_mods,
        candidates,
        context.max_var_mods,
    )

    rows: List[VariantTaskRow] = []
    for sites, total_mods, delta_mass in mod_combos:
        if context.require_var_mod and total_mods == 0:
            continue
        mh_plus = base_mass + delta_mass
        if mh_plus < context.min_mass or mh_plus > context.max_mass:
            continue
        mass_bin10 = int(mh_plus * 10.0)
        sites_text = format_var_mod_sites(sites)
        rows.append(
            (
                seq_id,
                mh_plus,
                prev_aa,
                next_aa,
                sites_text,
                total_mods,
                mass_bin10,
                sites.copy(),
            )
        )
    return rows


def init_variant_worker(context: VariantContext) -> None:
    global _VARIANT_WORKER_CONTEXT
    _VARIANT_WORKER_CONTEXT = context


def enumerate_variant_chunk(tasks: Sequence[SequenceTask]) -> List[VariantTaskRow]:
    if _VARIANT_WORKER_CONTEXT is None:
        raise RuntimeError("Variant worker context not initialized.")
    rows: List[VariantTaskRow] = []
    for task in tasks:
        rows.extend(enumerate_sequence_variants(task, _VARIANT_WORKER_CONTEXT))
    return rows


__all__ = [
    "PROTON_MASS",
    "HYDROGEN_MONO",
    "OXYGEN_MONO",
    "SequenceTask",
    "VariantTaskRow",
    "VariantContext",
    "build_aa_masses",
    "apply_set_masses",
    "get_static_mods",
    "build_var_mod_candidates",
    "enumerate_var_mods",
    "format_var_mod_sites",
    "compute_peptide_mass",
    "resolve_thread_count",
    "enumerate_sequence_variants",
    "init_variant_worker",
    "enumerate_variant_chunk",
]

