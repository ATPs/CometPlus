"""Regression tests for peptide variants with multiple protein occurrences."""

from __future__ import annotations

import unittest

from utils.peptide_index_tables import build_tables
from utils.pyLoadParameters import CometParams, EnzymeDefinition, VarMod


def make_var_mod(
    index: int,
    mass: float,
    residues: str,
    term_distance: int = -1,
    which_term: int = 0,
) -> VarMod:
    return VarMod(
        index=index,
        mass=mass,
        residues=residues,
        binary_mod=0,
        min_per_pep=0,
        max_per_pep=1,
        term_distance=term_distance,
        which_term=which_term,
        require_this_mod=0,
        neutral_loss1=0.0,
        neutral_loss2=0.0,
    )


def make_params(
    variable_mods=None,
    clip_nterm_methionine: int = 0,
    add_nterm_protein: float = 0.0,
    add_cterm_protein: float = 0.0,
) -> CometParams:
    trypsin = EnzymeDefinition(
        number=1,
        name="Trypsin",
        offset=1,
        break_aa="KR",
        no_break_aa="P",
        raw_line="1. Trypsin 1 KR P",
    )
    no_enzyme = EnzymeDefinition(
        number=0,
        name="Cut_everywhere",
        offset=0,
        break_aa="-",
        no_break_aa="-",
        raw_line="0. Cut_everywhere 0 - -",
    )
    return CometParams(
        params_path="test.params",
        version="test",
        raw_params={},
        parsed={
            "mass_type_parent": 1,
            "search_enzyme2_number": 0,
            "num_enzyme_termini": 2,
            "allowed_missed_cleavage": 2,
            "peptide_length_range": (1, 60),
            "digest_mass_range": (0.0, 10000.0),
            "clip_nterm_methionine": clip_nterm_methionine,
            "max_variable_mods_in_peptide": 3,
            "require_variable_mod": 0,
            "add_Nterm_protein": add_nterm_protein,
            "add_Cterm_protein": add_cterm_protein,
        },
        variable_mods=variable_mods or [],
        enzymes={"search": trypsin, "search2": no_enzyme, "sample": trypsin},
        enzyme_table=[no_enzyme, trypsin],
    )


def build_for(proteins, params, threads: int = 1):
    return build_tables(
        params,
        protein_sequences=proteins,
        threads=threads,
    )


def variants_for(tables, peptide: str):
    return [row for row in tables["peptide_variant"] if row[2] == peptide]


def variant_keys(tables, peptide: str):
    return {(row[6], row[9]) for row in variants_for(tables, peptide)}


class PeptideOccurrenceContextTests(unittest.TestCase):
    peptide = "EVIDENCEK"
    nterm_only = "EVIDENCEKAAAK"
    cterm_only = "AAAKEVIDENCEK"
    internal = "AAAKEVIDENCEKAAAK"

    def test_later_protein_nterm_occurrence_adds_terminal_variant_once(self):
        params = make_params([make_var_mod(1, 42.010565, "n", 0, 0)])
        tables = build_for([self.internal, self.nterm_only], params)

        self.assertEqual(variant_keys(tables, self.peptide), {("", ""), ("9:1", "")})
        self.assertEqual(len(variants_for(tables, self.peptide)), 2)

    def test_reversing_protein_order_preserves_variant_key_set(self):
        params = make_params([make_var_mod(1, 42.010565, "n", 0, 0)])
        forward = build_for([self.internal, self.nterm_only], params)
        reversed_tables = build_for([self.nterm_only, self.internal], params)

        self.assertEqual(
            variant_keys(forward, self.peptide),
            variant_keys(reversed_tables, self.peptide),
        )

    def test_internal_only_has_no_protein_nterm_variant(self):
        params = make_params([make_var_mod(1, 42.010565, "n", 0, 0)])
        tables = build_for([self.internal], params)

        self.assertEqual(variant_keys(tables, self.peptide), {("", "")})

    def test_clipped_initiator_methionine_is_protein_nterm(self):
        params = make_params(
            [make_var_mod(1, 42.010565, "n", 0, 0)],
            clip_nterm_methionine=1,
        )
        tables = build_for(["M" + self.nterm_only], params)

        self.assertIn(("9:1", ""), variant_keys(tables, self.peptide))

    def test_later_protein_cterm_occurrence_adds_terminal_variant_once(self):
        params = make_params([make_var_mod(1, 17.002740, "c", 0, 1)])
        tables = build_for([self.internal, self.cterm_only], params)

        self.assertEqual(variant_keys(tables, self.peptide), {("", ""), ("10:1", "")})
        self.assertEqual(len(variants_for(tables, self.peptide)), 2)

    def test_terminal_candidates_are_not_combined_across_contexts(self):
        params = make_params(
            [
                make_var_mod(1, 42.010565, "n", 0, 0),
                make_var_mod(2, 17.002740, "c", 0, 1),
            ]
        )
        split = build_for([self.nterm_only, self.cterm_only], params)
        split_keys = variant_keys(split, self.peptide)
        self.assertEqual(split_keys, {("", ""), ("9:1", ""), ("10:2", "")})
        self.assertNotIn(("9:1;10:2", ""), split_keys)

        whole = build_for([self.nterm_only, self.cterm_only, self.peptide], params)
        self.assertIn(("9:1;10:2", ""), variant_keys(whole, self.peptide))

    def test_positive_protein_term_distance_uses_union_of_real_contexts(self):
        params = make_params([make_var_mod(1, 12.0, "V", 1, 0)])
        tables = build_for([self.internal, self.nterm_only], params)

        self.assertEqual(variant_keys(tables, self.peptide), {("", ""), ("1:1", "")})
        internal_only = build_for([self.internal], params)
        self.assertEqual(variant_keys(internal_only, self.peptide), {("", "")})

    def test_protein_terminal_fixed_mods_affect_sites_and_mass(self):
        params = make_params(
            add_nterm_protein=10.0,
            add_cterm_protein=20.0,
        )
        tables = build_for(
            [self.internal, self.nterm_only, self.cterm_only, self.peptide],
            params,
        )
        rows = {row[9]: row for row in variants_for(tables, self.peptide)}
        self.assertEqual(
            set(rows),
            {
                "",
                "9:add_Nterm_protein",
                "10:add_Cterm_protein",
                "9:add_Nterm_protein;10:add_Cterm_protein",
            },
        )
        base = rows[""][3]
        self.assertAlmostEqual(rows["9:add_Nterm_protein"][3] - base, 10.0)
        self.assertAlmostEqual(rows["10:add_Cterm_protein"][3] - base, 20.0)
        self.assertAlmostEqual(
            rows["9:add_Nterm_protein;10:add_Cterm_protein"][3] - base,
            30.0,
        )

    def test_equivalent_contexts_deduplicate_variants_but_keep_locations(self):
        params = make_params([make_var_mod(1, 42.010565, "n", 0, 0)])
        tables = build_for([self.nterm_only, self.nterm_only], params)
        peptide_id = next(row[0] for row in tables["peptide_sequence"] if row[1] == self.peptide)
        locations = [row for row in tables["peptide_protein_location"] if row[1] == peptide_id]

        self.assertEqual(len(locations), 2)
        self.assertEqual(len(variants_for(tables, self.peptide)), 2)
        peptide_variant_ids = {row[0] for row in variants_for(tables, self.peptide)}
        peptide_mod_rows = [
            row for row in tables["peptide_variant_mod"] if row[0] in peptide_variant_ids
        ]
        self.assertEqual(len(peptide_mod_rows), 1)

    def test_unimod_disabled_keeps_distinct_internal_variant_keys(self):
        params = make_params(
            [
                make_var_mod(1, 42.010565, "n", 0, 0),
                make_var_mod(2, 17.002740, "c", 0, 1),
            ]
        )
        tables = build_for([self.peptide, self.peptide], params)
        rows = variants_for(tables, self.peptide)

        self.assertEqual({row[6] for row in rows}, {"", "9:1", "10:2", "9:1;10:2"})
        self.assertTrue(all(row[7] is None for row in rows))

    def test_threaded_enumeration_is_deterministic(self):
        params = make_params(
            [
                make_var_mod(1, 42.010565, "n", 0, 0),
                make_var_mod(2, 17.002740, "c", 0, 1),
            ]
        )
        proteins = [self.internal, self.nterm_only, self.cterm_only, self.peptide]
        one_thread = build_for(proteins, params, threads=1)
        two_threads = build_for(proteins, params, threads=2)

        self.assertEqual(one_thread["peptide_variant"], two_threads["peptide_variant"])
        self.assertEqual(one_thread["peptide_variant_mod"], two_threads["peptide_variant_mod"])

    def test_context_independent_mods_match_single_occurrence_variant_rows(self):
        params = make_params([make_var_mod(1, 15.994915, "M", -1, 0)])
        peptide = "MEVIDENCEK"
        one_occurrence = build_for([peptide + "AAAK"], params)
        many_occurrences = build_for([peptide + "AAAK", "AAAK" + peptide + "AAAK"], params)

        self.assertEqual(
            [row[2:] for row in variants_for(one_occurrence, peptide)],
            [row[2:] for row in variants_for(many_occurrences, peptide)],
        )


if __name__ == "__main__":
    unittest.main()
