**pyPeptideIndex.py**
This document explains how `pyPeptideIndex.py` works and describes every output file and column.

**Overview**
`pyPeptideIndex.py` builds a Comet-like peptide index from a FASTA database and a `comet.params` file. It digests protein sequences using Comet enzyme rules, enumerates variable modifications, computes peptide masses, and writes multiple TSV files (no headers). Progress is printed to stderr to help debug long runs.

**How It Works (High Level)**
1. Parse `comet.params` using `pyLoadParameters.py`, following Comet parsing rules (comments at `#`, `variable_modNN` requires 8 fields, and `[COMET_ENZYME_INFO]` is parsed at EOF).
2. Determine input proteins from `--database` FASTA files, `--protein` inline sequences, or `database_name` in the params file.
3. Normalize sequences to uppercase and keep only letters and `*`. `*` is treated as a stop and splits the protein into segments for digestion.
4. Digest each segment using the search enzyme (and optional second enzyme) with `num_enzyme_termini`, `allowed_missed_cleavage`, and `peptide_length_range`. If `clip_nterm_methionine=1`, also generate peptides starting at position 2 when the protein starts with `M`.
5. Record unique peptide sequences across all proteins. The "primary" protein for a peptide is the earliest occurrence (smallest FASTA file offset).
6. Compute base peptide mass (MH+) using `mass_type_parent` (mono or average), `set_` masses, static `add_` mods, and terminal static mods.
7. Enumerate variable mods using `variable_modNN` rules, `max_variable_mods_in_peptide`, and `require_variable_mod`. Each variant is filtered by `digest_mass_range`.
8. Write TSV tables. Sorting is stable and reproducible (variants are sorted by mass, sequence, and mod sites).

**Command Line Options**
1. `--params` / `-P`: Path to `comet.params`. Used for enzyme rules, digestion filters, and mods.
2. `--database` / `-D`: One or more FASTA files (repeatable). Overrides `database_name`.
3. `--protein`: Inline protein sequences (repeatable). Synthetic headers and offsets are used.
4. `--prefix` / `-N`: Output file prefix (can include a directory).
5. `--max-record`: Maximum number of proteins to process. Default is all proteins. Intended for test/debug.
6. `--use-protein-name`: Use the protein name from the FASTA header as `protein_id` instead of a numeric counter. The protein name is the first whitespace-delimited token in the FASTA header.

**Variable Mod Enumeration Details**
1. A `variable_modNN` entry creates a mod with 1-based index `NN` (01..15).
2. Candidate positions are:
   - Residue positions `0..len(seq)-1` where the residue matches the `residues` field.
   - Peptide N-terminus position `len(seq)` if `residues` includes `n`.
   - Peptide C-terminus position `len(seq)+1` if `residues` includes `c`.
3. `term_distance` and `which_term` filter candidate positions:
   - `which_term=0` protein N, `1` protein C, `2` peptide N, `3` peptide C.
4. `binary=1` enforces at most one site per peptide for that mod.
5. `min_per_pep` and `max_per_pep` are enforced per mod.
6. `require_this_mod > 0` forces at least one site for that mod.
7. `require_this_mod < 0` (exclusive sets) is not enforced by this script.

**Output Files (TSV, No Headers)**
All output files are TSV with no header row. Columns are listed in order. `run_id` is always `1` in the current implementation and can be used as a join key.

**`<prefix>.index_run.tsv`**
One row describing the index build.
1. `run_id`: Integer, currently always `1`.
2. `comet_version`: Version string read from `# comet_version` in the params file.
3. `params_path`: The params file path.
4. `database_label`: Space-separated FASTA paths or `<inline_sequences>`.
5. `digest_mass_min`: Minimum MH+ from `digest_mass_range`.
6. `digest_mass_max`: Maximum MH+ from `digest_mass_range`.
7. `peptide_len_min`: Minimum length from `peptide_length_range`.
8. `peptide_len_max`: Maximum length from `peptide_length_range`.
9. `created_at_utc`: UTC timestamp in ISO 8601, suffixed with `Z`.
10. `params_json`: JSON snapshot of parsed params, variable mods, and enzymes.

**`<prefix>.static_mod.tsv`**
Static modifications applied to all peptides.
1. `run_id`
2. `residue`: Single-letter residue or `-` for terminal mods.
3. `delta_mass`: Mass delta applied.
4. `site`: One of `residue`, `N-term`, `C-term`, `protein N-term`, `protein C-term`.

**`<prefix>.variable_mod.tsv`**
Variable modification definitions.
1. `run_id`
2. `mod_index`: 1-based `variable_modNN` index.
3. `residues`: Residues string from params (may include `n` and `c`).
4. `delta_mass`: Variable mod mass.
5. `binary_mod`: `1` if at most one site per peptide, else `0`.
6. `min_per_pep`: Minimum count if this mod is present.
7. `max_per_pep`: Maximum count for this mod.
8. `term_distance`: Distance filter, `-1` disables.
9. `which_term`: `0` protein N, `1` protein C, `2` peptide N, `3` peptide C.
10. `require_this_mod`: `>0` means required; `<0` indicates an exclusive set in Comet.
11. `neutral_loss1`: First neutral loss mass.
12. `neutral_loss2`: Second neutral loss mass (0 if not provided).

**`<prefix>.protein.tsv`**
Protein table from the input FASTA.
1. `run_id`
2. `protein_id`: 1-based ID assigned in input order, or the FASTA protein name if `--use-protein-name` is set.
3. `fasta_offset`: Byte offset of the header line (`>`) in the FASTA file. Synthetic for `--protein`.
4. `header`: FASTA header line without the leading `>`, with tabs/newlines normalized to spaces.

**`<prefix>.peptide_sequence.tsv`**
Unique peptide sequences across all proteins.
1. `run_id`
2. `peptide_sequence_id`: 1-based ID, assigned after sorting sequences.
3. `sequence`: Peptide sequence string.
4. `length`: Peptide length.
5. `primary_protein_id`: Protein ID of the earliest occurrence (smallest offset). If `--use-protein-name` is set, this is the FASTA protein name.

**`<prefix>.peptide_sequence_protein.tsv`**
Many-to-many mapping between peptides and proteins.
1. `run_id`
2. `peptide_sequence_id`
3. `protein_id` (FASTA protein name if `--use-protein-name` is set)

**`<prefix>.peptide_variant.tsv`**
Each peptide sequence expanded into variable-modified variants.
1. `run_id`
2. `variant_id`: 1-based ID assigned during enumeration.
3. `peptide_sequence_id`
4. `mh_plus`: Calculated MH+ mass (includes static and variable mods).
5. `prev_aa`: Residue before the peptide in the primary protein, or `-` if N-terminus or after `*`.
6. `next_aa`: Residue after the peptide in the primary protein, or `-` if C-terminus or before `*`.
   - If a peptide occurs in multiple contexts with different flanks, only one flank pair is stored.
   - The stored flanks come from the "primary" occurrence with the smallest FASTA file offset (earliest protein).
   - This matches Comet’s peptide index tie-breaker, which keeps flanks from the lowest protein index when duplicates exist.
7. `var_mod_sites`: `pos:mod_index` pairs separated by `;`. `pos` is 0-based residue index. `pos=len(seq)` is peptide N-term, `pos=len(seq)+1` is peptide C-term.
8. `var_mod_count`: Total number of variable mods in this variant.
9. `mass_bin10`: `int(mh_plus * 10.0)` for binning.

**`<prefix>.peptide_variant_mod.tsv`**
One row per modified site in each variant.
1. `run_id`
2. `variant_id`
3. `position`: Same encoding as `var_mod_sites` (0-based residue, N-term `len`, C-term `len+1`).
4. `mod_index`: 1-based variable mod index.

**Notes and Practical Usage**
1. `require_variable_mod=1` removes unmodified variants entirely.
2. Peptides containing residues with no mass (e.g., B, J, X, Z) are kept in `peptide_sequence.tsv` but produce no variants.
3. Use `--max-record` for quick validation on large databases; output will be incomplete by design.
4. Progress logs are printed every 1000 proteins and every 1000 sequences during variant enumeration.
5. `X` means unknown residue. In the current implementation `X` has mass `0.0`, so any peptide containing `X` is excluded from `peptide_variant.tsv` and `peptide_variant_mod.tsv` (no calculable mass), but the sequence still appears in `peptide_sequence.tsv` and `peptide_sequence_protein.tsv`.
6. When `--use-protein-name` is set, the protein identifier columns become strings from the FASTA header (first token). Ensure those names are unique if you rely on `protein_id` as a key.
