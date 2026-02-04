# Peptide index -> PostgreSQL table plan

## Scope

This note captures how Comet builds the peptide index (`create_peptide_index`) and maps that workflow into a set of relational tables you can populate from Python. It focuses on:

- how `comet.params` is parsed (exact file + flow)

- where peptide enumeration happens

- how the `.idx` file is structured so we can design equivalent tables

- practical TSV/CSV export guidance for Postgres `COPY`

## Quick glossary (database terms)

- `PK`: Primary key. A column (or set of columns) that uniquely identifies each row in a table.

- `FK`: Foreign key. A column that references a primary key in another table, used to link tables.

- `run_id`: A synthetic identifier for one index-generation run; use it to group all rows created from the same FASTA + params.

- `sequence_id`, `variant_id`, `protein_id`: Synthetic identifiers for rows in their respective tables.

## Source map (files to read)

- `Comet.cpp` - `LoadParameters()` parses `comet.params` and feeds `CometSearchManager::SetParam()`.

- `CometSearch/CometSearchManager.cpp` - `InitializeStaticParams()` transfers parsed params into `g_staticParams` and sets defaults; `DoSearch()` dispatches `create_peptide_index` to `CometPeptideIndex::WritePeptideIndex()`.

- `CometSearch/CometSearch.cpp` - FASTA traversal + digestion + variable-mod enumeration; this is where `g_pvDBIndex` entries are created.

- `CometSearch/CometPeptideIndex.cpp` - dedup, protein list construction, mass-sorted peptide list, and binary `.idx` writing.

- `CometSearch/CometDataInternal.h` - data structures (`DBIndex`, `IndexProteinStruct`) and key constants (`MAX_PEPTIDE_LEN`, `MAX_PEPTIDE_LEN_P2`, `VMODS`, `WIDTH_REFERENCE`).

- `CometSearch/Common.h` - binning macros; relevant if you want identical mass-bin rounding.

Note: this repo does not include a `comet.params` file. You will need the params file you actually used to build the index.

## How `comet.params` is parsed (exact behavior)

Parsing happens in `Comet.cpp::LoadParameters()`.

### High-level flow

1. **Version check (first pass).** Scans file for a line starting with `# comet_version` and validates via `ICometSearchManager::IsValidCometVersion()`. If no valid version is found, Comet exits with an error.

2. **Main parameters block.** Rewinds and reads line-by-line until `[COMET_ENZYME_INFO]` is reached.

3. **Enzyme block.** Continues reading from `[COMET_ENZYME_INFO]` to EOF, capturing enzyme definitions and selecting the entries that match `search_enzyme_number`, `search_enzyme2_number`, and `sample_enzyme_number`.

4. **Post-parse validation.** If `output_percolatorfile` was never seen, it flags the params file as "outdated" and errors. It also errors if any of the enzyme numbers point to a missing entry.

### Line parsing rules

- Lines are stripped of comments by truncating at the first `#`.

- Only lines containing `=` are parsed as parameters.

- The parameter name is the first token on the left-hand side (parsed with `sscanf(..., "%s")`).

- The parameter value is the raw text after `=`; for most params it is parsed via `sscanf` (whitespace-delimited), so **values with embedded spaces are not supported** unless explicitly handled.

- Unknown parameter names produce a warning and are ignored.

### Special parameter handling

- **String paths with trimming:** `database_name`, `peff_obo`, `spectral_library_name` call a `trim_whitespace()` helper first.

- **Ranges:** `peptide_length_range`, `digest_mass_range`, `scan_range`, etc. are parsed as two numbers.

- **Vector lists:** `mass_offsets` and `precursor_NL_ions` are parsed as whitespace-delimited numbers, non-negative only, then sorted.

- **Variable mods:** `variable_mod01` .. `variable_mod15` must each have **exactly 8 fields**. The parser expects:

  1) `dVarModMass` (double)

  2) `szVarModChar` (AA string)

  3) `iBinaryMod` (int)

  4) `min|max` count spec (string, either `N` or `min,max`)

  5) `iVarModTermDistance` (int)

  6) `iWhichTerm` (int)

  7) `iRequireThisMod` (int)

  8) `neutral_loss` (string, either `x` or `x,y`)

  Parsing details:

  - Field 4 is parsed into `iMinNumVarModAAPerMod` and `iMaxNumVarModAAPerMod`. If no comma is present, `min=0` and `max=value`.

  - Field 8 is parsed into `dNeutralLoss` and optional `dNeutralLoss2`.

  - The raw text of the parameter line (minus trailing newline) is stored in the params map.

### Where parsed values go

- `LoadParameters()` calls `ICometSearchManager::SetParam(...)`, which stores typed values in the `CometSearchManager` parameter map.

- `CometSearchManager::InitializeStaticParams()` later maps those stored values into `g_staticParams` (the global runtime configuration). This is where defaults and derived flags (like `bVarModSearch`, `bNoEnzymeSelected`, mass ranges, etc.) are set.

## Peptide index generation (create_peptide_index)

The indexing path is triggered by:

- command line `-j` **or** `create_peptide_index=1` in params

- `CometSearchManager::DoSearch()` calls `CometPeptideIndex::WritePeptideIndex()`

### FASTA walk + protein names

- `CometSearch::RunSearch()` reads the FASTA file.

- For each `>` header:

  - `dbe.lProteinFilePosition = ftell(fp)` is captured.

  - The protein name (header, truncated to `WIDTH_REFERENCE-1`) is stored into `g_pvProteinNames` when index creation is enabled.

  - The protein sequence is read and uppercased (A-Z); `*` is preserved as a stop-codon marker.

### Peptide enumeration

Peptide enumeration happens inside `CometSearch::SearchProtein()` (called via thread pool in `RunSearch()`). Key steps:

- Apply digestion rules: `CheckEnzymeTermini()` (and helper start/end checks) uses enzyme settings from `g_staticParams.enzymeInformation`.

- Apply length constraint: `peptide_length_range`.

- Apply mass constraint (`digest_mass_range`) **only at specific points**:

  - **Unmodified peptides:** included if mass is within range **and** `require_variable_mod` is not set.

  - **Variable-mod peptides:** included if mass is within range **and** at least one variable mod is present.

  - When creating an index, peptides are first collected by length and digestion; var mods can shift mass later.

- Variable-mod enumeration uses `VariableModSearch()` + `MergeVarMods()`; when a specific mod placement is found, `DBIndex` is created with `pcVarModSites` encoding the sites.

- For `create_peptide_index`, the **protein-variable-mod filter file** is not read (`ReadProteinVarModFilterFile()` is called only in the fragment index / classic search path), so protein-specific mod restrictions are ignored in this mode.

### DBIndex entries (in memory)

Each peptide variant is stored as a `DBIndex`:

- `szPeptide`: peptide sequence

- `cPrevAA` / `cNextAA`: flanking residues (`'-'` at termini)

- `pcVarModSites`: length `len+2`, each entry is `0` (no mod) or **1-based** mod index

  - position `0..len-1` = residues

  - position `len` = N-term mod

  - position `len+1` = C-term mod

- `dPepMass`: **MH+** mass including static + variable mods

- `lIndexProteinFilePosition`: initially a FASTA header file offset

### Dedup + protein list (WritePeptideIndex)

`CometPeptideIndex::WritePeptideIndex()` performs:

1. **Sort by peptide + mod state + protein offset.**

2. **Build `pvProteinsListLocal`:**

   - For each **unique peptide sequence** (mod state ignored), collect all protein offsets.

   - The list is sorted and de-duplicated.

   - Each `DBIndex.lIndexProteinFilePosition` is replaced by the index into this list.

3. **Deduplicate** peptide variants by sequence + mod state (`operator==`).

4. **Sort by mass** and write peptides with a 0.1 Da mass index.

## Suggested PostgreSQL table design

Below is a relational schema that mirrors the peptide index data while staying query-friendly. It assumes you want one import per "index run".

### 1) Index run metadata

`index_run`

- `run_id` (PK)

- `comet_version` (text)

- `params_path` (text)

- `database_path` (text)

- `digest_mass_low`, `digest_mass_high` (double)

- `peptide_len_min`, `peptide_len_max` (int)

- `created_at` (timestamp)

- `params_json` (jsonb, optional full parsed params)

### 2) Mod definitions

`static_mod`

- `run_id` (FK)

- `residue` (char)

- `delta_mass` (double)

- `site` (text: residue/N-term/C-term/protein-term)

`variable_mod`

- `run_id` (FK)

- `mod_index` (int, 1-based, matches `pcVarModSites` encoding)

- `residue_set` (text)

- `delta_mass` (double)

- `binary_mod` (int)

- `min_per_pep`, `max_per_pep` (int)

- `term_distance` (int)

- `which_term` (int)

- `require_this_mod` (int)

- `neutral_loss1`, `neutral_loss2` (double)

### 3) Proteins

`protein`

- `run_id` (FK)

- `protein_id` (PK)

- `file_offset` (bigint) - from FASTA `ftell()` at header

- `header` (text)

### 4) Peptide sequences (unique by sequence)

`peptide_sequence`

- `run_id` (FK)

- `sequence_id` (PK)

- `sequence` (text)

- `length` (int)

- `primary_protein_id` (FK, optional) - the lowest `file_offset` protein (matches Comet's tie-breaker for flanks)

`peptide_sequence_protein`

- `run_id` (FK)

- `sequence_id` (FK)

- `protein_id` (FK)

### 5) Peptide variants (sequence + mod state)

`peptide_variant`

- `run_id` (FK)

- `variant_id` (PK)

- `sequence_id` (FK)

- `mh_plus_mass` (double)

- `prev_aa` (char)

- `next_aa` (char)

- `var_mod_sites` (bytea or text) - positions 0..len+1, 1-based mod ids

- `var_mod_count` (int)

- `mass_bin10` (int) - `floor(mh_plus_mass * 10.0)`

Optional normalized mod table for faster mod queries:

`peptide_variant_mod`

- `run_id` (FK)

- `variant_id` (FK)

- `position` (int) - 0..len+1

- `mod_index` (int)

### Indexing suggestions

- `peptide_variant(run_id, mass_bin10)` btree index (fast range queries in 0.1 Da bins).

- `peptide_variant(run_id, mh_plus_mass)` btree index (for exact mass windows).

- `peptide_sequence(sequence)` unique index.

- `peptide_sequence_protein(sequence_id)` and `(protein_id)` for joins.

### Example rows (TSV-style)

Columns are shown in the same order listed above. Values are illustrative only.

`index_run`

```
1	2023.01	/data/params/comet.params	/data/db.fasta	500.0	5000.0	7	35	2026-01-28T12:00:00Z	{}
```

`static_mod`

```
1	C	57.02146	residue
```

`variable_mod`

```
1	1	M	15.9949	0	0	3	-1	0	0	0.0	0.0
```

`protein`

```
1	10	123456	sp|P12345|PROT_HUMAN Example protein
```

`peptide_sequence`

```
1	100	ACDEFGHIK	9	10
```

`peptide_sequence_protein`

```
1	100	10
```

`peptide_variant`

```
1	1000	100	1000.1234	K	R	3:1;8:2	2	10001
```

`peptide_variant_mod`

```
1	1000	3	1
```

### Notes on fidelity

- Comet stores **one** `prev_aa`/`next_aa` per peptide variant, derived from the lowest protein file offset (after sorting); storing `primary_protein_id` allows you to replicate this behavior.

- `pcVarModSites` uses 1-based mod indices; keep this identical to Comet so that `variable_mod.mod_index` matches.

- For mass binning and search windows, Comet uses truncation via `int(dMass * 10.0)` for indexing; if you need identical bin behavior, keep the same formula.

## TSV/CSV output recommendations

For Postgres `COPY`, TSV is the simplest.

### General export rules

- Use **TSV** with `\t` delimiter and `\n` row terminators.

- Use `\N` for NULL (Postgres default with `COPY`).

- Avoid quoting unless required; keep fields ASCII/UTF-8 and strip newlines from protein headers.

- If you must keep raw headers, output a sanitized column (e.g., replace `\t` with space, trim). Keep the raw header in a separate JSON or escaped column if needed.

### Suggested encoding for `var_mod_sites`

Pick one consistent format:

- **Text list:** `pos:modId;pos:modId` (empty string for no mods)

- **Array literal:** `{pos,modId,...}` (pairs in a flattened array)

- **Bytea:** store raw bytes for `len+2` positions (fast + compact, but less human readable)

If you use the normalized `peptide_variant_mod` table, you can avoid this column entirely.

### Deterministic ordering

To make repeated runs reproducible:

- Sort proteins by `file_offset`.

- Sort sequences lexicographically.

- Sort variants by `(mh_plus_mass, sequence, var_mod_sites)`.

### Example COPY setup

- `COPY protein FROM STDIN WITH (FORMAT csv, DELIMITER E'\t', NULL '\\N', HEADER false);`

- `COPY peptide_variant_mod FROM STDIN WITH (FORMAT csv, DELIMITER E'\t', NULL '\\N');`

## Implementation checklist (Python)

- Parse `comet.params` with the same tokenization rules as `LoadParameters()`.

- Apply `InitializeStaticParams()`-equivalent defaults (especially `digest_mass_range`, `peptide_length_range`, enzyme definitions, static + variable mods).

- Digest FASTA (respect `clip_nterm_methionine`, enzyme cleavage rules, missed cleavages).

- Enumerate variable mods using Comet's rules (`VariableModSearch()` + `MergeVarMods()` behavior). Start with a faithful implementation; optimize later.

- Build:

  - a unique peptide sequence table

  - a peptide->protein mapping table

  - peptide variants with mod sites + MH+ mass

- Dedup by `(sequence, var_mod_sites, mh_plus_mass)` before export.
