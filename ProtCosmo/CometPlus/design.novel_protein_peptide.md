# CometPlus Novel Protein/Peptide Design (Detailed)

## 1. Scope and Added Options
CometPlus adds eight options in this design scope:

1. `--novel_protein <file>`
2. `--novel_peptide <file>`
3. `--scan <file>`
4. `--scan_numbers <list>`
5. `--output-folder <dir>`
6. `--output_internal_novel_peptide <file_internal_novel_peptide>`
7. `--internal_novel_peptide <file_internal_novel_peptide>`
8. `--stop-after-saving-novel-peptide`

Goal:

1. Build novel candidates from novel inputs.
2. Subtract known peptides under `equal_I_and_L` policy.
3. Search known + retained novel space.
4. Reduce scanned spectra by explicit scan constraints and novel-aware precursor-mass plausibility.

## 2. Option Matrix (Exact Behavior)
| Option | Input type | Parser behavior | Can combine with `-i/-j` |
|---|---|---|---|
| `--novel_protein` | FASTA | digested using current Comet params | No |
| `--novel_peptide` | FASTA or tokenized text | auto-detect FASTA by `>` line; otherwise token parser | No |
| `--scan` | text file of scan integers | delimiters: comma or any whitespace | No |
| `--scan_numbers` | inline scan list | same integer parser as `--scan` | No |
| `--output-folder` | directory path | output root; created recursively if missing | Yes |
| `--output_internal_novel_peptide` | output file path | write internal TSV; no-dir path resolves to output-folder | No |
| `--internal_novel_peptide` | TSV file path | load internal TSV and skip subtraction | Yes |
| `--stop-after-saving-novel-peptide` | flag | requires `--output_internal_novel_peptide` | No |

Novel option combination:

1. `--novel_protein` and `--novel_peptide` can be used together in one run.
2. They can also be combined with `--scan` and/or `--scan_numbers`.

Hard constraints:

1. Any of `--novel_*` / `--scan*` requires at least one spectrum input file, except stop-after mode (`--stop-after-saving-novel-peptide`).
2. `--novel_*` and `--scan*` are rejected with `-i` or `-j`.
3. Novel mode requires known DB from `--database` or `database_name` in params.
4. Known DB inputs must be all FASTA or all `.idx`.
5. If known DB is `.idx`, all `.idx` must be same index type.
6. `--output_internal_novel_peptide` requires at least one of `--novel_protein` / `--novel_peptide`.
7. `--internal_novel_peptide` is mutually exclusive with `--novel_protein` and `--novel_peptide`.
8. `--stop-after-saving-novel-peptide` requires `--output_internal_novel_peptide`.
9. If `--name` and `--output-folder` are both set, `--name` must not contain path separators.
10. Spectrum input can be omitted only in stop-after mode.

## 3. `--novel_peptide`: FASTA vs Tokenized Text

### 3.1 Auto-detection rule
`--novel_peptide <file>` is parsed as FASTA if the file contains at least one non-empty line whose first non-space char is `>`.

Important:

1. If any `>` line exists, the whole file is treated as FASTA mode.
2. Do not mix FASTA and tokenized-list style in one file.

### 3.2 FASTA mode
Behavior:

1. Header lines start with `>`.
2. Non-header non-empty lines are concatenated as sequence for that entry.
3. Sequence is normalized by keeping alphabetic chars only and converting to uppercase.
4. Duplicates are removed after normalization.

Example FASTA input:

```text
>pep_1
PEPTIDEK
>pep_2
acdm[+16]efg
>pep_3
PEP
TIDE
```

Parsed peptides:

1. `PEPTIDEK`
2. `ACDMEFG`
3. `PEPTIDE`

### 3.3 Tokenized text mode
Definition: plain text peptide tokens split by delimiters.

Supported delimiters:

1. `,`
2. space
3. tab
4. newline (`\n`)
5. CRLF (`\r\n`)

Normalization per token:

1. Keep alphabetic chars only (`A-Z`/`a-z`).
2. Convert to uppercase.
3. Empty-after-normalization tokens are discarded.
4. Duplicates are removed.

Example tokenized input:

```text
PEPTIDEK,pepTIDEL ACD[+57]EFG
ILPEPTIDE<TAB>MPEPTIDE
K.PEPTIDE.R
12345
```

Parsed peptides:

1. `PEPTIDEK`
2. `PEPTIDEL`
3. `ACDEFG`
4. `ILPEPTIDE`
5. `MPEPTIDE`
6. `KPEPTIDER`

Notes:

1. `K.PEPTIDE.R` is not interpreted as flanking-AA syntax; dots are removed, letters remain (`KPEPTIDER`).
2. `12345` becomes empty after normalization and is dropped.
3. Inline comments are not supported in this parser.

### 3.4 Practical recommendation for tokenized text
Use one clean peptide token per entry (or comma-separated list), containing only amino-acid letters. Avoid flanking residues, punctuation, or annotation columns.

## 4. `--scan <file>` File Format
`--scan` reads a text file and parses positive integers with the same parser as `--scan_numbers`.

Supported separators:

1. comma
2. any whitespace (space/tab/newline/CRLF)

Set semantics:

1. Duplicates are deduplicated.
2. Order is irrelevant.

Valid example file:

```text
1001,1002,1003
1004 1005
1006<TAB>1007
1008
```

Equivalent parsed set:

`{1001,1002,1003,1004,1005,1006,1007,1008}`

Invalid examples:

1. `1001;1002` (semicolon is not a delimiter)
2. `1001-1005` (range syntax unsupported)
3. `1001 # comment` (comment text becomes invalid token)
4. `0` or negative numbers (must be `> 0`)

## 5. `--scan_numbers <list>` Accepted List Styles
`--scan_numbers` uses exactly the same integer parser as `--scan`.

Valid CLI examples:

```bash
--scan_numbers 1001,1002,1003
--scan_numbers "1001 1002 1003"
--scan_numbers $'1001\t1002\n1003'
```

Invalid CLI example:

```bash
--scan_numbers 1001;1002
```

Answering format questions directly:

1. `xx;xx` -> not supported.
2. `xx,xx` -> supported.
3. `"xx xx"` -> supported (quoted whitespace-separated list).

## 6. Explicit Scan Set Composition
When both `--scan_numbers` and `--scan` are provided:

1. Parse both to integer sets.
2. Union the sets.
3. Later intersect with requested scan range (`-F/-L` or params `scan_range`) during per-spectrum filtering.

Example:

1. `--scan_numbers 1999,2346,2505,5001`
2. `--scan scan_file.txt` where file adds `2617`
3. `-F2000 -L5000`

Effective explicit scan candidates become:

`{2346,2505,2617}`

(1999 removed by start boundary 2000, 5001 removed by end boundary 5000).

## 7. Novel-aware Scan Reduction (Detailed)

### 7.1 When prefiltering is enabled
Prefilter runs if either condition is true:

1. Novel mode active (`--novel_protein` or `--novel_peptide`).
2. Explicit scan filter active (`--scan` or `--scan_numbers`).

Each input is converted to a temporary filtered MGF and search runs on that subset.

### 7.2 Keep/drop logic per spectrum
A spectrum is kept only if all required checks pass:

```text
keep = (scan_number > 0)
   AND in_requested_scan_range
   AND (not use_explicit_set OR scan in explicit_set)
   AND (not use_novel_mass_filter OR passes_novel_mass_plausibility)
```

Where:

1. `in_requested_scan_range` comes from `-F/-L` or `scan_range` param.
2. `use_novel_mass_filter` is true only in novel mode.

### 7.3 What `passes_novel_mass_plausibility` means
In novel mode, CometPlus first computes retained novel peptides after subtraction vs known DB, then collects their peptide masses (`vNovelMasses`).

A spectrum passes novel plausibility if at least one candidate precursor charge/mass state can match at least one retained novel mass under current Comet tolerance model, including:

1. `digest_mass_range`
2. `peptide_mass_tolerance_lower/upper`
3. `peptide_mass_units` (amu/mmu/ppm)
4. `precursor_tolerance_type`
5. `isotope_error`
6. `mass_offsets`
7. `precursor_charge`, `override_charge`
8. `min_precursor_charge`, `max_precursor_charge`
9. `correct_mass`

And spectrum reading uses `ms_level` filter (`MS2` default, or MS1/MS3 if configured).

### 7.4 Consequences
1. Novel-aware reduction is precursor-mass plausibility filtering, not fragment-level confirmation.
2. Some non-novel spectra may remain if precursor mass overlaps novel mass windows.
3. If no novel peptides remain after subtraction, novel mass set is empty and all spectra fail novel plausibility (effectively zero kept scans), with warning logged.

### 7.5 Combined behavior examples
Example A: only explicit scan filter

1. Command includes `--scan_numbers 1001,1002,1003`
2. No novel options.
3. Kept scans are explicit-set scans that also pass `-F/-L` range.

Example B: only novel options

1. Command includes `--novel_peptide novel.txt`
2. No explicit scan set.
3. Kept scans are those passing range + novel mass plausibility.

Example C: both explicit and novel

1. Command includes `--novel_protein novel.fasta --scan scan_ids.txt`
2. Kept scans must satisfy explicit-set membership and novel mass plausibility simultaneously.

## 8. Known/Novel Subtraction Semantics

### 8.1 Known DB handling
1. Known `.idx` files are reused directly (no mutation).
2. Known FASTA can be indexed temporarily for extraction when needed.
3. Known DB type mixing (FASTA + `.idx`) is rejected.

### 8.2 Novel candidate sources
1. `--novel_protein`: digest protein FASTA using current Comet settings and extract peptide candidates from the generated peptide index.
2. `--novel_peptide`: parse peptide input into peptide candidates (FASTA/text parsing only; no digestion step).

### 8.2.1 Mixed `--novel_protein` + `--novel_peptide` flow
1. CometPlus first adds peptide candidates from `--novel_protein`.
2. Then it adds peptide candidates from `--novel_peptide`.
3. Both sources are merged into one peptide-level map keyed by normalized peptide identity (`equal_I_and_L` aware), so overlap across the two sources is deduplicated before subtraction.
4. Known-db subtraction runs once on the merged set, yielding one retained novel peptide set.

### 8.3 Subtraction key
Comparison uses normalized peptide key with `equal_I_and_L` policy:

1. `equal_I_and_L=1`: `I` and `L` treated as equivalent.
2. `equal_I_and_L=0`: exact identity.

Retained novel set is candidates not found in known set under this policy.

### 8.4 Novel name generation and duplicate protein headers
Name generation for novel entries:

1. `--novel_protein` FASTA headers are used only in the intermediate digestion/indexing stage, not as final novel scoring entry names.
2. `--novel_peptide` FASTA headers are also not preserved as final search protein names.
3. After merged-set subtraction, retained novel peptides are written to one temporary FASTA with synthetic headers:
   `>COMETPLUS_NOVEL_<n>`.
4. The retained novel scoring database (after subtraction) is the same `COMETPLUS_NOVEL_<n>`-header FASTA (or its `.idx` form, when known DB input is `.idx`).
5. Result: every retained novel peptide from either `--novel_protein` or `--novel_peptide` is represented as `COMETPLUS_NOVEL_<n>` on the novel side.

Duplicate protein header behavior in `--novel_protein` input:

1. CometPlus does not perform header-name deduplication or auto-renaming in `--novel_protein` FASTA.
2. Duplicate FASTA headers are accepted and processed; indexing/search identity is tracked by internal positions/protein-set references, not by requiring globally unique header strings.
3. During novel-mode aggregation, candidates are collapsed at peptide level (normalized peptide key), so overlap from duplicated protein entries mainly affects peptide multiplicity before deduplication, not final retained peptide uniqueness.

## 9. Logging and Error Patterns (Useful for Users)
Representative runtime messages:

1. `Error - no peptide entries were parsed from --novel_peptide input.`
2. `Error - invalid scan token "..." in --scan_numbers/--scan.`
3. `Error - scan number out of range ...`
4. `Error - --novel_protein/--novel_peptide/--scan/--scan_numbers cannot be used with -i or -j.`
5. `Warning - no novel peptides remain after subtraction ...`
6. `Scan prefilter: "<input>" -> <N> scans retained.`

## 10. Recommended CLI Examples

### 10.1 Novel peptide from tokenized text + explicit inline scans
```bash
cometplus \
  --database known.fasta \
  --novel_peptide novel_peptides.txt \
  --scan_numbers "2346 2505 2617" \
  -F2000 -L5000 \
  sample.mzML
```

### 10.2 Novel protein + scan file
```bash
cometplus \
  --database known.fasta \
  --novel_protein novel_proteins.fasta \
  --scan scan_ids.txt \
  sample.mzML
```

### 10.3 Combine `--scan` and `--scan_numbers`
```bash
cometplus \
  --database known.fasta \
  --scan scan_ids.txt \
  --scan_numbers 1001,1002,1003 \
  sample.mzML
```

The two sources are unioned first, then intersected by range constraints and optional novel mass plausibility.

## 11. Output Routing And Conflict Precheck

### 11.1 `--output-folder` routing
1. Default output folder is current directory (`.`).
2. CometPlus rewrites each input basename to `<output-folder>/<input-stem>`.
3. If `--name` is used (single-input mode), basename is `<output-folder>/<name>`.
4. Result: output files are no longer written to input spectrum directories.

### 11.2 Planned output conflict detection
Before search starts, CometPlus computes all planned outputs for each input:

1. `sqt`
2. `txt` (with configured text extension)
3. `pep.xml`
4. `mzid`
5. `pin`
6. decoy-side files when `decoy_search=2`
7. internal TSV target when `--output_internal_novel_peptide` is set

Then two checks run:

1. Internal collisions inside this invocation (same path planned multiple times).
2. Existing files already on disk.

Any collision prints a conflict list and exits before search/prefilter.

## 12. Internal Novel Peptide TSV (Export And Reuse)

### 12.1 Export: `--output_internal_novel_peptide`
This file stores retained novel peptides after known-db subtraction.

Format is fixed TSV with header:

```text
peptide	peptide_id	protein_id
```

Column semantics:

1. `peptide`: normalized peptide amino-acid sequence.
2. `peptide_id`: `COMETPLUS_NOVEL_<n>` for newly generated records.
3. `protein_id`:
   - novel_protein source: source protein IDs, semicolon-separated.
   - novel_peptide source: normalized peptide sequence itself.
   - mixed source: merged semicolon list.

Path behavior:

1. If path includes directory, missing directories are created recursively.
2. If path has no directory component, file is written under `--output-folder`.

### 12.2 Stop-after mode
`--stop-after-saving-novel-peptide`:

1. Requires `--output_internal_novel_peptide`.
2. Exits with code `0` right after TSV is written.
3. Skips subsequent spectrum prefilter and search stages.
4. Can be used without spectrum inputs.

### 12.3 Reuse: `--internal_novel_peptide`
1. Input must be the TSV format produced above.
2. Required columns: `peptide`, `peptide_id`, `protein_id`.
3. Imported records preserve `peptide_id` values from file (no renumbering).
4. This mode is mutually exclusive with `--novel_protein` / `--novel_peptide`.
5. Known-vs-novel subtraction is skipped.
6. Workflow continues from novel mass collection + prefilter + search.

## 13. Novel-Only Spectrum Output Filtering

In novel mode, all output writers apply query-level gate:

1. `txt`
2. `sqt`
3. `pepXML`
4. `mzid`
5. `pin`

Rule:

1. A spectrum is retained only if at least one target-side printable PSM maps to `COMETPLUS_NOVEL_...`.
2. Otherwise the whole spectrum is omitted from all outputs.
3. With `decoy_search=2`, once a spectrum passes target-side gate, target and decoy sides are both kept for that spectrum.

Non-novel mode behavior is unchanged.

## 14. Temporary File Directory Policy

To avoid system tmp usage:

1. CometPlus temporary artifacts are created under `--output-folder` root:
   - merged FASTA
   - temporary params
   - temporary index files
   - filtered MGF
   - novel scoring FASTA
2. CometSearch-side temporary inflated files (`*.mzML.gz`, etc.) also land under output basename directory, which now points to `--output-folder`.

## 15. Runtime Logging And Performance Interpretation

CometPlus logs timestamped stage messages with elapsed seconds for major steps:

1. known peptide extraction
2. novel candidate assembly
3. subtraction
4. internal TSV import/export
5. novel mass calculation
6. per-input scan prefilter and total prefilter time
7. search launch summary

Observed thread behavior (`--thread 1` vs `--thread 20`) can be close in novel workflows because large portions of runtime are outside main search threads (index generation subprocesses, I/O, single-thread prefilter). Thread override is now forwarded to temporary index-generation subprocesses to reduce this gap.
