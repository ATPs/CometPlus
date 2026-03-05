# pyComet

`ProtCosmo/pyComet` provides Python tools to process Comet-related files:

1. Convert PIN/percolator files across text/gzip/parquet formats and normalize peptide modifications to UniMod IDs.
2. Validate/resolve active modifications in `comet.params` against UniMod.
3. Build Comet-like peptide index TSV tables from FASTA + `comet.params`.
4. Build `peptide/peptide_id/protein_id` internal novel peptide TSV from digested database and/or peptide inputs.

## Files In This Folder

- `pyCometPinConverter.py`: PIN/percolator format conversion + peptide/protein normalization.
- `pyUnimodMatch.py`: active Comet modification to UniMod ID matching report.
- `pyPeptideIndex.py`: peptide index table generator.
- `pyCometProteinDigest.py`: digest database and/or peptide inputs to internal novel peptide TSV.
- `pyCometPinConverter.md`: detailed notes for `pyCometPinConverter.py`.
- `pyPeptideIndex.md`: detailed schema/behavior for `pyPeptideIndex.py`.
- `pyCometProteinDigest.md`: detailed usage notes for `pyCometProteinDigest.py`.
- `resource/common_unimod.csv`: default UniMod mapping table used by all tools.
- `utils/`: shared modules (`pyLoadParameters`, digestion/enumeration helpers, UniMod matcher, converter helpers).

## Requirements

- Python 3
- Optional but needed for parquet conversion in `pyCometPinConverter.py`: `pandas`, `pyarrow`
- Optional for nicer progress in `pyPeptideIndex.py`: `tqdm`

Run from repo root:

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/<script>.py --help
```

## `pyCometPinConverter.py`

Convert among `.pin`, `.pin.gz`, `.pin.parquet`, `.pin.parquet.gz` (also accepts `.tsv`/`.parquet` variants).
CLI:

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py \
  (--input <file> | --input-str <peptide>) \
  [--input-percolator] [--output <file>] \
  [--proteins-keep N] [--unimod <csv_or_tsv>] [--unimod-ppm 10]
```

Behavior summary:

1. PIN text input (`.pin`/`.pin.gz`): `Peptide` flanks are removed and `[mass]` tags are mapped to `[U:id]`; `Proteins` are merged with commas and truncated by `--proteins-keep` (default `1`, `0` keeps all).
2. PIN parquet input: `Peptide` text is kept unchanged, and `Proteins` are normalized with `--proteins-keep`.
3. Percolator mode (`--input-percolator`): non-parquet input is parsed as TSV regardless of extension; only `peptide` and `proteinIds` values are normalized; column names are preserved.
4. Extra columns are preserved.

Default output path:

1. If input is text (`.pin`, `.pin.gz`, or percolator TSV) and `--output` is omitted, output defaults to parquet (`<input>.parquet`).
2. If input is parquet (`.parquet`/`.parquet.gz`), `--output` is required.

Examples:

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py --input x.pin
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py --input x.pin.gz --output x.pin.parquet.gz
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py --input x.pin.parquet --output x.pin
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py --input xx.decoy.psms.tsv --input-percolator
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometPinConverter.py --input-str '-.n[42.0106]M[15.9949]LQFLLEVNK.S'
```

Example conversion result:

- Input: `R.LHWLVM[15.9949]RK.N`
- Output: `LHWLVM[U:35]RK`

## `pyUnimodMatch.py`

Check active modifications in `comet.params` and resolve them to UniMod IDs.

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyUnimodMatch.py \
  [-P comet.params] [--unimod ProtCosmo/pyComet/resource/common_unimod.csv] [--ppm 10]
```

Behavior summary:

1. Reads active variable mods and non-zero active fixed mods (`add_*`, terminal mods).
2. Matches by residue compatibility and mass ppm tolerance.
3. Prints variable/fixed mapping results.
4. Exits with code `1` when any active mod cannot be matched.
5. If `--unimod` is `None` or missing, prints a message and exits successfully.

## `pyPeptideIndex.py`

Generate Comet-like peptide index tables from FASTA plus `comet.params`.

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyPeptideIndex.py \
  [-P comet.params] \
  [-D db1.fasta -D db2.fasta | --protein SEQ1 --protein SEQ2] \
  [-N out/comet] [--max-record N] [--use-protein-name] [--thread N] \
  [--unimod ProtCosmo/pyComet/resource/common_unimod.csv] [--unimod-ppm 10]
```

It writes these TSV files (header line starts with `#`):

1. `<prefix>.index_run.tsv`
2. `<prefix>.static_mod.tsv`
3. `<prefix>.variable_mod.tsv`
4. `<prefix>.protein.tsv`
5. `<prefix>.peptide_sequence.tsv`
6. `<prefix>.peptide_sequence_protein.tsv`
7. `<prefix>.peptide_protein_location.tsv`
8. `<prefix>.peptide_variant.tsv`
9. `<prefix>.peptide_variant_mod.tsv`

Notes:

1. If no `-D/--database` or `--protein` is given, it falls back to `database_name` in `comet.params`.
2. `--thread 0` uses all CPUs; `--thread 1` disables multiprocessing.
3. If a valid UniMod file is provided and any active mod is unmatched, indexing exits early.
4. If UniMod is missing or set to `None`, indexing continues without UniMod IDs.

See `pyPeptideIndex.md` for full column-by-column schema and algorithm details.

## `pyCometProteinDigest.py`

Generate a 3-column internal novel peptide TSV from digested FASTA and/or peptide inputs:

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometProteinDigest.py \
  [-P comet.params] [-D db1.fasta -D db2.fasta] [--max-record N] \
  [--peptide <item> ...] [--peptide-exclude <exclude.tsv|exclude.sqlite>] \
  -O output.tsv
```

Output header and column order are fixed:

```text
peptide	peptide_id	protein_id
```

Behavior summary:

1. `peptide_id` is assigned as `COMETPLUS_NOVEL_<n>` in first-seen order.
2. `--database` peptides are generated by Comet digestion rules from `comet.params`.
3. `--peptide` accepts sequence tokens, token files, or FASTA files; `--database` and `--peptide` can be used together.
4. `--peptide-exclude` is optional and supports TSV or SQLite inputs based on `pep_seq`.
5. If one peptide maps to multiple proteins, `protein_id` uses `;`-joined IDs.

See `pyCometProteinDigest.md` for full details and examples.
