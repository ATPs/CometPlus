# CometPlus

CometPlus is a thin CLI layer on top of CometSearch that keeps upstream-compatible search behavior while adding quality-of-life options (notably repeated `--database`) used by ProtCosmo workflows.

Supported spectrum inputs:
- `mzXML`
- `mzML`
- `mzMLb` (requires HDF5-enabled MSToolkit build)
- `mgf`
- `mzXML.gz`
- `mzML.gz`
- `mgf.gz`
- `ms2` variants (`ms2`, `cms2`, `bms2`)
- Thermo RAW (platform/toolchain dependent)

## Build Requirements

- Linux/macOS toolchain with C++14 support
- Existing repo subprojects (`MSToolkit`, `AScorePro`, `CometSearch`)
- For fully static `mzMLb` binaries: static HDF5 build with C and C++ libraries

For static Linux builds, this project uses:
- `x86_64-conda-linux-gnu-g++`
- `x86_64-conda-linux-gnu-gcc`

Reproducible static build notes and validation commands are documented in:
- `ProtCosmo/CometPlus/build.md`

## Build Matrix

### 1) Default build
```bash
make -C ProtCosmo/CometPlus
```

### 2) Static build (no mzMLb/HDF5 requirement)
```bash
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static WITH_MZMLB=0
```

### 3) Static build with mzMLb enabled
```bash
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp
```

Optional explicit overrides:
```bash
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/abs/path/to/hdf5-static-cxx \
  HDF5_LIB_C=/abs/path/to/hdf5-static-cxx/lib/libhdf5.a \
  HDF5_LIB_CPP=/abs/path/to/hdf5-static-cxx/lib/libhdf5_cpp.a
```

## HDF5 Prerequisite for mzMLb

Environment-specific note:
- tested mzMLb-enabled static prefix: `/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp`
- do not use `/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp` for mzMLb runtime in this environment (observed: `required filter 'deflate' is not registered`)

`WITH_MZMLB=1` requires both of these files:
- `include/H5Cpp.h`
- `lib/libhdf5_cpp.a`

It also requires HDF5 to have deflate/zlib filter support enabled at build time. If deflate support is missing, mzMLb reads can fail at runtime with:
- `required filter 'deflate' is not registered`

If your HDF5 install is C-only, build a static C++ HDF5 install (example):
```bash
cd /data/p/hdf5/hdf5-1.14.4-3
mkdir -p /data/p/hdf5/hdf5-1.14.4-3_build_static_cpp

CC=x86_64-conda-linux-gnu-gcc \
CXX=x86_64-conda-linux-gnu-g++ \
./configure \
  --prefix=/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp \
  --enable-cxx --enable-hl \
  --enable-static --disable-shared \
  --disable-tools --disable-tests \
  --with-zlib=/data/p/anaconda3 \
  --with-szlib=no

make -j"$(nproc)"
make install
```

## Verify Static Binary

```bash
file ProtCosmo/CometPlus/cometplus
ldd ProtCosmo/CometPlus/cometplus
```

Expected for static target:
- `statically linked`
- `not a dynamic executable`
- binary size is typically much larger than dynamic builds (for example, around `10-15 MB` vs around `4 MB`)

For `WITH_MZMLB=0` (no mzMLb), confirm no HDF5 link symbols:
```bash
nm -A ProtCosmo/CometPlus/cometplus 2>/dev/null | rg -i "\\bH5[A-Za-z0-9_]*\\b|hdf5"
```
Expected output: no matches.

For `WITH_MZMLB=1` (mzMLb enabled), the same command should show `H5...` symbols.

## Gzip Input Notes

- `mzXML.gz`, `mzML.gz`, and `mgf.gz` are handled by CometPlus/CometSearch using temporary inflated files:
  - `.mzXML.gz` -> temp `.mzXML`
  - `.mzML.gz` -> temp `.mzML`
  - `.mgf.gz` -> temp `.mgf`
- Temporary file behavior:
  - temporary file directory is derived from `--output-folder` (default current directory),
  - the original input path remains the original `.gz` path for reporting/metadata,
  - temporary file is removed automatically at normal completion and handled error exits.
- `ms2.gz` is not supported in this milestone.

## Usage

### Help
```bash
./ProtCosmo/CometPlus/cometplus --help
```

### Full Help And Generic Param Overrides
```bash
./ProtCosmo/CometPlus/cometplus --params /path/to/comet.params --help-full
```

`--help-full` prints a long per-key guide. For each `--<param_key>`, it shows:
- the current/default example value parsed from the params file
- the inline `# ...` guidance comment from that params line when present
- fallback built-in examples/comments when the params file is not readable

Generic override syntax:
- `--<param_key> <value>`
- `--<param_key>=<value>`

Rules:
- `param_key` must exist in the loaded params file before `[COMET_ENZYME_INFO]`.
- Multi-token values must be quoted as one shell argument.
- Dedicated overrides are not accepted as generic keys:
  - `database_name` -> use `--database`
  - `num_threads` -> use `--thread`
  - `scan_range` -> use `--first-scan`/`--last-scan`
  - `spectrum_batch_size` -> use `-B`

Examples:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --decoy_search 2 \
  --fragment_bin_tol 0.02 \
  --digest_mass_range "600 5000" \
  --variable_mod01 "15.994915 M 0 3 -1 0 0  0.0 # Oxidation (M)" \
  /path/to/input.mzML
```

### Single database search
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/proteins.fasta \
  --name run1 \
  /path/to/input.mzML
```

### Multiple FASTA databases (CometPlus merges internally)
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/target.fasta \
  --database /path/to/decoy.fasta \
  --name run_target_decoy \
  /path/to/input.mzMLb
```

### Multiple `.idx` databases
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/db1.idx \
  --database /path/to/db2.idx \
  --name run_multi_idx \
  /path/to/input.mgf
```

## Novel Inputs and Scan Subset Options

This section documents design and usage for:
- `--novel_protein <file>`
- `--novel_peptide <file>`
- `--scan <file>`
- `--scan_numbers <list>`
- `--output-folder <dir>`
- `--output_internal_novel_peptide <file>`
- `--internal_novel_peptide <file>`
- `--stop-after-saving-novel-peptide`
- `--keep-tmp`
- `--run-comet-each`

These options are implemented as additive orchestration in CometPlus and preserve normal Comet behavior when not used.

### Option Semantics

- `--novel_protein <file>`
  - Input must be FASTA.
  - Proteins are digested with the active Comet settings to generate peptide candidates.
  - Candidate peptides are merged with `--novel_peptide` candidates (if provided), then one subtraction pass is applied against known DB peptides.

- `--novel_peptide <file>`
  - Input supports either FASTA or tokenized text.
  - FASTA mode is auto-detected if any non-empty trimmed line begins with `>`.
  - Tokenized mode accepts delimiters: comma, space, tab, newline.
  - In mixed mode (`--novel_protein` + `--novel_peptide`), both sources share one deduplicated candidate map before subtraction.

- `--output-folder <dir>`
  - Default: current directory (`.`).
  - All normal search outputs are routed here (instead of input spectrum directories).
  - CometPlus temporary artifacts used by novel/scan workflows are also created under this directory.
  - In novel mode with multiple spectrum inputs:
    - default path: filtered MGF shards are merged and searched once,
    - `--run-comet-each` path (peptide `.idx` known DB only): filtered shards are grouped into balanced merged task-MGF files, then searched in child processes and merged as one final pin.
  - Default merged output base is `<output-folder>/cometplus_novel_merged` unless `--name` is provided.

- `--output_internal_novel_peptide <file>`
  - Exports subtraction result as detailed internal TSV with columns:
    `peptide`, `peptide_id`, `protein_id`, `peptide_with_mod`, `charge`, `mz`, `mz_window_min`, `mz_window_max`.
  - Requires at least one of `--novel_protein` / `--novel_peptide`.
  - If `<file>` contains no directory component, file is created under `--output-folder`.
  - If `<file>` includes a directory path, missing directories are created recursively.

- `--internal_novel_peptide <file>`
  - Reuses an exported internal TSV and skips novel-vs-known subtraction.
  - Accepts both legacy 3-column TSV (`peptide`, `peptide_id`, `protein_id`) and the new detailed TSV format.
  - When detailed TSV is provided, precomputed masses are imported directly for prefilter fast-path.
  - Mutually exclusive with `--novel_protein` and `--novel_peptide`.
  - If `<file>` is relative, it is resolved from current working directory (not `--output-folder`).

- `--stop-after-saving-novel-peptide`
  - Must be used with `--output_internal_novel_peptide`.
  - Workflow exits successfully after saving internal TSV and does not perform spectrum prefilter/search.

- `--keep-tmp`
  - Keeps CometPlus temporary artifacts on exit for debugging.
  - By default (without this flag), temporary artifacts are removed on exit.

- `--run-comet-each`
  - Optional acceleration mode for novel multi-input runs with known peptide `.idx` databases (`-j` type).
  - Parent process still runs novel workflow + prefilter, then builds grouped merged task-MGF files from filtered shards.
  - One child `cometplus` process is launched per grouped task-MGF, each producing a task `.pin`; final output is one merged `.pin`.
  - Only pin output is guaranteed in this mode (other output formats are disabled for task child jobs).
  - If prerequisites are not met, CometPlus logs a warning and falls back to the normal merged-MGF single-search path.
  - Grouping and thread policy:
    - total threads `T`: use `--thread` when `>0`, else params `num_threads` when `>0`, else error.
    - target grouped task count `N`: `N = max(1, T / 4)` (integer division), then cap by filtered shard count (`N <= shard_count`).
    - if `N == shard_count`, no grouping merge is needed (each shard is one task).
    - when grouping is needed, larger shards are assigned first to the current smallest group to balance merged task-MGF sizes.
    - `T` is split across `N` tasks as evenly as possible (difference at most 1 thread per task).
    - child process concurrency equals grouped task count `N`.

- `--scan <file>`
  - File-based explicit scan list.
  - Tokens use the same delimiters as `--scan_numbers`.

- `--scan_numbers <list>`
  - Inline explicit scan list, for example: `1001,1002,1003`.

### End-to-End Design Behavior

When any novel option or explicit scan option is used, CometPlus runs this flow:

1. Resolve known DB input from repeatable `--database` values (or `database_name` in params if CLI does not provide one).
2. Resolve output root (`--output-folder`, default `.`), create directory if needed, and route output basenames to this directory.
3. Pre-check planned output collisions:
   - internal collisions within this run (same target path planned multiple times),
   - pre-existing files on disk.
   Any collision causes early error and exit.
4. Build novel candidates:
   - from `--novel_protein` digestion and/or `--novel_peptide` parsing, then merge and deduplicate by normalized peptide identity,
   - or load from `--internal_novel_peptide` and skip subtraction (detailed TSV can also provide precomputed masses).
5. Parse known peptide universe for subtraction when using fresh novel inputs:
   - known `.idx`: read directly, no rebuild and no mutation.
   - known FASTA: temporary peptide index generation is used for extraction.
6. Normalize peptide identity for subtraction using `equal_I_and_L`:
   - `equal_I_and_L=1`: `I` and `L` are equivalent.
   - `equal_I_and_L=0`: `I` and `L` are distinct.
7. Remove novel candidates already present in known DB(s) (skipped when `--internal_novel_peptide` is used).
8. Assign/retain novel IDs:
   - fresh subtraction records are named `COMETPLUS_NOVEL_<n>`,
   - imported internal records preserve input `peptide_id`.
9. Materialize retained novel peptides as temporary scoring DB input and compute retained novel masses.
10. Optionally write detailed internal TSV (`--output_internal_novel_peptide`) with columns:
   - `peptide`: normalized peptide sequence,
   - `peptide_id`: `COMETPLUS_NOVEL_<n>` or imported ID,
   - `protein_id`: semicolon-separated source list (protein IDs for protein source; peptide sequence for peptide source),
   - `peptide_with_mod`: variable-mod mass-annotated peptide string,
   - `charge`: charge state used for window expansion,
   - `mz`: theoretical precursor m/z for (`peptide_with_mod`, `charge`),
   - `mz_window_min` / `mz_window_max`: tolerance-only m/z window bounds.
   Isotope-error and mass-offset expansion are still applied at runtime prefilter.
11. If `--stop-after-saving-novel-peptide` is set, exit after step 10.
12. Parse explicit scans from `--scan_numbers` and `--scan`, union and deduplicate.
13. Intersect explicit scans with any scan-range constraint (`-F/-L` or `--first-scan/--last-scan`).
14. Filter spectra into temporary MGF files (parallel workers based on `--thread` or `num_threads`):
   - by explicit scan set if provided,
   - and, in novel mode, by precursor-mass plausibility against retained novel peptide masses.
15. If novel mode has more than one spectrum input:
   - default path: merge filtered MGFs into one temporary MGF while preserving each spectrum `TITLE=`.
   - `--run-comet-each` path (peptide `.idx` known DB only): regroup filtered shards into `N` balanced merged task-MGF files (`N = max(1, total_threads / 4)`, capped by shard count).
16. Search:
   - default merged novel path: one search input (one output set),
   - `--run-comet-each` path: one child search per grouped task-MGF and merged final `.pin`,
   - otherwise: one search input per filtered file (legacy behavior).
17. In novel mode, output files (`txt/sqt/pepXML/mzid/pin`) keep only spectra that pass target-side printable PSM gate:
   - no printable novel target PSM => drop spectrum,
   - printable target PSM all novel-only => keep spectrum,
   - printable target PSM mixed novel+known => keep only if every best-`Xcorr` tie is novel-only (no known assignment in best-score ties).
   With `decoy_search=2`, target and decoy records are retained/dropped together at spectrum level.
18. In merged novel mode, Percolator `SpecId` keeps `<base>_<scan>_<charge>_<rank>` format, and `base` is read from each spectrum `TITLE=` prefix. Input basenames should therefore be unique across merged inputs.

### Input Validation and Error Conditions

- Scan tokens must be positive integers (`1..INT_MAX`); malformed or out-of-range tokens fail fast.
- `--novel_peptide` with no parsed peptide entries fails fast.
- `--novel_protein/--novel_peptide/--scan/--scan_numbers` (and internal TSV export) cannot be used with index-creation modes `-i` or `-j`.
- Novel mode requires a known database source (`--database` or `database_name` in params).
- At least one spectrum input file is required for novel/scan-subset runs, except `--stop-after-saving-novel-peptide`.
- `--output_internal_novel_peptide` requires `--novel_protein` and/or `--novel_peptide`.
- `--internal_novel_peptide` cannot be combined with `--novel_protein` or `--novel_peptide`.
- `--stop-after-saving-novel-peptide` requires `--output_internal_novel_peptide`.
- When `--name` and `--output-folder` are both specified, `--name` must be a base name (no path separators).
- `--name` with multiple spectrum inputs is accepted only for merged multi-input novel mode; otherwise it remains single-input only.
- `--run-comet-each` is effective only when all prerequisites are met (novel mode + multi-input + known peptide `.idx`); otherwise it logs warning and falls back.
- Effective `--run-comet-each` requires positive total threads from `--thread` or params `num_threads`.

### Usage Examples

Only explicit scans from inline list:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.fasta \
  --scan_numbers 1001,1002,1003 \
  /path/to/input.mzML
```

Explicit scans from file + scan-range intersection:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.fasta \
  --scan /path/to/scan_ids.txt \
  --first-scan 2000 \
  --last-scan 5000 \
  /path/to/input.mgf
```

Novel peptide text with explicit scan subset:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.idx \
  --novel_peptide /path/to/novel_peptides.txt \
  --scan_numbers 2104,2456,3001 \
  /path/to/input.mzML.gz
```

Combined novel protein + novel peptide + dual scan sources:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known_target.idx \
  --database /path/to/known_decoy.idx \
  --output-folder /path/to/out \
  --novel_protein /path/to/novel_proteins.fasta \
  --novel_peptide /path/to/novel_peptides.fasta \
  --scan /path/to/scan_ids.txt \
  --scan_numbers 3001,3002,3003 \
  --first-scan 2500 \
  --last-scan 6000 \
  --name run_novel_scan_subset \
  /path/to/input.mzMLb
```

Multi-input novel mode (parallel prefilter + merged single search):
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.idx \
  --output-folder /path/to/out \
  --novel_protein /path/to/novel_proteins.fasta \
  --name merged_run \
  /path/to/input_a.mzMLb /path/to/input_b.mzMLb
```

Multi-input novel peptide-idx mode with grouped-task cometplus and merged pin:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known_a.idx \
  --database /path/to/known_b.idx \
  --internal_novel_peptide /path/to/internal_novel.tsv \
  --run-comet-each \
  --thread 32 \
  --output-folder /path/to/out \
  --name merged_run \
  /path/to/input_a.mzMLb /path/to/input_b.mzMLb
```

Export internal novel TSV and stop before search:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.idx \
  --output-folder /path/to/out \
  --novel_protein /path/to/novel_proteins.fasta \
  --output_internal_novel_peptide internal_novel.tsv \
  --stop-after-saving-novel-peptide
```

Reuse internal novel TSV and run search directly:
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.idx \
  --output-folder /path/to/out \
  --internal_novel_peptide /path/to/out/internal_novel.tsv \
  /path/to/input.mzMLb
```

### Related Design Docs

- `ProtCosmo/CometPlus/design.novel_protein_peptide.md`
- `ProtCosmo/CometPlus/design.multiple.input.db.md`

## Troubleshooting

- `WITH_MZMLB=1 requires static HDF5 C++ library`:
  your `HDF5_DIR` does not contain `libhdf5_cpp.a`; rebuild/install HDF5 with `--enable-cxx` and static libs.

- Build succeeds but `mzMLb` input fails at runtime:
  confirm MSToolkit was rebuilt with HDF5 (`WITH_MZMLB=1` on CometPlus build invocation) and that input extension is `.mzMLb`.

- Link errors from HDF5 dependency chain:
  add extra link flags with `HDF5_EXTRA_LIBS` (default includes `-ldl -lz -lm -lpthread`).

- Mixed database types error:
  use all FASTA paths or all `.idx` paths in one invocation; do not mix them.

## Reproducible Parity Workflow

A full all-scan parity harness for `mzML` vs `mzMLb` vs `mgf` is available at:
- `ProtCosmo/CometPlus/test_mzmlb/run_mzmlb_comparison.py`

Study design and artifact documentation:
- `ProtCosmo/CometPlus/test_mzmlb/README.md`

## Terminology

- `Peptide string`:
  Comet-style peptide token that may include flanking residues and modifications, for example `K.AC[57.0215]DEFGHIK.R`.
- `Flanking residues`:
  The residues outside the matched peptide (`K.` prefix and `.R` suffix in the example above).
- `Core peptide sequence`:
  The middle peptide segment between flanks, including modification annotations if present (`AC[57.0215]DEFGHIK` in the example above).
- `Core peptide (I/L-equivalent)`:
  Core peptide after treating `I` and `L` as the same symbol for equivalence testing.
- `Top hit (Top-1)`:
  Best-scoring peptide-spectrum match for one scan.
- `PSM`:
  Peptide-spectrum match entry; one scan can have multiple ranked PSMs.
- `Label`:
  Target/decoy assignment in `.pin` (commonly `1` for target and `-1` for decoy).
