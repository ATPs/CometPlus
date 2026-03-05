# Comet + CometPlus

Forked from [UWPR/Comet](https://github.com/UWPR/Comet).

This repository contains:
- **Comet**: the core MS/MS search engine (`comet.exe`).
- **CometPlus**: an optional CLI layer in [ProtCosmo/CometPlus](ProtCosmo/CometPlus) that preserves Comet behavior and adds workflow-oriented features.

## Quick Start

### Build Comet (root project)
```bash
make
```
Produces `comet.exe` in the repository root.

### Build CometPlus
```bash
make -C ProtCosmo/CometPlus
```
Produces:
- `ProtCosmo/CometPlus/cometplus`

### Show CometPlus help
```bash
./ProtCosmo/CometPlus/cometplus --help
./ProtCosmo/CometPlus/cometplus --params /path/to/comet.params --help-full
```

## CometPlus Capability Highlights

CometPlus now includes:
- Repeatable `--database` with multi-FASTA or multi-`.idx` search in one invocation.
- Multiple input spectrum files in one invocation.
- `--input_files <list_file>` for large spectrum input lists (one path per line).
- Generic CLI override for normal Comet params: `--<param_key> <value>` / `--<param_key>=<value>`.
- Dedicated overrides for common keys:
  - `database_name` -> `--database`
  - `num_threads` -> `--thread`
  - `scan_range` -> `--first-scan` / `--last-scan`
  - `spectrum_batch_size` -> `-B`
- Novel workflow options:
  - `--novel_protein`, `--novel_peptide`
  - `--output_internal_novel_peptide`, `--internal_novel_peptide`
  - `--stop-after-saving-novel-peptide`
  - `--run-comet-each`
- Explicit scan subset options: `--scan`, `--scan_numbers`.
- Output routing and temp-artifact control:
  - `--output-folder`
  - `--keep-tmp`
- Gzip spectrum input support: `.mzXML.gz`, `.mzML.gz`, `.mgf.gz`.
- `mzMLb` input support (when built with `WITH_MZMLB=1` and static HDF5 C/C++ libs).
- mzMLb process-parallel prefilter worker mode using the same executable (`cometplus --job <job.tsv>`).

## CometPlus Input / Database Support

### Supported spectrum formats
- `mzXML`
- `mzML`
- `mzMLb` (requires `WITH_MZMLB=1`)
- `mzXML.gz`
- `mzML.gz`
- `mgf`
- `mgf.gz`
- `ms2`, `cms2`, `bms2`
- Thermo RAW (platform/toolchain dependent)

Notes:
- gzip support is currently `mzXML.gz`, `mzML.gz`, `mgf.gz`.
- `ms2.gz` is currently **not supported**.

### Database inputs
- `--database` is repeatable.
- You can pass multiple FASTA files, or multiple `.idx` files.
- Mixed types in one run are rejected (FASTA + `.idx` is not allowed).
- For multi-`.idx`, all `.idx` files must be the same index type.
- If `--database` is not passed, CometPlus can fall back to `database_name` in params.

## Common CometPlus Commands

### Single input + single database
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/proteins.fasta \
  --name run1 \
  /path/to/input.mzML
```

### Generic parameter overrides
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/proteins.fasta \
  --decoy_search 2 \
  --fragment_bin_tol 0.02 \
  --digest_mass_range "600 5000" \
  /path/to/input.mzML
```

### Multiple FASTA databases
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/target.fasta \
  --database /path/to/decoy.fasta \
  /path/to/input.mzMLb
```

### Multiple `.idx` databases
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/db1.idx \
  --database /path/to/db2.idx \
  /path/to/input.mgf
```

### Large input list via `--input_files`
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/db1.idx \
  --database /path/to/db2.idx \
  --input_files /path/to/input_list.txt
```

`--input_files` rules:
- one spectrum input path per line
- empty lines and lines whose first non-space character is `#` are ignored
- relative paths are resolved from current working directory
- cannot be combined with positional spectrum inputs
- cannot be passed more than once
- list file with no valid input paths is rejected

## Novel + Scan Workflow

Novel/scan options:
- `--novel_protein <file>`: novel protein FASTA (digested under active Comet settings).
- `--novel_peptide <file>`: FASTA or tokenized peptide text.
- `--scan <file>`: scan list file (comma/space/tab/newline delimiters).
- `--scan_numbers <list>`: inline scan list.
- `--output_internal_novel_peptide <file>`: export detailed internal TSV.
- `--internal_novel_peptide <file>`: reuse internal TSV and skip subtraction.
- `--stop-after-saving-novel-peptide`: stop after internal TSV export.
- `--run-comet-each`: grouped child-search path for multi-input novel mode on peptide `.idx`.
- `--output-folder <dir>` and `--keep-tmp`.

High-level behavior in novel mode:
1. Resolve known DB (`--database` or params fallback).
2. Build novel candidates (`--novel_protein` and/or `--novel_peptide`, or load `--internal_novel_peptide`).
3. Subtract known peptides from novel candidates (unless internal TSV reuse mode).
4. Build temporary novel scoring DB and mass set.
5. Apply explicit scan filtering (`--scan` / `--scan_numbers`) and scan-range intersection.
6. Prefilter spectra, then search.
7. In multi-input novel mode:
   - default path: merge filtered MGF shards and run one merged search;
   - `--run-comet-each` path: regroup filtered shards, run child searches, and merge to one final `.pin`.

### Example: combined novel + scan filtering
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known_target.idx \
  --database /path/to/known_decoy.idx \
  --output-folder /path/to/out \
  --novel_protein /path/to/novel_proteins.fasta \
  --novel_peptide /path/to/novel_peptides.txt \
  --scan /path/to/scan_ids.txt \
  --scan_numbers 3001,3002,3003 \
  --first-scan 2500 \
  --last-scan 6000 \
  --name run_novel_scan_subset \
  /path/to/input.mzML.gz
```

### Example: export internal novel TSV and stop
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known.idx \
  --output-folder /path/to/out \
  --novel_protein /path/to/novel_proteins.fasta \
  --output_internal_novel_peptide internal_novel.tsv \
  --stop-after-saving-novel-peptide
```

### Example: reuse internal TSV + run-comet-each
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known_a.idx \
  --database /path/to/known_b.idx \
  --internal_novel_peptide /path/to/out/internal_novel.tsv \
  --run-comet-each \
  --thread 32 \
  --output-folder /path/to/out \
  --name merged_run \
  /path/to/input_a.mzMLb /path/to/input_b.mzMLb
```

Important constraints:
- `--input_files` cannot be combined with positional spectrum inputs and cannot be repeated.
- `--novel_*`, `--scan*`, and internal-novel export options cannot be used with `-i` or `-j`.
- Novel mode requires known DB from `--database` or `database_name` in params.
- At least one spectrum input is required for novel/scan mode, except stop-after mode.
- `--output_internal_novel_peptide` requires `--novel_protein` and/or `--novel_peptide`.
- `--internal_novel_peptide` is mutually exclusive with `--novel_protein` / `--novel_peptide`.
- `--stop-after-saving-novel-peptide` requires `--output_internal_novel_peptide`.
- `--name` is valid for one input, or for merged multi-input novel mode.
- `--run-comet-each` is effective only when:
  - novel mode is active,
  - multiple spectrum inputs are provided,
  - known databases are all `.idx`,
  - `.idx` type is peptide index (`-j`),
  - total threads from `--thread` or params `num_threads` is positive.

## CometPlus Build Matrix

### Default build
```bash
make -C ProtCosmo/CometPlus
```

### Static build (no mzMLb)
```bash
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static WITH_MZMLB=0
```

### Static build with mzMLb
```bash
make -C MSToolkit clean
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/abs/path/to/hdf5-static-cxx
```

When switching between `WITH_MZMLB=0` and `WITH_MZMLB=1`, run `make -C MSToolkit clean` first to avoid stale objects.

For reproducible `mzMLb`/HDF5 details and troubleshooting, see [ProtCosmo/CometPlus/build.md](ProtCosmo/CometPlus/build.md).

## Core Comet (Upstream-Compatible) Notes

<img src="https://uwpr.github.io/Comet/images/cometlogo_1_small.png" align="right">

Comet is an open source tandem mass spectrometry (MS/MS) sequence database search tool written primarily in C/C++.

Project site (release notes and parameter docs): [https://uwpr.github.io/Comet/](https://uwpr.github.io/Comet/)

### Build on Linux/macOS
```bash
make
```

### Build on Windows (Visual Studio 2022, v143)
1. Install [MSFileReader from Thermo Fisher Scientific](https://uwpr.github.io/Comet/notes/20220228_rawfile.html) for RAW support.
2. Open `Comet.sln`.
3. Select `Release` + `x64`.
4. Build the `Comet` project.

## Additional Documentation

- CometPlus overview and usage: [ProtCosmo/CometPlus/README.md](ProtCosmo/CometPlus/README.md)
- CometPlus reproducible build notes: [ProtCosmo/CometPlus/build.md](ProtCosmo/CometPlus/build.md)
- Novel protein/peptide design: [ProtCosmo/CometPlus/design.novel_protein_peptide.md](ProtCosmo/CometPlus/design.novel_protein_peptide.md)
- Novel protein/peptide implementation plan: [ProtCosmo/CometPlus/plan.novel_protein_peptide.md](ProtCosmo/CometPlus/plan.novel_protein_peptide.md)

## Integrated Libraries

Comet integrates:
- Mike Hoopmann's [MSToolkit](https://github.com/mhoopmann/mstoolkit) for spectrum I/O.
- Matthew Belmonte's C implementation of the [Twiddle algorithm](https://www.netlib.org/toms-2014-06-10/382) for modification permutation generation.
- C++ port of Gygi Lab's [AScorePro](https://github.com/gygilab/MPToolkit/) for PTM localization.
