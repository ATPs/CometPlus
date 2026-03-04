# Comet + CometPlus

Forked from [UWPR/Comet](https://github.com/UWPR/Comet).

This repository contains:
- **Comet**: the core MS/MS search engine (`comet.exe`).
- **CometPlus**: an optional CLI wrapper in [ProtCosmo/CometPlus](ProtCosmo/CometPlus) that keeps Comet-compatible behavior and adds workflow-oriented capabilities:
  - repeated `--database` (multi-database search),
  - multi-input search in one run,
  - novel peptide/protein search-space augmentation,
  - explicit scan subset filtering,
  - gzip input support (`.mzML.gz`, `.mzXML.gz`, `.mgf.gz`),
  - `mzMLb` support when built with HDF5-enabled MSToolkit.

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
Produces `ProtCosmo/CometPlus/cometplus`.

### Show CometPlus help
```bash
./ProtCosmo/CometPlus/cometplus --help
```

## CometPlus Input/Database Support

### Supported spectrum formats
- `mzXML`
- `mzML`
- `mzMLb` (requires `WITH_MZMLB=1` build with static HDF5 C/C++ libs)
- `mzXML.gz`
- `mzML.gz`
- `mgf`
- `mgf.gz`
- `ms2` family: `ms2`, `cms2`, `bms2`
- Thermo RAW (platform/toolchain dependent)

Notes:
- gzip support in this milestone is `mzXML.gz`, `mzML.gz`, `mgf.gz`.
- `ms2.gz` is currently **not supported**.

### Database inputs
- `--database` is repeatable.
- You can pass multiple FASTA files, or multiple `.idx` files.
- Mixed types in one run are rejected (FASTA + `.idx` together is not allowed).
- For multi-`.idx`, all `.idx` files must be the same index type.
- If `--database` is not passed, CometPlus can fall back to `database_name` from params.

## Common CometPlus Commands

### Single input + single database
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/proteins.fasta \
  --name run1 \
  /path/to/input.mzML
```

### Multiple input files in one run
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/proteins.fasta \
  /path/to/file1.mzML /path/to/file2.mgf /path/to/file3.mzML.gz
```

### Multiple FASTA databases
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/target.fasta \
  --database /path/to/decoy.fasta \
  --name run_multi_fasta \
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

## Novel Protein/Peptide + Scan Subset Workflow

CometPlus adds:
- `--novel_protein <file>`: FASTA, digested with active Comet settings.
- `--novel_peptide <file>`: FASTA or tokenized text (auto-detected by `>` headers).
- `--scan <file>`: scan IDs from file (comma/whitespace delimited).
- `--scan_numbers <list>`: inline scan IDs (comma/whitespace delimited).
- `--keep-tmp`: keep temporary artifacts on exit (debugging helper).

High-level behavior when novel/scan options are used:
1. Resolve known DB from `--database` (or params fallback).
2. Build known peptide universe.
3. Build novel candidates from `--novel_protein` and/or `--novel_peptide`.
4. Subtract known peptides from novel candidates (respecting `equal_I_and_L`).
5. Add retained novel set to scoring DB.
6. Union scan IDs from `--scan` and `--scan_numbers`, then apply scan range intersection if set.
7. Prefilter spectra and run normal Comet search on filtered inputs.

Example (combined novel + scan filtering):
```bash
./ProtCosmo/CometPlus/cometplus \
  --params /path/to/comet.params \
  --database /path/to/known_target.idx \
  --database /path/to/known_decoy.idx \
  --novel_protein /path/to/novel_proteins.fasta \
  --novel_peptide /path/to/novel_peptides.txt \
  --scan /path/to/scan_ids.txt \
  --scan_numbers 3001,3002,3003 \
  --first-scan 2500 \
  --last-scan 6000 \
  --name run_novel_scan_subset \
  /path/to/input.mzML.gz
```

Important constraints:
- `--novel_*` and `--scan*` cannot be used with `-i` or `-j` index-creation modes.
- Novel mode requires known DB from `--database` or `database_name` in params.
- At least one spectrum input file is required when novel/scan subset options are used.
- `-N`/`--name` is intended for single-input runs.

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
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/abs/path/to/hdf5-static-cxx
```

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
- Novel protein/peptide design: [ProtCosmo/CometPlus/design.novel_protein_peptide.md](ProtCosmo/CometPlus/design.novel_protein_peptide.md)
- Novel protein/peptide implementation plan: [ProtCosmo/CometPlus/plan.novel_protein_peptide.md](ProtCosmo/CometPlus/plan.novel_protein_peptide.md)

## Integrated Libraries

Comet integrates:
- Mike Hoopmann's [MSToolkit](https://github.com/mhoopmann/mstoolkit) for spectrum I/O.
- Matthew Belmonte's C implementation of the [Twiddle algorithm](https://www.netlib.org/toms-2014-06-10/382) for modification permutation generation.
- C++ port of Gygi Lab's [AScorePro](https://github.com/gygilab/MPToolkit/) for PTM localization.
