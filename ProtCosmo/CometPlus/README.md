# CometPlus

CometPlus is a thin CLI layer on top of CometSearch that keeps upstream-compatible search behavior while adding quality-of-life options (notably repeated `--database`) used by ProtCosmo workflows.

Supported spectrum inputs:
- `mzXML`
- `mzML`
- `mzMLb` (requires HDF5-enabled MSToolkit build)
- `mgf`
- `mgf.gz` (inflated to a temporary `.mgf` in the output directory during search)
- `ms2` variants (`ms2`, `cms2`, `bms2`)
- Thermo RAW (platform/toolchain dependent)

## Build Requirements

- Linux/macOS toolchain with C++14 support
- Existing repo subprojects (`MSToolkit`, `AScorePro`, `CometSearch`)
- For fully static `mzMLb` binaries: static HDF5 build with C and C++ libraries

For static Linux builds, this project uses:
- `x86_64-conda-linux-gnu-g++`
- `x86_64-conda-linux-gnu-gcc`

## Build Matrix

### 1) Default build
```bash
make -C ProtCosmo/CometPlus
```

### 2) Static build (no mzMLb/HDF5 requirement)
```bash
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static
```

### 3) Static build with mzMLb enabled
```bash
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/abs/path/to/hdf5-static-cxx
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

## Gzip Input Notes

- `mzXML.gz` and `mzML.gz` are read through the native MSToolkit gzip path.
- `mgf.gz` is handled by CometPlus/CometSearch using a temporary `.mgf` file:
  - temporary file directory is derived from the output base path (`-N/--name` or input basename),
  - the original input path remains `.mgf.gz` for reporting/metadata,
  - temporary file is removed automatically at normal completion and handled error exits.
- `ms2.gz` is not supported in this milestone.

## Usage

### Help
```bash
./ProtCosmo/CometPlus/cometplus --help
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
