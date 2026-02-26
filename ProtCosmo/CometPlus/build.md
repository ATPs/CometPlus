# CometPlus Reproducible Build Guide (Static, WITH_MZMLB=0/1)

This guide documents how to build `cometplus` reproducibly in this environment, including how `/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp` was created and how to prepare a working HDF5 prefix for mzMLb support.

## Scope

- Repository root: `/data/p/comet/Comet`
- Binary target: `/data/p/comet/Comet/ProtCosmo/CometPlus/cometplus`
- Static build target: `ldd` prints `not a dynamic executable`
- Profiles covered:
  - `WITH_MZMLB=0` (no mzMLb support, no HDF5 linkage)
  - `WITH_MZMLB=1` (mzMLb enabled, static HDF5 linkage)

## Toolchain And Build Environment

Use these commands to record your local build context before building:

```bash
date
uname -a
make --version | head -n 1
ld --version | head -n 1

command -v x86_64-conda-linux-gnu-gcc
command -v x86_64-conda-linux-gnu-g++
x86_64-conda-linux-gnu-gcc --version | head -n 1
x86_64-conda-linux-gnu-g++ --version | head -n 1

command -v /usr/bin/gcc
command -v /usr/bin/g++
/usr/bin/gcc --version | head -n 1
/usr/bin/g++ --version | head -n 1
```

If `x86_64-conda-linux-gnu-gcc/g++` are not found in `PATH`, use:

```bash
export PATH=/data/p/anaconda3/bin:$PATH
```

Observed in this environment:

- `x86_64-conda-linux-gnu-gcc (conda-forge gcc 13.2.0-13) 13.2.0`
- `x86_64-conda-linux-gnu-g++ (conda-forge gcc 13.2.0-13) 13.2.0`
- `/usr/bin/gcc (GCC) 11.5.0 20240719`
- `/usr/bin/g++ (GCC) 11.5.0 20240719`
- `GNU Make 4.3`

## Required Libraries And Link Notes

CometPlus `WITH_MZMLB=1` static link path uses:

- Static HDF5 C archive: `libhdf5.a`
- Static HDF5 C++ archive: `libhdf5_cpp.a`
- Extra link libs (from `ProtCosmo/CometPlus/Makefile`): `-ldl -lz -lm -lpthread`

Default variables in `ProtCosmo/CometPlus/Makefile`:

```make
WITH_MZMLB ?= 0
HDF5_DIR ?=
HDF5_LIB_C ?= $(HDF5_DIR)/lib/libhdf5.a
HDF5_LIB_CPP ?= $(HDF5_DIR)/lib/libhdf5_cpp.a
HDF5_EXTRA_LIBS ?= -ldl -lz -lm -lpthread
```

If your HDF5 build has additional dependencies (for example `-lsz`), add them through `HDF5_EXTRA_LIBS`.

## HDF5 Source And Prefixes Under `/data/p/hdf5`

List available HDF5 directories:

```bash
ls -d /data/p/hdf5/*
```

Relevant paths used here:

- Source: `/data/p/hdf5/hdf5-1.14.4-3`
- Historical prefix: `/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp`
- Recommended prefix for mzMLb runtime in this environment: `/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp`

If source is missing, prepare it from tarball:

```bash
cd /data/p/hdf5
test -d hdf5-1.14.4-3 || tar -xzf hdf5-1.14.4-3.tar.gz
```

Check zlib headers that each build path typically uses:

```bash
rg -n "#define ZLIB_VERSION" /data/p/anaconda3/include/zlib.h /usr/include/zlib.h
```

Observed in this environment:

- `/data/p/anaconda3/include/zlib.h`: `1.3.1`
- `/usr/include/zlib.h`: `1.2.11`

Verify required files in a prefix:

```bash
HDF5_PREFIX=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp
ls -la \
  "$HDF5_PREFIX/include/H5Cpp.h" \
  "$HDF5_PREFIX/lib/libhdf5.a" \
  "$HDF5_PREFIX/lib/libhdf5_cpp.a"
```

## How `/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp` Was Created

This is the exact historical configure recipe from its `config.log`.

```bash
cd /data/p/hdf5
mkdir -p hdf5-1.14.4-3_build_static_cpp_work
cd hdf5-1.14.4-3_build_static_cpp_work

CC=x86_64-conda-linux-gnu-gcc \
CXX=x86_64-conda-linux-gnu-g++ \
/data/p/hdf5/hdf5-1.14.4-3/configure \
  --prefix=/data/p/hdf5/hdf5-1.14.4-3_build_static_cpp \
  --enable-cxx --enable-hl \
  --enable-static --disable-shared \
  --disable-tools --disable-tests \
  --with-zlib=/data/p/anaconda3 \
  --with-szlib=no

make -j"$(nproc)"
make install
```

Verification:

```bash
grep -n "C Compiler\\|C++ Compiler\\|I/O filters (external)" \
  /data/p/hdf5/hdf5-1.14.4-3_build_static_cpp/lib/libhdf5.settings
```

Observed issue in this environment:

- `I/O filters (external)` was empty in this prefix.
- mzMLb runtime hit: `required filter 'deflate' is not registered`.
- This path was built with conda cross toolchain (`x86_64-conda-linux-gnu-gcc/g++`) and `--with-zlib=/data/p/anaconda3`.

Conclusion: keep this section for reproducibility/history, but do not use this prefix for current mzMLb runtime validation.

## Recommended: Rebuild Working HDF5 Prefix For mzMLb

This is the known working approach in this environment.

```bash
cd /data/p/hdf5
mkdir -p hdf5-1.14.4-3_build_native_cpp_work
cd hdf5-1.14.4-3_build_native_cpp_work

CC=/usr/bin/gcc \
CXX=/usr/bin/g++ \
/data/p/hdf5/hdf5-1.14.4-3/configure \
  --prefix=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp \
  --enable-cxx --enable-hl \
  --enable-static --disable-shared \
  --disable-tools --disable-tests \
  --with-szlib=no

make -j"$(nproc)"
make install
```

Verification:

```bash
grep -n "C Compiler\\|C++ Compiler\\|I/O filters (external)" \
  /data/p/hdf5/hdf5-1.14.4-3_build_native_cpp/lib/libhdf5.settings
```

Expected for working mzMLb support:

- C compiler: `/usr/bin/gcc`
- C++ compiler: `/usr/bin/g++`
- `I/O filters (external): deflate(zlib)`

## Build CometPlus: Static Without mzMLb (`WITH_MZMLB=0`)

From repo root:

```bash
cd /data/p/comet/Comet
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static WITH_MZMLB=0
```

Verify:

```bash
file ProtCosmo/CometPlus/cometplus
ldd ProtCosmo/CometPlus/cometplus
nm -A ProtCosmo/CometPlus/cometplus 2>/dev/null | rg -i "\\bH5[A-Za-z0-9_]*\\b|hdf5"
```

Expected:

- static binary (`statically linked`, `not a dynamic executable`)
- no `H5`/`hdf5` symbols

## Build CometPlus: Static With mzMLb (`WITH_MZMLB=1`)

From repo root, using the working HDF5 prefix:

```bash
cd /data/p/comet/Comet
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp
```

Equivalent explicit archive invocation:

```bash
cd /data/p/comet/Comet
make -C ProtCosmo/CometPlus clean
make -C ProtCosmo/CometPlus static \
  WITH_MZMLB=1 \
  HDF5_DIR=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp \
  HDF5_LIB_C=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp/lib/libhdf5.a \
  HDF5_LIB_CPP=/data/p/hdf5/hdf5-1.14.4-3_build_native_cpp/lib/libhdf5_cpp.a \
  HDF5_EXTRA_LIBS="-ldl -lz -lm -lpthread"
```

Verify:

```bash
ls -lh ProtCosmo/CometPlus/cometplus
file ProtCosmo/CometPlus/cometplus
ldd ProtCosmo/CometPlus/cometplus || true
nm -A ProtCosmo/CometPlus/cometplus 2>/dev/null | rg -i "\\bH5[A-Za-z0-9_]*\\b|hdf5" | head
strings ProtCosmo/CometPlus/cometplus | rg -n "mzMLb|Supported input formats include"
```

Expected:

- static binary (`statically linked`, `not a dynamic executable`)
- `H5...` symbols present
- `mzMLb` appears in strings/help text
- typical size increase vs dynamic build in this environment: dynamic ~`8.7M`, static ~`16M`

## Smoke Validation Commands

### Direct `.mzML.gz` validation (Numpress+zlib file)

```bash
python3 ProtCosmo/CometPlus/test/cometplus_mzmlgz_direct_validation.py \
  --run-tag 20260209_static_full_v1
```

### mzMLb smoke run (scan 1000-1200)

```bash
BIN=/data/p/comet/Comet/ProtCosmo/CometPlus/cometplus
PARAMS=/XCLabServer002_fastIO/cometplus_tst/run_all_full/comet.params.runtime
PEP=/XCLabServer002_fastIO/cometplus_tst/run_all_full/idx/combined.pep.idx
FRAG=/XCLabServer002_fastIO/cometplus_tst/run_all_full/idx/combined.frag.idx
IN=/XCLabServer002_fastIO/cometplus_test_mzmlb/data/1554877.mzMLb

"$BIN" --params "$PARAMS" --thread 12 --name static_mzmlb_smoke_pepidx \
  --database "$PEP" --first-scan 1000 --last-scan 1200 "$IN"

"$BIN" --params "$PARAMS" --thread 12 --name static_mzmlb_smoke_fragidx \
  --database "$FRAG" --first-scan 1000 --last-scan 1200 "$IN"
```

## Troubleshooting

- Error: `WITH_MZMLB=1 requires static HDF5 C++ library`
  - Cause: `HDF5_DIR` missing `lib/libhdf5_cpp.a`.
  - Fix: rebuild HDF5 with `--enable-cxx --enable-static --disable-shared`.

- Runtime: `required filter 'deflate' is not registered`
  - Cause: HDF5 prefix lacks usable deflate(zlib) filter.
  - Fix: use/rebuild prefix where `libhdf5.settings` shows `I/O filters (external): deflate(zlib)`.

- Link errors when enabling mzMLb
  - Cause: missing transitive system libs from static HDF5 stack.
  - Fix: add libs through `HDF5_EXTRA_LIBS` (baseline `-ldl -lz -lm -lpthread`; add more only if linker asks).

- Binary unexpectedly small (for example ~4 MB)
  - Cause: dynamic build or non-static path used.
  - Fix: run `make -C ProtCosmo/CometPlus static ...` and recheck with `file` + `ldd`.
