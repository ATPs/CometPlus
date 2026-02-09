# MSToolkit Guidelines

## Scope
`MSToolkit/` provides mass-spectrometry reader and parser libraries consumed by Comet.
Main outputs are `libmstoolkit.a` and `libmstoolkitextern.a`, built from:
- `src/MSToolkit/` (core toolkit objects)
- `src/mzParser/` (format parsing)
- `extern/` (bundled expat and zlib sources)

## Build and Clean
- Build this module: `make all` (or just `make`).
- Build internal components only: `make mst`.
- Build bundled externals only: `make extern`.
- Clean intermediates: `make clean`.
- Full cleanup/extract reset: `make realclean`.

## Dependency and Vendor Notes
- Comet primarily consumes headers under `MSToolkit/include/` (`MSReader.h`, `Spectrum.h`, `MSObject.h`, `mzParser.h`).
- `extern/expat-2.2.9` and `extern/zlib-1.2.11` are vendored third-party code. Avoid local style rewrites there unless you are intentionally updating vendor code.
- Some files under `include/` are generated/copied during `make` (`expat*.h`, `zlib.h`, `zconf.h`); prefer changing source/vendor inputs instead of editing generated outputs directly.

## Coding Expectations
- Keep platform behavior consistent: Linux builds exclude `RAWReader.cpp` while Windows uses Visual Studio paths.
- Preserve existing compile flags and macros unless you are deliberately changing ABI or parser behavior.

## Validation
- After edits, run `make` here and then rebuild the repo root (`../make`) to verify integration with `CometSearch` and `Comet.cpp`.
