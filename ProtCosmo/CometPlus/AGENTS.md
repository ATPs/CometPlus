# CometPlus Agent Guide

## Mission
Maintain `cometplus` as an additive layer over upstream-style Comet behavior:

- Preserve legacy single-db behavior.
- Add modern CLI ergonomics.
- Support multi-db combined semantics.
- Keep core `CometSearch` edits minimal and surgical.

## Non-Negotiables
1. Do not break legacy short options or existing single-db search behavior.
2. Do not change root `Makefile` unless explicitly requested.
3. Multi-idx mode must reject mixed index types and mixed FASTA/`.idx` lists.
4. Multi-idx mode must keep deterministic protein-set mapping and duplicate-merging behavior.
5. Multi-idx + mzIdentML must fail fast with a clear message.

## Ownership Boundaries
### Preferred location for new work
- `ProtCosmo/CometPlus/CometPlus/CometPlus.cpp`
- `ProtCosmo/CometPlus/CometPlus/`
- `ProtCosmo/CometPlus/Makefile`
- `ProtCosmo/CometPlus/*.md`
each new c++ code file best not more than 1000 lines. include instructions about the major aim of each file in folder `ProtCosmo/CometPlus/CometPlus/`

### Allowed minimal hooks in core engine
- `CometSearch/CometPlusMultiDB.{h,cpp}`
- `CometSearch/CometSearchManager.cpp`
- `CometSearch/CometSearch.cpp`
- `CometSearch/CometFragmentIndex.cpp`
- `CometSearch/CometMassSpecUtils.cpp`
- `CometSearch/Makefile`


Avoid broad refactors outside these files for CometPlus features.

## Multi-DB Invariants
1. `database_name_list` is newline-delimited and preserves CLI order.
2. `cometplus_multi_db_mode` values currently used:
   - `multi_fasta`
   - `multi_idx`
3. In multi-idx mode:
   - all dbs are `.idx`
   - all `.idx` headers are compatibility-checked
   - all protein-set references are remapped to global set IDs
   - duplicate peptides across files union protein sets

## CLI Contract
`cometplus` must continue supporting:

- Legacy short options: `-D -P -N -F -L -p -q -i -j` (and `-B` pass-through).
- Long options:
  - `--database` (repeatable)
  - `--thread`
  - `--params`
  - `--name`
  - `--first-scan`
  - `--last-scan`

`--thread` should only override `num_threads`; core thread semantics remain in `CometSearchManager`.

## Novel Protein/Peptide + Scan Milestone Rules
1. Keep new option orchestration in `ProtCosmo/CometPlus/CometPlus/CometPlus.cpp` first.
2. Known `.idx` inputs must be reused; do not regenerate or mutate known `.idx` files.
3. `--novel_peptide` must accept both FASTA and tokenized text input.
4. Explicit scan filters must support comma/whitespace/newline delimiters and intersect with `-F/-L`.
5. Preserve existing behavior for runs that do not use new options.
6. Allow only minimal hooks in:
- `CometSearch/CometSearchManager.cpp`
- `CometSearch/CometPreprocess.cpp`
- `CometSearch/CometData.h`
- `CometSearch/CometPreprocess.h`
7. Require validation to include:
- legacy run parity,
- known `.idx` non-regeneration path,
- filtered-scan-only output behavior.

## Output Expectations
For multi-idx mode in this milestone:

- Supported: `txt`, `sqt`, `pep.xml`, `pin`.
- Unsupported: `mzid` (must fail pre-search).

## Validation Workflow
At minimum before handing off:

1. `make -C ProtCosmo/CometPlus` passes.
2. `./ProtCosmo/CometPlus/cometplus --help` includes modern options.
3. Negative CLI checks:
   - mixed FASTA + `.idx` rejected
   - mixed `.idx` types rejected
4. If datasets available, run parity checks against baseline merged-db runs.
5. Verify known `.idx` mode does not trigger index creation.
6. Verify `--novel_protein/--novel_peptide` + `--scan/--scan_numbers` combinations.
7. Verify `equal_I_and_L` subtraction behavior in both states.

## Practical Notes
- Multi-idx protein-name output does not rely on a single open `.idx` handle.
- Keep deterministic ordering for tie-breaks by earliest global protein ordinal.
- If adding new multi-idx output formats, route through the global protein registry helpers first.

## mzMLb Milestone Rules
1. Native `.mzMLb` support requires MSToolkit compiled with `MZP_HDF`.
2. Full-static `cometplus` builds for `.mzMLb` require an external static HDF5 install.
3. Prefer CometPlus-local changes first:
   - `ProtCosmo/CometPlus/CometPlus/CometPlus.cpp`
   - `ProtCosmo/CometPlus/Makefile`
   - `ProtCosmo/CometPlus/*.md`
4. Minimal allowed hooks for `.mzMLb` compatibility:
   - `CometSearch/CometSearchManager.cpp` for input extension gating.
   - `CometSearch/CometWriteMzIdentML.cpp` for file format metadata mapping.
   - `MSToolkit/Makefile` for HDF5 compile/archive control.
5. Keep root `Makefile` unchanged unless explicitly requested.
6. Keep CLI semantics compatible with upstream Comet behavior.
7. `.mgf.gz` support is implemented via temporary inflate in the output directory:
   - keep original `.mgf.gz` input path for metadata/reporting,
   - use temporary `.mgf` only for spectrum reading,
   - remove temporary files at normal completion/handled failure paths.
8. Defer remaining gzip-expansion work (`.ms2.gz`, broader strategy) unless explicitly reopened.

# build cometplus
check this file [build.md](ProtCosmo/CometPlus/build.md) 