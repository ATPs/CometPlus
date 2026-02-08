# AScorePro Guidelines

## Scope
`AScorePro/` builds `libascorepro.a`, which provides phosphorylation site-scoring utilities used by `CometSearch`.
Implementation sources are in this directory and public interfaces are in `include/`.

## Build and Clean
- Build library: `make` (default target `all`).
- Clean outputs: `make clean`.
- Main output artifact: `libascorepro.a` plus objects under `obj/`.

## Dependency Boundaries
- `CometSearch` includes `AScoreOptions.h`, `AScoreFactory.h`, `AScoreDllInterface.h`, and related API headers from `include/`.
- Keep factory/interface behavior backward-compatible unless you also update `CometSearch` call sites in the same change.

## Coding Expectations
- Preserve current lightweight library structure: algorithms in `*.cpp`, interfaces in `include/*.h`.
- Keep floating-point/scoring behavior changes explicit and documented in commit/PR notes because they can affect search ranking outputs.

## Validation
- Run `make` in this folder.
- Rebuild from repo root with `make` and run a representative workflow that exercises post-analysis/AScore paths.
