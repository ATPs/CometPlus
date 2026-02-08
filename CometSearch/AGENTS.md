# CometSearch Guidelines

## Scope
`CometSearch/` builds `libcometsearch.a`, which is the main engine library used by top-level `Comet.cpp`.
Core orchestration lives in `CometSearchManager.cpp`, with search logic in `CometSearch.cpp` and preprocessing/post-analysis in `CometPreprocess.cpp` and `CometPostAnalysis.cpp`.

## Build and Clean
- From repo root: `make` (preferred end-to-end build).
- From this folder: `make` to build `libcometsearch.a`.
- Clean this module only: `make clean`.

## Dependency Boundaries
- Public includes expected by this module come from `../MSToolkit/include` and `../AScorePro/include`.
- If you add a new cross-module include, verify `CometSearch/Makefile` include paths and dependency rules still cover it.
- Keep Comet interfaces stable when changing `CometInterfaces.h` or `CometData*.h` because they are used by the top-level executable.

## Coding Expectations
- Follow `../CometCodingStyleGuidelines.txt` (Allman braces, spaces-only indentation, `//` comments).
- Keep headers self-contained and include only what is required.
- Prefer existing logging/status paths (`CometStatus`, `logout`, `logerr`) rather than adding ad hoc output paths.

## Validation
- Rebuild from repo root with `make`.
- Run a representative search workflow to verify preprocessing, search, and output writers still behave correctly.
