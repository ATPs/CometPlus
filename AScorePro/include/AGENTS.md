# AScorePro Public Headers Guidelines

## Scope
This folder defines the public AScore API consumed by `CometSearch`.
Key integration headers include:
- `AScoreOptions.h`
- `AScoreFactory.h`
- `AScoreDllInterface.h`
- `AScoreOutput.h`
- `AScorePeptideBuilder.h`

## API Change Rules
- Prefer non-breaking changes; many headers are transitively included by `CometDataInternal.h`.
- If changing structs/classes used across module boundaries, update corresponding implementation files and dependent `CometSearch` code together.
- Keep ownership/lifetime expectations clear in signatures and comments.

## Header Hygiene
- Keep dependencies minimal to reduce rebuild fan-out.
- Avoid introducing unnecessary STL-heavy or platform-specific includes in shared interfaces.
- Maintain existing naming/style conventions used across AScore headers.

## Validation
- Rebuild `AScorePro`, then rebuild `CometSearch` and top-level `Comet` to ensure API compatibility.
