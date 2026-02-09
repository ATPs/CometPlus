# MSToolkit Public Headers Guidelines

## Scope
This folder contains public headers used by Comet and other modules.
Directly used entry points from Comet include:
- `MSReader.h`
- `Spectrum.h`
- `MSObject.h`
- `mzParser.h`

## API Stability
- Treat these headers as a public interface surface.
- Prefer additive changes over breaking signature/type changes.
- If a breaking change is required, update all dependent call sites in `CometSearch/` in the same change.

## Include Hygiene
- Keep includes minimal and order them deterministically.
- Avoid introducing platform-specific dependencies into shared headers unless guarded and documented.
- Preserve compatibility macros already used by this module.

## Validation
- Rebuild `MSToolkit` and `CometSearch`, then run top-level `make` to verify header-level compatibility.
