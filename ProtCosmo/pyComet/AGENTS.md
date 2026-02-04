# pyComet Guidelines

## Purpose
Python utilities that mirror Comet behavior for prototyping and data prep. Keep outputs consistent with the C++ engine where feasible.

## Module Map
- `pyLoadParameters.py`: parses `comet.params` and mirrors `Comet.cpp::LoadParameters()` behavior.
- `pyPeptideIndex.py`: digests FASTA and emits peptide index tables; depends on `pyLoadParameters`.
- Add new Python modules in this folder and keep them importable without external dependencies unless explicitly approved.

## Run
- `python3 pyPeptideIndex.py --help`

## Coding Style
- Python 3, standard library only by default.
- 4-space indentation, CRLF line endings, no trailing whitespace.
- Prefer type hints and dataclasses where they clarify data flow.
- If you change parsing or digestion logic, cross-check the equivalent code in `Comet.cpp` or `CometSearch/` and call out intentional deviations.

## Testing
No automated tests. Validate changes by running `pyPeptideIndex.py` on a small FASTA and a known `comet.params`, then sanity-check the generated TSVs. If you add tests, keep them adjacent to the module and document how to run them.
