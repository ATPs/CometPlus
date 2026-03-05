# ProtCosmo

ProtCosmo is a toolkit folder for proteogenomics research, focused on finding novel peptides from MS/MS data.

check [ProtCosmo](https://github.com/ATPs/ProtCosmo) page.

Here is tools that will be used by ProtCosmo

## Main Components

- [CometPlus/](CometPlus/): a CLI layer on top of CometSearch for search workflows, including repeated `--database`, novel peptide/protein inputs, and scan-subset options.
- [pyComet/](pyComet/): Python utilities for Comet-related file processing, including PIN/percolator conversion, UniMod matching, and peptide index table generation.

## Typical Workflow

1. Build/select known databases and novel candidate inputs.
2. Run searches with `CometPlus`.
3. Post-process PSMs with Percolator and related tools.
4. Use `pyComet` utilities to convert/process outputs and support downstream analysis.

## Notes

- ProtCosmo currently includes CometPlus and pyComet.
- Additional tools may be added under `ProtCosmo/` as workflows expand.
