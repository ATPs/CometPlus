# Comet indexing notes

## Flags and dispatch
- `-i` sets `create_fragment_index=1` and `create_peptide_index=0` (Comet.cpp).
- `-j` sets `create_peptide_index=1` and `create_fragment_index=0` (Comet.cpp).
- `CometSearchManager::DoSearch` routes `-j` to `CometPeptideIndex::WritePeptideIndex` and `-i` to `CometFragmentIndex::WriteFIPlainPeptideIndex` (CometSearch/CometSearchManager.cpp).
- If the database name ends with `.idx`, the first line decides index type: `Comet peptide index` or `Comet fragment ion index` (CometSearch/CometSearchManager.cpp).
- If `.idx` is requested but missing and the FASTA exists, Comet flips on fragment index creation (CometSearch/CometSearchManager.cpp).

## Common types and rounding
- File offsets are stored as `comet_fileoffset_t` (POSIX `off64_t`, Windows `__int64`); the on-disk byte width is `clSizeCometFileOffset = sizeof(comet_fileoffset_t)` (CometSearch/Common.h, CometSearch/CometFragmentIndex.cpp).
- Counts use a mix of `size_t`, `unsigned long`, and `uint64_t` depending on file type; this makes the binary layout platform-size and endianness dependent (no explicit byte order conversion).
- Fragment binning uses `BIN(dMass) = (int)(dMass * dInverseBinWidth + dOneMinusBinOffset)` where `dInverseBinWidth = 1 / fragment_bin_tol` and `dOneMinusBinOffset = 1 - fragment_bin_offset`. The cast to `int` truncates toward zero, so for positive masses this is effectively a floor after offset (CometSearch/Common.h, CometSearch/CometSearchManager.cpp).
- Precursor binning uses `BINPREC(dMass) = (int)(dMass / dMS1BinSize)`, also truncating toward zero (CometSearch/Common.h).

## Fragment ion index mode (-i)

### Generation flow
1. `WriteFIPlainPeptideIndex` selects the output `.idx` name. If the input already ends in `.idx`, it temporarily strips the extension to read the FASTA (CometSearch/CometFragmentIndex.cpp).
2. `CometSearch::RunSearch` runs in create-index mode and walks the FASTA. For fragment indexing it collects unmodified peptides by length and enzyme rules, skipping the lower mass bound because later mods can increase mass (CometSearch/CometSearch.cpp).
3. Protein names are captured while reading the FASTA when create-index mode is on (CometSearch/CometSearch.cpp).
4. Peptides are de-duplicated, then `g_pvProteinsList` is built to map each unique peptide to all protein offsets; peptides are then sorted by mass (CometSearch/CometFragmentIndex.cpp).
5. Mod permutation tables are generated for fragment indexing using the first `FRAGINDEX_VMODS` variable mods and the combination limits (CometSearch/CometFragmentIndex.cpp).
6. The file is written and the run exits if there are no input spectra (CometSearch/CometSearchManager.cpp).

### File content (on disk, with types)
Header lines (ASCII text):
- `Comet fragment ion index plain peptides.  Comet version ...`
- `InputDB`, `MassRange`, `LengthRange`, `MassType`
- `Enzyme` and `Enzyme2`
- `NumPeptides`
- `StaticMod` (A-Z plus N/C term peptide and protein)
- `VariableMod` (first `FRAGINDEX_VMODS` entries only)
- `ProteinModList` and `RequireVariableMod`

Binary sections (order):
- Protein names: each stored as fixed-width `WIDTH_REFERENCE` bytes (raw chars).
- Peptides block:
  - `size_t tNumPeptides`
  - For each peptide: `int iLen`, `char[iLen] sequence`, `char cPrevAA`, `char cNextAA`, `double dPepMass`, `unsigned short siVarModProteinFilter`, `comet_fileoffset_t lIndexProteinFilePosition`.
- Protein list (`g_pvProteinsList`): `comet_fileoffset_t tSize`, then for each row: `size_t tNumProteinOffsets`, then `comet_fileoffset_t` offsets (each offset points into the protein-name region).
- Mod permutations:
  - `unsigned long ulSizeModSeqs`, `unsigned long ulSizevRawPeptides`, `unsigned long ulModNumSize`
  - `int[ulSizeModSeqs] MOD_SEQ_MOD_NUM_START`, `int[ulSizeModSeqs] MOD_SEQ_MOD_NUM_CNT`, `int[ulSizevRawPeptides] PEPTIDE_MOD_SEQ_IDXS`
  - For each modifiable sequence: `int len`, `char[len]`
  - For each mod-number entry: `int modStringLen`, `char[modStringLen]` (mod IDs as stored by `CometModificationsPermuter`).
- Footer: three `comet_fileoffset_t` values pointing to the peptide block, protein list block, and permutation block (CometSearch/CometFragmentIndex.cpp).

### Runtime fragment index (in memory, not on disk)
- `ReadPlainPeptideIndex` loads the header, peptides, protein list, and permutation tables into memory (CometSearch/CometFragmentIndex.cpp).
- `CreateFragmentIndex` expands each plain peptide into modified variants (`g_vFragmentPeptides`) and fills `g_iFragmentIndex` bins (`unsigned int**`) with indices into `g_vFragmentPeptides` (CometSearch/CometFragmentIndex.cpp).
- Each fragment mass bin stores `unsigned int` peptide indices; counts are tracked in `g_iCountFragmentIndex` (`unsigned int*`).
- If `fragindex_skipreadprecursors=0` and input spectra exist, precursor bins are read first and used to skip peptide masses not seen in the data (`g_bIndexPrecursors[BIN(dCalcPepMass)]`) (CometSearch/CometSearchManager.cpp, CometSearch/CometPreprocess.cpp).

### Rounding/bins used in fragment indexing
- Fragment masses are accepted only if `dBion` or `dYion` are strictly between `fragindex_min_fragmentmass` and `fragindex_max_fragmentmass`; bin index is `BIN(dIonMass)` (CometSearch/CometFragmentIndex.cpp).
- Precursor filtering bins a tolerance window around each MS/MS precursor using `BIN(dMassLow)`/`BIN(dMassHigh)` where `dMassLow/High` already include tolerance and isotope offsets; bins are inclusive of the start/end indices (CometSearch/CometPreprocess.cpp).

## Peptide index mode (-j)

### Generation flow
1. `WritePeptideIndex` runs `CometSearch::RunSearch` with `bCreatePeptideIndex` to enumerate peptides (including variable-modified forms) within the configured mass range (CometSearch/CometPeptideIndex.cpp, CometSearch/CometSearch.cpp).
2. Peptides are de-duplicated by sequence + mod state; `pvProteinsListLocal` maps each unique peptide sequence to all proteins (CometSearch/CometPeptideIndex.cpp).
3. Peptides are sorted by mass and written with a 0.1 Da mass index for fast lookup (CometSearch/CometPeptideIndex.cpp).

### File content (on disk, with types)
Header lines (ASCII text):
- `Comet peptide index database.  Comet version ...`
- `InputDB`, `MassRange`, `LengthRange`, `MassType`
- `Enzyme` and `Enzyme2`
- `NumPeptides`
- `StaticMod` (A-Z plus N/C term peptide and protein)
- `VariableMod` (all `VMODS` entries)

Binary sections (order):
- Protein names: fixed-width `WIDTH_REFERENCE` bytes (raw chars).
- Protein list: `comet_fileoffset_t tSize`, then rows of `size_t tNumProteinOffsets` followed by `comet_fileoffset_t` offsets.
- Peptides block: `int iLen`, `char[iLen] sequence`, `char cPrevAA`, `char cNextAA`, `unsigned char cNumMods`, then `cNumMods` pairs of (`unsigned char position`, `char modId`), then `double dPepMass`, `comet_fileoffset_t lIndexProteinFilePosition` (CometSearch/CometPeptideIndex.cpp).
  - `position` ranges over `0..len+1` where `len` is N-term and `len+1` is C-term; `modId` is 1-based into `varModList`.
- Footer: `int minMass`, `int maxMass`, `uint64_t numPeptides`, `comet_fileoffset_t massIndex[]` at 0.1 Da resolution, then `comet_fileoffset_t lEndOfPeptides`, and `comet_fileoffset_t clProteinsFilePos` (CometSearch/CometPeptideIndex.cpp).

### Rounding/mass index behavior
- The mass index uses `iMass10 = (int)(dPepMass * 10.0)` (truncation toward zero) when writing. `lIndex[iMass10]` stores the file offset of the first peptide at or above that 0.1 Da bin (CometSearch/CometPeptideIndex.cpp).
- `minMass` is written as `(int)dPeptideMassLow` and `maxMass` as `(int)dPeptideMassHigh` (truncation).
- During search, the index lookup window uses `iStart10 = (int)(dMinMass * 10.0 - 0.5)` and `iEnd10 = (int)(dMaxMass * 10.0 + 0.5)` before truncation, which approximates rounding to the nearest 0.1 Da (CometSearch/CometSearch.cpp).

## Key differences
- Fragment index `.idx` stores plain peptides plus mod-permutation tables; the fragment index bins are built in memory each run and depend on `fragindex_*` settings and optional precursor filtering.
- Peptide index `.idx` stores fully enumerated peptides (including variable mod states) plus a 0.1 Da mass index so searches can jump directly into the peptide list.
- Fragment index mode limits variable mods to `FRAGINDEX_VMODS` and `FRAGINDEX_MAX_MODS_PER_MOD`; peptide index mode uses full `VMODS` and writes explicit mod sites per peptide.
- Both `.idx` formats use native sizes and endianness for binary sections (`size_t`, `unsigned long`, `comet_fileoffset_t`), so files are not guaranteed portable across platforms/architectures.
