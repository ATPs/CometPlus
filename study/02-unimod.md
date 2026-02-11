# Comet + Unimod: conclusions

## Scope of Unimod support
- Unimod is used in two places: (1) PEFF parsing when `peff_format` is Unimod, and (2) mzIdentML output where some common mass/residue pairs are mapped to Unimod IDs. Outside these, modifications are handled as mass deltas, not Unimod IDs.

## PEFF Unimod ingestion (search-time)
- Config: `peff_format` controls PEFF mode; Unimod is selected by `peff_format = 2` (and internal modes include `4 = Unimod only`). `peff_obo` is required; if `peff_format` is set and `peff_obo` is empty, Comet errors out early. (`Comet.cpp`, `CometSearchManager.cpp`, `CometSearch.cpp`)
- OBO parsing: `ReadOBO` reads `[Term]` blocks, captures `id:` plus `xref: delta_mono_mass` / `xref: delta_avge_mass` (Unimod) or `xref: DiffMono` / `DiffAvg` (PSI-MOD). Entries are kept only when mono mass is non-zero, sorted, and then accessed by binary search. (`CometSearch.cpp`)
- PEFF header parsing: when Unimod is selected, Comet looks for `\ModResUnimod=` in PEFF headers. It expects entries like `(pos[,pos]|UNIMOD:XXXX|...)`, converts positions from 1-based to 0-based, and maps the Unimod ID to mass via the OBO table. Invalid positions are skipped with warnings; missing closing parentheses are treated as errors. (`CometSearch.cpp`)
- Stored representation: each PEFF mod becomes a `PeffModStruct` with `iPosition`, `dMassDiffMono`, `dMassDiffAvg`, and the original mod ID string (e.g., `UNIMOD:1`). (`CometDataInternal.h`, `CometSearch.cpp`)

## How Unimod/PEFF mods are applied in the search
- PEFF mods are only considered when the peptide is *not* a PEFF variant (AA substitution/insertion/deletion); variants and PEFF mods are not applied together. (`CometSearch.cpp`)
- For a candidate peptide, PEFF mods are grouped by position. Mass screening uses **only one PEFF mod at a time**: `WithinMassTolerancePeff` checks whether *any single* PEFF mod makes the peptide mass fall within tolerance. No combinations across multiple PEFF positions are attempted. (`CometSearch.cpp`)
- During full permutation (`MergeVarMods`), each PEFF mod is tried individually on top of any variable mods. PEFF mods are encoded as negative values in `piVarModSites`, pointing back to `PeffModStruct` indices. (`CometSearch.cpp`)
- Fragment mass calculations add PEFF mass deltas using the stored mono mass (`dMassDiffMono`). Average mass is stored but not used for scoring. (`CometSearch.cpp`)

## Output/annotation behavior
- Text output: when PEFF is enabled, the "mod string with PEFF" prints `UNIMOD:XXXX` (from `pszMod`) at the modified residue; encoded mod list uses `P` entries with the mass delta. (`CometWriteTxt.cpp`, `CometDataInternal.h`)
- mzIdentML: Comet always registers the UNIMOD CV in the header. It assigns specific Unimod IDs **heuristically** based on mass/residue matches (0.01 Da tolerance) for a fixed list of common mods; otherwise it uses PSI-MS "unknown modification." This mapping applies to both static and variable mods defined in params, not to arbitrary Unimod IDs unless they match the hard-coded mass/residue rules. (`CometWriteMzIdentML.cpp`)

## Unimod IDs hard-coded in Comet (mzIdentML mapping)
These are the "common" Unimod IDs Comet will emit in mzIdentML when a mod mass and residue match its hard-coded rules (0.01 Da tolerance). (`CometWriteMzIdentML.cpp`)
- UNIMOD:35 Oxidation (+15.994915) on DKNPFYRMCWHGUEILQSTV
- UNIMOD:21 Phospho (+79.966331) on TSYDHCRK
- UNIMOD:1 Acetyl (+42.010565) on nKCSTYHR
- UNIMOD:3 Biotin (+226.077598) on nK
- UNIMOD:4 Carbamidomethyl (+57.021464) on nCKHDESTYHM
- UNIMOD:5 Carbamyl (+43.005814) on nKRCMSTY
- UNIMOD:6 Carboxymethyl (+58.005479) on nCKWU
- UNIMOD:7 Deamidated (+0.984016) on QNRF
- UNIMOD:27 Glu->pyro-Glu (-18.010565) on E
- UNIMOD:23 Dehydrated (-18.010565) on NQSTYDC
- UNIMOD:28 Gln->pyro-Glu (-17.026549) on Q
- UNIMOD:385 Ammonia-loss (-17.026549) on TSCN
- UNIMOD:34 Methyl (+14.01565) on ncCHKNQRILEDST
- UNIMOD:121 GG (+114.042927) on nKSTC
- UNIMOD:122 Formyl (+27.994915) on nKST
- UNIMOD:275 Nitrosyl (+28.990164) on C
- UNIMOD:737 TMT6plex (+229.162932) on nKHST
- UNIMOD:738 TMT2plex (+225.155833) on nKHST
- UNIMOD:739 TMT (+224.152478) on nKHST
- UNIMOD:2016 TMTpro 16plex (+304.207146) on nKHST

## extra
- UNIMOD:737 TMT10plex (+229.162932) on nKHST
- UNIMOD:2016 TMTpro 18plex (+304.207146) on nKHST
- UNIMOD:2017 TMTpro_zero (+295.189592) on nKHSTY
- UNIMOD:2015 shTMT (+235.176741) on nK
- UNIMOD:2050 shTMTpro (+313.231019) on nK
- UNIMOD:2122 Label:13C(6)15N(2)+TMT6plex (+237.177131) on K
- UNIMOD:2123 Label:13C(6)15N(2)+TMTpro (+312.221344) on K
- UNIMOD:984 cysTMT (+299.166748) on C
- UNIMOD:985 cysTMT6plex (+304.177202) on C
- UNIMOD:1341 iodoTMT (+324.216141) on CDEHK
- UNIMOD:1342 iodoTMT6plex (+329.226595) on CDEHK
- UNIMOD:765 Met-loss (-131.040485) on nM
- UNIMOD:766 Met-loss+Acetyl (-89.02992) on nM
- UNIMOD:36 Dimethyl (+28.0313) on nKR
- UNIMOD:37 Trimethyl (+42.04695) on nKR
- UNIMOD:425 Dioxidation (+31.989829) on MCWY
- UNIMOD:345 Trioxidation (+47.984744) on CWFY
- UNIMOD:24 Propionamide (+71.037114) on nCK
- UNIMOD:40 Sulfo (+79.956815) on STY
- UNIMOD:41 Hex (+162.052824) on NST
- UNIMOD:43 HexNAc (+203.079373) on NST
- UNIMOD:214 iTRAQ4plex (+144.102063) on nK
- UNIMOD:730 iTRAQ8plex (+304.20536) on nK
- UNIMOD:731 iTRAQ8plex:13C(6)15N(2) (+304.19904) on nK
- UNIMOD:747 Malonyl (+86.000394) on CSK



## Practical implications
- Unimod is **not** a general mod-definition mechanism in Comet; for search-time behavior, it is only used to interpret PEFF `ModResUnimod` entries via an OBO file.
- PEFF Unimod mods are treated as **single-mod alternatives**, not combinatorial multi-PEFF modifications on the same peptide.
- If you need Unimod IDs preserved in output, PEFF + Unimod OBO is the direct path; otherwise mzIdentML will only label mods that match the built-in mass/residue heuristics.
