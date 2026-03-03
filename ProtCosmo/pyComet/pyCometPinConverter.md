## pyCometPinConverter convert pin file format

## aim

develop a tool to convert comet pin files to unimod format, and output format can be .pin.parquet.

this tool can also convert .pin.parquet to .pin format
```
comet pin files first fiew lines looks like
SpecId	Label	ScanNr	ExpMass	CalcMass	lnrSp	deltLCn	deltCn	lnExpect	Xcorr	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzN	enzC	enzInt	lnNumSP	dM	absdM	Peptide	Proteins
1554451_10_2_1	-1	10	1099.628222	1099.625856	0.000000	1.000000	0.333333	6.906755	0.246000	5.000000	0.0556	1099.628222	10	0	1	0	0	0	0	1	1	0	6.553933	0.000002	0.000002	R.AGSSLLWLPR.S	DECOY_245585	DECOY_305321	DECOY_308119
1554451_10_2_2	1	10	1099.628222	1098.624082	0.693147	1.000000	0.298781	6.906755	0.164000	2.000000	0.0714	1099.628222	8	0	1	0	0	0	0	1	1	1	6.553933	0.000914	0.000914	R.LHWLVM[15.9949]RK.N	204486	249065	314100
```

use ProtCosmo/pyComet/resource/common_unimod.csv

change the Peptide column modifications [15.9949] to [U:35] and remove the N and C terminal amino acids. for example `R.LHWLVM[15.9949]RK.N` change to `LHWLVMRK`

if it is `-.n[42.0106]M[15.9949]LQFLLEVNK.S` here n[42.0106] is the n terminal modification of peptides, change to `[U:1]-M[U:35]LQFLLEVNK`

For the Pin file, if there is more than one Proteins let's add option --proteins-keep default 1 
if more than one proteins, when converting, join the Proteins with ','

use argparse --input --output
also support argparse `--input-str` for single peptide conversion (no input/output file required)
also support argparse `--input-percolator`:
- input is treated as percolator output
- if input file suffix is not `.parquet`/`.parquet.gz`, parse as TSV even if suffix is not `.pin`
- `peptide` column is treated like Comet `Peptide` and gets the same unimod conversion
- `proteinIds` column is treated like Comet `Proteins` and gets `--proteins-keep` normalization
- **do not change column names** in percolator mode (only update values in existing columns)

if input can ends with .pin .pin.gz .pin.parquet .pin.parquet.gz
output can also be these format.
based on the file extention name, do format converstion.

example:
`python pyCometPinConverter.py --input-str "-.n[42.0106]M[15.9949]LQFLLEVNK.S"`
print:
`[U:1]-M[U:35]LQFLLEVNK`

default output format is .pin.parquet

if pin has extra columns not listed in the schema table below, keep them in output.
for those extra columns, parquet dtype infer rule is:
- int -> int32
- float -> float32
- others -> string

--input --output is a filename with path
if --output is not provided, then use the same name as --input with different extention name

may use functions fomr ProtCosmo/pyComet/pyUnimodMatch.py and may need other options


## module structure

- `ProtCosmo/pyComet/pyCometPinConverter.py`:
  CLI + format conversion orchestration + parquet/pin IO.
- `ProtCosmo/pyComet/utils/pin_converter_utils.py`:
  reusable unimod loader + peptide mass-tag to `[U:id]` conversion helpers.
