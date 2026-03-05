**pyCometProteinDigest.py**

`pyCometProteinDigest.py` generates a 3-column internal novel peptide TSV:

```text
peptide	peptide_id	protein_id
```

`peptide_id` is assigned as `COMETPLUS_NOVEL_<n>` in first-seen order.

**CLI**

```bash
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometProteinDigest.py \
  [-P comet.params] \
  [-D db1.fasta -D db2.fasta ...] [--max-record N] \
  [--peptide <item> ...] \
  [--peptide-exclude <exclude.tsv|exclude.sqlite>] \
  -O output.tsv
```

**Parameters**

1. `-D, --database`: FASTA database file(s), digested with Comet settings from `--params`.
2. `--max-record`: max protein count to read from `--database`.
3. `-P, --params`: Comet params file (default `comet.params`).
4. `--peptide`: can be:
   - one peptide sequence,
   - comma/whitespace separated peptide list,
   - a text file (one or more tokens per line, lines starting with `#` ignored),
   - a FASTA file.
   `--peptide` and `--database` can be used together.
5. `--peptide-exclude` (optional):
   - TSV: use column `pep_seq` (header optional), or first column when no header is present.
   - SQLite: file must contain exactly one table with column `pep_seq`.
6. `-O, --output`: output TSV path.

**protein_id rules**

1. From `--database`: first token of FASTA header.
2. From `--peptide` token/text input: peptide sequence itself.
3. From `--peptide` FASTA input: FASTA header text.
4. If one peptide maps to multiple proteins, `protein_id` joins them with `;`.

**Behavior notes**

1. Peptides are normalized to uppercase alphabetic tokens.
2. Exclusion (`--peptide-exclude`) applies to both digested and `--peptide` inputs.
3. Output always includes header and keeps first-seen peptide order.

**Examples**

```bash
# digest databases only
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometProteinDigest.py \
  -P comet.params \
  -D known1.fasta -D known2.fasta \
  -O novel_input.tsv

# peptide tokens only
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometProteinDigest.py \
  --peptide "PEPTIDE,ACDEFGHIK" \
  -O novel_input.tsv

# mix database + peptide file, and exclude known peptides
/data/p/anaconda3/bin/python ProtCosmo/pyComet/pyCometProteinDigest.py \
  -P comet.params -D known.fasta \
  --peptide novel_peptides.fasta \
  --peptide-exclude known.peptide_sequence.sqlite \
  -O novel_input.tsv
```
