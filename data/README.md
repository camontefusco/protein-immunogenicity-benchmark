# Data

Expected benchmark input format:

```csv
sequence,label
SYFPEITHI,1
GILGFVFTL,0
```

For the first coding pass, use `toy_peptides.csv` only to verify that the software works.

For the real benchmark, candidate sources to evaluate carefully include:

- IEDB epitope assay exports
- curated peptide immunogenicity datasets from publications
- internal or manually curated ADA/immunogenicity sequence examples, if shareable

Do not mix assay types casually. Keep labels, species, HLA context, assay modality, and curation choices explicit.

