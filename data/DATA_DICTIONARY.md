# Data dictionary

## Raw source

`data/raw/tcell_table_export_1769046013.csv` is the 161-column IEDB export used to initiate curation. Its checksum and Drive provenance are recorded beside the file.

## Curated assay table

| Field | Meaning |
|---|---|
| `peptide` | Normalized uppercase peptide sequence |
| `length` | Peptide length |
| `label` | Binary qualitative assay outcome |
| `assay_method` | IEDB assay platform |
| `readout` | Functional readout metadata |
| `outcome_raw` | Original qualitative outcome text |
| `n_tested`, `n_positive` | Available response counts |
| `response_freq_pct` | Recorded response frequency |
| `virus`, `virus_species`, `virus_strain` | Viral source metadata |
| `antigen`, `protein` | Antigen and source-protein metadata |
| `earliest_pub_year`, `latest_pub_year` | Publication chronology used for temporal audits |

## Peptide-level table

`peptide_level_hiconf_with_year.csv` contains peptides with unambiguous qualitative labels. `label_hiconf` is the high-confidence binary target; conflicted peptides are excluded from this target. `pos_rate` and evidence-count fields are descriptive and should not be confused with prospective individual response probabilities.

## Context table

`context_dataset_and_splits.csv` is a separate peptide-context task. It contains viral, antigen, protein, assay, host, evidence-count, target-construction, and split-assignment fields. Its rows and target are not interchangeable with the peptide-only benchmark.
