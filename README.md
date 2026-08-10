# Protein Immunogenicity Benchmark

Reproducibility repository for sequence-based ranking of peptide IFN-gamma assay outcomes and the revision of manuscript `CBAC-D-26-03155`.

## Repository status

The original repository contained a small portfolio starter built around amino-acid composition and a ten-row toy dataset. That starter package remains available in `src/`, `configs/`, `data/toy_peptides.csv`, and `tests/`.

The manuscript-revision evidence is under [`revision/`](revision/README.md). It contains:

- fixed peptide-level validation assignments;
- executable analysis scripts and publication notebooks;
- corrected independently tuned model comparisons;
- grouped-bootstrap confidence intervals;
- context-target sensitivity analyses;
- consolidated publication tables and figures;
- analysis governance, provenance, and an explicit limitations register.

## Headline corrected results

Best sequence-only PR-AUC by validation design:

| Validation design | Best model | PR-AUC |
|---|---:|---:|
| Exact peptide-disjoint | Stacked ensemble | 0.479 |
| Edit-distance <=2 cluster split | Random forest | 0.404 |
| Temporal peptide-disjoint | XGBoost | 0.151 |

For the majority-label context target, sequence-plus-context PR-AUC was 0.771, 0.755, and 0.487, respectively. The paired uplift over the matched sequence model was +0.236, +0.254, and +0.047. The direction persisted when tied outcomes were excluded and when analysis was restricted to completely consistent labels.

These are database-target prediction results. They do not establish a causal biological effect of context variables.

## Important limitations

- No independent external validation dataset is available.
- Temporal testing is peptide-disjoint but not simultaneously edit-distance-cluster-purged.
- Edit distance <=2 addresses one- and two-substitution similarity but is not CD-HIT percentage-identity clustering.
- The ESM-2 supporting analysis uses a small frozen model, not fine-tuning or a comprehensive PLM benchmark.
- The original raw 161-column IEDB export is unavailable, so curation cannot be reproduced from the original download.
- Reliable HLA restriction is absent, preventing defensible HLA-anchor validation.

See [`revision/governance/LIMITATIONS_REGISTER.md`](revision/governance/LIMITATIONS_REGISTER.md) for the controlled wording.

## Reproducing the publication layer

The notebooks in `revision/notebooks/` are executed records that read the registered aggregate artifacts and regenerate publication tables and figures. Model-training scripts are retained separately in `revision/scripts/`.

This repository intentionally excludes submission correspondence, author-identifying attachments, superseded/quarantined results, serialized model binaries, embedding arrays, and row-level test predictions.

