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
- biological, assay, and individual-feature context ablations;
- consolidated publication tables and figures;
- analysis governance, provenance, and an explicit limitations register.

## Headline corrected results

Best sequence-only PR-AUC by validation design:

| Validation design | Best model | PR-AUC |
|---|---:|---:|
| Exact peptide-disjoint | Stacked ensemble | 0.479 |
| Edit-distance <=2 cluster split | Random forest | 0.404 |
| Temporal peptide-disjoint | XGBoost | 0.151 |

For the majority-label context target, the cross-sectional sequence-plus-context PR-AUC was 0.771 and 0.755 for exact-peptide and edit-distance-2 designs. The original temporal estimate of 0.487 is superseded by the strict-censor rerun: 0.363 versus 0.324 for the matched sequence-only model, with paired uplift +0.039 (95% CI 0.016 to 0.061). The temporal result is therefore a small, assay-metadata-driven uplift rather than evidence of broad temporal generalization.

The principal sequence benchmark contains 9,668 valid high-confidence peptides (9.2% positive). The peptide-context analysis is a separate 20,785-row task (30.0% positive) with a different target, unit, prevalence, and split assignments. Absolute PR-AUC values are not comparable across these tasks; only the paired context uplift against its matched sequence-only model is interpreted.

These are database-target prediction results. They do not establish a causal biological effect of context variables.

## Important limitations

- No independent external validation dataset is available.
- Temporal testing is peptide-disjoint but not simultaneously edit-distance-cluster-purged.
- Edit distance <=2 addresses one- and two-substitution similarity but is not CD-HIT percentage-identity clustering.
- The ESM-2 supporting analysis uses a small frozen model, not fine-tuning or a comprehensive PLM benchmark.
- The original raw 161-column IEDB export is included under `revision/data/raw/` through Git LFS. Its SHA-256 is recorded in `revision/data/raw/README.md`; the same file is archived in the project Drive folder.
- Reliable HLA restriction is absent, preventing defensible HLA-anchor validation.
- One frozen outer test set is used per design; grouped-bootstrap intervals are conditional on those sets.
- Sequence tuning used one three-fold grouped CV run, modest grids, and one seed; context and ESM-2 downstream classifiers were fixed.
- Efficiency records omit peak memory and complete stacked-model inference time and bundle size.

See [`revision/governance/LIMITATIONS_REGISTER.md`](revision/governance/LIMITATIONS_REGISTER.md) for the controlled wording.

## Provenance and manuscript crosswalk

The raw IEDB export (31,629 rows and 161 columns) is the source table. The curated assay-level and peptide-level tables under `revision/data/curated/` are downstream products of the documented curation and aggregation steps. The fixed splits, rerun scripts, notebooks, aggregate tables, and figures are the auditable analysis layer. Manuscript and response-letter DOCX files are intentionally kept outside this code repository; the repository contains the evidence needed to verify their reported methods and results.

The Drive folder used during recovery is the project archive: <https://drive.google.com/drive/folders/1l1aj9JQ1F1E7hSVDe3eb38_1Puuhzmg6>.

## Reproducing the publication layer

The notebooks in `revision/notebooks/` are executed records that read the registered aggregate artifacts and regenerate publication tables and figures. Model-training scripts are retained separately in `revision/scripts/`.

This repository intentionally excludes submission correspondence, author-identifying attachments, superseded/quarantined results, serialized model binaries, embedding arrays, and row-level test predictions.
