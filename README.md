# Protein Immunogenicity Benchmark

Reproducible analysis of peptide-level IFN-gamma assay outcomes from an IEDB export. This repository contains the data audit, leakage-aware validation designs, model reruns, context sensitivity analyses, calibration checks, frozen ESM-2 comparison, publication tables, and figures used for the manuscript revision.

## Project layout

```text
data/          raw IEDB export, curated tables, toy example, and fixed split assignments
configs/       validation-design configuration
scripts/       executable audits, splits, training, uncertainty, calibration, and sensitivities
notebooks/     six executed publication notebooks
results/       validated metrics, forensic audits, manuscript tables, and figures
governance/    analysis rules, reviewer crosswalk, quality reviews, and limitations
provenance/    runtime metadata and publication dependency specification
src/           reusable baseline package retained from the starter project
tests/         unit tests
```

The repository uses a regular project layout; there is no separate revision-only directory.

## Data and provenance

The source file is `data/raw/tcell_table_export_1769046013.csv`, stored with Git LFS. It contains 31,629 IEDB records and 161 columns. Its SHA-256 and Drive provenance are recorded in `data/raw/README.md`.

The downstream assay-level, peptide-level, and peptide-context tables under `data/curated/` are derived products. Fixed peptide assignments under `data/splits/` are reused across model comparisons.

## Validated analysis scope

The primary sequence task is peptide-level ranking. The context task is separate, with a different unit, target, prevalence, and split structure. Absolute PR-AUC values must not be compared between these tasks; context is interpreted through paired uplift against its matched sequence-only model.

Corrected headline results are exact-peptide PR-AUC 0.479, edit-distance-2 cluster PR-AUC 0.404, and strict temporal sequence-plus-context PR-AUC 0.363 versus 0.324 sequence-only (paired uplift 0.039; 95% CI 0.016 to 0.061).

These are retrospective database-target results. They do not establish external generalization, HLA-specific recognition, vaccine efficacy, protection, causal biological effects, or individual clinical immunogenicity.

## Reproduce or audit

1. Install `provenance/requirements-publication.txt`.
2. Review `governance/ANALYSIS_GOVERNANCE.md` and `governance/LIMITATIONS_REGISTER.md`.
3. Inspect `data/splits/` and result provenance files.
4. Run scripts in numerical order when a full rerun is required; scripts 17 to 20 contain forensic audits.
5. Open the six notebooks in `notebooks/` for the executed publication workflow.
6. Compare outputs with `results/CORRECTED_RESULTS_SUMMARY.md` and `results/publication/tables/`.

All publication notebooks use repository-relative paths and do not require a Google Drive mount.

## Manuscript crosswalk

`results/publication/figures/` contains the four manuscript figures and graphical abstract. `results/publication/tables/` contains manuscript-facing aggregate tables. Manuscript DOCX files and the response letter remain in the paper project because they contain submission formatting and correspondence; this repository contains the computational evidence needed to audit their reported methods and results.

## Historical Drive notebooks

The Drive archive contains earlier Colab notebooks from exploratory development. Several use obsolete paths, superseded splits, or earlier claims, so they are not mixed into the validated workflow. Their names, Drive IDs, and status are recorded in `governance/DRIVE_NOTEBOOK_INVENTORY.md`. The six notebooks in `notebooks/` are authoritative for the publication layer.

## Limitations

- No independent external validation dataset is available.
- Temporal testing is peptide-disjoint but not simultaneously cluster-purged.
- Edit distance <=2 is not identical to CD-HIT percentage-identity clustering.
- ESM-2 is a small frozen model, not fine-tuning or a comprehensive PLM benchmark.
- Reliable HLA restriction is absent, preventing defensible HLA-anchor validation.
- One frozen outer test set is used per design; grouped-bootstrap intervals are conditional on those sets.
- Exact reproduction remains dependent on documented curation rules and software environment.

See `governance/LIMITATIONS_REGISTER.md` for controlled wording.
