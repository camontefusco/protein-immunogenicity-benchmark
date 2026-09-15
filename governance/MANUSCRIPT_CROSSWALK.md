# Manuscript crosswalk

This file records which repository artifacts support the manuscript claims and flags draft statements that must not be copied into the final manuscript without checking the corrected results.

| Manuscript topic | Repository evidence | Interpretation |
|---|---|---|
| Dataset construction | `data/raw/`, `data/curated/`, `scripts/01_audit_dataset.py` | 31,629 raw IEDB rows were reduced to 31,502 curated assay records and then aggregated into peptide and context tasks. |
| Sequence discrimination | `results/publication/figures/figure_01_sequence_models.png`, `results/publication/tables/table_02_sequence_performance.csv` | Report design-specific PR-AUC and uncertainty; do not present a single pooled score as general immunogenicity performance. |
| Context sensitivity | `results/publication/figures/figure_02_context_sensitivity.png`, `results/publication/tables/table_03_context_sensitivity.csv`, `results/context_sensitivity/` | Use paired uplift and state that the context task has a different unit and target. |
| Calibration | `results/publication/figures/figure_03_calibration.png`, `results/publication/tables/table_04_calibration.csv` | Treat probabilities as decision-support scores, not individual response probabilities. |
| ESM-2 and motifs | `results/publication/figures/figure_04_esm2_motifs.png`, `results/publication/tables/table_05_esm2_comparison.csv` | Describe this as a small frozen-model/supporting experiment and exploratory association analysis. |
| Temporal sensitivity | `results/forensic_sensitivities/`, `scripts/17_forensic_sensitivities.py`, `scripts/18_strict_context_sensitivity.py` | The strict-censor temporal context result is 0.363 versus 0.324 sequence-only, not the earlier 0.487 estimate. |
| Reviewer audit trail | `governance/reviewer_analysis_matrix.md`, `governance/QUALITY_REVIEW_V5_1_RELEASE_2026-08-10.md` | Use these files to map reviewer requests to analyses and controlled limitations. |

## Draft claims requiring correction

Older manuscript drafts report four or seven-figure layouts and context values such as 0.534 to 0.756. Those values belong to superseded analysis states and are not the authoritative repository results. The current publication layer has four figures and uses the corrected validation designs, context constructions, strict temporal sensitivity, and explicit limitations.

The repository does not claim external validation, HLA-specific recognition, vaccine efficacy, protection, causal biological effects, or individual clinical immunogenicity.
