# Reviewer-to-analysis matrix

Status values are `planned`, `running`, `complete`, or `manuscript-only`.

| ID | Reviewer concern | Required evidence/change | Status |
|---|---|---|---|
| R1.1 | Novelty versus prior IFN-gamma predictors | Reframe as a leakage-safe benchmarking and screening-utility study; compare explicitly with prior evaluation designs and state the new knowledge generated. | manuscript-only |
| R1.2 | Context uplift may encode database/assay bias | Pre-specified single-feature and grouped ablations; sequence-only, biological-context-only, assay-context-only, evidence-only, and exclusion models. | complete |
| R1.3 | IEDB-only internal validation | Attempt temporal and source/laboratory holdouts; describe them as robustness analyses, not independent external validation; explicitly retain the external-validation limitation. | complete |
| R1.4 | Similar peptides may cross splits | Cluster peptides by sequence identity and evaluate cluster-disjoint splits at documented thresholds; quantify exact and near-neighbor overlap. | complete |
| R1.5 | Neural comparison omits pretrained models | Add at least one frozen pretrained protein-language-model embedding baseline if computationally feasible; moderate all deep-learning conclusions regardless. | complete |
| R1.6 | Small model differences lack uncertainty | Paired stratified bootstrap confidence intervals for PR-AUC, ROC-AUC, Top-K precision and differences on identical test predictions. | complete |
| R1.7 | Add sequencing-ML references | Evaluate PMID 34730875 and PMID 33848577 for relevance and cite accurately. | manuscript-only |
| R1.8 | Hyperparameter search insufficiently described | Save search spaces, selection metric, folds/groups, iteration counts, random seeds, early stopping and final parameters in machine-readable artifacts. | complete |
| R1.9 | Class-imbalance handling unclear | Record class weights/resampling/threshold rules per model; use threshold-free PR-AUC as primary and keep threshold optimization nested. | complete |
| R1.10 | Motif interpretation descriptive | Test enrichment/position and physicochemical summaries; map cautiously to known HLA anchors only where allele information supports it. | complete |
| R1.11 | Translational claims overstate assay prediction | Replace vaccine-efficacy/clinical-immunogenicity implications with IFN-gamma assay-outcome prioritization language. | manuscript-only |
| R1.12 | Discussion repetitive | Condense repeated sequence-only limitation statements. | manuscript-only |
| R1.13 | “Sequence-only ceiling” undefined | Remove the phrase or define it operationally as performance under the specified representation, data and validation design—not a biological upper bound. | manuscript-only |
| R1.14 | Duplicate reference | Deduplicate Patronov and Doytchinova and audit all references. | manuscript-only |
| R1.15 | Calibration metrics incomplete | Report Brier score and ECE with an explicit binning definition, alongside reliability curves. | complete |
| R1.16 | Computational efficiency absent | Record wall-clock training/inference time, hardware, peak memory where available, serialized model size and parameter count where applicable. | complete |
| R1.17 | Modern protein-language-model baselines requested | Same analysis as R1.5, reported separately from models trained from scratch. | complete |

## Scientific decision rules

- No historical point estimate will be retained solely because it appears in the submitted manuscript.
- A revised number must be traceable to saved predictions from a named split and configuration.
- Context variables measured after, or mathematically derived from, the assay outcome will be excluded from deployable prediction models and may only appear in clearly labeled bias-diagnostic analyses.
- Temporal or database-source holdouts will not be described as external validation unless the evaluation dataset is genuinely independent of the development database and curation process.
