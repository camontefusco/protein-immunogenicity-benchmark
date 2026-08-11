# Forensic metric and leakage audit

Date: 2026-08-11  
Scope: every reported discrimination, early-retrieval, calibration, context, sensitivity, motif, and efficiency result in revision v5.1, plus targeted reruns where the original metric could be biased.

## Executive decision

The numerical values in the saved prediction files are internally consistent: sequence discrimination, Brier, and Top-K metrics recompute to numerical precision. There is no evidence that metrics were fabricated or calculated from the wrong prediction column.

The current v5.1 manuscript nevertheless should not be submitted unchanged. One material temporal-label problem and several interpretation problems were found:

1. The temporal context training labels include assay evidence published in or after the 2021 test period. This is label look-ahead, not merely ordinary dataset shift. Strict time censoring lowers temporal context PR-AUC from 0.487 to 0.363 and the matched sequence-only value from 0.440 to 0.324. The paired uplift remains positive (0.039; grouped-bootstrap 95% CI 0.016-0.061).
2. Exact-split Top-20 precision of 0.90 is real for that test set but is predominantly a near-neighbour result: 90% of the top 20 have a same-length training peptide with at least 80% identity, and 75% have at least 90% identity. It must not be presented as broad novelty generalization.
3. Cross-sectional context metrics are strongly compatible with study/library proxy learning. Biological context performs well cross-sectionally but fails temporally after strict censoring. The uplift is associative and data-source-specific.
4. Temporal sequence probabilities remain worse than a constant-prevalence predictor by Brier score even after Platt scaling (Brier skill score -0.047). Temporal scores may be discussed as weak rankings, not as reliable probabilities.
5. Efficiency fields omit feature-extraction time; stack timing and size have different scopes. The current efficiency table is not an apples-to-apples computational benchmark.
6. The recovered assay-level table and recovered notebook reproduce the 9,672-row high-confidence peptide table exactly. The manuscript's statement that the historical aggregation rule cannot be reconstructed is incorrect. Upstream reconstruction from the original 161-column IEDB export remains impossible because that export is missing.

## Metric-by-metric assessment

### Dataset counts and labels

- Assay-level counts (31,502 assays, 17,336 peptides, 10,233 positive and 21,269 negative records) are consistent with the recovered curated assay table.
- Peptide benchmark counts (9,672 labeled peptides; 9,668 canonical sequences; 890 positives before sequence validation) are consistent.
- The high-confidence label is reproducible from the recovered assay table: sum tested and positive responses by peptide; assign missing if assay labels conflict; otherwise assign positive when the pooled positive fraction is at least 0.5, negative when it is zero, and missing otherwise. The reconstruction has zero label mismatches with the recovered peptide table.
- “High confidence” should be treated as a historical table name, not a guarantee of strong evidence per peptide. 1,973 labeled peptides have fewer than five tested responses, 3,654 have fewer than ten, and 9,153 are supported by at most one reference.
- A minimum-total-tested sensitivity (at least 10) changes prevalence and absolute PR-AUC but preserves the main pattern: strong cross-sectional ranking and weak temporal transfer. Because prevalence differs, raw PR-AUC values from this sensitivity must not be compared directly with the principal values.

### Sequence PR-AUC and ROC-AUC

- All saved PR-AUC and ROC-AUC values recompute exactly from the test predictions.
- Exact split: stacking PR-AUC 0.479 and Random Forest 0.472 are plausible but optimistic for novel-sequence use because exact identity alone is excluded; close sequence neighbours and shared-study libraries remain.
- Edit-distance-component split: Random Forest PR-AUC 0.404 is the most defensible cross-sectional sequence result among the principal designs. It prevents connected one/two-edit neighbours but does not prevent shared publications or all percentage-identity similarities.
- Temporal split: best PR-AUC 0.151 is only modestly above prevalence (0.116). Strictly removing training peptides whose evidence extends into 2021 or later gives XGBoost PR-AUC 0.148. Therefore the weak temporal sequence conclusion is stable.
- No model family consistently dominates. Paired differences are small relative to split uncertainty, so architectural-superiority language is unwarranted.

### Top-K precision and enrichment

- The calculations are correct, but K=20 is only 20 observations and has high sampling uncertainty.
- Exact Random Forest Precision@20=0.90 is near-neighbour-driven: 18/20 top-ranked peptides have at least 80% same-length identity to a training peptide and 15/20 have at least 90% identity.
- Component Random Forest Precision@20=0.65 occurs without any top-20 peptide reaching 80% same-length identity, but study/library dependence remains.
- Temporal early retrieval is not useful: the original best reported Precision@20=0.10 is below prevalence enrichment; under strict temporal censoring XGBoost Precision@20=0.05.
- Recommendation: retain Top-K only as secondary descriptive evidence, add uncertainty or the underlying numerator, and explicitly distinguish close-neighbour triage from novel-family generalization.

### Bootstrap confidence intervals and paired differences

- Grouped resampling is aligned with each frozen split: peptide groups for exact/temporal and edit-distance components for the component design.
- Paired differences use shared bootstrap replicates and are methodologically appropriate.
- Intervals are conditional on one selected outer split. They do not quantify variability across alternative split realizations, time cutoffs, studies, or curation choices.
- The component outer split was selected from 512 candidates to match test size and prevalence without model performance. This is acceptable but produces a deliberately balanced test and should not be confused with repeated external validation.

### Context PR-AUC and ablations

- Cross-sectional context metrics recompute from the saved predictions. They are not directly comparable with the 9,668-peptide benchmark because the row unit, target, prevalence, and split population differ.
- Exact/component context PR-AUC values around 0.75-0.77 are plausible as retrospective repository prediction, but too high for a biological-generalization claim. Almost all test categories are seen in training and publications overlap heavily; virus/protein/strain fields can identify study libraries and collection practices.
- The original temporal context construction has material label look-ahead: 2,745 of 16,723 temporal training rows have supporting peptide evidence extending to 2021 or later.
- Strict temporal rerun, majority construction: sequence-only PR-AUC 0.324; sequence+assay 0.386; sequence+all safe context 0.363; biological-only 0.286; sequence+biological 0.281. Thus assay context, not stable biological context, supplies the temporal uplift.
- Strict temporal target sensitivities preserve the direction of paired all-context uplift:
  - majority including ties: +0.039 (95% CI 0.016-0.061);
  - ties excluded: +0.037 (0.013-0.060);
  - consistent repeated labels only: +0.043 (0.017-0.068).
- Recommendation: replace all temporal context values in the manuscript and response with strict-censor values; keep cross-sectional context results but label them study-structure-sensitive retrospective associations.

### Study/publication sensitivity

- Shared PMID is not an input feature, so this is not direct feature leakage. It is a dependence/confounding problem: 98.0% of exact and 95.8% of component test peptides share at least one PMID with training; temporal overlap is 64.2%.
- A single PMID-component-disjoint sensitivity reduces Random Forest PR-AUC to 0.141. A combined PMID- and edit-distance-component-disjoint sensitivity gives 0.207. These are one-split diagnostics with very large components, not definitive replacement estimates.
- Recommendation: report this sensitivity as evidence that study-disjoint external evaluation is needed. Do not replace the main frozen designs with one post hoc study split.

### Calibration, Brier score, and ECE

- Calibration computations reproduce. Platt scaling leaves ranking metrics unchanged, as expected.
- Exact calibrated Random Forest: Brier 0.070 versus null 0.084; Brier skill +0.156; equal-width ECE 0.019; equal-frequency ECE 0.029.
- Component calibrated Random Forest: Brier 0.076 versus null 0.084; Brier skill +0.085; equal-width ECE 0.024; equal-frequency ECE 0.034.
- Temporal calibrated Random Forest: Brier 0.108 versus null 0.103; Brier skill -0.047; equal-width ECE 0.060; equal-frequency ECE 0.058.
- ECE alone makes temporal calibration look more favorable than Brier skill. The probabilities are not good enough for transported threshold decisions.
- Context probabilities are especially poor temporally: strict all-context Brier 0.313 versus a null Brier near 0.207. PR-AUC uplift does not imply calibrated probability improvement.
- Recommendation: add null Brier/Brier skill, state that calibration only adds value cross-sectionally, and remove the claim that temporal calibrated scores support resource-allocation thresholds.

### Frozen ESM-2 baseline

- The design is a small frozen encoder plus a trained regularized downstream classifier. It is not fine-tuning and not a broad PLM benchmark.
- The reported PR-AUC values (0.217 exact, 0.194 component, 0.147 temporal) do not look implausibly high and do not improve on the strongest classical model.
- No direct leakage mechanism was identified in the described embedding workflow. A strict temporal downstream refit would be cleaner but is unlikely to change the scientific conclusion because only 89 sequence-training peptides are removed and the existing temporal result is already weak.
- Recommendation: retain as a bounded representation check; do not infer that PLMs or deep learning are generally inferior.

### Motif counts and FDR

- Motif selection uses training data and Fisher tests use test data, avoiding direct reuse of the same labels for selection and testing.
- The Fisher/BH analysis treats peptide rows as independent. Peptides within study libraries and sequence components can be correlated, so the quoted FDR is not cluster-adjusted and may be anti-conservative.
- The absence of temporal hits is consistent with instability, but counts of 5 and 9 must remain exploratory and cannot support HLA-anchor or mechanistic claims.
- Recommendation: describe these as peptide-level exploratory Fisher/BH associations, or rerun with study/component-aware permutation before emphasizing FDR control.

### Efficiency and serialized size

- Base-model fit timing starts after TF-IDF vectorization; inference timing excludes test transformation. These are estimator-only timings, not complete pipeline timings.
- Stack fit time includes out-of-fold/base fitting, but stack inference is missing and serialized size contains only the meta-learner. Rows therefore have different measurement boundaries.
- Recommendation: either rerun end-to-end wall time and serialize complete fitted pipelines on the same platform, or relabel the table as incomplete diagnostic measurements and remove comparative efficiency conclusions. Peak memory remains unmeasured.

## Leakage taxonomy

| Issue | Classification | Severity | Required action |
|---|---|---:|---|
| Exact peptide overlap | Prevented in all principal tests | None | Retain integrity check |
| One/two-edit overlap | Prevented only in component design | Expected design difference | Preserve explicit scope |
| Future evidence in temporal sequence labels | Small (89 training peptides); sensitivity stable | Low for ranking conclusion | Report strict sensitivity |
| Future evidence in temporal context labels | Material label look-ahead | High | Replace temporal context results |
| Shared publications/study libraries | Dependence and proxy confounding, not direct column leakage | High for generalization claims | Add sensitivity and limitation |
| Outcome-derived context counts | Excluded from principal model | Prevented | Retain exclusion audit |
| Category fitting on test | Not detected; encoders fit on development data | None detected | Retain pipeline checks |
| Hyperparameter selection on test | Not detected | None detected | Retain frozen-test provenance |
| Motif hypothesis selection/testing reuse | Prevented by train/test separation | Low | Add cluster-dependence caveat |

## Required manuscript corrections before submission

1. Replace temporal context metrics with strict-censor values and revise Figure 2, supplementary context tables, abstract, Results, Discussion, response letter, highlights, and graphical abstract.
2. Correct the recovered-data provenance: downstream high-confidence aggregation is exactly reproducible from the recovered assay table/notebook; only original-export download and upstream curation remain unavailable.
3. Add study/publication dependence and the post hoc study-disjoint sensitivity, without presenting one split as external validation.
4. Add the near-neighbour composition of exact Top-20 and lower its rhetorical prominence.
5. Add null Brier and Brier skill; state that temporal and context probabilities are not reliable for absolute thresholding.
6. Relabel or rerun efficiency metrics with consistent end-to-end measurement boundaries.
7. Qualify motif FDR as peptide-level and not cluster-adjusted.
8. Keep all existing external-validation, HLA, edit-distance-versus-CD-HIT, single-split, modest-CV, and small-frozen-ESM-2 limitations.

## Bottom line

The project still supports a defensible paper, but the defensible contribution is an audit of how validation design, time, study structure, and target construction change retrospective IFN-gamma ranking. It does not support a high-performance general immunogenicity predictor. The strongest trustworthy claims are: (i) arithmetic and frozen-split implementation are reproducible downstream; (ii) sequence ranking degrades under near-neighbour separation and becomes weak temporally; (iii) strict temporal context provides a small ranking uplift driven mainly by assay metadata; and (iv) neither absolute probabilities nor external transportability have been established.
