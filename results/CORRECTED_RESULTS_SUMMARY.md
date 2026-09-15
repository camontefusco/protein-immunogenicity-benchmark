# Corrected results summary

## Sequence models

Hyperparameters were selected independently within each validation design's
development partition using three-fold StratifiedGroupKFold and PR-AUC. Test
peptide overlap with the corresponding development partition is zero.

| Design | ElasticNet logistic | Random Forest | XGBoost | Stack |
|---|---:|---:|---:|---:|
| Exact-peptide | 0.361 | 0.472 | 0.409 | 0.479 |
| Edit-distance-2 cluster | 0.200 | 0.404 | 0.204 | 0.372 |
| Temporal | 0.122 | 0.145 | 0.151 | 0.147 |

PR-AUC 95% bootstrap intervals use peptides for exact/temporal designs and complete
edit-distance components for the cluster design:

- Exact stack: 0.479 (0.406 to 0.552).
- Cluster Random Forest: 0.404 (0.273 to 0.546).
- Temporal XGBoost: 0.151 (0.126 to 0.183).

Random Forest minus stack is -0.007 (95% CI -0.019 to 0.005) for exact,
+0.032 (-0.005 to 0.064) for cluster, and -0.003 (-0.005 to -0.0004) for temporal.
No practically meaningful stacking advantage is established. All temporal models
show weak absolute performance.

## Context models

Point estimates are unchanged because the original fit/test partitions were valid.
Uncertainty now uses a paired hierarchical bootstrap, carrying all context rows for
each sampled peptide or edit-distance component.

- Exact sequence plus all safe context: 0.771 (0.746 to 0.795); difference versus
  sequence-only +0.236 (0.205 to 0.266).
- Cluster sequence plus all safe context: 0.755 (0.725 to 0.785); difference +0.254
  (0.218 to 0.293).
- Temporal sequence plus all safe context: 0.487 (0.456 to 0.519); difference +0.047
  (0.022 to 0.072).
- Temporal context without sequence versus sequence-only: +0.018 (-0.015 to 0.051).
- Temporal biological metadata alone performs below sequence-only by -0.053
  (-0.080 to -0.024).

These results support strong cross-sectional database-conditional context signal but
limited temporal transfer. They do not demonstrate causal biological determinants
or external generalization.
# Post-release forensic audit (2026-08-11)

Revision v5.1 temporal context values are superseded by the strict-censor results in `forensic_sensitivities/`. The corrected majority-label temporal context PR-AUC is 0.363 versus 0.324 for the matched sequence-only model (paired uplift 0.039; grouped-bootstrap 95% CI 0.016-0.061). See `revision/governance/FORENSIC_METRIC_AND_LEAKAGE_AUDIT_2026-08-11.md` for the complete metric, calibration, near-neighbour, study-dependence, and provenance review.
