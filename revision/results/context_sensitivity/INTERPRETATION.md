# Context-label construction sensitivity

Three target definitions were evaluated on unchanged split assignments and identical
feature preprocessing:

1. majority labels with exact 0.5 ties assigned negative (20,785 rows);
2. majority labels after excluding 662 ties (20,123 rows);
3. only context pairs with completely consistent assay labels, excluding all 1,087
   mixed-label rows (19,698 rows).

## Sequence plus all safe context PR-AUC

| Construction | Exact | Cluster | Temporal |
|---|---:|---:|---:|
| Majority, all rows | 0.771 | 0.755 | 0.487 |
| Majority, ties excluded | 0.788 | 0.778 | 0.527 |
| Consistent labels only | 0.791 | 0.782 | 0.528 |

## Context uplift over the matched sequence-only model

| Construction | Exact difference (95% CI) | Cluster difference (95% CI) | Temporal difference (95% CI) |
|---|---:|---:|---:|
| Majority, all rows | +0.236 (0.205 to 0.266) | +0.254 (0.218 to 0.293) | +0.047 (0.022 to 0.072) |
| Majority, ties excluded | +0.232 (0.201 to 0.263) | +0.259 (0.221 to 0.298) | +0.060 (0.033 to 0.085) |
| Consistent labels only | +0.235 (0.203 to 0.267) | +0.272 (0.234 to 0.312) | +0.066 (0.035 to 0.094) |

The main qualitative conclusion is not produced by assigning tied outcomes to the
negative class. Removing ties or all mixed-label contexts modestly increases both
sequence-only and context-model PR-AUC, consistent with a cleaner and easier target.
The cross-sectional context uplift remains large, whereas temporal discrimination
remains substantially below the cross-sectional result.

For context without sequence, the temporal difference versus sequence-only is not
clear under the majority constructions: +0.018 (95% CI -0.015 to 0.051) with all
rows and +0.022 (-0.013 to 0.057) after removing ties. In the consistent-label subset
it is modestly positive at +0.039 (0.002 to 0.075). This subset-dependent result
reinforces the need to present context performance as target- and database-dependent.

This sensitivity analysis does not convert the context association into evidence of
causal biological determinants or external generalization. Virus, strain, antigen
and protein categories can encode historical database and study-composition
structure. The manuscript should use the fully consistent-label analysis as a
sensitivity result and retain the majority construction as the primary context task
only if its operational rationale is explained explicitly.
