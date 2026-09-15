# Analysis governance

## Evidence classes

- **Historical:** submitted files and recovered notebooks. Immutable; may explain
  provenance but cannot establish reproducibility.
- **Validated:** clean rerun satisfying every publication gate.
- **Quarantined:** an output with known contamination, incorrect uncertainty unit,
  incomplete provenance or another material defect.
- **Exploratory:** useful diagnostic evidence not designated as confirmatory.

## Validation designs

### Cluster-disjoint primary analysis

Training and test observations share no connected component defined by Levenshtein
edit distance less than or equal to two. Hyperparameter selection and stacking must
occur entirely inside the cluster-training partition.

### Exact-peptide sensitivity analysis

Identical peptides are disjoint. Hyperparameter selection must occur entirely inside
the exact-training partition. Near-sequence overlap is measured and reported.

### Temporal sensitivity analysis

The test set contains peptides first reported on or after the frozen cutoff.
Hyperparameter selection must occur inside the pre-cutoff training partition.
Near-sequence cluster overlap is measured and reported. This is internal temporal
validation, not external validation.

## Statistical unit

- Peptide-level tasks: resample peptides for ordinary exact/temporal estimates and
  edit-distance components for cluster-disjoint estimates.
- Context tasks: never resample context rows independently. Resample peptide or
  edit-distance clusters, carrying all associated context rows together.
- Model differences use paired resamples of the same sampling units.

## Model selection

- PR-AUC is the sole primary selection criterion.
- Search spaces are declared before execution.
- Test data do not influence parameters, split seeds, thresholds, feature selection
  or stopping decisions.
- A stacker uses group-safe out-of-fold predictions only.

## Interpretation constraints

- IEDB-only evaluations are internal validation.
- Context prediction is database-conditional association, not causal biology.
- IFN-gamma assay labels are not vaccine efficacy, protection or clinical
  immunogenicity.
- Protein-language-model conclusions apply only to the evaluated checkpoint,
  pooling and classifier.
- N-gram associations are not HLA anchors without peptide-level HLA evidence.
