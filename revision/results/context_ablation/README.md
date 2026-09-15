# Validated context ablation

This directory promotes the recovered, frozen-split context ablations into the
canonical project. `metrics.csv` contains grouped feature-set comparisons;
`individual_feature_ablation.csv` contains feature-only and
sequence-plus-one-feature comparisons. The original row-level bootstrap was not
used. Manuscript uncertainty comes from the corrected grouped files, which
resample peptides for exact/temporal designs and edit-distance-2 components for
the component design.

Principal models exclude `n_assays_context`, `positive_fraction`, and other
outcome-derived evidence variables. The context analysis uses a fixed balanced
logistic regression (`C=1.0`) for all feature sets, so it measures differences
under a common classifier rather than separately optimized maxima.

The self-contained input is
`02_data/curated/context_dataset_and_splits.csv`. This is a recovered curated
analysis table; it does not restore the missing original 161-column IEDB export
or the unrecovered upstream curation path.
