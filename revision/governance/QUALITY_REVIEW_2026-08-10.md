# Quality review - 2026-08-10

## Release decision

The recovery analysis is scientifically informative but is not yet cleared for
manuscript use. Two defects require rerun.

## Q-001: cross-design hyperparameter contamination

`10_tuned_stacked_comparison.py` selected parameters on the complete cluster-training
partition and reused them for exact and temporal tests. The selection data include
1,520 exact-test peptides and 1,538 temporal-test peptides. Cluster-test overlap is
zero.

- Valid: tuned cluster comparison.
- Quarantined: tuned exact and tuned temporal comparison and their intervals.
- Correction: tune independently inside each design's development partition, or use
  a genuinely pre-specified configuration.

## Q-002: incorrect context bootstrap unit

Context confidence intervals resampled context rows independently, although several
rows can share a peptide and edit-distance component.

- Valid: context point estimates.
- Quarantined: context confidence intervals and paired-difference intervals.
- Correction: paired hierarchical bootstrap by peptide for exact/temporal designs
  and by edit-distance component for the cluster design.

## Additional scope constraints

- The temporal peptide test has 110 edit-distance components spanning train/test.
- The temporal context test has 423 spanning components.
- No locked unified software environment exists yet.
- The original raw 161-column IEDB export and deterministic curation pipeline remain
  unavailable.
- Executed publication notebooks have not yet been generated.
