# Quality review of revision v4

Date: 2026-08-10

## Release decision

**Cleared as a resubmission draft, subject to author verification of journal
administrative requirements and final figure-file upload.** No remaining
manuscript statement requires a principal-model rerun.

## Resolved v3 findings

- The primary sequence cohort is now defined as 9,672 recovered
  high-confidence peptides (890 positive; 9.2%), with four invalid sequences
  excluded and 9,668 modeled.
- The peptide-context task is explicitly defined as a separate 20,785-row,
  30.0%-positive task with a different target, unit, prevalence, and splits.
- The text states that absolute PR-AUC values cannot be compared between the
  sequence and context tasks; context uplift is paired only within the context
  task.
- The stale legacy Figure 1 cross-reference was removed.
- Table 1 and Supplementary Tables S1-S7 were rebuilt into one reconciled
  submission workbook matching the manuscript captions.
- Individual-feature values are explicitly described as descriptive rather
  than independent importance estimates.
- Reviewer statuses for external validation, near-neighbour clustering, PLM
  scope, HLA interpretation, and incomplete efficiency accounting now say that
  the concern was addressed to the extent feasible with the limitation
  retained.

## Verification

- Abstract length: 215 words.
- Clean manuscript: 12 rendered pages, visually inspected.
- Cumulative tracked manuscript: 18 rendered pages, 58 insertions and 64
  deletions, visually inspected.
- Reviewer response: 5 rendered pages, visually inspected.
- Submission workbook: Table 1 and Tables S1-S7 rendered and visually inspected;
  no formula-error tokens detected.
- No stale submitted headline values (`0.864`, `0.756`) remain.
- No claim of external validation, causal context effect, HLA anchor mechanism,
  vaccine efficacy, clinical immunogenicity, or comprehensive PLM superiority
  is made.

## Limitations intentionally retained

- No independent external dataset.
- One frozen outer test set per design.
- Modest three-fold grouped CV grids and one seed.
- Temporal split is peptide-disjoint but not component-purged.
- Edit distance <=2 is not CD-HIT percentage-identity clustering.
- Small frozen ESM-2 rather than fine-tuning or a comprehensive PLM benchmark.
- Missing original 161-column IEDB export and unreproducible upstream curation.
- Missing reliable HLA restriction.
- Incomplete peak-memory and complete stacked-model efficiency accounting.
