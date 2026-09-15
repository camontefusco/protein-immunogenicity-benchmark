# Quality review — v5 release (2026-08-10)

## Release decision

The v5 revision package is internally consistent and suitable for author review before journal resubmission. The rerun supports a leakage-aware internal benchmark, not a clinical, causal, or externally validated predictor.

## Checks completed

- The clean manuscript contains no revision markup; the tracked copy retains cumulative insertions and deletions relative to the submitted manuscript.
- All 29 bibliography entries are cited and no undefined numerical citations remain.
- Obsolete submitted headline claims (9.3-fold enrichment and PR-AUC 0.534 to 0.756) were removed from the manuscript and submission documents.
- Exact peptide-disjoint, edit-distance <=2 component-disjoint, and temporal peptide-disjoint designs are distinguished consistently.
- Model tuning, one-hot handling, grouped bootstrap uncertainty, calibration, frozen ESM-2 scope, and Random Forest Figure 4 comparator are described explicitly.
- The context sensitivity analysis compares majority labels, exclusion of 662 ties, and consistent-label-only rows.
- Outcome-derived evidence variables are excluded from principal context models. Retained evidence rows are explicitly labelled retrospective bias diagnostics, not deployable models.
- The manuscript, response, highlights, cover letter, title page, figures, graphical abstract, and eight-sheet supplementary workbook were rendered or inspected for layout.

## Interpretation supported by the rerun

- Sequence ranking performance depends strongly on the evaluation design: best PR-AUC was 0.479 exact, 0.404 component-disjoint, and 0.151 temporal.
- Safe context variables were associated with positive peptide-context outcomes across three label constructions, but temporal uplift was smaller.
- The small frozen ESM-2 baseline did not outperform the strongest classical sequence model under the frozen tests.
- Motif analysis is exploratory and cannot be interpreted as HLA-anchor validation.

## Limitations retained

- No independent external dataset.
- Only one frozen outer test set per design; grouped bootstrap intervals are conditional on those test sets.
- One three-fold grouped CV run, modest hyperparameter grids, and one seed; downstream context and ESM-2 classifiers were fixed.
- Temporal tests are peptide-disjoint but not simultaneously edit-distance-component-purged.
- Edit distance <=2 addresses one/two substitutions but is not CD-HIT percentage-identity clustering.
- The ESM-2 experiment uses a small frozen model, not fine-tuning or a comprehensive PLM benchmark.
- The original raw 161-column IEDB export is unavailable, preventing reconstruction of upstream curation from the original download.
- Reliable HLA restriction is absent.
- Efficiency records omit peak memory and complete stacked inference/serialized size.

## Remaining author action

Read the tracked manuscript and response once for preferred tone, then replace the repository URL with a DOI/accession if an archival release is created before resubmission.
