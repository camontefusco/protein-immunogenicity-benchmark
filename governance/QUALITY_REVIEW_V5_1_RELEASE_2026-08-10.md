# Quality review — v5.1 release (2026-08-10)

## Decision

The narrow v5.1 correction pass clears the four findings from the independent v5 audit. No model rerun or numerical result changed.

## Corrections verified

- The abstract identifies the exact peptide-disjoint design without the ambiguous adverb “exactly.”
- The manuscript and reviewer response distinguish held-out statistical motif enrichment from biological validation. They state that HLA-anchor validation is not defensible and that a post hoc physicochemical interpretation was not retained without prespecified independent validation.
- Figure 4 now reports “No motifs passed FDR < 0.05 temporally,” rather than claiming temporal absence.
- The manuscript, Highlights, and Cover Letter consistently state that the small frozen ESM-2 baseline did not outperform the strongest classical model in any design.

## QA checks

- Clean manuscript: no revision markup.
- Tracked manuscript: cumulative insertions and deletions retained.
- Bibliography: all 29 references cited; no undefined citations.
- Obsolete submitted headline claims remain absent.
- Revised manuscript, response, Highlights, Cover Letter, and Title Page rendered without clipping or overlap.
- Corrected Figure 4 and graphical abstract inspected visually.
- Supplementary workbook is numerically unchanged from the verified v5 workbook.
- Repository tests: 4 passed.

## Retained limitations

No independent external dataset; one frozen outer test set per design; one three-fold grouped CV run with modest grids and one seed; temporal test not simultaneously component-purged; edit distance <=2 is not CD-HIT identity clustering; small frozen ESM-2 only; unavailable original raw export; absent reliable HLA restriction; heterogeneous retrospective assay labels.

The package is suitable for final author reading before journal resubmission.
