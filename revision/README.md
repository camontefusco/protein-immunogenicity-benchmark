# Manuscript revision evidence layer

This directory is a publication-safe projection of the canonical local revision project. It separates validated primary analyses from supporting or exploratory work and excludes manuscript administration files and large model artifacts.

## Directory map

- `governance/`: analysis rules, quality review, reviewer-to-analysis mapping, and mandatory limitations.
- `configs/`: validation-design configuration.
- `scripts/`: deterministic audit, splitting, training, uncertainty, calibration, ESM-2, motif, and sensitivity scripts.
- `notebooks/`: six executed publication notebooks that consume registered artifacts and produce consolidated outputs.
- `data/splits/`: fixed peptide-level assignments used by the corrected analyses.
- `results/corrected_sequence/`: independently tuned per-design model metrics and grouped-bootstrap uncertainty.
- `results/context_sensitivity/`: majority-label, tie-excluded, and consistent-label sensitivity results.
- `results/publication/`: manuscript-facing aggregate tables and figures.
- `provenance/`: observed runtime metadata and a portable publication dependency specification.

## Evidence hierarchy

1. Corrected per-design tuning and context sensitivity are the primary validated results.
2. Calibration and small frozen ESM-2 analyses are supporting experiments and must not be presented as comprehensive benchmarks.
3. Motif results are exploratory associations and cannot support HLA-anchor or mechanistic claims.
4. Historical results identified as invalid or superseded are not published here.

## Data boundary

`peptide_split_assignments_v1.csv` contains peptide-level assignments derived from the locally recovered curated dataset. The original raw 161-column IEDB download is unavailable. Consequently, this repository supports reproduction from the curated analysis state, but not independent reconstruction of the original curation from the raw download.

## Recommended execution order

1. Review `governance/ANALYSIS_GOVERNANCE.md` and `governance/LIMITATIONS_REGISTER.md`.
2. Inspect the fixed assignments in `data/splits/`.
3. Run or audit the relevant scripts in numerical order after supplying the recovered curated inputs through the documented environment variables.
4. Review the six notebooks in `notebooks/` as executed publication records. They retain the canonical local-project layout assumptions; the consolidated CSV and PNG outputs are supplied under `results/publication/`.
5. Reconcile outputs against `results/CORRECTED_RESULTS_SUMMARY.md`.

Because the original raw IEDB export and some supporting intermediate artifacts cannot be redistributed, this projection is not a claim of raw-download-to-paper reproducibility.
