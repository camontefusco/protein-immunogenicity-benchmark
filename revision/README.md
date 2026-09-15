# Manuscript revision evidence layer

This directory is a publication-safe projection of the canonical local revision project. It separates validated primary analyses from supporting or exploratory work and excludes manuscript administration files and large model artifacts.

## Directory map

- `governance/`: analysis rules, quality review, reviewer-to-analysis mapping, and mandatory limitations.
- `configs/`: validation-design configuration.
- `scripts/`: deterministic audit, splitting, training, uncertainty, calibration, ESM-2, motif, and sensitivity scripts.
- `notebooks/`: six executed publication notebooks that consume registered artifacts and produce consolidated outputs.
- `data/splits/`: fixed peptide-level assignments used by the corrected analyses.
- `data/curated/`: recovered assay-level, peptide-level, and peptide-context tables used by the downstream reruns.
- `results/corrected_sequence/`: independently tuned per-design model metrics and grouped-bootstrap uncertainty.
- `results/context_ablation/`: grouped biological/assay ablations and individual-feature comparisons with corrected grouped uncertainty.
- `results/context_sensitivity/`: majority-label, tie-excluded, and consistent-label sensitivity results.
- `results/publication/`: manuscript-facing aggregate tables and figures.
- `provenance/`: observed runtime metadata and a portable publication dependency specification.

## Evidence hierarchy

1. Corrected per-design tuning and context sensitivity are the primary validated results.
2. Calibration and small frozen ESM-2 analyses are supporting experiments and must not be presented as comprehensive benchmarks.
3. Motif results are exploratory associations and cannot support HLA-anchor or mechanistic claims.
4. Historical results identified as invalid or superseded are not published here.

## Data boundary

`peptide_split_assignments_v1.csv` contains peptide-level assignments derived from the recovered curated dataset. `data/curated/assay_level_with_year.csv` and the recovered aggregation rule reproduce `peptide_level_hiconf_with_year.csv` exactly; `context_dataset_and_splits.csv` supplies the context-analysis input. The original 161-column IEDB export is supplied under `data/raw/` through Git LFS. Consequently, the source table and downstream analysis state are available for audit; exact reproduction of curation still depends on the documented field selection, normalization, and software environment.

## Recommended execution order

1. Review `governance/ANALYSIS_GOVERNANCE.md`, `governance/LIMITATIONS_REGISTER.md`, and `data/raw/README.md`.
2. Inspect the fixed assignments in `data/splits/`.
3. Run or audit the relevant scripts in numerical order. Scripts 17-20 reproduce the strict temporal, evidence-count, study-disjoint, and extended-calibration audits.
4. Review the six notebooks in `notebooks/` as executed publication records. They retain the canonical local-project layout assumptions; the consolidated CSV and PNG outputs are supplied under `results/publication/`.
5. Reconcile outputs against `results/CORRECTED_RESULTS_SUMMARY.md`.

Some supporting intermediate artifacts, serialized models, and manuscript administration files remain excluded. The supplied raw export plus curated tables and scripts support an auditable source-to-results workflow, subject to the stated environment and curation limitations.
