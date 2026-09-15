# Protein Immunogenicity Benchmark

Reproducible analysis of peptide-level IFN-gamma assay outcomes from an IEDB export. This repository contains the data audit, leakage-aware validation designs, model reruns, context sensitivity analyses, calibration checks, frozen ESM-2 comparison, publication tables, and figures used for the manuscript revision.

## Project layout

```text
data/          raw IEDB export, curated tables, toy example, and fixed split assignments
configs/       validation-design configuration
scripts/       executable audits, splits, training, uncertainty, calibration, and sensitivities
notebooks/     six executed publication notebooks
apps/          interactive Streamlit app provenance and deployment notes
results/       validated metrics, forensic audits, manuscript tables, and figures
governance/    analysis rules, reviewer crosswalk, quality reviews, and limitations
provenance/    runtime metadata and publication dependency specification
src/           reusable baseline package retained from the starter project
tests/         unit tests
CITATION.cff   citation metadata
```

The repository uses a regular project layout; there is no separate revision-only directory.

## Data and provenance

The source file is `data/raw/tcell_table_export_1769046013.csv`, stored with Git LFS. It contains 31,629 IEDB records and 161 columns. Its SHA-256 and Drive provenance are recorded in `data/raw/README.md`.

The downstream assay-level, peptide-level, and peptide-context tables under `data/curated/` are derived products. Fixed peptide assignments under `data/splits/` are reused across model comparisons.

The field meanings and target boundaries are summarized in [data/DATA_DICTIONARY.md](data/DATA_DICTIONARY.md).

## Validated analysis scope

The primary sequence task is peptide-level ranking. The context task is separate, with a different unit, target, prevalence, and split structure. Absolute PR-AUC values must not be compared between these tasks; context is interpreted through paired uplift against its matched sequence-only model.

Corrected headline results are exact-peptide PR-AUC 0.479, edit-distance-2 cluster PR-AUC 0.404, and strict temporal sequence-plus-context PR-AUC 0.363 versus 0.324 sequence-only (paired uplift 0.039; 95% CI 0.016 to 0.061).

The main visual results are linked here:

- [Figure 1 sequence-model discrimination](results/publication/figures/figure_01_sequence_models.png): sequence-only performance across exact, cluster, and temporal validation designs.
- [Figure 2 context sensitivity](results/publication/figures/figure_02_context_sensitivity.png): context-target constructions and paired uplift against sequence-only models.
- [Figure 3 calibration](results/publication/figures/figure_03_calibration.png): group-safe reliability and calibration summaries.
- [Figure 4 ESM-2 and motif analysis](results/publication/figures/figure_04_esm2_motifs.png): the small frozen protein-language-model comparison and exploratory motif enrichment.
- [Graphical abstract](results/publication/figures/graphical_abstract_v5.png): overview of the benchmark workflow and its interpretation.

The corresponding machine-readable tables are [sequence performance](results/publication/tables/table_02_sequence_performance.csv), [context sensitivity](results/publication/tables/table_03_context_sensitivity.csv), [calibration](results/publication/tables/table_04_calibration.csv), and [ESM-2 comparison](results/publication/tables/table_05_esm2_comparison.csv). The complete corrected summary is in [CORRECTED_RESULTS_SUMMARY.md](results/CORRECTED_RESULTS_SUMMARY.md).

These are retrospective database-target results. They do not establish external generalization, HLA-specific recognition, vaccine efficacy, protection, causal biological effects, or individual clinical immunogenicity.

## Conclusion

The benchmark shows that peptide sequence contains reproducible but bounded signal for ranking IFN-gamma assay outcomes. Performance decreases when validation becomes more stringent, especially for edit-distance-2 and temporal separation. Adding assay and biological context improves ranking in the separate peptide-context task, but the strict temporal uplift is small and is driven mainly by assay metadata. The practical conclusion is therefore methodological: leakage-safe splits, explicit target definitions, early-retrieval metrics, calibration checks, and context sensitivity analyses are necessary before interpreting public IFN-gamma prediction results. The repository supports auditing and prioritization research, not a validated clinical or vaccine-response predictor.

## Use in real-world problems

The results are most useful as a pre-experimental triage and audit workflow:

- **Screening large candidate lists:** use sequence-only scores to rank candidates for a limited experimental shortlist, then report Precision@K and enrichment at the actual screening budget rather than relying only on ROC-AUC.
- **Choosing validation designs:** use exact-peptide, edit-distance-2, and temporal splits to distinguish performance on repeated-like candidates, near-neighbour candidates, and later evidence. A model should not be deployed based on a random peptide split alone.
- **Using assay metadata responsibly:** context features can improve retrospective ranking when the relevant metadata will genuinely be available at prediction time. Metadata that encode future assay evidence must be excluded from prospective use.
- **Planning experiments:** calibrated outputs can support relative prioritization and threshold selection for follow-up assays, but they should not be interpreted as the probability that a peptide will respond in a particular person.
- **Auditing database-derived models:** the raw export, curated tables, split assignments, provenance files, and forensic analyses allow teams to check label construction, duplicate handling, temporal overlap, and sensitivity to tied or inconsistent context labels.
- **Designing the next model:** reliable HLA restriction, prospective external data, study-disjoint evaluation, and larger frozen protein-language-model comparisons are the main additions needed before making stronger translational claims.

In a real deployment, the workflow should be treated as a ranking aid alongside antigen-processing, HLA-binding, safety, conservation, and laboratory evidence. It should not be used as a stand-alone decision rule for vaccine composition, patient treatment, or claims of protective immunity.

## Reproduce or audit

1. Install `provenance/requirements-publication.txt`.
2. Review `governance/ANALYSIS_GOVERNANCE.md` and `governance/LIMITATIONS_REGISTER.md`.
3. Inspect `data/splits/` and result provenance files.
4. Run scripts in numerical order when a full rerun is required; scripts 17 to 20 contain forensic audits.
5. Open the six notebooks in `notebooks/` for the executed publication workflow.
6. Compare outputs with `results/CORRECTED_RESULTS_SUMMARY.md` and `results/publication/tables/`.

For a fast integrity check without rerunning analyses, run `python scripts/audit_repository.py`. It verifies required artifacts, the raw-export checksum, and the validation configuration.

All publication notebooks use repository-relative paths and do not require a Google Drive mount.

## Interactive application

The Drive archive also contains an `iedb_streamlit_app` with `app.py`, `requirements.txt`, and a short README. It is an interactive exploratory interface for browsing or scoring IEDB-derived records, not the source of the validated publication metrics. The repository records its provenance under `apps/iedb_streamlit/README.md`; the batch scripts and publication notebooks remain authoritative for the manuscript results.

## Manuscript crosswalk

`results/publication/figures/` contains the four manuscript figures and graphical abstract. `results/publication/tables/` contains manuscript-facing aggregate tables. Manuscript DOCX files and the response letter remain in the paper project because they contain submission formatting and correspondence; this repository contains the computational evidence needed to audit their reported methods and results.

The manuscript-to-repository mapping, including the corrected values that supersede older draft prose, is documented in [MANUSCRIPT_CROSSWALK.md](governance/MANUSCRIPT_CROSSWALK.md).

## Historical Drive notebooks

The Drive archive contains earlier Colab notebooks from exploratory development. Several use obsolete paths, superseded splits, or earlier claims, so they are not mixed into the validated workflow. Their names, Drive IDs, and status are recorded in `governance/DRIVE_NOTEBOOK_INVENTORY.md`. The six notebooks in `notebooks/` are authoritative for the publication layer.

## Limitations

- No independent external validation dataset is available.
- Temporal testing is peptide-disjoint but not simultaneously cluster-purged.
- Edit distance <=2 is not identical to CD-HIT percentage-identity clustering.
- ESM-2 is a small frozen model, not fine-tuning or a comprehensive PLM benchmark.
- Reliable HLA restriction is absent, preventing defensible HLA-anchor validation.
- One frozen outer test set is used per design; grouped-bootstrap intervals are conditional on those sets.
- Exact reproduction remains dependent on documented curation rules and software environment.

See `governance/LIMITATIONS_REGISTER.md` for controlled wording.
