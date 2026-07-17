# Protein Immunogenicity Benchmark

Portfolio project for sequence-to-function modeling of immunogenic peptide/protein candidates.

This repo is designed to become a compact, reproducible benchmark that maps directly to biomolecular ML roles:

- peptide/protein sequence featurization
- baseline sequence-to-function prediction
- protein language model embeddings
- calibration and uncertainty checks
- scientific error analysis

## First Milestone

Build a clean classical baseline before adding large protein language models.

The initial pipeline:

1. Load a peptide-level dataset with columns `sequence` and `label`.
2. Validate amino acid sequences.
3. Convert each sequence into interpretable amino acid composition features.
4. Train a scikit-learn classifier.
5. Report AUROC, average precision, accuracy, and calibration error.

## Why This Project Exists

The Bayer role asks for evidence of:

- sequence-to-function modeling
- protein language models
- biomolecular interaction and therapeutic candidate workflows
- scientific rigor
- collaborative software engineering practices

This project is intentionally structured to show those signals in code, documentation, and analysis.

## Repository Layout

```text
protein-immunogenicity-benchmark/
  configs/              Experiment configs
  data/                 Local datasets or download notes
  notebooks/            Exploratory analyses
  reports/              Figures, tables, and written findings
  src/pib/              Reusable Python package
  tests/                Lightweight regression tests
```

## Quick Start

Create a virtual environment and install the package:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Run tests:

```bash
pytest
```

Run the starter baseline on the toy dataset:

```bash
python -m pib.train_baseline --data data/toy_peptides.csv
```

## Three-Month Target

By the end of this project, the repo should include:

- classical peptide descriptors
- ESM2 or ProtT5 embeddings
- at least two model families
- calibration and uncertainty evaluation
- biologically grounded error analysis
- reproducible CLI workflows
- a polished technical write-up

