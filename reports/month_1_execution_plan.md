# Month 1 Execution Plan

## Week 1: Classical Baseline

Goal: make the project runnable, testable, and scientifically honest.

Coding tasks:

1. Create the package layout.
2. Add sequence validation.
3. Add amino acid composition features.
4. Train a logistic regression baseline.
5. Report AUROC, average precision, accuracy, and calibration error.

Scientific tasks:

1. Define the prediction problem carefully.
2. Decide what counts as an immunogenic positive label.
3. Track assay type, species, peptide length, and HLA context.
4. Write down likely sources of label noise and dataset leakage.

## Week 2: Real Dataset

Goal: replace the toy data with a defensible curated dataset.

Candidate source:

- IEDB assay exports, filtered into a clear binary task.

Minimum metadata to preserve:

- source database
- assay type
- host species
- HLA allele or restriction, if available
- peptide length
- source antigen or protein
- publication or accession reference

## Week 3: Protein Language Model Embeddings

Goal: add one modern representation-learning baseline.

Start with ESM2 embeddings through Hugging Face or the ESM package.

Compare:

- amino acid composition baseline
- ESM2 frozen embeddings plus logistic regression
- ESM2 frozen embeddings plus gradient boosting or shallow neural net

## Week 4: Error Analysis

Goal: show scientific maturity, not just metrics.

Analyze:

- false positives and false negatives
- performance by peptide length
- performance by source antigen family
- calibration curves
- whether near-duplicate sequences leak across folds

