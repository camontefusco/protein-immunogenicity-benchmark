# Mandatory limitations register

The following limitations must appear explicitly in the revised manuscript and the
response to reviewers. They are not optional discussion points.

## L-001 - No external validation

All development and evaluation data originate from IEDB. Exact, cluster-disjoint
and temporal evaluations are internal validation or sensitivity analyses. None is
an independent external dataset generated under separate curation and laboratory
processes.

## L-002 - Temporal validation is not cluster-purged

The temporal peptide test is peptide-disjoint but has 110 edit-distance-2 components
spanning training and test. The temporal context test has 423 spanning components.
Temporal results must not be described as simultaneously sequence-cluster-disjoint.

## L-003 - Edit distance is not CD-HIT identity clustering

Connected components at Levenshtein edit distance less than or equal to two directly
test the reviewer's concern about one- or two-residue variants. This is not identical
to CD-HIT or another global percentage-sequence-identity threshold, particularly
across peptides of different lengths.

## L-004 - Limited protein-language-model comparison

The modern baseline uses the small frozen ESM-2 8M checkpoint, mean residue pooling
and logistic regression. It does not test larger checkpoints, fine-tuning,
task-specific representation learning, HLA-aware architectures or a comprehensive
PLM search.

## L-005 - Original raw export unavailable

The frozen assay and peptide tables are checksum-registered, but the original
161-column IEDB export and a deterministic raw-export-to-frozen-table curation run
have not been recovered. Data curation therefore cannot yet be reproduced from the
original download.

## L-006 - HLA-anchor validation unavailable

Reliable peptide-level HLA restriction is absent from the frozen modeling table.
N-gram coefficients and held-out enrichment are statistical associations and cannot
establish HLA anchors or mechanistic residues.
