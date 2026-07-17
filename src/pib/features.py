"""Feature extraction for peptide and protein sequences."""

from __future__ import annotations

import pandas as pd

from pib.validation import CANONICAL_AMINO_ACIDS, validate_sequences


def amino_acid_composition(sequences: list[str]) -> pd.DataFrame:
    """Compute normalized amino acid composition features for each sequence."""
    normalized = validate_sequences(sequences)
    residues = sorted(CANONICAL_AMINO_ACIDS)
    rows: list[dict[str, float]] = []

    for sequence in normalized:
        length = len(sequence)
        row = {f"aa_frac_{residue}": sequence.count(residue) / length for residue in residues}
        row["length"] = float(length)
        rows.append(row)

    return pd.DataFrame(rows)

