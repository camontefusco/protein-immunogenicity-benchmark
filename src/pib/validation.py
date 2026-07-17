"""Sequence validation helpers."""

CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def normalize_sequence(sequence: str) -> str:
    """Return an uppercase sequence without surrounding whitespace."""
    return sequence.strip().upper()


def is_valid_protein_sequence(sequence: str) -> bool:
    """Check whether a sequence contains only canonical amino acids."""
    normalized = normalize_sequence(sequence)
    return bool(normalized) and all(residue in CANONICAL_AMINO_ACIDS for residue in normalized)


def validate_sequences(sequences: list[str]) -> list[str]:
    """Normalize sequences and raise a clear error if any sequence is invalid."""
    normalized = [normalize_sequence(sequence) for sequence in sequences]
    invalid = [sequence for sequence in normalized if not is_valid_protein_sequence(sequence)]
    if invalid:
        examples = ", ".join(invalid[:5])
        raise ValueError(f"Invalid protein sequences found: {examples}")
    return normalized

