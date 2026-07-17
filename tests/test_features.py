import pandas as pd
import pytest

from pib.features import amino_acid_composition


def test_amino_acid_composition_returns_expected_columns() -> None:
    features = amino_acid_composition(["ACD"])

    assert isinstance(features, pd.DataFrame)
    assert features.loc[0, "aa_frac_A"] == pytest.approx(1 / 3)
    assert features.loc[0, "aa_frac_C"] == pytest.approx(1 / 3)
    assert features.loc[0, "aa_frac_D"] == pytest.approx(1 / 3)
    assert features.loc[0, "length"] == 3.0


def test_amino_acid_composition_rejects_invalid_sequences() -> None:
    with pytest.raises(ValueError, match="Invalid protein sequences"):
        amino_acid_composition(["PEPTIDE*"])

