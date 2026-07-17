"""Train a starter amino-acid-composition baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from pib.features import amino_acid_composition
from pib.metrics import binary_classification_report


def load_dataset(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Load a sequence classification dataset."""
    data = pd.read_csv(path)
    required_columns = {"sequence", "label"}
    missing = required_columns.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    features = amino_acid_composition(data["sequence"].tolist())
    labels = data["label"].astype(int)
    return features, labels


def train_and_evaluate(data_path: Path) -> dict[str, float]:
    """Evaluate a logistic regression baseline with stratified cross-validation."""
    x, y = load_dataset(data_path)
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=1000, class_weight="balanced"),
    )
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
    probabilities = cross_val_predict(model, x, y, cv=cv, method="predict_proba")[:, 1]
    return binary_classification_report(y.to_numpy(), probabilities)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="CSV with sequence and label columns.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metrics = train_and_evaluate(args.data)
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

