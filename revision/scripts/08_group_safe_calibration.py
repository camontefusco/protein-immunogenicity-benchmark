from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline


SEED = 42
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "splits" / "peptide_split_assignments.csv"
OUT = ROOT / "artifacts" / "calibration"
MODEL_DIR = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)


def base_model() -> Pipeline:
    return Pipeline(
        [
            (
                "tfidf",
                TfidfVectorizer(
                    analyzer="char", ngram_range=(2, 4), lowercase=False,
                    min_df=2, sublinear_tf=True, dtype=np.float32,
                ),
            ),
            (
                "classifier",
                RandomForestClassifier(
                    n_estimators=500, min_samples_leaf=2, max_features="sqrt",
                    class_weight="balanced_subsample", random_state=SEED, n_jobs=-1,
                ),
            ),
        ]
    )


def reliability(y: np.ndarray, probability: np.ndarray, bins: int = 10) -> pd.DataFrame:
    edges = np.linspace(0, 1, bins + 1)
    membership = np.clip(np.digitize(probability, edges[1:-1], right=True), 0, bins - 1)
    rows = []
    for idx in range(bins):
        mask = membership == idx
        rows.append(
            {
                "bin": idx,
                "lower": edges[idx],
                "upper": edges[idx + 1],
                "n": int(mask.sum()),
                "mean_probability": float(probability[mask].mean()) if mask.any() else np.nan,
                "observed_fraction": float(y[mask].mean()) if mask.any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


def ece(table: pd.DataFrame) -> float:
    populated = table[table["n"] > 0]
    return float(
        (
            populated["n"] / populated["n"].sum()
            * (populated["observed_fraction"] - populated["mean_probability"]).abs()
        ).sum()
    )


data = pd.read_csv(SOURCE)
metrics, predictions, reliability_tables, split_details = [], [], [], []
for split_name, split_column in {
    "exact": "exact_split", "cluster": "cluster_split", "temporal": "temporal_split"
}.items():
    development = data[data[split_column] == "train"].copy()
    test = data[data[split_column] == "test"].copy()

    # Keep edit-distance clusters intact between model fitting and calibration.
    splitter = GroupShuffleSplit(n_splits=256, test_size=0.20, random_state=SEED)
    candidates = []
    for candidate_id, (fit_idx, calibration_idx) in enumerate(
        splitter.split(development, development["label_hiconf"], development["edit_distance_2_cluster"])
    ):
        calibration = development.iloc[calibration_idx]
        score = abs(len(calibration) / len(development) - 0.20) + abs(
            calibration["label_hiconf"].mean() - development["label_hiconf"].mean()
        )
        candidates.append((score, candidate_id, fit_idx, calibration_idx))
    _, selected_candidate, fit_idx, calibration_idx = min(candidates, key=lambda value: value[:2])
    fit, calibration = development.iloc[fit_idx], development.iloc[calibration_idx]

    model = base_model()
    model.fit(fit["peptide"], fit["label_hiconf"])
    calibration_probability = model.predict_proba(calibration["peptide"])[:, 1]
    test_probability = model.predict_proba(test["peptide"])[:, 1]

    eps = 1e-6
    calibration_logit = np.log(np.clip(calibration_probability, eps, 1 - eps) / np.clip(1 - calibration_probability, eps, 1 - eps))
    test_logit = np.log(np.clip(test_probability, eps, 1 - eps) / np.clip(1 - test_probability, eps, 1 - eps))
    platt = LogisticRegression(C=1e6, max_iter=2000, random_state=SEED, solver="lbfgs")
    platt.fit(calibration_logit.reshape(-1, 1), calibration["label_hiconf"])
    calibrated_probability = platt.predict_proba(test_logit.reshape(-1, 1))[:, 1]

    joblib.dump(model, MODEL_DIR / f"{split_name}__base_random_forest.joblib", compress=3)
    joblib.dump(platt, MODEL_DIR / f"{split_name}__platt_calibrator.joblib", compress=3)

    for status, probability in (("uncalibrated", test_probability), ("platt", calibrated_probability)):
        table = reliability(test["label_hiconf"].to_numpy(dtype=int), probability)
        table.insert(0, "split", split_name)
        table.insert(1, "calibration", status)
        reliability_tables.append(table)
        metrics.append(
            {
                "split": split_name,
                "calibration": status,
                "n_fit": len(fit),
                "n_calibration": len(calibration),
                "n_test": len(test),
                "roc_auc": float(roc_auc_score(test["label_hiconf"], probability)),
                "pr_auc": float(average_precision_score(test["label_hiconf"], probability)),
                "brier_score": float(brier_score_loss(test["label_hiconf"], probability)),
                "ece_10_equal_width": ece(table),
            }
        )
        result = test[["peptide", "label_hiconf", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
        result.insert(0, "split", split_name)
        result.insert(1, "calibration", status)
        result["probability"] = probability
        predictions.append(result)

    split_details.append(
        {
            "split": split_name,
            "selected_candidate": selected_candidate,
            "fit_calibration_cluster_overlap": len(
                set(fit["edit_distance_2_cluster"]) & set(calibration["edit_distance_2_cluster"])
            ),
            "fit_test_peptide_overlap": len(set(fit["peptide"]) & set(test["peptide"])),
            "calibration_test_peptide_overlap": len(set(calibration["peptide"]) & set(test["peptide"])),
        }
    )

pd.DataFrame(metrics).to_csv(OUT / "metrics.csv", index=False)
pd.concat(predictions, ignore_index=True).to_csv(OUT / "test_predictions.csv", index=False)
pd.concat(reliability_tables, ignore_index=True).to_csv(OUT / "reliability_bins.csv", index=False)
with (OUT / "split_integrity.json").open("w", encoding="utf-8") as stream:
    json.dump(split_details, stream, indent=2, sort_keys=True)
print(pd.DataFrame(metrics).to_string(index=False))
print(pd.DataFrame(split_details).to_string(index=False))
