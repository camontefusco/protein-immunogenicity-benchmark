from __future__ import annotations

import json
import os
import platform
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.base import clone
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline


SEED = 42
ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "artifacts" / "splits" / "peptide_split_assignments.csv"
OUT = ROOT / "artifacts" / "sequence_baselines"
MODEL_DIR = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)


def vectorizer() -> TfidfVectorizer:
    return TfidfVectorizer(
        analyzer="char",
        ngram_range=(2, 4),
        lowercase=False,
        min_df=2,
        sublinear_tf=True,
        dtype=np.float32,
    )


MODELS = {
    "dummy_prior": DummyClassifier(strategy="prior"),
    "logistic_l2": LogisticRegression(
        C=1.0,
        class_weight="balanced",
        max_iter=3000,
        random_state=SEED,
        solver="liblinear",
    ),
    "sgd_elasticnet": SGDClassifier(
        loss="log_loss",
        penalty="elasticnet",
        alpha=1e-5,
        l1_ratio=0.15,
        class_weight="balanced",
        max_iter=3000,
        tol=1e-4,
        random_state=SEED,
    ),
    "random_forest": RandomForestClassifier(
        n_estimators=500,
        max_features="sqrt",
        min_samples_leaf=2,
        class_weight="balanced_subsample",
        random_state=SEED,
        n_jobs=-1,
    ),
    "extra_trees": ExtraTreesClassifier(
        n_estimators=500,
        max_features="sqrt",
        min_samples_leaf=2,
        class_weight="balanced",
        random_state=SEED,
        n_jobs=-1,
    ),
}


def expected_calibration_error(y: np.ndarray, p: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    bucket = np.clip(np.digitize(p, edges[1:-1], right=True), 0, bins - 1)
    result = 0.0
    for idx in range(bins):
        mask = bucket == idx
        if mask.any():
            result += mask.mean() * abs(y[mask].mean() - p[mask].mean())
    return float(result)


def topk(y: np.ndarray, p: np.ndarray, k: int) -> tuple[float, float]:
    k = min(k, len(y))
    order = np.argsort(-p, kind="stable")[:k]
    precision = float(y[order].mean())
    return precision, precision / float(y.mean())


data = pd.read_csv(SPLITS)
split_columns = {"exact": "exact_split", "cluster": "cluster_split", "temporal": "temporal_split"}
metrics_rows: list[dict] = []
predictions: list[pd.DataFrame] = []
efficiency_rows: list[dict] = []

for split_name, split_column in split_columns.items():
    train = data[data[split_column] == "train"].copy()
    test = data[data[split_column] == "test"].copy()
    x_train, y_train = train["peptide"], train["label_hiconf"].to_numpy(dtype=int)
    x_test, y_test = test["peptide"], test["label_hiconf"].to_numpy(dtype=int)

    for model_name, estimator in MODELS.items():
        pipe = Pipeline([("tfidf", vectorizer()), ("model", clone(estimator))])
        start = time.perf_counter()
        pipe.fit(x_train, y_train)
        train_seconds = time.perf_counter() - start

        start = time.perf_counter()
        probability = pipe.predict_proba(x_test)[:, 1]
        inference_seconds = time.perf_counter() - start

        model_path = MODEL_DIR / f"{split_name}__{model_name}.joblib"
        joblib.dump(pipe, model_path, compress=3)

        row = {
            "split": split_name,
            "model": model_name,
            "n_train": len(train),
            "n_test": len(test),
            "test_positive_rate": float(y_test.mean()),
            "roc_auc": float(roc_auc_score(y_test, probability)),
            "pr_auc": float(average_precision_score(y_test, probability)),
            "brier_score": float(brier_score_loss(y_test, probability)),
            "ece_10_equal_width": expected_calibration_error(y_test, probability),
        }
        for k in (20, 50, 100, 200, 500):
            row[f"precision_at_{k}"], row[f"enrichment_at_{k}"] = topk(y_test, probability, k)
        metrics_rows.append(row)

        frame = test[["peptide", "label_hiconf", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
        frame.insert(0, "split", split_name)
        frame.insert(1, "model", model_name)
        frame["probability"] = probability
        predictions.append(frame)

        efficiency_rows.append(
            {
                "split": split_name,
                "model": model_name,
                "train_seconds": train_seconds,
                "inference_seconds": inference_seconds,
                "inference_microseconds_per_peptide": inference_seconds * 1e6 / len(test),
                "serialized_model_bytes": model_path.stat().st_size,
            }
        )
        print(split_name, model_name, f"PR-AUC={row['pr_auc']:.4f}", f"ROC-AUC={row['roc_auc']:.4f}")

pd.DataFrame(metrics_rows).to_csv(OUT / "metrics.csv", index=False)
pd.concat(predictions, ignore_index=True).to_csv(OUT / "test_predictions.csv", index=False)
pd.DataFrame(efficiency_rows).to_csv(OUT / "efficiency.csv", index=False)

environment = {
    "python": sys.version,
    "platform": platform.platform(),
    "processor": platform.processor(),
    "logical_cpu_count": os.cpu_count(),
    "numpy": np.__version__,
    "pandas": pd.__version__,
    "scikit_learn": sklearn.__version__,
    "seed": SEED,
    "threading_note": "Tree ensembles used n_jobs=-1; timings are wall-clock measurements on this host.",
}
with (OUT / "environment.json").open("w", encoding="utf-8") as stream:
    json.dump(environment, stream, indent=2, sort_keys=True)
