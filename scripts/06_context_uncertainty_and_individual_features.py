from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler


SEED = 42
N_BOOT = 2000
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "context_ablation"
CONTEXT = OUT / "context_dataset_and_splits.csv"
PREDICTIONS = OUT / "test_predictions.csv"
BIOLOGICAL = ["virus_species", "virus_strain", "antigen", "protein"]
ASSAY = ["assay_method", "readout", "host"]
FEATURES = BIOLOGICAL + ASSAY + ["length", "log_evidence_count"]


def metric(y: np.ndarray, probability: np.ndarray, name: str) -> float:
    if name == "pr_auc":
        return float(average_precision_score(y, probability))
    return float(roc_auc_score(y, probability))


def pipeline(feature: str, include_sequence: bool) -> Pipeline:
    transformers = []
    if include_sequence:
        transformers.append(
            (
                "sequence",
                TfidfVectorizer(
                    analyzer="char", ngram_range=(2, 4), lowercase=False,
                    min_df=2, sublinear_tf=True, dtype=np.float32,
                ),
                "peptide",
            )
        )
    if feature in BIOLOGICAL + ASSAY:
        transformers.append(("feature", OneHotEncoder(handle_unknown="ignore", min_frequency=2), [feature]))
    else:
        transformers.append(
            (
                "feature",
                Pipeline(
                    [
                        ("log1p", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
                        ("scale", StandardScaler()),
                    ]
                ),
                [feature],
            )
        )
    return Pipeline(
        [
            ("features", ColumnTransformer(transformers)),
            ("classifier", LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=SEED, solver="liblinear")),
        ]
    )


# Paired, outcome-stratified uncertainty for the pre-specified main context models.
pred = pd.read_csv(PREDICTIONS)
ci_rows, diff_rows = [], []
for split_name, split_frame in pred.groupby("split", sort=True):
    wide = split_frame.pivot(index="row_id", columns="model", values="probability")
    labels = split_frame.drop_duplicates("row_id").set_index("row_id")["y"].loc[wide.index]
    y = labels.to_numpy(dtype=int)
    positive, negative = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    rng = np.random.default_rng(SEED)
    samples = [
        np.concatenate([rng.choice(positive, len(positive), replace=True), rng.choice(negative, len(negative), replace=True)])
        for _ in range(N_BOOT)
    ]
    for model_name in wide.columns:
        probability = wide[model_name].to_numpy(dtype=float)
        for metric_name in ("pr_auc", "roc_auc"):
            values = np.array([metric(y[idx], probability[idx], metric_name) for idx in samples])
            ci_rows.append(
                {
                    "split": split_name,
                    "model": model_name,
                    "metric": metric_name,
                    "point_estimate": metric(y, probability, metric_name),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                }
            )
    sequence_probability = wide["sequence_only"].to_numpy(dtype=float)
    for model_name in sorted(set(wide.columns) - {"sequence_only"}):
        probability = wide[model_name].to_numpy(dtype=float)
        for metric_name in ("pr_auc", "roc_auc"):
            values = np.array(
                [metric(y[idx], probability[idx], metric_name) - metric(y[idx], sequence_probability[idx], metric_name) for idx in samples]
            )
            diff_rows.append(
                {
                    "split": split_name,
                    "model_a": model_name,
                    "model_b": "sequence_only",
                    "metric": metric_name,
                    "difference_a_minus_b": metric(y, probability, metric_name) - metric(y, sequence_probability, metric_name),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                }
            )

pd.DataFrame(ci_rows).to_csv(OUT / "bootstrap_confidence_intervals.csv", index=False)
pd.DataFrame(diff_rows).to_csv(OUT / "paired_differences_vs_sequence.csv", index=False)

# Individual feature contribution: feature alone and sequence plus exactly one
# feature, using the already frozen context-row split assignments.
context = pd.read_csv(CONTEXT)
individual_rows = []
for split_name in ("exact", "cluster", "temporal"):
    split_column = f"{split_name}_split"
    train = context[context[split_column] == "train"]
    test = context[context[split_column] == "test"]
    y_train, y_test = train["y"].to_numpy(dtype=int), test["y"].to_numpy(dtype=int)
    for feature in FEATURES:
        for include_sequence in (False, True):
            model = pipeline(feature, include_sequence)
            columns = (["peptide"] if include_sequence else []) + [feature]
            start = time.perf_counter()
            model.fit(train[columns], y_train)
            probability = model.predict_proba(test[columns])[:, 1]
            elapsed = time.perf_counter() - start
            individual_rows.append(
                {
                    "split": split_name,
                    "feature": feature,
                    "configuration": "sequence_plus_feature" if include_sequence else "feature_only",
                    "pr_auc": float(average_precision_score(y_test, probability)),
                    "roc_auc": float(roc_auc_score(y_test, probability)),
                    "fit_and_predict_seconds": elapsed,
                    "n_train": len(train),
                    "n_test": len(test),
                }
            )
            print(split_name, feature, include_sequence, f"PR-AUC={individual_rows[-1]['pr_auc']:.4f}")

pd.DataFrame(individual_rows).to_csv(OUT / "individual_feature_ablation.csv", index=False)
print("\nMain PR-AUC intervals")
print(pd.DataFrame(ci_rows).query("metric == 'pr_auc'").to_string(index=False))
