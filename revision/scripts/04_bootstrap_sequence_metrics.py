from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


SEED = 42
N_BOOT = 2000
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "sequence_baselines" / "test_predictions.csv"
OUT = ROOT / "artifacts" / "sequence_baselines"


def metric(y: np.ndarray, probability: np.ndarray, name: str) -> float:
    if name == "pr_auc":
        return float(average_precision_score(y, probability))
    if name == "roc_auc":
        return float(roc_auc_score(y, probability))
    if name.startswith("precision_at_"):
        k = min(int(name.rsplit("_", 1)[1]), len(y))
        return float(y[np.argsort(-probability, kind="stable")[:k]].mean())
    raise KeyError(name)


pred = pd.read_csv(SOURCE)
metric_names = ["pr_auc", "roc_auc", "precision_at_20", "precision_at_50", "precision_at_100"]
ci_rows = []
diff_rows = []

for split_name, split_frame in pred.groupby("split", sort=True):
    wide = split_frame.pivot(index="peptide", columns="model", values="probability")
    labels = split_frame.drop_duplicates("peptide").set_index("peptide")["label_hiconf"].loc[wide.index]
    y = labels.to_numpy(dtype=int)
    positive = np.flatnonzero(y == 1)
    negative = np.flatnonzero(y == 0)
    rng = np.random.default_rng(SEED)
    samples = [
        np.concatenate(
            [rng.choice(positive, len(positive), replace=True), rng.choice(negative, len(negative), replace=True)]
        )
        for _ in range(N_BOOT)
    ]

    for model_name in wide.columns:
        p = wide[model_name].to_numpy(dtype=float)
        for metric_name in metric_names:
            point = metric(y, p, metric_name)
            values = np.array([metric(y[idx], p[idx], metric_name) for idx in samples])
            ci_rows.append(
                {
                    "split": split_name,
                    "model": model_name,
                    "metric": metric_name,
                    "point_estimate": point,
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                    "bootstrap_design": "stratified by binary outcome",
                }
            )

    reference = "random_forest"
    p_reference = wide[reference].to_numpy(dtype=float)
    for comparator in sorted(set(wide.columns) - {reference}):
        p_comparator = wide[comparator].to_numpy(dtype=float)
        for metric_name in metric_names:
            point = metric(y, p_reference, metric_name) - metric(y, p_comparator, metric_name)
            values = np.array(
                [metric(y[idx], p_reference[idx], metric_name) - metric(y[idx], p_comparator[idx], metric_name) for idx in samples]
            )
            diff_rows.append(
                {
                    "split": split_name,
                    "model_a": reference,
                    "model_b": comparator,
                    "metric": metric_name,
                    "difference_a_minus_b": point,
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                }
            )

pd.DataFrame(ci_rows).to_csv(OUT / "bootstrap_confidence_intervals.csv", index=False)
pd.DataFrame(diff_rows).to_csv(OUT / "paired_bootstrap_differences.csv", index=False)
print(pd.DataFrame(ci_rows).query("metric == 'pr_auc'").to_string(index=False))
