from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


SEED = 42
N_BOOT = 2000
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "05_results" / "validated" / "corrected_per_design_tuning"
pred = pd.read_csv(OUT / "test_predictions.csv")


def score(y: np.ndarray, probability: np.ndarray, metric: str) -> float:
    if metric == "pr_auc":
        return float(average_precision_score(y, probability))
    if metric == "roc_auc":
        return float(roc_auc_score(y, probability))
    k = min(20, len(y))
    return float(y[np.argsort(-probability, kind="stable")[:k]].mean())


ci_rows, difference_rows, integrity_rows = [], [], []
for split_name, frame in pred.groupby("split", sort=True):
    wide = frame.pivot(index="peptide", columns="model", values="probability").sort_index()
    metadata = (
        frame.drop_duplicates("peptide")
        .set_index("peptide")
        .loc[wide.index, ["label_hiconf", "edit_distance_2_cluster"]]
    )
    y = metadata["label_hiconf"].to_numpy(dtype=int)
    unit_values = (
        metadata["edit_distance_2_cluster"].astype(str).to_numpy()
        if split_name == "cluster"
        else wide.index.astype(str).to_numpy()
    )
    unit_name = "edit_distance_2_cluster" if split_name == "cluster" else "peptide"
    unique_units = np.unique(unit_values)
    unit_indices = {unit: np.flatnonzero(unit_values == unit) for unit in unique_units}
    rng = np.random.default_rng(SEED)
    samples = []
    rejected = 0
    while len(samples) < N_BOOT:
        sampled_units = rng.choice(unique_units, size=len(unique_units), replace=True)
        indices = np.concatenate([unit_indices[unit] for unit in sampled_units])
        if np.unique(y[indices]).size < 2:
            rejected += 1
            continue
        samples.append(indices)
    integrity_rows.append(
        {
            "split": split_name, "resampling_unit": unit_name,
            "unique_units": len(unique_units), "replicates": N_BOOT,
            "rejected_single_class_replicates": rejected,
        }
    )
    for model_name in wide.columns:
        probability = wide[model_name].to_numpy(dtype=float)
        for metric in ("pr_auc", "roc_auc", "precision_at_20"):
            values = np.array([score(y[idx], probability[idx], metric) for idx in samples])
            ci_rows.append(
                {
                    "split": split_name, "model": model_name, "metric": metric,
                    "point_estimate": score(y, probability, metric),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                    "resampling_unit": unit_name,
                }
            )
    reference = wide["random_forest"].to_numpy(dtype=float)
    for comparator_name in sorted(set(wide.columns) - {"random_forest"}):
        comparator = wide[comparator_name].to_numpy(dtype=float)
        for metric in ("pr_auc", "roc_auc", "precision_at_20"):
            values = np.array(
                [score(y[idx], reference[idx], metric) - score(y[idx], comparator[idx], metric) for idx in samples]
            )
            difference_rows.append(
                {
                    "split": split_name, "model_a": "random_forest", "model_b": comparator_name,
                    "metric": metric,
                    "difference_a_minus_b": score(y, reference, metric) - score(y, comparator, metric),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                    "resampling_unit": unit_name,
                }
            )

pd.DataFrame(ci_rows).to_csv(OUT / "bootstrap_confidence_intervals.csv", index=False)
pd.DataFrame(difference_rows).to_csv(OUT / "paired_differences_vs_random_forest.csv", index=False)
pd.DataFrame(integrity_rows).to_csv(OUT / "bootstrap_integrity.csv", index=False)
print(pd.DataFrame(ci_rows).query("metric == 'pr_auc'").to_string(index=False))
print("\nRandom Forest minus stack")
print(pd.DataFrame(difference_rows).query("metric == 'pr_auc' and model_b == 'stacked_ensemble'").to_string(index=False))
