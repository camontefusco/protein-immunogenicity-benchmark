from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


SEED = 42
N_BOOT = 2000
ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT.parent / "04_reproducible_analysis" / "artifacts" / "context_ablation"
SOURCE = RECOVERY / "test_predictions.csv"
OUT = ROOT / "05_results" / "validated" / "corrected_context_uncertainty"
OUT.mkdir(parents=True, exist_ok=True)


def metric(y: np.ndarray, probability: np.ndarray, name: str) -> float:
    if name == "pr_auc":
        return float(average_precision_score(y, probability))
    if name == "roc_auc":
        return float(roc_auc_score(y, probability))
    raise KeyError(name)


predictions = pd.read_csv(SOURCE)
ci_rows: list[dict] = []
difference_rows: list[dict] = []
integrity_rows: list[dict] = []

for split_name, split_frame in predictions.groupby("split", sort=True):
    wide = split_frame.pivot(index="row_id", columns="model", values="probability").sort_index()
    metadata = (
        split_frame.drop_duplicates("row_id")
        .set_index("row_id")
        .loc[wide.index, ["peptide", "y", "edit_distance_2_cluster"]]
    )
    y = metadata["y"].to_numpy(dtype=int)
    unit_column = "edit_distance_2_cluster" if split_name == "cluster" else "peptide"
    unit_values = metadata[unit_column].astype(str).to_numpy()
    unique_units = np.unique(unit_values)
    unit_indices = {unit: np.flatnonzero(unit_values == unit) for unit in unique_units}
    rng = np.random.default_rng(SEED)

    bootstrap_indices = []
    rejected = 0
    while len(bootstrap_indices) < N_BOOT:
        sampled_units = rng.choice(unique_units, size=len(unique_units), replace=True)
        indices = np.concatenate([unit_indices[unit] for unit in sampled_units])
        if np.unique(y[indices]).size < 2:
            rejected += 1
            continue
        bootstrap_indices.append(indices)

    integrity_rows.append(
        {
            "split": split_name,
            "resampling_unit": unit_column,
            "test_rows": len(y),
            "unique_resampling_units": len(unique_units),
            "bootstrap_replicates": N_BOOT,
            "rejected_single_class_replicates": rejected,
        }
    )

    for model_name in wide.columns:
        probability = wide[model_name].to_numpy(dtype=float)
        for metric_name in ("pr_auc", "roc_auc"):
            values = np.array([metric(y[idx], probability[idx], metric_name) for idx in bootstrap_indices])
            ci_rows.append(
                {
                    "split": split_name,
                    "model": model_name,
                    "metric": metric_name,
                    "point_estimate": metric(y, probability, metric_name),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                    "resampling_unit": unit_column,
                }
            )

    sequence_probability = wide["sequence_only"].to_numpy(dtype=float)
    for model_name in sorted(set(wide.columns) - {"sequence_only"}):
        probability = wide[model_name].to_numpy(dtype=float)
        for metric_name in ("pr_auc", "roc_auc"):
            values = np.array(
                [
                    metric(y[idx], probability[idx], metric_name)
                    - metric(y[idx], sequence_probability[idx], metric_name)
                    for idx in bootstrap_indices
                ]
            )
            difference_rows.append(
                {
                    "split": split_name,
                    "model_a": model_name,
                    "model_b": "sequence_only",
                    "metric": metric_name,
                    "difference_a_minus_b": metric(y, probability, metric_name)
                    - metric(y, sequence_probability, metric_name),
                    "ci_low": float(np.quantile(values, 0.025)),
                    "ci_high": float(np.quantile(values, 0.975)),
                    "bootstrap_replicates": N_BOOT,
                    "resampling_unit": unit_column,
                }
            )

pd.DataFrame(ci_rows).to_csv(OUT / "bootstrap_confidence_intervals.csv", index=False)
pd.DataFrame(difference_rows).to_csv(OUT / "paired_differences_vs_sequence.csv", index=False)
pd.DataFrame(integrity_rows).to_csv(OUT / "bootstrap_integrity.csv", index=False)
with (OUT / "provenance.json").open("w", encoding="utf-8") as stream:
    json.dump(
        {
            "source_predictions": str(SOURCE),
            "seed": SEED,
            "replicates": N_BOOT,
            "design": "paired nonparametric cluster bootstrap carrying all context rows for each sampled unit",
            "exact_unit": "peptide",
            "cluster_unit": "edit_distance_2_cluster",
            "temporal_unit": "peptide",
        },
        stream,
        indent=2,
        sort_keys=True,
    )

summary = pd.DataFrame(ci_rows).query("metric == 'pr_auc'")
print(summary.to_string(index=False))
print("\nMain differences versus sequence")
print(
    pd.DataFrame(difference_rows)
    .query("metric == 'pr_auc' and model_a in ['sequence_plus_all_no_evidence', 'biological_only', 'all_context_no_sequence_no_evidence']")
    .to_string(index=False)
)
