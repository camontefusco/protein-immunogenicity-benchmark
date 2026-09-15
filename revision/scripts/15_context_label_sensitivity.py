from __future__ import annotations

import json
import time
from pathlib import Path

import joblib
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
ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "02_data" / "curated" / "context_dataset_and_splits.csv"
OUT = ROOT / "05_results" / "validated" / "context_label_sensitivity"
MODEL_DIR = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

BIOLOGICAL = ["virus_species", "virus_strain", "antigen", "protein"]
ASSAY = ["assay_method", "readout", "host"]
MODEL_SPECS = {
    "sequence_only": (True, [], []),
    "biological_only": (False, BIOLOGICAL, ["length"]),
    "assay_only": (False, ASSAY, []),
    "sequence_plus_biological": (True, BIOLOGICAL, ["length"]),
    "sequence_plus_assay": (True, ASSAY, []),
    "sequence_plus_all_no_evidence": (True, BIOLOGICAL + ASSAY, ["length"]),
    "all_context_no_sequence_no_evidence": (False, BIOLOGICAL + ASSAY, ["length"]),
}


def build_pipeline(include_sequence: bool, categorical: list[str], numeric: list[str]) -> Pipeline:
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
    if categorical:
        transformers.append(("categorical", OneHotEncoder(handle_unknown="ignore", min_frequency=2), categorical))
    if numeric:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    [
                        ("log1p", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric,
            )
        )
    return Pipeline(
        [
            ("features", ColumnTransformer(transformers, remainder="drop")),
            (
                "classifier",
                LogisticRegression(
                    C=1.0, class_weight="balanced", max_iter=3000,
                    random_state=SEED, solver="liblinear",
                ),
            ),
        ]
    )


def pr_auc(y: np.ndarray, probability: np.ndarray) -> float:
    return float(average_precision_score(y, probability))


context = pd.read_csv(SOURCE)
constructions = {
    "majority_all": np.ones(len(context), dtype=bool),
    "majority_excluding_ties": ~context["majority_tie"].astype(bool).to_numpy(),
    "consistent_labels_only": context["n_unique_labels"].eq(1).to_numpy(),
}

metrics_rows: list[dict] = []
predictions: list[pd.DataFrame] = []
construction_rows: list[dict] = []
for construction_name, inclusion_mask in constructions.items():
    subset = context.loc[inclusion_mask].copy()
    construction_rows.append(
        {
            "construction": construction_name,
            "rows": len(subset),
            "excluded_rows": len(context) - len(subset),
            "unique_peptides": subset["peptide"].nunique(),
            "positive_rate": float(subset["y"].mean()),
            "ties_remaining": int(subset["majority_tie"].sum()),
            "mixed_label_rows_remaining": int((subset["n_unique_labels"] > 1).sum()),
        }
    )
    for split_name in ("exact", "cluster", "temporal"):
        split_column = f"{split_name}_split"
        train = subset[subset[split_column] == "train"].copy()
        test = subset[subset[split_column] == "test"].copy()
        y_train = train["y"].to_numpy(dtype=int)
        y_test = test["y"].to_numpy(dtype=int)
        for model_name, (with_sequence, categorical, numeric) in MODEL_SPECS.items():
            feature_columns = (["peptide"] if with_sequence else []) + categorical + numeric
            model = build_pipeline(with_sequence, categorical, numeric)
            start = time.perf_counter()
            model.fit(train[feature_columns], y_train)
            probability = model.predict_proba(test[feature_columns])[:, 1]
            elapsed = time.perf_counter() - start
            model_path = MODEL_DIR / f"{construction_name}__{split_name}__{model_name}.joblib"
            joblib.dump(model, model_path, compress=3)
            metrics_rows.append(
                {
                    "construction": construction_name,
                    "split": split_name,
                    "model": model_name,
                    "n_train_rows": len(train),
                    "n_test_rows": len(test),
                    "n_train_peptides": train["peptide"].nunique(),
                    "n_test_peptides": test["peptide"].nunique(),
                    "test_positive_rate": float(y_test.mean()),
                    "roc_auc": float(roc_auc_score(y_test, probability)),
                    "pr_auc": pr_auc(y_test, probability),
                    "fit_and_predict_seconds": elapsed,
                }
            )
            result = test[["peptide", "y", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
            result.insert(0, "construction", construction_name)
            result.insert(1, "split", split_name)
            result.insert(2, "model", model_name)
            result["row_id"] = result.index
            result["probability"] = probability
            predictions.append(result)
            print(construction_name, split_name, model_name, f"PR-AUC={metrics_rows[-1]['pr_auc']:.4f}")

metrics = pd.DataFrame(metrics_rows)
prediction_table = pd.concat(predictions, ignore_index=True)
metrics.to_csv(OUT / "metrics.csv", index=False)
prediction_table.to_csv(OUT / "test_predictions.csv", index=False)
pd.DataFrame(construction_rows).to_csv(OUT / "construction_summary.csv", index=False)

# Paired hierarchical bootstrap for the main context uplift within each target
# construction and validation design.
interval_rows: list[dict] = []
difference_rows: list[dict] = []
integrity_rows: list[dict] = []
for (construction_name, split_name), frame in prediction_table.groupby(["construction", "split"], sort=True):
    wide = frame.pivot(index="row_id", columns="model", values="probability").sort_index()
    metadata = (
        frame.drop_duplicates("row_id")
        .set_index("row_id")
        .loc[wide.index, ["peptide", "y", "edit_distance_2_cluster"]]
    )
    y = metadata["y"].to_numpy(dtype=int)
    unit_column = "edit_distance_2_cluster" if split_name == "cluster" else "peptide"
    units = metadata[unit_column].astype(str).to_numpy()
    unique_units = np.unique(units)
    indices_by_unit = {unit: np.flatnonzero(units == unit) for unit in unique_units}
    rng = np.random.default_rng(SEED)
    bootstrap_indices = []
    rejected = 0
    while len(bootstrap_indices) < N_BOOT:
        sampled_units = rng.choice(unique_units, size=len(unique_units), replace=True)
        indices = np.concatenate([indices_by_unit[unit] for unit in sampled_units])
        if np.unique(y[indices]).size < 2:
            rejected += 1
            continue
        bootstrap_indices.append(indices)
    integrity_rows.append(
        {
            "construction": construction_name,
            "split": split_name,
            "resampling_unit": unit_column,
            "unique_units": len(unique_units),
            "replicates": N_BOOT,
            "rejected_single_class_replicates": rejected,
        }
    )
    for model_name in wide.columns:
        probability = wide[model_name].to_numpy(dtype=float)
        values = np.array([pr_auc(y[idx], probability[idx]) for idx in bootstrap_indices])
        interval_rows.append(
            {
                "construction": construction_name,
                "split": split_name,
                "model": model_name,
                "metric": "pr_auc",
                "point_estimate": pr_auc(y, probability),
                "ci_low": float(np.quantile(values, 0.025)),
                "ci_high": float(np.quantile(values, 0.975)),
                "resampling_unit": unit_column,
                "replicates": N_BOOT,
            }
        )
    sequence_probability = wide["sequence_only"].to_numpy(dtype=float)
    for model_name in sorted(set(wide.columns) - {"sequence_only"}):
        probability = wide[model_name].to_numpy(dtype=float)
        values = np.array(
            [pr_auc(y[idx], probability[idx]) - pr_auc(y[idx], sequence_probability[idx]) for idx in bootstrap_indices]
        )
        difference_rows.append(
            {
                "construction": construction_name,
                "split": split_name,
                "model_a": model_name,
                "model_b": "sequence_only",
                "metric": "pr_auc",
                "difference_a_minus_b": pr_auc(y, probability) - pr_auc(y, sequence_probability),
                "ci_low": float(np.quantile(values, 0.025)),
                "ci_high": float(np.quantile(values, 0.975)),
                "resampling_unit": unit_column,
                "replicates": N_BOOT,
            }
        )

pd.DataFrame(interval_rows).to_csv(OUT / "bootstrap_confidence_intervals.csv", index=False)
pd.DataFrame(difference_rows).to_csv(OUT / "paired_differences_vs_sequence.csv", index=False)
pd.DataFrame(integrity_rows).to_csv(OUT / "bootstrap_integrity.csv", index=False)
with (OUT / "provenance.json").open("w", encoding="utf-8") as stream:
    json.dump(
        {
            "source_context_dataset": str(SOURCE),
            "target_constructions": {
                "majority_all": "majority rule; exact 0.5 ties assigned negative",
                "majority_excluding_ties": "majority rule after excluding exact 0.5 ties",
                "consistent_labels_only": "retain only context rows with one unique assay label",
            },
            "splits": "unchanged frozen exact, cluster and temporal assignments",
            "bootstrap": "paired hierarchical bootstrap by peptide or edit-distance-2 component",
            "seed": SEED,
            "replicates": N_BOOT,
        },
        stream,
        indent=2,
        sort_keys=True,
    )

print("\nConstruction summary")
print(pd.DataFrame(construction_rows).to_string(index=False))
print("\nContext uplift versus sequence")
print(
    pd.DataFrame(difference_rows)
    .query("model_a == 'sequence_plus_all_no_evidence'")
    .to_string(index=False)
)
