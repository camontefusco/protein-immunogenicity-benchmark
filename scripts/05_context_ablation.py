from __future__ import annotations

import json
import os
import platform
import time
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder, StandardScaler


SEED = 42
SOURCE = Path(os.environ.get("PIB_ASSAY_DATA", "data/curated/assay_level_with_year.csv"))
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "generated" / "context_ablation"
MODEL_DIR = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

BIOLOGICAL = ["virus_species", "virus_strain", "antigen", "protein"]
ASSAY = ["assay_method", "readout", "host"]
SAFE_CONTEXT = ["peptide", *BIOLOGICAL, *ASSAY]
FORBIDDEN = ["label", "outcome_raw", "n_positive", "response_freq_pct", "pos_rate"]


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def deletion_keys(sequence: str):
    yield sequence
    for count in (1, 2):
        for positions in combinations(range(len(sequence)), count):
            blocked = set(positions)
            yield "".join(char for idx, char in enumerate(sequence) if idx not in blocked)


def distance_at_most_two(left: str, right: str) -> bool:
    if abs(len(left) - len(right)) > 2:
        return False
    previous = list(range(len(right) + 1))
    for i, char_left in enumerate(left, 1):
        current = [i] + [3] * len(right)
        low, high = max(1, i - 2), min(len(right), i + 2)
        row_min = 3
        for j in range(low, high + 1):
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + (char_left != right[j - 1]),
            )
            row_min = min(row_min, current[j])
        if row_min > 2:
            return False
        previous = current
    return previous[len(right)] <= 2


def make_clusters(sequences: list[str]) -> np.ndarray:
    union = UnionFind(len(sequences))
    deletion_index: dict[str, list[int]] = defaultdict(list)
    candidates: set[tuple[int, int]] = set()
    for idx, sequence in enumerate(sequences):
        for key in deletion_keys(sequence):
            for other in deletion_index[key]:
                candidates.add((other, idx))
            deletion_index[key].append(idx)
    for left, right in sorted(candidates):
        if distance_at_most_two(sequences[left], sequences[right]):
            union.union(left, right)
    roots = [union.find(idx) for idx in range(len(sequences))]
    numbering = {root: idx for idx, root in enumerate(sorted(set(roots)))}
    return np.array([numbering[root] for root in roots], dtype=int)


def choose_group_split(frame: pd.DataFrame, group_column: str) -> tuple[np.ndarray, np.ndarray, int]:
    target = frame["y"].mean()
    best = None
    for split_seed in range(512):
        splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=split_seed)
        train, test = next(splitter.split(frame, frame["y"], frame[group_column]))
        score = abs(len(test) / len(frame) - 0.20) + abs(frame.iloc[test]["y"].mean() - target)
        candidate = (score, split_seed, train, test)
        if best is None or candidate[:2] < best[:2]:
            best = candidate
    return best[2], best[3], best[1]


def ece(y: np.ndarray, probability: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    membership = np.clip(np.digitize(probability, edges[1:-1], right=True), 0, bins - 1)
    return float(
        sum(
            (membership == idx).mean()
            * abs(y[membership == idx].mean() - probability[membership == idx].mean())
            for idx in range(bins)
            if (membership == idx).any()
        )
    )


def topk(y: np.ndarray, probability: np.ndarray, k: int) -> tuple[float, float]:
    k = min(k, len(y))
    precision = float(y[np.argsort(-probability, kind="stable")[:k]].mean())
    return precision, precision / float(y.mean())


def build_pipeline(include_sequence: bool, categorical: list[str], numeric: list[str]) -> Pipeline:
    transformers = []
    if include_sequence:
        transformers.append(
            ("sequence", TfidfVectorizer(analyzer="char", ngram_range=(2, 4), lowercase=False, min_df=2, sublinear_tf=True, dtype=np.float32), "peptide")
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
            ("classifier", LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=SEED, solver="liblinear")),
        ]
    )


raw = pd.read_csv(SOURCE)
raw["peptide"] = raw["peptide"].astype(str).str.upper().str.strip()
alphabet = set("ACDEFGHIKLMNPQRSTVWY")
raw = raw[raw["peptide"].map(lambda value: set(value) <= alphabet)].copy()
for column in BIOLOGICAL + ASSAY:
    raw[column] = raw[column].fillna("__MISSING__").astype(str)

# Remove only rows identical across every exported field. Counts after this step are
# retained solely for the explicitly retrospective evidence-count ablations.
raw_deduplicated = raw.drop_duplicates().copy()
grouped = raw_deduplicated.groupby(SAFE_CONTEXT, as_index=False, dropna=False)
context = grouped.agg(
    n_assays_context=("label", "size"),
    positive_fraction=("label", "mean"),
    n_unique_labels=("label", "nunique"),
    earliest_pub_year=("earliest_pub_year", "min"),
)
context["y"] = (context["positive_fraction"] > 0.5).astype(int)
context["majority_tie"] = context["positive_fraction"].eq(0.5)
context["length"] = context["peptide"].str.len().astype(int)
context["log_evidence_count"] = context["n_assays_context"].astype(float)

unique_peptides = sorted(context["peptide"].unique())
clusters = make_clusters(unique_peptides)
cluster_map = dict(zip(unique_peptides, clusters))
context["edit_distance_2_cluster"] = context["peptide"].map(cluster_map).astype(int)

exact_train, exact_test, exact_seed = choose_group_split(context, "peptide")
cluster_train, cluster_test, cluster_seed = choose_group_split(context, "edit_distance_2_cluster")
year = pd.to_numeric(context["earliest_pub_year"], errors="coerce")
year_valid = context.index[year.notna()].to_numpy()
candidate_years = sorted(year.loc[year_valid].astype(int).unique())
temporal_cutoff = min(candidate_years, key=lambda value: abs((year.loc[year_valid] >= value).mean() - 0.20))
temporal_train = context.index[year.notna() & (year < temporal_cutoff)].to_numpy()
temporal_test = context.index[year.notna() & (year >= temporal_cutoff)].to_numpy()

splits = {
    "exact": (exact_train, exact_test),
    "cluster": (cluster_train, cluster_test),
    "temporal": (temporal_train, temporal_test),
}
context["exact_split"] = "excluded"
context.loc[exact_train, "exact_split"] = "train"
context.loc[exact_test, "exact_split"] = "test"
context["cluster_split"] = "excluded"
context.loc[cluster_train, "cluster_split"] = "train"
context.loc[cluster_test, "cluster_split"] = "test"
context["temporal_split"] = "excluded"
context.loc[temporal_train, "temporal_split"] = "train"
context.loc[temporal_test, "temporal_split"] = "test"
context.to_csv(OUT / "context_dataset_and_splits.csv", index=False)

ABLATIONS = {
    "sequence_only": (True, [], []),
    "biological_only": (False, BIOLOGICAL, ["length"]),
    "assay_only": (False, ASSAY, []),
    "evidence_only": (False, [], ["log_evidence_count"]),
    "sequence_plus_biological": (True, BIOLOGICAL, ["length"]),
    "sequence_plus_assay": (True, ASSAY, []),
    "sequence_plus_all_no_evidence": (True, BIOLOGICAL + ASSAY, ["length"]),
    "sequence_plus_all_with_evidence": (True, BIOLOGICAL + ASSAY, ["length", "log_evidence_count"]),
    "all_context_no_sequence_no_evidence": (False, BIOLOGICAL + ASSAY, ["length"]),
}

metrics, predictions, efficiency = [], [], []
for split_name, (train_indices, test_indices) in splits.items():
    train, test = context.loc[train_indices].copy(), context.loc[test_indices].copy()
    y_train, y_test = train["y"].to_numpy(dtype=int), test["y"].to_numpy(dtype=int)
    for model_name, (with_sequence, categorical, numeric) in ABLATIONS.items():
        model = build_pipeline(with_sequence, categorical, numeric)
        feature_columns = (["peptide"] if with_sequence else []) + categorical + numeric
        start = time.perf_counter()
        model.fit(train[feature_columns], y_train)
        train_seconds = time.perf_counter() - start
        start = time.perf_counter()
        probability = model.predict_proba(test[feature_columns])[:, 1]
        inference_seconds = time.perf_counter() - start
        model_path = MODEL_DIR / f"{split_name}__{model_name}.joblib"
        joblib.dump(model, model_path, compress=3)
        row = {
            "split": split_name,
            "model": model_name,
            "n_train_rows": len(train),
            "n_test_rows": len(test),
            "n_train_peptides": train["peptide"].nunique(),
            "n_test_peptides": test["peptide"].nunique(),
            "test_positive_rate": float(y_test.mean()),
            "roc_auc": float(roc_auc_score(y_test, probability)),
            "pr_auc": float(average_precision_score(y_test, probability)),
            "brier_score": float(brier_score_loss(y_test, probability)),
            "ece_10_equal_width": ece(y_test, probability),
        }
        for k in (20, 50, 100, 200, 500):
            row[f"precision_at_{k}"], row[f"enrichment_at_{k}"] = topk(y_test, probability, k)
        metrics.append(row)
        result = test[["peptide", "y", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
        result.insert(0, "split", split_name)
        result.insert(1, "model", model_name)
        result["row_id"] = result.index
        result["probability"] = probability
        predictions.append(result)
        efficiency.append(
            {
                "split": split_name,
                "model": model_name,
                "train_seconds": train_seconds,
                "inference_seconds": inference_seconds,
                "inference_microseconds_per_row": inference_seconds * 1e6 / len(test),
                "serialized_model_bytes": model_path.stat().st_size,
            }
        )
        print(split_name, model_name, f"PR-AUC={row['pr_auc']:.4f}")

pd.DataFrame(metrics).to_csv(OUT / "metrics.csv", index=False)
pd.concat(predictions, ignore_index=True).to_csv(OUT / "test_predictions.csv", index=False)
pd.DataFrame(efficiency).to_csv(OUT / "efficiency.csv", index=False)

summary = {
    "raw_rows": int(len(raw)),
    "deduplicated_rows": int(len(raw_deduplicated)),
    "context_rows": int(len(context)),
    "context_unique_peptides": int(context["peptide"].nunique()),
    "context_positive_rate": float(context["y"].mean()),
    "mixed_label_context_rows": int((context["n_unique_labels"] > 1).sum()),
    "majority_tie_rows_labeled_negative": int(context["majority_tie"].sum()),
    "cluster_count": int(len(set(clusters))),
    "largest_cluster": int(pd.Series(clusters).value_counts().max()),
    "exact_split_seed_selected_by_size_and_prevalence_only": int(exact_seed),
    "cluster_split_seed_selected_by_size_and_prevalence_only": int(cluster_seed),
    "cluster_overlap": int(
        len(set(context.loc[cluster_train, "edit_distance_2_cluster"]) & set(context.loc[cluster_test, "edit_distance_2_cluster"]))
    ),
    "temporal_cutoff_test_inclusive": int(temporal_cutoff),
    "forbidden_predictors": FORBIDDEN,
    "evidence_feature_warning": "n_assays_context is retrospective and appears only in explicitly labeled diagnostic ablations.",
    "software_platform": platform.platform(),
}
with (OUT / "summary.json").open("w", encoding="utf-8") as stream:
    json.dump(summary, stream, indent=2, sort_keys=True)
print(json.dumps(summary, indent=2, sort_keys=True))
