from __future__ import annotations

import json
import platform
import time
from itertools import product
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
import xgboost
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold


SEED = 42
ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "02_data" / "splits" / "peptide_split_assignments_v1.csv"
OUT = ROOT / "05_results" / "validated" / "corrected_per_design_tuning"
MODEL_DIR = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)


def vectorizer() -> TfidfVectorizer:
    return TfidfVectorizer(
        analyzer="char", ngram_range=(2, 4), lowercase=False,
        min_df=2, sublinear_tf=True, dtype=np.float32,
    )


LOGISTIC_GRID = [
    {"alpha": alpha, "l1_ratio": ratio}
    for alpha, ratio in product([1e-5, 1e-4, 1e-3], [0.15, 0.50, 0.85])
]
RF_GRID = [
    {"max_features": features, "min_samples_leaf": leaf}
    for features, leaf in product(["sqrt", "log2"], [1, 3])
]
XGB_GRID = [
    {"max_depth": depth, "learning_rate": rate, "min_child_weight": weight}
    for depth, rate, weight in product([3, 5], [0.03, 0.10], [1])
]


def estimator(name: str, params: dict, positive_weight: float):
    if name == "elastic_logistic":
        return SGDClassifier(
            loss="log_loss", penalty="elasticnet", alpha=params["alpha"],
            l1_ratio=params["l1_ratio"], class_weight="balanced",
            max_iter=3000, tol=1e-4, random_state=SEED,
        )
    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=500, max_features=params["max_features"],
            min_samples_leaf=params["min_samples_leaf"],
            class_weight="balanced_subsample", random_state=SEED, n_jobs=-1,
        )
    if name == "xgboost":
        return xgboost.XGBClassifier(
            n_estimators=500, max_depth=params["max_depth"],
            learning_rate=params["learning_rate"], min_child_weight=params["min_child_weight"],
            subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0,
            objective="binary:logistic", eval_metric="aucpr", tree_method="hist",
            scale_pos_weight=positive_weight, random_state=SEED, n_jobs=8,
        )
    raise KeyError(name)


def topk(y: np.ndarray, probability: np.ndarray, k: int) -> tuple[float, float]:
    k = min(k, len(y))
    precision = float(y[np.argsort(-probability, kind="stable")[:k]].mean())
    return precision, precision / float(y.mean())


data = pd.read_csv(SOURCE)
grids = {"elastic_logistic": LOGISTIC_GRID, "random_forest": RF_GRID, "xgboost": XGB_GRID}
metrics_rows, prediction_frames, efficiency_rows, search_rows, integrity_rows = [], [], [], [], []
selected_parameters = {}
split_columns = {"exact": "exact_split", "cluster": "cluster_split", "temporal": "temporal_split"}
for split_name, split_column in split_columns.items():
    train = data[data[split_column] == "train"].reset_index(drop=True)
    test = data[data[split_column] == "test"].reset_index(drop=True)
    y_train = train["label_hiconf"].to_numpy(dtype=int)
    y_test = test["label_hiconf"].to_numpy(dtype=int)
    inner = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED)
    inner_folds = list(inner.split(train["peptide"], y_train, train["edit_distance_2_cluster"]))

    # Hyperparameter selection is repeated independently inside this design's
    # development partition. The associated test set is never visible here.
    fold_matrices = []
    for inner_train, inner_valid in inner_folds:
        vec = vectorizer()
        x_inner_train = vec.fit_transform(train.loc[inner_train, "peptide"])
        x_inner_valid = vec.transform(train.loc[inner_valid, "peptide"])
        fold_matrices.append((inner_train, inner_valid, x_inner_train, x_inner_valid))

    best_params = {}
    for model_name, grid in grids.items():
        design_candidates = []
        for candidate_id, params in enumerate(grid):
            fold_scores = []
            start = time.perf_counter()
            for inner_train, inner_valid, x_inner_train, x_inner_valid in fold_matrices:
                inner_y = y_train[inner_train]
                valid_y = y_train[inner_valid]
                weight = float((inner_y == 0).sum() / (inner_y == 1).sum())
                candidate_model = estimator(model_name, params, weight)
                candidate_model.fit(x_inner_train, inner_y)
                fold_scores.append(
                    float(average_precision_score(valid_y, candidate_model.predict_proba(x_inner_valid)[:, 1]))
                )
            candidate_row = {
                "split": split_name,
                "model": model_name,
                "candidate_id": candidate_id,
                "parameters_json": json.dumps(params, sort_keys=True),
                "mean_group_cv_pr_auc": float(np.mean(fold_scores)),
                "std_group_cv_pr_auc": float(np.std(fold_scores, ddof=1)),
                "fold_pr_auc_json": json.dumps(fold_scores),
                "search_seconds": time.perf_counter() - start,
            }
            search_rows.append(candidate_row)
            design_candidates.append(candidate_row)
            print("search", split_name, model_name, params, np.mean(fold_scores))
        best_row = sorted(
            design_candidates,
            key=lambda row: (-row["mean_group_cv_pr_auc"], row["candidate_id"]),
        )[0]
        best_params[model_name] = json.loads(best_row["parameters_json"])
    selected_parameters[split_name] = best_params
    print("selected", split_name, best_params)

    integrity_rows.append(
        {
            "split": split_name,
            "development_test_peptide_overlap": len(set(train["peptide"]) & set(test["peptide"])),
            "development_test_cluster_overlap": len(
                set(train["edit_distance_2_cluster"]) & set(test["edit_distance_2_cluster"])
            ),
            "selection_rows": len(train),
            "test_rows": len(test),
        }
    )
    oof = np.zeros((len(train), len(grids)), dtype=float)
    test_probabilities = {}
    full_models = {}
    total_start = time.perf_counter()

    for model_position, model_name in enumerate(grids):
        for inner_train, inner_valid in inner_folds:
            vec = vectorizer()
            x_inner_train = vec.fit_transform(train.loc[inner_train, "peptide"])
            x_inner_valid = vec.transform(train.loc[inner_valid, "peptide"])
            inner_y = y_train[inner_train]
            weight = float((inner_y == 0).sum() / (inner_y == 1).sum())
            fold_model = estimator(model_name, best_params[model_name], weight)
            fold_model.fit(x_inner_train, inner_y)
            oof[inner_valid, model_position] = fold_model.predict_proba(x_inner_valid)[:, 1]

        vec = vectorizer()
        x_train = vec.fit_transform(train["peptide"])
        x_test = vec.transform(test["peptide"])
        weight = float((y_train == 0).sum() / (y_train == 1).sum())
        fitted = estimator(model_name, best_params[model_name], weight)
        start = time.perf_counter()
        fitted.fit(x_train, y_train)
        fit_seconds = time.perf_counter() - start
        start = time.perf_counter()
        probability = fitted.predict_proba(x_test)[:, 1]
        inference_seconds = time.perf_counter() - start
        test_probabilities[model_name] = probability
        full_models[model_name] = {"vectorizer": vec, "estimator": fitted}
        model_path = MODEL_DIR / f"{split_name}__{model_name}.joblib"
        joblib.dump(full_models[model_name], model_path, compress=3)
        efficiency_rows.append(
            {
                "split": split_name, "model": model_name,
                "full_fit_seconds": fit_seconds, "test_inference_seconds": inference_seconds,
                "serialized_model_bytes": model_path.stat().st_size,
            }
        )

    meta = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000, random_state=SEED)
    meta.fit(oof, y_train)
    stacked_matrix = np.column_stack([test_probabilities[name] for name in grids])
    test_probabilities["stacked_ensemble"] = meta.predict_proba(stacked_matrix)[:, 1]
    stack_path = MODEL_DIR / f"{split_name}__stack_meta.joblib"
    joblib.dump({"meta": meta, "base_order": list(grids)}, stack_path, compress=3)
    efficiency_rows.append(
        {
            "split": split_name, "model": "stacked_ensemble",
            "full_fit_seconds": time.perf_counter() - total_start,
            "test_inference_seconds": np.nan,
            "serialized_model_bytes": stack_path.stat().st_size,
        }
    )

    for model_name, probability in test_probabilities.items():
        row = {
            "split": split_name, "model": model_name,
            "n_train": len(train), "n_test": len(test),
            "positive_rate": float(y_test.mean()),
            "roc_auc": float(roc_auc_score(y_test, probability)),
            "pr_auc": float(average_precision_score(y_test, probability)),
            "brier_score": float(brier_score_loss(y_test, probability)),
        }
        for k in (20, 50, 100, 200, 500):
            row[f"precision_at_{k}"], row[f"enrichment_at_{k}"] = topk(y_test, probability, k)
        metrics_rows.append(row)
        frame = test[["peptide", "label_hiconf", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
        frame.insert(0, "split", split_name)
        frame.insert(1, "model", model_name)
        frame["probability"] = probability
        prediction_frames.append(frame)
        print("final", split_name, model_name, row["pr_auc"], row["roc_auc"])

pd.DataFrame(search_rows).to_csv(OUT / "hyperparameter_search.csv", index=False)
with (OUT / "selected_parameters.json").open("w", encoding="utf-8") as stream:
    json.dump(selected_parameters, stream, indent=2, sort_keys=True)
pd.DataFrame(integrity_rows).to_csv(OUT / "split_integrity.csv", index=False)
pd.DataFrame(metrics_rows).to_csv(OUT / "metrics.csv", index=False)
pd.concat(prediction_frames, ignore_index=True).to_csv(OUT / "test_predictions.csv", index=False)
pd.DataFrame(efficiency_rows).to_csv(OUT / "efficiency.csv", index=False)
with (OUT / "provenance.json").open("w", encoding="utf-8") as stream:
    json.dump(
        {
            "selection_dataset": "each validation design's own training partition",
            "selection_metric": "mean 3-fold StratifiedGroupKFold PR-AUC selected independently per design",
            "selection_groups": "edit-distance-2 connected components",
            "random_seed": SEED,
            "stacking": "3-fold group-safe out-of-fold base predictions; logistic meta-learner",
            "python_platform": platform.platform(),
            "scikit_learn": sklearn.__version__,
            "xgboost": xgboost.__version__,
        },
        stream,
        indent=2,
        sort_keys=True,
    )
