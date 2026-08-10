from __future__ import annotations

import json
import platform
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer


SEED = 42
MODEL_ID = "facebook/esm2_t6_8M_UR50D"
BATCH_SIZE = 128
ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "artifacts" / "splits" / "peptide_split_assignments.csv"
OUT = ROOT / "artifacts" / "esm2_frozen"
MODEL_OUT = OUT / "models"
OUT.mkdir(parents=True, exist_ok=True)
MODEL_OUT.mkdir(parents=True, exist_ok=True)
EMBEDDINGS = OUT / "esm2_t6_mean_embeddings.npz"

torch.manual_seed(SEED)
np.random.seed(SEED)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

data = pd.read_csv(SPLITS)
eligible = data[data["valid_sequence"]].copy().reset_index(drop=True)
sequences = eligible["peptide"].tolist()

embedding_seconds = 0.0
if EMBEDDINGS.exists():
    cached = np.load(EMBEDDINGS, allow_pickle=False)
    if not np.array_equal(cached["peptides"], np.asarray(sequences)):
        raise ValueError("Cached embedding peptide order does not match the frozen split file")
    matrix = cached["embeddings"]
else:
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModel.from_pretrained(MODEL_ID).to(device)
    model.eval()
    batches = []
    start = time.perf_counter()
    with torch.inference_mode():
        for begin in range(0, len(sequences), BATCH_SIZE):
            batch = sequences[begin : begin + BATCH_SIZE]
            tokens = tokenizer(batch, return_tensors="pt", padding=True, truncation=True)
            tokens = {name: value.to(device) for name, value in tokens.items()}
            hidden = model(**tokens).last_hidden_state
            mask = tokens["attention_mask"].bool()
            # Exclude BOS and EOS special tokens from residue mean pooling.
            residue_mask = mask.clone()
            residue_mask[:, 0] = False
            lengths = mask.sum(dim=1)
            residue_mask[torch.arange(len(batch), device=device), lengths - 1] = False
            pooled = (hidden * residue_mask.unsqueeze(-1)).sum(dim=1) / residue_mask.sum(dim=1, keepdim=True)
            batches.append(pooled.cpu().numpy().astype(np.float32))
            print(f"embedded {min(begin + BATCH_SIZE, len(sequences))}/{len(sequences)}")
    embedding_seconds = time.perf_counter() - start
    matrix = np.concatenate(batches, axis=0)
    np.savez_compressed(EMBEDDINGS, peptides=np.asarray(sequences), embeddings=matrix)

metrics, predictions, efficiency = [], [], []
for split_name, split_column in {
    "exact": "exact_split", "cluster": "cluster_split", "temporal": "temporal_split"
}.items():
    train_mask = eligible[split_column].eq("train").to_numpy()
    test_mask = eligible[split_column].eq("test").to_numpy()
    x_train, x_test = matrix[train_mask], matrix[test_mask]
    y_train = eligible.loc[train_mask, "label_hiconf"].to_numpy(dtype=int)
    y_test = eligible.loc[test_mask, "label_hiconf"].to_numpy(dtype=int)
    classifier = Pipeline(
        [
            ("scale", StandardScaler()),
            ("classifier", LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000, random_state=SEED, solver="liblinear")),
        ]
    )
    start = time.perf_counter()
    classifier.fit(x_train, y_train)
    train_seconds = time.perf_counter() - start
    start = time.perf_counter()
    probability = classifier.predict_proba(x_test)[:, 1]
    inference_seconds = time.perf_counter() - start
    model_path = MODEL_OUT / f"{split_name}__esm2_frozen_logistic.joblib"
    joblib.dump(classifier, model_path, compress=3)
    metrics.append(
        {
            "split": split_name,
            "model": "ESM2-8M frozen mean embedding + logistic regression",
            "n_train": int(train_mask.sum()),
            "n_test": int(test_mask.sum()),
            "positive_rate": float(y_test.mean()),
            "roc_auc": float(roc_auc_score(y_test, probability)),
            "pr_auc": float(average_precision_score(y_test, probability)),
            "brier_score": float(brier_score_loss(y_test, probability)),
        }
    )
    result = eligible.loc[test_mask, ["peptide", "label_hiconf", "edit_distance_2_cluster", "earliest_pub_year"]].copy()
    result.insert(0, "split", split_name)
    result["probability"] = probability
    predictions.append(result)
    efficiency.append(
        {
            "split": split_name,
            "embedding_seconds_all_peptides_this_run": embedding_seconds,
            "classifier_train_seconds": train_seconds,
            "classifier_inference_seconds": inference_seconds,
            "serialized_classifier_bytes": model_path.stat().st_size,
            "embedding_dimensions": int(matrix.shape[1]),
        }
    )
    print(split_name, metrics[-1])

pd.DataFrame(metrics).to_csv(OUT / "metrics.csv", index=False)
pd.concat(predictions, ignore_index=True).to_csv(OUT / "test_predictions.csv", index=False)
pd.DataFrame(efficiency).to_csv(OUT / "efficiency.csv", index=False)
with (OUT / "provenance.json").open("w", encoding="utf-8") as stream:
    json.dump(
        {
            "model_id": MODEL_ID,
            "representation": "mean of final hidden-state residue embeddings, excluding BOS/EOS/padding",
            "fine_tuned": False,
            "device": str(device),
            "torch": torch.__version__,
            "platform": platform.platform(),
            "seed": SEED,
            "batch_size": BATCH_SIZE,
        },
        stream,
        indent=2,
        sort_keys=True,
    )
