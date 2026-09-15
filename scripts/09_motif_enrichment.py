from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression


SEED = 42
TOP_PER_DIRECTION = 50
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "splits" / "peptide_split_assignments.csv"
OUT = ROOT / "artifacts" / "motif_analysis"
OUT.mkdir(parents=True, exist_ok=True)

HYDROPHOBIC = set("AILMFWVY")
POSITIVE_CHARGE = set("KRH")
NEGATIVE_CHARGE = set("DE")
AROMATIC = set("FWY")
SMALL_POLAR = set("STNQ")


def fraction(motif: str, residues: set[str]) -> float:
    return sum(character in residues for character in motif) / len(motif)


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.empty_like(ranked, dtype=float)
    running = 1.0
    for index in range(len(ranked) - 1, -1, -1):
        candidate = ranked[index] * len(ranked) / (index + 1)
        running = min(running, candidate)
        adjusted[index] = running
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


data = pd.read_csv(SOURCE)
all_rows = []
summary = []
for split_name, split_column in {
    "exact": "exact_split", "cluster": "cluster_split", "temporal": "temporal_split"
}.items():
    train = data[data[split_column] == "train"]
    test = data[data[split_column] == "test"]
    vectorizer = TfidfVectorizer(
        analyzer="char", ngram_range=(2, 4), lowercase=False,
        min_df=2, sublinear_tf=True, dtype=np.float64,
    )
    matrix = vectorizer.fit_transform(train["peptide"])
    model = LogisticRegression(
        C=1.0, penalty="l1", solver="liblinear", class_weight="balanced",
        max_iter=3000, random_state=SEED,
    )
    model.fit(matrix, train["label_hiconf"])
    motifs = vectorizer.get_feature_names_out()
    coefficients = model.coef_[0]
    nonzero = np.flatnonzero(coefficients)
    positive = nonzero[np.argsort(coefficients[nonzero])[-TOP_PER_DIRECTION:]]
    negative = nonzero[np.argsort(coefficients[nonzero])[:TOP_PER_DIRECTION]]
    selected = np.concatenate([positive, negative])

    y = test["label_hiconf"].to_numpy(dtype=int)
    peptide = test["peptide"].astype(str).to_numpy()
    split_rows = []
    for feature_index in selected:
        motif = motifs[feature_index]
        present = np.array([motif in sequence for sequence in peptide])
        positive_present = int(((y == 1) & present).sum())
        positive_absent = int(((y == 1) & ~present).sum())
        negative_present = int(((y == 0) & present).sum())
        negative_absent = int(((y == 0) & ~present).sum())
        odds_ratio, p_value = fisher_exact(
            [[positive_present, positive_absent], [negative_present, negative_absent]],
            alternative="two-sided",
        )
        split_rows.append(
            {
                "split": split_name,
                "motif": motif,
                "coefficient": float(coefficients[feature_index]),
                "coefficient_direction": "positive" if coefficients[feature_index] > 0 else "negative",
                "test_positive_present": positive_present,
                "test_positive_absent": positive_absent,
                "test_negative_present": negative_present,
                "test_negative_absent": negative_absent,
                "test_odds_ratio": float(odds_ratio),
                "test_fisher_p": float(p_value),
                "frac_hydrophobic": fraction(motif, HYDROPHOBIC),
                "frac_positive_charge": fraction(motif, POSITIVE_CHARGE),
                "frac_negative_charge": fraction(motif, NEGATIVE_CHARGE),
                "frac_aromatic": fraction(motif, AROMATIC),
                "frac_small_polar": fraction(motif, SMALL_POLAR),
            }
        )
    split_frame = pd.DataFrame(split_rows)
    split_frame["test_fisher_fdr"] = benjamini_hochberg(split_frame["test_fisher_p"].to_numpy())
    split_frame["direction_consistent_with_test_odds"] = np.where(
        split_frame["coefficient_direction"].eq("positive"),
        split_frame["test_odds_ratio"] > 1,
        split_frame["test_odds_ratio"] < 1,
    )
    all_rows.append(split_frame)
    summary.append(
        {
            "split": split_name,
            "vocabulary_size": int(len(motifs)),
            "nonzero_coefficients": int(len(nonzero)),
            "selected_motifs": int(len(split_frame)),
            "heldout_fdr_below_0_05": int((split_frame["test_fisher_fdr"] < 0.05).sum()),
            "direction_consistent_count": int(split_frame["direction_consistent_with_test_odds"].sum()),
        }
    )

pd.concat(all_rows, ignore_index=True).to_csv(OUT / "heldout_motif_enrichment.csv", index=False)
pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
with (OUT / "interpretation_constraints.json").open("w", encoding="utf-8") as stream:
    json.dump(
        {
            "selection": "Top 50 positive and top 50 negative nonzero L1-logistic coefficients selected using training data only.",
            "heldout_test": "Two-sided Fisher exact tests for motif presence versus label; Benjamini-Hochberg correction within split.",
            "hla_constraint": "Reliable peptide-level HLA restriction is absent from the frozen modeling table. No motif is labeled an HLA anchor from these data alone.",
            "mechanistic_constraint": "N-gram association does not establish residue-level mechanism or causal IFN-gamma induction.",
        },
        stream,
        indent=2,
        sort_keys=True,
    )
print(pd.DataFrame(summary).to_string(index=False))
print("\nExact split strongest held-out motifs")
combined = pd.concat(all_rows, ignore_index=True)
print(combined[combined["split"] == "exact"].sort_values("test_fisher_p").head(15).to_string(index=False))
