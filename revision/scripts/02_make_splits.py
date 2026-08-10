from __future__ import annotations

import json
import os
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, train_test_split


SEED = 42
SOURCE = Path(os.environ.get("PIB_PEPTIDE_DATA", "data/curated/peptide_level_hiconf_with_year.csv"))
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "generated" / "splits"
OUT.mkdir(parents=True, exist_ok=True)


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def deletion_keys(sequence: str, max_deletions: int = 2):
    yield sequence
    for n_delete in range(1, max_deletions + 1):
        for positions in combinations(range(len(sequence)), n_delete):
            blocked = set(positions)
            yield "".join(char for idx, char in enumerate(sequence) if idx not in blocked)


def bounded_edit_distance(a: str, b: str, limit: int = 2) -> int:
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, start=1):
        current = [i] + [limit + 1] * len(b)
        lo = max(1, i - limit)
        hi = min(len(b), i + limit)
        row_min = limit + 1
        for j in range(lo, hi + 1):
            current[j] = min(
                previous[j] + 1,
                current[j - 1] + 1,
                previous[j - 1] + (char_a != b[j - 1]),
            )
            row_min = min(row_min, current[j])
        if row_min > limit:
            return limit + 1
        previous = current
    return previous[len(b)]


data = pd.read_csv(SOURCE)
data["peptide"] = data["peptide"].astype(str).str.upper().str.strip()
alphabet = set("ACDEFGHIKLMNPQRSTVWY")
data["valid_sequence"] = data["peptide"].map(lambda value: set(value) <= alphabet)
data["valid_year"] = pd.to_numeric(data["earliest_pub_year"], errors="coerce")

# Invalid sequences remain documented but are excluded from every model split.
eligible = data.index[data["valid_sequence"]].to_numpy()
train_idx, test_idx = train_test_split(
    eligible,
    test_size=0.20,
    random_state=SEED,
    stratify=data.loc[eligible, "label_hiconf"],
)
data["exact_split"] = "excluded"
data.loc[train_idx, "exact_split"] = "train"
data.loc[test_idx, "exact_split"] = "test"

# Cluster sequences connected by edit distance <= 2. This directly addresses the
# reviewer's one- or two-residue substitution concern and additionally catches
# one- or two-residue indels. Components are transitive by design.
seqs = data.loc[eligible, "peptide"].tolist()
uf = UnionFind(len(seqs))
index: dict[str, list[int]] = defaultdict(list)
candidate_pairs: set[tuple[int, int]] = set()
for local_idx, sequence in enumerate(seqs):
    for key in deletion_keys(sequence):
        for other in index[key]:
            candidate_pairs.add((other, local_idx))
        index[key].append(local_idx)

for left, right in sorted(candidate_pairs):
    if bounded_edit_distance(seqs[left], seqs[right], 2) <= 2:
        uf.union(left, right)

roots = [uf.find(idx) for idx in range(len(seqs))]
root_to_cluster = {root: number for number, root in enumerate(sorted(set(roots)))}
clusters = np.array([root_to_cluster[root] for root in roots])
data["edit_distance_2_cluster"] = -1
data.loc[eligible, "edit_distance_2_cluster"] = clusters

# Search deterministic group splits for the closest match to 20% test size and
# global class prevalence, never selecting a split based on predictive performance.
best = None
target_rate = data.loc[eligible, "label_hiconf"].mean()
for split_seed in range(512):
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=split_seed)
    local_train, local_test = next(
        splitter.split(np.zeros(len(eligible)), data.loc[eligible, "label_hiconf"], clusters)
    )
    rate = data.loc[eligible[local_test], "label_hiconf"].mean()
    size = len(local_test) / len(eligible)
    score = abs(rate - target_rate) + abs(size - 0.20)
    candidate = (score, split_seed, local_train, local_test)
    if best is None or candidate[:2] < best[:2]:
        best = candidate

_, selected_seed, local_train, local_test = best
data["cluster_split"] = "excluded"
data.loc[eligible[local_train], "cluster_split"] = "train"
data.loc[eligible[local_test], "cluster_split"] = "test"

# Temporal robustness analysis: train on peptides first published before the
# smallest year that gives approximately the most recent 20% as test. Peptides
# without a year are excluded rather than assigned retrospectively.
year_eligible = data.index[data["valid_sequence"] & data["valid_year"].notna()]
candidate_years = sorted(data.loc[year_eligible, "valid_year"].astype(int).unique())
cutoff = min(
    candidate_years,
    key=lambda year: abs((data.loc[year_eligible, "valid_year"] >= year).mean() - 0.20),
)
data["temporal_split"] = "excluded"
data.loc[year_eligible[data.loc[year_eligible, "valid_year"] < cutoff], "temporal_split"] = "train"
data.loc[year_eligible[data.loc[year_eligible, "valid_year"] >= cutoff], "temporal_split"] = "test"

data.to_csv(OUT / "peptide_split_assignments.csv", index=False)

def split_summary(column: str):
    return {
        name: {
            "n": int(len(group)),
            "positive_n": int(group["label_hiconf"].sum()),
            "positive_rate": float(group["label_hiconf"].mean()),
            "unique_peptides": int(group["peptide"].nunique()),
            "clusters": int(group.loc[group["edit_distance_2_cluster"] >= 0, "edit_distance_2_cluster"].nunique()),
        }
        for name, group in data.groupby(column)
    }

cluster_train = set(data.loc[data["cluster_split"] == "train", "edit_distance_2_cluster"])
cluster_test = set(data.loc[data["cluster_split"] == "test", "edit_distance_2_cluster"])
summary = {
    "seed": SEED,
    "cluster_definition": "transitive connected components at Levenshtein edit distance <= 2",
    "cluster_split_search_seeds": 512,
    "cluster_split_selected_seed": int(selected_seed),
    "cluster_count": int(len(set(clusters))),
    "largest_cluster": int(pd.Series(clusters).value_counts().max()),
    "cluster_overlap_train_test": int(len(cluster_train & cluster_test)),
    "temporal_cutoff_year_test_inclusive": int(cutoff),
    "exact": split_summary("exact_split"),
    "cluster": split_summary("cluster_split"),
    "temporal": split_summary("temporal_split"),
}
with (OUT / "split_summary.json").open("w", encoding="utf-8") as stream:
    json.dump(summary, stream, indent=2, sort_keys=True)
print(json.dumps(summary, indent=2, sort_keys=True))
