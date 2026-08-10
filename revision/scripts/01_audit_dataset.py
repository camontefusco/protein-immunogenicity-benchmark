from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path(os.environ.get("PIB_CURATED_DATA_DIR", "data/curated"))
OUT = Path(__file__).resolve().parents[1] / "results" / "generated" / "dataset_audit"
OUT.mkdir(parents=True, exist_ok=True)

assay = pd.read_csv(SOURCE / "assay_level_with_year.csv")
peptide = pd.read_csv(SOURCE / "peptide_level_hiconf_with_year.csv")

required_assay = {
    "peptide", "label", "assay_method", "virus_species", "protein",
    "earliest_pub_year", "latest_pub_year", "n_references",
}
required_peptide = {"peptide", "label_hiconf", "pos_rate", "total_tested"}
missing = {
    "assay": sorted(required_assay - set(assay.columns)),
    "peptide": sorted(required_peptide - set(peptide.columns)),
}
if missing["assay"] or missing["peptide"]:
    raise ValueError(f"Missing required columns: {missing}")

assay["peptide"] = assay["peptide"].astype(str).str.upper().str.strip()
peptide["peptide"] = peptide["peptide"].astype(str).str.upper().str.strip()

valid_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
invalid_assay = assay["peptide"].map(lambda s: not set(s) <= valid_alphabet)
invalid_peptide = peptide["peptide"].map(lambda s: not set(s) <= valid_alphabet)

label_by_peptide = assay.groupby("peptide")["label"].agg(["count", "nunique", "mean"])
year_counts = (
    peptide.assign(year=pd.to_numeric(peptide["earliest_pub_year"], errors="coerce"))
    .groupby("year", dropna=False)
    .agg(peptides=("peptide", "size"), positives=("label_hiconf", "sum"))
    .reset_index()
)
year_counts["positive_rate"] = year_counts["positives"] / year_counts["peptides"]
year_counts.to_csv(OUT / "peptide_counts_by_earliest_year.csv", index=False)

summary = {
    "assay_rows": int(len(assay)),
    "assay_unique_peptides": int(assay["peptide"].nunique()),
    "assay_positive_rate": float(assay["label"].mean()),
    "assay_exact_duplicate_rows": int(assay.duplicated().sum()),
    "assay_invalid_sequence_rows": int(invalid_assay.sum()),
    "assay_peptides_with_conflicting_labels": int((label_by_peptide["nunique"] > 1).sum()),
    "hiconf_rows": int(len(peptide)),
    "hiconf_unique_peptides": int(peptide["peptide"].nunique()),
    "hiconf_positive_rate": float(peptide["label_hiconf"].mean()),
    "hiconf_exact_duplicate_rows": int(peptide.duplicated().sum()),
    "hiconf_duplicate_peptide_rows": int(peptide.duplicated("peptide").sum()),
    "hiconf_invalid_sequence_rows": int(invalid_peptide.sum()),
    "hiconf_length_min": int(peptide["peptide"].str.len().min()),
    "hiconf_length_max": int(peptide["peptide"].str.len().max()),
    "hiconf_earliest_year_missing": int(peptide["earliest_pub_year"].isna().sum()),
    "hiconf_latest_year_missing": int(peptide["latest_pub_year"].isna().sum()),
    "outcome_derived_columns_excluded_from_predictors": [
        "label", "outcome_raw", "n_positive", "response_freq_pct",
        "total_positive", "pos_rate", "label_hiconf",
    ],
    "post_observation_or_provenance_columns_restricted_to_diagnostics": [
        "latest_pub_year", "n_references", "pmids", "total_tested",
    ],
}

with (OUT / "summary.json").open("w", encoding="utf-8") as stream:
    json.dump(summary, stream, indent=2, sort_keys=True)

print(json.dumps(summary, indent=2, sort_keys=True))
