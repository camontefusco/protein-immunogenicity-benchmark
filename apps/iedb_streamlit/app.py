import streamlit as st
import pandas as pd
import numpy as np

from pathlib import Path

from sklearn.pipeline import Pipeline
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, mean_squared_error, mean_absolute_error, r2_score

# ============================================================
# Streamlit App — IEDB IFN-γ Peptide Screening (Baseline)
# ============================================================
# Audience: vaccine design / epitope screening teams (non-ML specialists)
# Purpose: quickly rank input peptides by predicted IFN-γ positivity likelihood
#          using simple, reproducible baseline models trained on IEDB frozen data.
# ============================================================

st.set_page_config(page_title="IEDB IFN-γ Peptide Screening (Baseline)", layout="wide")

# -----------------------------
# Helpers
# -----------------------------
DATA_CANDIDATES = [
    Path("data/curated"),
    Path("data/frozen"),
    Path("."),                  # local
    Path("/mnt/data"),          # chat environment
    Path("/content"),           # colab
    Path("/content/drive/MyDrive"),
]

def find_file(fname: str) -> Path:
    for d in DATA_CANDIDATES:
        p = d / fname
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find {fname}. Put it in data/frozen/ or the working directory.")

def group_split(df: pd.DataFrame, group_col: str, test_size: float = 0.2, random_state: int = 42):
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(gss.split(df, groups=df[group_col]))
    return df.iloc[train_idx].copy(), df.iloc[test_idx].copy()

@st.cache_data
def load_data():
    assay = pd.read_csv(find_file("assay_level_with_year.csv"))
    rates = assay.groupby("peptide", as_index=False).agg(
        total_tested=("n_tested", "sum"),
        total_positive=("n_positive", "sum"),
    )
    rates["pos_rate"] = rates["total_positive"] / rates["total_tested"].replace(0, np.nan)
    hiconf = pd.read_csv(find_file("peptide_level_hiconf_with_year.csv"))
    return assay, rates, hiconf

@st.cache_resource
def train_models(hiconf: pd.DataFrame, rates: pd.DataFrame, random_state: int = 42):
    # --- Classification on HI-CONF (label_hiconf)
    cls_df = hiconf.dropna(subset=["peptide", "label_hiconf"]).copy()
    cls_train, cls_test = group_split(cls_df, group_col="peptide", test_size=0.2, random_state=random_state)

    X_train = cls_train["peptide"].astype(str)
    y_train = cls_train["label_hiconf"].astype(int)
    X_test  = cls_test["peptide"].astype(str)
    y_test  = cls_test["label_hiconf"].astype(int)

    vec = TfidfVectorizer(analyzer="char", ngram_range=(2,4), lowercase=False)

    dummy_cls = Pipeline([
        ("vec", vec),
        ("clf", DummyClassifier(strategy="most_frequent")),
    ])
    logreg = Pipeline([
        ("vec", vec),
        ("clf", LogisticRegression(max_iter=3000, class_weight="balanced", random_state=random_state)),
    ])

    dummy_cls.fit(X_train, y_train)
    logreg.fit(X_train, y_train)

    # probabilities for metrics
    yhat_dummy = dummy_cls.predict_proba(X_test)[:,1]
    yhat_lr = logreg.predict_proba(X_test)[:,1]

    cls_metrics = {
        "Dummy": {
            "ROC_AUC": float(roc_auc_score(y_test, yhat_dummy)),
            "PR_AUC": float(average_precision_score(y_test, yhat_dummy)),
            "Accuracy": float(accuracy_score(y_test, (yhat_dummy>=0.5).astype(int))),
        },
        "LogReg_char_2_4": {
            "ROC_AUC": float(roc_auc_score(y_test, yhat_lr)),
            "PR_AUC": float(average_precision_score(y_test, yhat_lr)),
            "Accuracy": float(accuracy_score(y_test, (yhat_lr>=0.5).astype(int))),
        }
    }

    # --- Regression on peptide rates (pos_rate)
    reg_df = rates.dropna(subset=["peptide", "pos_rate"]).copy()
    reg_train, reg_test = group_split(reg_df, group_col="peptide", test_size=0.2, random_state=random_state)

    Xr_train = reg_train["peptide"].astype(str)
    yr_train = reg_train["pos_rate"].astype(float)
    Xr_test  = reg_test["peptide"].astype(str)
    yr_test  = reg_test["pos_rate"].astype(float)

    dummy_reg = Pipeline([("vec", vec), ("reg", DummyRegressor(strategy="mean"))])
    ridge = Pipeline([("vec", vec), ("reg", Ridge(alpha=1.0, random_state=random_state))])

    dummy_reg.fit(Xr_train, yr_train)
    ridge.fit(Xr_train, yr_train)

    pred_dummy = dummy_reg.predict(Xr_test)
    pred_ridge = ridge.predict(Xr_test)

    rmse_dummy = float(mean_squared_error(yr_test, pred_dummy) ** 0.5)
    rmse_ridge = float(mean_squared_error(yr_test, pred_ridge) ** 0.5)

    reg_metrics = {
        "Dummy_mean": {
            "RMSE": rmse_dummy,
            "MAE": float(mean_absolute_error(yr_test, pred_dummy)),
            "R2": float(r2_score(yr_test, pred_dummy)),
        },
        "Ridge_char_2_4": {
            "RMSE": rmse_ridge,
            "MAE": float(mean_absolute_error(yr_test, pred_ridge)),
            "R2": float(r2_score(yr_test, pred_ridge)),
        },
    }

    return {
        "models": {"dummy_cls": dummy_cls, "logreg": logreg, "dummy_reg": dummy_reg, "ridge": ridge},
        "metrics": {"classification": cls_metrics, "regression": reg_metrics},
        "splits": {"cls_test": cls_test, "reg_test": reg_test},
    }

def sanitize_peptides(raw: str):
    # Accept newline/comma/space separated peptides
    toks = [t.strip().upper() for t in raw.replace(",", "\n").splitlines()]
    toks = [t for t in toks if t]
    # Basic amino acid alphabet (20 AA + optional X)
    aa = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")  # allow common extras; we will warn
    cleaned, bad = [], []
    for t in toks:
        if all(ch in aa for ch in t):
            cleaned.append(t)
        else:
            bad.append(t)
    return cleaned, bad

# -----------------------------
# Load
# -----------------------------
assay_level, peptide_rates, peptide_hiconf = load_data()
bundle = train_models(peptide_hiconf, peptide_rates)
models = bundle["models"]
metrics = bundle["metrics"]

# -----------------------------
# UI
# -----------------------------
st.title("IEDB IFN-gamma Peptide Screening Baseline App")
st.caption("Exploratory baseline models trained on repository-curated human viral T-cell IFN-gamma assay data. This interface is not the publication benchmark or a clinical predictor.")

with st.expander("What this app does (in plain language)", expanded=False):
    st.markdown(
        """
- You paste peptide sequences (short amino-acid strings).
- The app uses **simple baseline models** to estimate:
  - **P(IFN-γ positive)** (classification on *high-confidence* peptides)
  - **Expected response rate** (regression on peptide-level pos_rate)
- Output helps you **rank** which peptides are *more likely* to trigger an IFN-γ response, so you can test fewer candidates first.
        """
    )

colA, colB = st.columns([1,1])

with colA:
    st.subheader("1) Paste peptide sequences")
    raw = st.text_area(
        "One peptide per line (or comma-separated). Example: AADLDDFSKQLQQSM",
        height=180
    )
    peptides, bad = sanitize_peptides(raw)
    if bad:
        st.warning(f"{len(bad)} sequences contain non-amino-acid characters and will be ignored: {bad[:5]}{'...' if len(bad)>5 else ''}")
    st.write(f"Valid peptides: **{len(peptides)}**")

with colB:
    st.subheader("2) Optional context filters (for exploration only)")
    st.caption("These filters do NOT change the baseline model (sequence-only), but help you explore how often a peptide appears in IEDB.")
    virus = st.selectbox("Virus", options=["(any)"] + sorted(assay_level["virus"].dropna().unique().tolist())[:200])
    assay_method = st.selectbox("Assay method", options=["(any)"] + sorted(assay_level["assay_method"].dropna().unique().tolist()))
    show_matches = st.checkbox("Show matching assay records for pasted peptides", value=True)

st.divider()

left, right = st.columns([1,1])

with left:
    st.subheader("Baseline performance snapshot")
    st.markdown("**Classification (HI-CONF, label_hiconf)**")
    st.dataframe(pd.DataFrame(metrics["classification"]).T.reset_index().rename(columns={"index":"Model"}), use_container_width=True)
    st.markdown("**Regression (pos_rate)**")
    st.dataframe(pd.DataFrame(metrics["regression"]).T.reset_index().rename(columns={"index":"Model"}), use_container_width=True)

with right:
    st.subheader("Dataset snapshot")
    st.write(f"Assay-level rows: **{len(assay_level):,}**")
    st.write(f"Unique peptides: **{assay_level['peptide'].nunique():,}**")
    st.write(f"HI-CONF peptides: **{peptide_hiconf['peptide'].nunique():,}**")
    st.write(f"Peptide rates rows: **{len(peptide_rates):,}**")
    st.caption("Tip: the HI-CONF subset is more conservative (drops peptides with conflicting labels).")

st.divider()

st.subheader("3) Predictions")
if not peptides:
    st.info("Paste at least 1 peptide sequence to get predictions.")
else:
    X = pd.Series(peptides, name="peptide")

    p_ifng = models["logreg"].predict_proba(X)[:,1]
    pred_rate = models["ridge"].predict(X)

    out = pd.DataFrame({
        "peptide": peptides,
        "pred_P_IFNg_positive": p_ifng,
        "pred_pos_rate": np.clip(pred_rate, 0, 1),
        "length": [len(p) for p in peptides],
    }).sort_values("pred_P_IFNg_positive", ascending=False)

    st.dataframe(out, use_container_width=True)

    st.download_button(
        "Download predictions (CSV)",
        data=out.to_csv(index=False).encode("utf-8"),
        file_name="iedb_ifng_predictions.csv",
        mime="text/csv",
    )

    st.caption("Interpretation: higher **pred_P_IFNg_positive** suggests higher likelihood of a positive IFN-γ readout; **pred_pos_rate** approximates expected response frequency across cohorts/studies.")

    if show_matches:
        st.subheader("4) Where do these peptides appear in IEDB?")
        df = assay_level[assay_level["peptide"].isin(peptides)].copy()
        if virus != "(any)":
            df = df[df["virus"] == virus]
        if assay_method != "(any)":
            df = df[df["assay_method"] == assay_method]

        if df.empty:
            st.write("No matching assay-level records found for the pasted peptides under the selected filters.")
        else:
            cols = ["peptide","virus","assay_method","label","n_tested","n_positive","response_freq_pct","protein","virus_strain"]
            cols = [c for c in cols if c in df.columns]
            st.dataframe(df[cols].sort_values(["peptide","virus"]), use_container_width=True)

st.divider()
st.subheader("Notes & next steps")
st.markdown(
"""
- This is a **baseline** app: sequence-only, fast, easy to reproduce.
- Next iterations can:
  - add **virus/protein context** features
  - replace char-ngrams with **CNN / LSTM / Transformer** peptide encoders
  - support **calibration** & decision thresholds for lab triage
"""
)
