from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "publication" / "tables"
FIGURES = ROOT / "results" / "publication" / "figures"
SUBMISSION = FIGURES
FIGURES.mkdir(parents=True, exist_ok=True)
SUBMISSION.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

# Corrected Figure 4: the classical comparator is explicitly Random Forest.
esm = pd.read_csv(TABLES / "table_05_esm2_comparison.csv").set_index("split").loc[["exact", "cluster", "temporal"]]
motif = pd.read_csv(TABLES / "table_05b_motif_summary.csv").set_index("split").loc[["exact", "cluster", "temporal"]]
x = np.arange(3)
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"wspace": 0.42})
w = 0.35
axs[0].bar(x - w / 2, esm["rf_pr_auc"], w, label="Random Forest", color="#2C7FB8")
axs[0].bar(x + w / 2, esm["esm2_frozen_pr_auc"], w, label="Frozen ESM-2 (8M)", color="#F28E2B")
axs[0].set_xticks(x, ["exact", "component", "temporal"])
axs[0].set_ylabel("PR-AUC")
axs[0].set_title("Frozen ESM-2 versus Random Forest")
axs[0].legend(frameon=False, fontsize=9)
axs[0].spines[["top", "right"]].set_visible(False)
axs[1].bar(x, motif["heldout_fdr_below_0_05"], color="#756BB1", width=0.65)
axs[1].set_xticks(x, ["exact", "component", "temporal"])
axs[1].set_ylabel("Held-out motifs with FDR < 0.05")
axs[1].set_title("No motifs passed FDR < 0.05 temporally")
axs[1].spines[["top", "right"]].set_visible(False)
fig.savefig(FIGURES / "figure_04_esm2_motifs.png", dpi=300, bbox_inches="tight")
fig.savefig(SUBMISSION / "figure_04_esm2_motifs.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Text-exact graphical abstract built as a native scientific diagram.
fig, ax = plt.subplots(figsize=(16, 9))
ax.set_xlim(0, 16); ax.set_ylim(0, 9); ax.axis("off")
ax.text(8, 8.35, "Leakage-aware benchmarking of viral T-cell IFN-gamma assay outcomes",
        ha="center", va="center", fontsize=23, fontweight="bold", color="#173F5F")

def box(x0, y0, width, height, title, body, color):
    patch = FancyBboxPatch((x0, y0), width, height, boxstyle="round,pad=0.03,rounding_size=0.18",
                           linewidth=1.4, edgecolor="#34495E", facecolor=color)
    ax.add_patch(patch)
    ax.text(x0 + 0.28, y0 + height - 0.30, title, ha="left", va="top", fontsize=14.5, fontweight="bold")
    ax.text(x0 + 0.28, y0 + height - 0.94, body, ha="left", va="top", fontsize=10.7, linespacing=1.30)

box(0.65, 4.75, 7.05, 2.45, "1  Recovered IEDB analysis state",
    "31,502 assay records; 17,336 unique peptides\n"
    "Original 161-column export unavailable\n"
    "Analysis begins from recovered curated tables", "#DCEAF7")
box(8.30, 4.75, 7.05, 2.45, "2  Leakage-aware sequence benchmark",
    "9,668 valid high-confidence peptides; 9.2% positive\n"
    "PR-AUC: 0.479 exact | 0.404 component | 0.151 temporal\n"
    "Validation design strongly changes apparent performance", "#E4F2DE")
box(0.65, 1.55, 7.05, 2.45, "3  Separate peptide-context task",
    "20,785 peptide-context rows\n"
    "Safe-context uplift: +0.236 exact | +0.254 component\n"
    "+0.047 temporal; persisted across three label constructions", "#FFF0C9")
box(8.30, 1.55, 7.05, 2.45, "4  Cautious interpretation",
    "Context associations are not causal; task PR-AUCs are not comparable\n"
    "No external dataset; temporal test not component-purged\n"
    "Small frozen ESM-2; missing raw export and HLA restriction\n"
    "Independent experimental validation is required", "#E9DFF2")
ax.annotate("", xy=(8.18, 5.97), xytext=(7.82, 5.97), arrowprops=dict(arrowstyle="-|>", lw=2, color="#173F5F"))
ax.annotate("", xy=(8.18, 2.77), xytext=(7.82, 2.77), arrowprops=dict(arrowstyle="-|>", lw=2, color="#173F5F"))
fig.savefig(SUBMISSION / "graphical_abstract_v5.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Wrote corrected Figure 4 and graphical abstract")
