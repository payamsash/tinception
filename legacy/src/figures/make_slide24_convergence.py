#!/usr/bin/env python
"""
Slide 24 figure for the TRI talk — "the same axis appears in resting-state EEG".

Two columns, three rows. Left: structural MRI (Tinception). Right: resting-state
EEG (TIDE). The point of the figure is the *parallel*, not a numerical comparison
— the two modalities are on different scales and are deliberately given their own
x-axes and their own units.

  row 1   how patients are distributed along each modality's deviation axis,
          split into the group that deviates from normative expectation and the
          group that looks normative. Areas are proportional to group size.
  row 2   median hearing in each group.
  row 3   the dose-response: proportion of patients in the deviant group across
          hearing-loss quartiles.

MRI numbers are computed here from the result files. EEG numbers are transcribed
from ~/Tinnorm (see EEG_* below for the per-value provenance) and the row-1 curves
are normal approximations from the reported mean and SD, which is stated on the
figure.

Usage
  venv/bin/python src/figures/make_slide24_convergence.py
  venv/bin/python src/figures/make_slide24_convergence.py --outdir /some/where
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats

TINCEPTION = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
DEFAULT_OUT = TINCEPTION / "plots" / "talk_tri"

ASSIGNMENTS = TINCEPTION / "biotypes_gmm" / "primary" / "biotype_assignments_primary_(excl._london_2).csv"
NORM_DIR    = TINCEPTION / "subcortical_roi" / "norm_models_main"

# the 11 tier-2/3 convergent regions the mixture model was fitted on
REGIONS = ["Accessory-Basal-nucleus-lh", "Basal-nucleus-lh", "Lateral-nucleus-lh",
           "GC-ML-DG-body-rh", "presubiculum-head-rh",
           "MDm-lh", "MDm-rh", "MDl-lh", "CeM-rh", "CeM-lh", "PuA-lh"]

# ── EEG constants, from ~/Tinnorm ────────────────────────────────────────────
# FN_SUBGROUP_ANALYSIS.md §4b   discriminant-axis projection, mean and SD
EEG_DEV   = dict(n=183, mu=+1.658, sd=3.255, label="detected by the model")
EEG_NORM  = dict(n=85,  mu=-0.473, sd=3.247, label="missed — looks like a control")
# ANALYSIS_RESULTS.md, misclassification raincloud   median HF-PTA, dB HL
EEG_HEAR  = {"deviant": 48.8, "normative": 38.6}
# FN_SUBGROUP_ANALYSIS.md SA3   FN rate by HF-PTA quartile -> detected rate
EEG_QUART = [100 - r for r in (37, 43, 31, 15)]
EEG_SIL   = 0.072          # FN_SUBGROUP_ANALYSIS.md §4a, k-means in SHAP-PCA space
MRI_SIL   = 0.121          # biotypes_gmm/BIOTYPING_RESULTS.md, primary run

# colourblind-safe, high contrast on a projector; matches the slide 2 figures
DEVIANT   = "#C2482D"
NORMATIVE = "#2C6E9B"
INK       = "#1A1A1A"
MUTED     = "#6B7280"
GRID      = "#E3E6EA"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,          # editable text in Illustrator / PowerPoint
    "ps.fonttype": 42,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


# ─────────────────────────────────────────────────────────────────────────────
# data
# ─────────────────────────────────────────────────────────────────────────────
def clean_region_names(cols: pd.Index) -> pd.Index:
    """Normalise region names to 'X-lh' / 'X-rh' — same rule as 13_biotype_gmm.py."""
    return (cols.str.removeprefix("Volume_of_")
                .str.replace(r"^Left-(.+)$",  r"\1-lh", regex=True)
                .str.replace(r"^Right-(.+)$", r"\1-rh", regex=True))


def load_mri():
    """Per-patient mean normative deviation across the 11 clustering regions."""
    asg = pd.read_csv(ASSIGNMENTS)
    frames = []
    for atlas in ("Amygdala", "Hippo", "Thalamic-nuclei"):
        z = pd.read_csv(NORM_DIR / atlas / "results" / "Z_test.csv")
        z.columns = clean_region_names(pd.Index(z.columns))
        z = z.rename(columns={"subject_ids": "tinception_id", "subject-ids": "tinception_id"})
        frames.append(z.set_index("tinception_id")[[c for c in REGIONS if c in z.columns]])
    Z = pd.concat(frames, axis=1)
    missing = [r for r in REGIONS if r not in Z.columns]
    if missing:
        raise ValueError(f"regions not found in the Z files: {missing}")

    df = asg.set_index("tinception_id").join(Z, how="inner")
    if len(df) != len(asg):
        raise ValueError(f"joined {len(df)} of {len(asg)} patients — check the ID column")
    df["meanZ"] = df[REGIONS].mean(axis=1)
    df["quartile"] = pd.qcut(df["PTA"], 4, labels=False)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# panels
# ─────────────────────────────────────────────────────────────────────────────
def panel_distribution(ax, groups, xlabel, note):
    """groups: list of (values_or_(mu,sd), n, colour, label). Areas ∝ n."""
    lo = min(g["lo"] for g in groups)
    hi = max(g["hi"] for g in groups)
    pad = 0.06 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, 600)

    for g in groups:
        dens = g["pdf"](grid) * g["n"]
        ax.fill_between(grid, dens, color=g["colour"], alpha=0.42, linewidth=0)
        ax.plot(grid, dens, color=g["colour"], linewidth=2.0)

    ax.set_xlabel(xlabel, fontsize=8.5, color=MUTED, labelpad=6)
    ax.set_ylabel("patients", fontsize=8.5, color=MUTED)
    ax.set_yticks([])
    ax.tick_params(labelsize=8.5, length=3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.text(0.5, -0.33, note, transform=ax.transAxes, ha="center", va="top",
            fontsize=7.2, color=MUTED, style="italic")


def panel_hearing(ax, deviant, normative, xlabel):
    vals = [normative, deviant]
    cols = [NORMATIVE, DEVIANT]
    ypos = [1, 0]
    ax.barh(ypos, vals, height=0.52, color=cols, alpha=0.85, linewidth=0)
    for y, v in zip(ypos, vals):
        ax.text(v + 0.9, y, f"{v:.1f}", va="center", ha="left",
                fontsize=9.5, fontweight="bold", color=INK)
    ax.set_yticks(ypos)
    ax.set_yticklabels(["looks normative", "deviates"], fontsize=8.5)
    ax.set_xlim(0, max(vals) * 1.30)
    ax.set_xlabel(xlabel, fontsize=8.5, color=MUTED, labelpad=4)
    ax.tick_params(labelsize=8.5, length=3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def panel_quartiles(ax, rates, ylabel):
    x = np.arange(4)
    ax.bar(x, rates, width=0.62, color=DEVIANT, alpha=0.85, linewidth=0)
    for xi, r in zip(x, rates):
        ax.text(xi, r + 2.5, f"{r:.0f}%", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=INK)
    ax.plot(x, rates, color=INK, linewidth=1.1, marker="o", markersize=3.4,
            markerfacecolor="white", markeredgecolor=INK, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels(["Q1\nbest", "Q2", "Q3", "Q4\nworst"], fontsize=8.5)
    ax.set_ylim(0, max(rates) * 1.32)
    ax.set_ylabel(ylabel, fontsize=8.5, color=MUTED)
    ax.set_xlabel("hearing-loss quartile", fontsize=8.5, color=MUTED, labelpad=2)
    ax.tick_params(labelsize=8.5, length=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


# ─────────────────────────────────────────────────────────────────────────────
def build(outdir: Path):
    mri = load_mri()
    dev = mri[mri.Biotype == 0]
    nrm = mri[mri.Biotype == 1]
    mri_rates = [100 * (g.Biotype == 0).mean() for _, g in mri.groupby("quartile")]
    mri_hear = {"deviant": dev.PTA.median(), "normative": nrm.PTA.median()}

    fig = plt.figure(figsize=(11.4, 9.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.20, 0.62, 1.00],
                          hspace=0.95, wspace=0.26,
                          left=0.075, right=0.975, top=0.845, bottom=0.150)

    # ── row 1 ────────────────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    kde_d, kde_n = stats.gaussian_kde(dev.meanZ), stats.gaussian_kde(nrm.meanZ)
    span = (mri.meanZ.min(), mri.meanZ.max())
    panel_distribution(
        ax,
        [dict(pdf=kde_n, n=len(nrm), colour=NORMATIVE, lo=span[0], hi=span[1]),
         dict(pdf=kde_d, n=len(dev), colour=DEVIANT,   lo=span[0], hi=span[1])],
        "mean normative deviation across 11 subcortical nuclei  (Z)",
        "observed distribution of all 270 patients",
    )
    ax_mri_dist = ax

    ax = fig.add_subplot(gs[0, 1])
    lo = min(g["mu"] - 3.2 * g["sd"] for g in (EEG_DEV, EEG_NORM))
    hi = max(g["mu"] + 3.2 * g["sd"] for g in (EEG_DEV, EEG_NORM))
    panel_distribution(
        ax,
        [dict(pdf=(lambda x, g=g: stats.norm.pdf(x, g["mu"], g["sd"])),
              n=g["n"], colour=c, lo=lo, hi=hi)
         for g, c in ((EEG_NORM, NORMATIVE), (EEG_DEV, DEVIANT))],
        "position on the control-to-tinnitus discriminant axis",
        "normal approximation from the reported mean ± SD (n = 268)",
    )
    ax_eeg_dist = ax

    # ── row 2 ────────────────────────────────────────────────────────────────
    panel_hearing(fig.add_subplot(gs[1, 0]),
                  mri_hear["deviant"], mri_hear["normative"],
                  "median PTA  (dB HL)")
    panel_hearing(fig.add_subplot(gs[1, 1]),
                  EEG_HEAR["deviant"], EEG_HEAR["normative"],
                  "median high-frequency PTA  (dB HL)")

    # ── row 3 ────────────────────────────────────────────────────────────────
    panel_quartiles(fig.add_subplot(gs[2, 0]), mri_rates,
                    "% of patients whose\nbrains deviate")
    panel_quartiles(fig.add_subplot(gs[2, 1]), EEG_QUART,
                    "% of patients whose\nbrains deviate")

    # ── column headers ───────────────────────────────────────────────────────
    for x, title, sub, sil in (
        (0.285, "Structural MRI", "Tinception · 9 sites · 270 patients · unsupervised clustering", MRI_SIL),
        (0.760, "Resting-state EEG", "TIDE · 7 sites · 268 patients · supervised classifier", EEG_SIL),
    ):
        fig.text(x, 0.898, title, ha="center", fontsize=13.5, fontweight="bold", color=INK)
        fig.text(x, 0.874, sub, ha="center", fontsize=8.6, color=MUTED)
        fig.text(x, 0.856, f"cluster silhouette {sil:.3f} — a gradient, not a boundary",
                 ha="center", fontsize=8.2, color=MUTED, style="italic")

    # ── title, legend, footer ────────────────────────────────────────────────
    fig.text(0.5, 0.975, "The same axis, in two cohorts and two modalities",
             ha="center", fontsize=17, fontweight="bold", color=INK)
    fig.text(0.5, 0.945,
             "Patients who deviate from normative expectation are the patients with more hearing loss — "
             "and the split is a gradient, not a category",
             ha="center", fontsize=10, color=MUTED)

    fig.legend(handles=[Patch(facecolor=DEVIANT, alpha=0.85,
                              label="deviates from normative expectation"),
                        Patch(facecolor=NORMATIVE, alpha=0.85,
                              label="looks like a control")],
               loc="lower center", bbox_to_anchor=(0.5, 0.055), ncol=2,
               frameon=False, fontsize=9.5, handlelength=1.5, columnspacing=2.4)

    fig.text(0.5, 0.012,
             "Separate cohorts, separate modalities, separate paradigms — the two x-axes are in different "
             "units and are not numerically comparable.",
             ha="center", fontsize=7.6, color=MUTED, style="italic")

    # divider between the two columns
    fig.add_artist(plt.Line2D([0.523, 0.523], [0.115, 0.855], transform=fig.transFigure,
                              color=GRID, linewidth=1.2))

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"slide24_mri_eeg_convergence.{ext}",
                    dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"MRI  deviant n={len(dev)} ({100*len(dev)/len(mri):.1f}%)  "
          f"median PTA {mri_hear['deviant']:.1f} vs {mri_hear['normative']:.1f} dB")
    print(f"MRI  deviant-group rate by PTA quartile: "
          + " · ".join(f"{r:.0f}%" for r in mri_rates))
    print(f"EEG  deviant-group rate by HF-PTA quartile: "
          + " · ".join(f"{r:.0f}%" for r in EEG_QUART))
    print(f"\nwrote {outdir / 'slide24_mri_eeg_convergence.pdf'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    build(ap.parse_args().outdir)


if __name__ == "__main__":
    main()
