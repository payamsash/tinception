#!/usr/bin/env python
"""
Slide 15 figure for the TRI talk — "44 regions, or 3?".

Two panels sharing an x-axis in Cohen's d.

  A  the three regions that survive the Bayesian criteria, each shown three
     ways: the clinical cohort alone, UK Biobank alone, and the hierarchical
     pooled effect with its 95% HDI. This makes the mechanism visible — UKB is
     tight and excludes zero, the clinical cohort is wide, and the pooled
     posterior is wider than UKB alone because the model pays for
     between-dataset heterogeneity. That is why the HDI crosses zero.

  B  every region that passes the nominal meta-analytic criterion (44 of them),
     as pooled effects with 95% HDIs, sorted. None excludes zero. The three from
     panel A are highlighted.

Everything is read from the Bayesian result table; the per-dataset confidence
intervals are rebuilt with the same Cohen's-d standard error the analysis used
(`cohend_se` in 10_compare_datasets_bayesian.py) and the same group sizes.

Usage
  venv/bin/python src/figures/make_slide15_forest.py
  venv/bin/python src/figures/make_slide15_forest.py --outdir /some/where
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

TINCEPTION = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
DEFAULT_OUT = TINCEPTION / "plots" / "talk_tri"
BAYES_CSV = TINCEPTION / "subcortical_roi" / "bayesian" / "bayesian_dmin_0.05.csv"

# group sizes and the equivalence bound, from 10_compare_datasets_bayesian.py:14-16
N_MAIN_TI, N_MAIN_CO = 285, 284
N_UKB_TI,  N_UKB_CO  = 2180, 2176
D_MIN = 0.05
HDI_PROB = 0.95        # 10_compare_datasets_bayesian.py:121 — az.hdi(..., hdi_prob=0.95)

PRETTY = {
    "Accessory-Basal-nucleus-lh": "Accessory-basal amygdala  (L)",
    "Medial-nucleus-rh":          "Medial amygdala  (R)",
    "MDm-rh":                     "MDm thalamus  (R)",
}

MAIN   = "#2C6E9B"
UKB    = "#1B7F5E"
POOLED = "#C2482D"
INK    = "#1A1A1A"
MUTED  = "#6B7280"
GRID   = "#E3E6EA"
FAINT  = "#B9BEC5"

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


def cohend_se(d, n1, n2):
    """Standard error of Cohen's d — identical to the analysis script."""
    return np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))


def load():
    df = pd.read_csv(BAYES_CSV)
    df["hdi_excludes_zero"] = ~((df.theta_hdi_lo < 0) & (df.theta_hdi_hi > 0))
    df["se_main"] = cohend_se(df.d_main, N_MAIN_TI, N_MAIN_CO)
    df["se_ukb"]  = cohend_se(df.d_ukb,  N_UKB_TI,  N_UKB_CO)
    bayes2 = df[df.tier_bayes.str.startswith("Tier 2")].copy()
    meta2  = df[df.tier_meta.str.startswith("Tier 2")].copy()
    if len(bayes2) != 3:
        raise ValueError(f"expected 3 Bayesian tier-2 regions, found {len(bayes2)}")
    missing = set(PRETTY) - set(bayes2.brain_label)
    if missing:
        raise ValueError(f"tier-2 regions changed; no display name for {missing}")
    return df, bayes2, meta2


def rope(ax):
    ax.axvspan(-D_MIN, D_MIN, color=GRID, alpha=0.75, zorder=0, linewidth=0)
    ax.axvline(0, color=MUTED, linewidth=1.0, zorder=1)


def panel_a(ax, bayes2):
    bayes2 = bayes2.sort_values("theta_mean", ascending=True).reset_index(drop=True)
    rope(ax)
    yticks, ylabels = [], []

    for i, r in bayes2.iterrows():
        base = i * 3.4
        rows = [
            (base + 0.85, r.d_main, 1.96 * r.se_main, MAIN, "o", 6.5),
            (base + 0.00, r.d_ukb,  1.96 * r.se_ukb,  UKB,  "s", 6.0),
        ]
        for y, d, half, colour, marker, size in rows:
            ax.plot([d - half, d + half], [y, y], color=colour, linewidth=2.0,
                    solid_capstyle="round", zorder=3)
            ax.plot(d, y, marker=marker, markersize=size, color=colour, zorder=4)

        # pooled posterior — heavier, with a diamond
        y = base - 0.85
        ax.plot([r.theta_hdi_lo, r.theta_hdi_hi], [y, y], color=POOLED,
                linewidth=3.4, solid_capstyle="round", zorder=3)
        ax.plot(r.theta_mean, y, marker="D", markersize=8.5, color=POOLED, zorder=4)

        yticks.append(base)
        ylabels.append(f"{PRETTY[r.brain_label]}\n"
                       f"P(|θ|>{D_MIN}) = {r.p_exceeds_dmin:.2f}   "
                       f"P(same sign) = {r.p_same_sign:.2f}")
        if i < len(bayes2) - 1:
            ax.axhline(base + 1.70, color=GRID, linewidth=0.9, zorder=0)

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=9.2, linespacing=1.7,
                       ha="right", multialignment="right")
    ax.set_ylim(-2.0, (len(bayes2) - 1) * 3.4 + 1.75)
    ax.set_xlabel("Cohen's d   (positive = larger in tinnitus)", fontsize=9, color=MUTED)
    ax.set_title("The 3 regions that survive", fontsize=12.5, fontweight="bold",
                 color=INK, pad=12, loc="left")
    ax.tick_params(labelsize=9, length=3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)


def panel_b(ax, meta2, bayes2):
    meta2 = meta2.sort_values("theta_mean").reset_index(drop=True)
    highlight = set(bayes2.brain_label)
    rope(ax)

    for i, r in meta2.iterrows():
        hit = r.brain_label in highlight
        colour = POOLED if hit else FAINT
        ax.plot([r.theta_hdi_lo, r.theta_hdi_hi], [i, i], color=colour,
                linewidth=2.6 if hit else 1.3, solid_capstyle="round",
                zorder=4 if hit else 2)
        ax.plot(r.theta_mean, i, marker="D" if hit else "o",
                markersize=6.5 if hit else 3.0, color=colour, zorder=5 if hit else 3)
    top = len(meta2) - 1
    ax.annotate("the 3 in the left panel",
                xy=(meta2.theta_hdi_hi.iloc[top] + 0.005, top - 1),
                xytext=(0.34, top - 7.5),
                fontsize=8.6, color=POOLED, fontweight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=POOLED, linewidth=1.1,
                                connectionstyle="arc3,rad=-0.25"))

    ax.set_yticks([])
    ax.set_ylim(-1.5, len(meta2) + 0.5)
    ax.set_xlabel("pooled θ with 95% HDI", fontsize=9, color=MUTED)
    ax.set_title(f"All {len(meta2)} regions the nominal criterion calls \"consistent\"",
                 fontsize=12.5, fontweight="bold", color=INK, pad=12, loc="left")
    ax.tick_params(labelsize=9, length=3)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.text(0.5, -0.105,
            f"Not one of the {len(meta2)} excludes zero — nor does any of the 97 regions tested.",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=8.6, color=INK, style="italic")


def build(outdir: Path):
    df, bayes2, meta2 = load()

    fig = plt.figure(figsize=(13.2, 7.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.0], wspace=0.16,
                          left=0.215, right=0.965, top=0.795, bottom=0.135)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    panel_a(ax_a, bayes2)
    panel_b(ax_b, meta2, bayes2)
    for ax in (ax_a, ax_b):             # same scale in both panels — they are both Cohen's d
        ax.set_xlim(-0.16, 0.42)

    fig.text(0.5, 0.955, "Forty-four regions, or three?",
             ha="center", fontsize=18, fontweight="bold", color=INK)
    fig.text(0.5, 0.910,
             "The same 97 regions, two criteria. Nominal significance in both datasets plus a "
             "direction match gives 44; a hierarchical model that\nprices in between-dataset "
             "heterogeneity gives 3 — and even those three have posteriors that cross zero.",
             ha="center", fontsize=9.6, color=MUTED, linespacing=1.55)

    fig.legend(handles=[
        Line2D([], [], color=MAIN,   marker="o", markersize=6.5, linewidth=2.0,
               label=f"clinical cohort, 95% CI  (n = {N_MAIN_TI}/{N_MAIN_CO})"),
        Line2D([], [], color=UKB,    marker="s", markersize=6.0, linewidth=2.0,
               label=f"UK Biobank, 95% CI  (n = {N_UKB_TI:,}/{N_UKB_CO:,})"),
        Line2D([], [], color=POOLED, marker="D", markersize=8.0, linewidth=3.4,
               label=f"pooled θ, {int(HDI_PROB*100)}% HDI  (hierarchical)"),
        Line2D([], [], color=GRID,   linewidth=9,
               label=f"|d| < {D_MIN}  (equivalence bound)"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.008), ncol=4, frameon=False,
        fontsize=8.6, handlelength=2.0, columnspacing=2.2)

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"slide15_bayesian_forest.{ext}",
                    dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"regions tested            : {len(df)}")
    print(f"same direction            : {int(df.same_sign.sum())}")
    print(f"nominal 'consistent'      : {len(meta2)}")
    print(f"Bayesian tier 2           : {len(bayes2)}")
    print(f"95% HDIs excluding zero   : {int(df.hdi_excludes_zero.sum())} of {len(df)}")
    print(f"of the {len(meta2)}: P(|θ|>{D_MIN})>=0.80 in {int((meta2.p_exceeds_dmin>=0.80).sum())}, "
          f"P(same sign)>=0.90 in {int((meta2.p_same_sign>=0.90).sum())}  <- the binding criterion")
    print(f"\nwrote {outdir / 'slide15_bayesian_forest.pdf'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    build(ap.parse_args().outdir)


if __name__ == "__main__":
    main()
