#!/usr/bin/env python
"""
Two slides for the Tinnitus Quest grant video.

These are NOT conference figures. The audience is Tinnitus Quest's executive board —
three patients and two researchers — watching on a laptop, possibly a phone. So:
one sentence of takeaway per slide, huge type, no axis jargon, no statistics on the
face of the slide, and both cohorts side by side so the convergence is the point.

  slide 1   the two "subtypes" overlap almost completely, in both cohorts
  slide 2   what is actually there: a hearing gradient, not a boundary

Data comes from make_slide24_convergence.py so there is exactly one source of truth.

Usage
  venv/bin/python src/figures/make_video_slides.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_slide24_convergence import (           # noqa: E402
    load_mri, EEG_DEV, EEG_NORM, EEG_QUART, DEFAULT_OUT,
)

DEVIANT   = "#C2482D"
NORMATIVE = "#2C6E9B"
INK       = "#1A1A1A"
MUTED     = "#5B6068"
GRID      = "#E3E6EA"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.color": INK,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

FIGSIZE = (13.33, 7.5)      # 16:9


def _frame(ax, title):
    ax.set_title(title, fontsize=20, fontweight="bold", color=INK, pad=16)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.spines["bottom"].set_linewidth(1.4)
    ax.tick_params(labelsize=15, length=4, width=1.2)


def headline(fig, big, small):
    fig.text(0.5, 0.935, big, ha="center", fontsize=31, fontweight="bold", color=INK)
    fig.text(0.5, 0.868, small, ha="center", fontsize=17, color=MUTED)


def slide_overlap(mri, outdir):
    """Two overlapping distributions per cohort — the core visual claim."""
    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(1, 2, wspace=0.13, left=0.055, right=0.965,
                          top=0.755, bottom=0.20)
    dev = mri[mri.Biotype == 0]
    nrm = mri[mri.Biotype == 1]

    # left — MRI, real observed distributions
    ax = fig.add_subplot(gs[0, 0])
    grid = np.linspace(mri.meanZ.min() - 0.15, mri.meanZ.max() + 0.15, 600)
    for sub, colour in ((nrm, NORMATIVE), (dev, DEVIANT)):
        y = stats.gaussian_kde(sub.meanZ)(grid) * len(sub)
        ax.fill_between(grid, y, color=colour, alpha=0.45, linewidth=0)
        ax.plot(grid, y, color=colour, linewidth=3.0)
    _frame(ax, "Brain scans   ·   270 patients")
    ax.set_yticks([])

    # right — EEG, normal approximation from the reported mean and SD
    ax = fig.add_subplot(gs[0, 1])
    lo = min(g["mu"] - 3.3 * g["sd"] for g in (EEG_DEV, EEG_NORM))
    hi = max(g["mu"] + 3.3 * g["sd"] for g in (EEG_DEV, EEG_NORM))
    grid = np.linspace(lo, hi, 600)
    for g, colour in ((EEG_NORM, NORMATIVE), (EEG_DEV, DEVIANT)):
        y = stats.norm.pdf(grid, g["mu"], g["sd"]) * g["n"]
        ax.fill_between(grid, y, color=colour, alpha=0.45, linewidth=0)
        ax.plot(grid, y, color=colour, linewidth=3.0)
    _frame(ax, "EEG   ·   268 patients")
    ax.set_yticks([])

    for ax in fig.axes:
        ax.set_xticks([])
        ax.set_xlabel("how far the brain is from normal",
                      fontsize=15, color=MUTED, labelpad=10)

    headline(fig,
             "The two “subtypes” I found overlap almost completely",
             "Two separate patient groups, two different brain measures — the same answer both times")

    fig.legend(handles=[
        Line2D([], [], color=DEVIANT,   linewidth=9, alpha=0.75,
               label="the group whose brains look different"),
        Line2D([], [], color=NORMATIVE, linewidth=9, alpha=0.75,
               label="the group whose brains look normal"),
    ], loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=2, frameon=False,
        fontsize=16, handlelength=1.6, columnspacing=3.0)

    fig.text(0.5, 0.015, "If these were genuinely different kinds of tinnitus, "
                         "the two curves would sit apart.",
             ha="center", fontsize=14, color=MUTED, style="italic")
    return fig


def slide_hearing(mri, outdir):
    """The dose-response that does hold, in both cohorts."""
    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(1, 2, wspace=0.20, left=0.075, right=0.965,
                          top=0.755, bottom=0.20)
    mri_rates = [100 * (g.Biotype == 0).mean() for _, g in mri.groupby("quartile")]

    for col, (rates, title) in enumerate((
        (mri_rates, "Brain scans   ·   270 patients"),
        (EEG_QUART, "EEG   ·   268 patients"),
    )):
        ax = fig.add_subplot(gs[0, col])
        x = np.arange(4)
        ax.bar(x, rates, width=0.66, color=DEVIANT, alpha=0.88, linewidth=0)
        for xi, r in zip(x, rates):
            ax.text(xi, r + max(rates) * 0.035, f"{r:.0f}%", ha="center",
                    va="bottom", fontsize=21, fontweight="bold", color=INK)
        ax.set_xticks(x)
        ax.set_xticklabels(["best\nhearing", "", "", "worst\nhearing"], fontsize=15)
        ax.set_ylim(0, max(rates) * 1.30)
        ax.set_yticks([])
        _frame(ax, title)
        if col == 0:
            ax.set_ylabel("share of patients whose\nbrain looks different",
                          fontsize=15, color=MUTED, labelpad=14)

    headline(fig,
             "What is real is hearing — and it comes as a slope",
             "The worse the hearing, the more likely the brain looks different. No line where one group stops and another starts.")

    fig.text(0.5, 0.045, "Tinnitus distress, by contrast, predicts almost nothing.",
             ha="center", fontsize=16, color=INK)
    return fig


def main():
    outdir = DEFAULT_OUT
    outdir.mkdir(parents=True, exist_ok=True)
    mri = load_mri()

    for name, fn in (("video_slide1_overlap", slide_overlap),
                     ("video_slide2_hearing", slide_hearing)):
        fig = fn(mri, outdir)
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{name}.{ext}", dpi=200,
                        bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"wrote {outdir / (name + '.pdf')}")


if __name__ == "__main__":
    main()
