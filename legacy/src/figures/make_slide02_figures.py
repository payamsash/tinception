#!/usr/bin/env python
"""
Slide 2 figures for the TRI talk — "the tinnitus literature doesn't converge".

Two panels, generated independently:

  A  literature landscape   same region, opposite reported directions, small N.
                            Driven by slide02_literature.csv, which YOU fill in
                            from your own reference list. Any row still marked
                            verified=no forces a DRAFT watermark onto the figure.

  B  power curves           why that disagreement is expected. Computed from
                            first principles (noncentral t), no external data,
                            using the effect sizes actually observed in this
                            project.

Usage
  venv/bin/python src/figures/make_slide02_figures.py
  venv/bin/python src/figures/make_slide02_figures.py --outdir /some/where
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = Path("/Volumes/Extreme_SSD/payam_data/Tinception/plots/talk_tri")
LIT_CSV = HERE / "slide02_literature.csv"

# colourblind-safe, high contrast on a projector
INCREASE = "#C2482D"
DECREASE = "#2C6E9B"
NULL     = "#9AA0A6"
INK      = "#1A1A1A"
MUTED    = "#6B7280"
GRID     = "#E3E6EA"
ACCENT   = "#1B7F5E"

DIR_COLOR = {"increase": INCREASE, "decrease": DECREASE, "null": NULL}
DIR_LABEL = {
    "increase": "Increase in tinnitus",
    "decrease": "Decrease in tinnitus",
    "null":     "No difference",
}

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
# helpers
# ─────────────────────────────────────────────────────────────────────────────
def read_literature(path: Path):
    """Read the study table, skipping # comments and blank lines."""
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    lines = [l for l in path.read_text(encoding="utf8").splitlines()
             if l.strip() and not l.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines))

    out = []
    for r in rows:
        d = (r.get("direction") or "").strip().lower()
        if d not in DIR_COLOR:
            raise ValueError(f"Bad direction {d!r} for study {r.get('study')!r}; "
                             f"expected one of {sorted(DIR_COLOR)}")
        if not (r.get("n_ti") or "").strip() or not (r.get("n_co") or "").strip():
            raise ValueError(f"Missing sample size for {r.get('study')!r} / {r.get('region')!r}; "
                             "leave the row commented out until N is confirmed")
        out.append({
            "source": (r.get("source") or "").strip(),
            "study": (r.get("study") or "").strip(),
            "year": int(r["year"]) if (r.get("year") or "").strip() else None,
            "region": (r.get("region") or "").strip(),
            "direction": d,
            "n_ti": int(r["n_ti"]),
            "n_co": int(r["n_co"]),
            "verified": (r.get("verified") or "no").strip().lower() in {"yes", "y", "true", "1"},
        })
    if not out:
        raise ValueError(f"No study rows found in {path}")
    return out


def watermark(fig, text="DRAFT — PLACEHOLDER DATA, NOT FOR PRESENTATION"):
    fig.text(0.5, 0.5, text, fontsize=26, color="#D32F2F", alpha=0.16,
             ha="center", va="center", rotation=24, weight="bold", zorder=1000)


def tidy(ax, spines=("top", "right")):
    for s in spines:
        ax.spines[s].set_visible(False)
    for s in ax.spines.values():
        s.set_linewidth(0.9)
    ax.tick_params(length=3, width=0.9, labelsize=10)


def save(fig, outdir: Path, stem: str):
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"  saved  {p}")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Panel A — literature landscape
# ─────────────────────────────────────────────────────────────────────────────
def figure_literature(rows, outdir: Path):
    """Dot matrix: one row per region, one dot per study, coloured by direction.

    Dots are ordered increase → decrease → null within a row, so a region
    reported in both directions is immediately visible as a two-colour row.
    """
    regions = {}
    for r in rows:
        regions.setdefault(r["region"], []).append(r)

    def is_conflicted(rs):
        return {"increase", "decrease"} <= {x["direction"] for x in rs}

    # conflicted regions first, then most-studied
    order = sorted(regions,
                   key=lambda k: (is_conflicted(regions[k]), len(regions[k])),
                   reverse=True)
    ypos = {reg: i for i, reg in enumerate(reversed(order))}
    max_studies = max(len(v) for v in regions.values())

    # size the plot area so dots sit ~1 inch apart regardless of how many
    # studies the busiest region has
    lab_in, tail_in = 3.0, 1.9
    dot_in = max(0.42, min(1.0, 9.0 / max_studies))
    figw = lab_in + dot_in * max_studies + tail_in
    figh = 0.60 * len(order) + 2.5
    fig, ax = plt.subplots(figsize=(figw, figh))
    fig.subplots_adjust(top=1 - 0.95 / figh, bottom=1.30 / figh,
                        left=lab_in / figw, right=1 - tail_in / figw)

    all_n = [min(r["n_ti"], r["n_co"]) for r in rows]
    nmin, nmax = min(all_n), max(all_n)

    def msize(n):
        if nmax == nmin:
            return 260.0
        return 110 + 420 * (n - nmin) / (nmax - nmin)

    rank = {"increase": 0, "decrease": 1, "null": 2}
    for reg in order:
        y = ypos[reg]
        if y % 2 == 0:
            ax.axhspan(y - 0.5, y + 0.5, color="#F7F9FB", zorder=0)
        rs = sorted(regions[reg], key=lambda r: (rank[r["direction"]],
                                                 -min(r["n_ti"], r["n_co"])))
        for i, r in enumerate(rs):
            ax.scatter(i, y, s=msize(min(r["n_ti"], r["n_co"])),
                       facecolor=DIR_COLOR[r["direction"]],
                       edgecolor="white", linewidth=1.8, alpha=0.95, zorder=3)
        if is_conflicted(rs):
            ax.text(max_studies - 0.35, y, "conflicting", fontsize=10,
                    color=MUTED, style="italic", va="center", ha="left",
                    clip_on=False)

    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels(list(ypos.keys()), fontsize=11.5)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlim(-0.6, max_studies - 0.4)
    ax.set_xticks([])
    ax.set_xlabel("each dot is one published finding   ·   dot size = sample size",
                  fontsize=10.5, color=MUTED, labelpad=12)
    tidy(ax, spines=("top", "right", "left", "bottom"))
    ax.tick_params(axis="y", length=0)

    n_conflict = sum(1 for reg in order if is_conflicted(regions[reg]))
    n_studies = len({(r["study"], r["year"]) for r in rows})
    fig.text(0.035, 1 - 0.36 / figh, "The same regions, reported in opposite directions",
             fontsize=16, weight="bold", color=INK, ha="left", va="bottom")
    fig.text(0.035, 1 - 0.66 / figh,
             f"{len(rows)} findings from {n_studies} studies across {len(order)} regions   ·   "
             f"{n_conflict} regions reported in both directions   ·   "
             f"median n = {int(np.median(all_n))} per group",
             fontsize=10.5, color=MUTED, ha="left", va="bottom")

    handles = [Line2D([], [], marker="o", linestyle="none", markersize=11,
                      markerfacecolor=DIR_COLOR[d], markeredgecolor="white",
                      markeredgewidth=1.4, label=DIR_LABEL[d])
               for d in ("increase", "decrease", "null")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.075),
              ncol=3, frameon=False, fontsize=10.5, handletextpad=0.6,
              columnspacing=1.8)

    if not all(r["verified"] for r in rows):
        watermark(fig)

    save(fig, outdir, "slide02a_literature_landscape")


# ─────────────────────────────────────────────────────────────────────────────
# Panel B — power
# ─────────────────────────────────────────────────────────────────────────────
def power_two_sample(n_per_group, d, alpha=0.05):
    """Exact power of a two-sided two-sample t-test (noncentral t)."""
    n = np.asarray(n_per_group, dtype=float)
    df = 2 * n - 2
    nc = d * np.sqrt(n / 2.0)
    crit = stats.t.ppf(1 - alpha / 2, df)
    with np.errstate(invalid="ignore"):
        pw = stats.nct.sf(crit, df, nc) + stats.nct.cdf(-crit, df, nc)
    # noncentral t loses precision at very large df; fall back to the normal
    # approximation there, where the two agree to <1e-4 anyway.
    approx = stats.norm.sf(crit - nc) + stats.norm.cdf(-crit - nc)
    return np.where(np.isfinite(pw), pw, approx)


def figure_power(outdir: Path):
    n = np.logspace(np.log10(8), np.log10(3000), 500)

    effects = [
        (0.35, INCREASE, "d = 0.35", "a large-ish effect"),
        (0.20, ACCENT,   "d = 0.20", "a small effect"),
        (0.10, DECREASE, "d = 0.10", "a very small effect"),
    ]

    fig, ax = plt.subplots(figsize=(11, 6.4))
    fig.subplots_adjust(top=0.80, left=0.09, right=0.97, bottom=0.13)

    ax.axhspan(0.8, 1.0, color="#F1F7F4", zorder=0)
    ax.axhline(0.8, color=ACCENT, linewidth=1.2, linestyle=(0, (5, 4)), zorder=2)
    ax.text(8.4, 0.815, "80% power", fontsize=10, color=ACCENT,
            weight="bold", va="bottom")

    handles = []
    for d, colour, lab, sub in effects:
        ax.plot(n, power_two_sample(n, d), color=colour, linewidth=2.8,
                zorder=4, solid_capstyle="round")
        handles.append(Line2D([], [], color=colour, linewidth=2.8,
                              label=f"{lab}   {sub}"))

    cohorts = [
        (30,   "Typical prior study\n~30 per group"),
        (285,  "This study, clinical cohort\n285 vs 285"),
        (2180, "UK Biobank\n2,180 vs 2,176"),
    ]
    for x, lab in cohorts:
        ax.axvline(x, color=MUTED, linewidth=1.0, linestyle=":", zorder=1)
        ax.text(x, 1.03, lab, fontsize=10, color=INK, ha="center",
                va="bottom", linespacing=1.4, transform=ax.get_xaxis_transform())

    # what a typical prior study could actually detect
    for d, colour, lab, _ in effects:
        p30 = power_two_sample(np.array([30.0]), d)[0]
        ax.scatter([30], [p30], s=78, facecolor=colour, edgecolor="white",
                   linewidth=1.7, zorder=6)
        ax.annotate(f"{p30*100:.0f}%", xy=(30, p30), xytext=(35, p30),
                    fontsize=11, color=colour, weight="bold", va="center")

    ax.set_xscale("log")
    ax.set_xlim(8, 3000)
    ax.set_ylim(0, 1.0)
    ax.set_xticks([10, 30, 100, 285, 1000, 2180])
    ax.set_xticklabels(["10", "30", "100", "285", "1,000", "2,180"])
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_yticklabels([f"{int(v*100)}%" for v in np.arange(0, 1.01, 0.2)])
    ax.set_xlabel("Participants per group  (log scale)", fontsize=11.5, labelpad=9)
    ax.set_ylabel("Probability of detecting the effect", fontsize=11.5, labelpad=9)
    ax.grid(color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    tidy(ax)

    ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=10.5,
              framealpha=0.96, edgecolor=GRID, borderpad=0.9,
              handlelength=1.8, handletextpad=0.9, labelspacing=0.7)

    fig.text(0.09, 0.945, "A literature of underpowered studies cannot converge",
             fontsize=16, weight="bold", color=INK, ha="left", va="bottom")
    fig.text(0.09, 0.915,
             "Two-sided two-sample t-test, α = 0.05, at effect sizes typical of structural brain differences in tinnitus",
             fontsize=10.5, color=MUTED, ha="left", va="bottom")

    save(fig, outdir, "slide02b_power_curves")

    # numbers for the speaker notes
    print("\n  power at n = 30 per group:")
    for d, _, lab, _ in effects:
        print(f"    {lab}: {power_two_sample(np.array([30.0]), d)[0]*100:5.1f}%")
    print("  n per group needed for 80% power:")
    grid = np.arange(5, 20001)
    for d, _, lab, _ in effects:
        pw = power_two_sample(grid, d)
        idx = np.argmax(pw >= 0.8)
        print(f"    {lab}: {grid[idx]:,}")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--csv", type=Path, default=LIT_CSV)
    args = ap.parse_args()

    print(f"output → {args.outdir}")
    rows = read_literature(args.csv)
    if not all(r["verified"] for r in rows):
        n_bad = sum(1 for r in rows if not r["verified"])
        print(f"  ⚠  {n_bad}/{len(rows)} rows are unverified — panel A gets a DRAFT watermark.")
        print(f"     Edit {args.csv}, then set verified=yes per row.")
    figure_literature(rows, args.outdir)
    figure_power(args.outdir)
    print("\ndone.")


if __name__ == "__main__":
    main()
