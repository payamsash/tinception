"""Human-readable QC report from qc_metrics.tsv / qc_exclusions.tsv.

Writes (into --qc-dir):
  qc_overview.xlsx   sheets: README, site_summary, all_subjects, review_list
  qc_report.html     self-contained page: per-site distributions (Euler, MRIQC IQMs) + tables
"""

from __future__ import annotations

import argparse
import base64
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tinception import config

# status encoding: colour + marker shape + legend label (never colour alone)
STATUS = {
    "include": {"color": "#9a9a94", "marker": "o", "size": 14, "label": "include"},
    "flag": {"color": "#e09a00", "marker": "^", "size": 46, "label": "flag (visual review)"},
    "exclude": {"color": "#d03b3b", "marker": "X", "size": 60, "label": "exclude"},
}
INK, MUTED, GRID = "#2b2b28", "#6b6b66", "#e6e6e1"
METRICS = [
    ("euler_total", "FreeSurfer Euler number (lh+rh)", "lower = more topological defects"),
    ("cjv", "MRIQC CJV", "higher = worse (noise / bias / motion)"),
    ("cnr", "MRIQC CNR", "lower = worse GM/WM contrast"),
    ("snr_total", "MRIQC SNR", "lower = worse"),
    ("efc", "MRIQC EFC", "higher = more ghosting / blurring"),
    ("qi_2", "MRIQC QI2", "higher = more artefacts in background"),
]
PHENO_COLS = ["participant_id", "site", "group", "group4", "sex", "age", "PTA4", "legacy_drop"]


def strip_plot(df: pd.DataFrame, col: str, title: str, note: str, sites: list[str]) -> str:
    fig, ax = plt.subplots(figsize=(11, 3.8), dpi=130)
    rng = np.random.default_rng(0)
    for i, site in enumerate(sites):
        d = df[df.site == site]
        med = d[col].median()
        if np.isfinite(med):
            ax.plot([i - 0.35, i + 0.35], [med, med], color=INK, lw=2, solid_capstyle="round", zorder=1)
        for status in ("include", "flag", "exclude"):
            s = d[d.qc_status == status]
            if s.empty:
                continue
            st = STATUS[status]
            ax.scatter(i + rng.uniform(-0.28, 0.28, len(s)), s[col], s=st["size"], c=st["color"],
                       marker=st["marker"], edgecolors="white", linewidths=0.8, zorder=2 if status == "include" else 3)
    ax.set_xticks(range(len(sites)), [f"{s}\n(n={int((df.site == s).sum())})" for s in sites], fontsize=8, color=INK,
                  rotation=30, ha="right", rotation_mode="anchor")
    ax.set_title(f"{title}  —  {note}", loc="left", fontsize=10, color=INK)
    ax.tick_params(axis="y", labelsize=8, colors=MUTED)
    ax.grid(axis="y", color=GRID, lw=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    handles = [plt.Line2D([], [], ls="", marker=v["marker"], color=v["color"], markersize=7, label=v["label"])
               for v in STATUS.values()]
    handles.append(plt.Line2D([], [], color=INK, lw=2, label="site median"))
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", type=Path, default=Path("qc/final"))
    ap.add_argument("--participants", type=Path, default=None)
    args = ap.parse_args()

    m = pd.read_csv(args.qc_dir / "qc_metrics.tsv", sep="\t", na_values="n/a")
    e = pd.read_csv(args.qc_dir / "qc_exclusions.tsv", sep="\t", keep_default_na=False)
    part = pd.read_csv(args.participants or config.paths().rawdata / "participants.tsv", sep="\t", na_values="n/a")
    m = m.drop(columns=[c for c in ("group", "legacy_drop") if c in m]).merge(
        part[PHENO_COLS], on=["participant_id", "site"], how="left", validate="one_to_one")
    m = m.merge(e[["participant_id", "qc_status", "manual_decision", "manual_note", "reviewer"]],
                on="participant_id", how="left", validate="one_to_one")
    sites = sorted(m.site.unique())

    summ = m.groupby("site").agg(
        n=("participant_id", "size"), include=("qc_status", lambda s: (s == "include").sum()),
        flag=("qc_status", lambda s: (s == "flag").sum()), exclude=("qc_status", lambda s: (s == "exclude").sum()),
        euler_median=("euler_total", "median"), euler_min=("euler_total", "min"),
        cjv_median=("cjv", "median"), cnr_median=("cnr", "median"), snr_median=("snr_total", "median"))
    summ.loc["ALL"] = [len(m), (m.qc_status == "include").sum(), (m.qc_status == "flag").sum(),
                       (m.qc_status == "exclude").sum(), m.euler_total.median(), m.euler_total.min(),
                       m.cjv.median(), m.cnr.median(), m.snr_total.median()]
    summ = summ.round(3).reset_index()

    cols = (PHENO_COLS[:7] + ["qc_status", "qc_auto", "qc_reasons", "euler_lh", "euler_rh", "euler_total", "euler_z",
            "cjv", "cnr", "snr_total", "efc", "fber", "qi_2", "n_iqm_outliers", "legacy_drop",
            "manual_decision", "manual_note", "reviewer"])
    allsub = m[[c for c in cols if c in m]].sort_values(["site", "participant_id"])
    review = allsub[allsub.qc_auto != "include"].sort_values(["qc_auto", "site", "participant_id"])

    readme = pd.DataFrame({"item": [
        "qc_status", "qc_auto", "flag", "exclude", "Euler", "MRIQC", "manual review"], "meaning": [
        "final status = manual_decision if filled in qc_exclusions.tsv, otherwise qc_auto",
        "automatic rule result (config/analysis.yaml -> qc)",
        "visual review needed: Euler or >=2 MRIQC metrics are site outliers (robust z beyond 3.5), or v1 DROPS subject",
        "recon-all failed, or Euler total < -434 (provisional floor)",
        "topological defects of the uncorrected surface; closer to 0 is better",
        "cjv/efc/qi_2: higher is worse; cnr/snr_total/fber: lower is worse (fber -1 = zeroed background, treated as missing)",
        "edit manual_decision (include/exclude), manual_note, reviewer in qc_exclusions.tsv (not in this workbook)"]})
    with pd.ExcelWriter(args.qc_dir / "qc_overview.xlsx") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        summ.to_excel(xw, sheet_name="site_summary", index=False)
        allsub.to_excel(xw, sheet_name="all_subjects", index=False)
        review.to_excel(xw, sheet_name="review_list", index=False)

    imgs = "".join(f'<img alt="{t}" src="data:image/png;base64,{strip_plot(m, c, t, n, sites)}">'
                   for c, t, n in METRICS if c in m)
    style = ("body{font-family:-apple-system,Helvetica,Arial,sans-serif;color:#2b2b28;background:#fcfcfb;"
             "max-width:1200px;margin:24px auto;padding:0 16px}img{width:100%;margin:8px 0}"
             "table{border-collapse:collapse;font-size:12px;margin:8px 0}td,th{border-bottom:1px solid #e6e6e1;"
             "padding:3px 8px;text-align:left}th{background:#f1f1ee}.k{color:#6b6b66;font-size:13px}")
    html = (f"<html><head><meta charset='utf-8'><title>Tinception QC</title><style>{style}</style></head><body>"
            f"<h1>Tinception structural QC</h1><p class='k'>{len(m)} subjects, {len(sites)} sites — "
            f"include {int((m.qc_status=='include').sum())}, flag {int((m.qc_status=='flag').sum())}, "
            f"exclude {int((m.qc_status=='exclude').sum())}. Marker shape and colour both encode status.</p>"
            f"<h2>Per-site summary</h2>{summ.to_html(index=False, na_rep='')}"
            f"<h2>Distributions</h2>{imgs}"
            f"<h2>Subjects to review ({len(review)})</h2>{review.to_html(index=False, na_rep='')}"
            f"</body></html>")
    (args.qc_dir / "qc_report.html").write_text(html)
    print(f"-> {args.qc_dir}/qc_overview.xlsx, qc_report.html ({len(m)} subjects, {len(review)} to review)")


if __name__ == "__main__":
    main()
