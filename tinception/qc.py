"""Image-quality control: collect metrics and apply site-wise outlier rules.

Metrics
-------
- FreeSurfer Euler number of the uncorrected surfaces (``orig.nofix``), from ``recon-all.log``.
  Strongest single predictor of unusable scans (Rosen et al. 2018, NeuroImage 169:407).
- MRIQC image-quality metrics (T1w): cjv, cnr, snr_total, efc, fber, qi_2.
- (later) CAT12 image-quality rating (IQR).

Rules (config/analysis.yaml -> qc)
----------------------------------
- exclude: recon-all error, or Euler total below an absolute floor.
- flag (visual review): Euler total or >= ``min_iqm_outliers`` IQMs beyond a robust site-wise z of
  ``mad_k`` in the "bad" direction. Robust z = (x - site median) / (1.4826 * site MAD).
Manual decisions in an existing qc_exclusions.tsv are preserved.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

EULER_RE = re.compile(r"orig\.nofix lheno\s*=\s*(-?\d+),\s*rheno\s*=\s*(-?\d+)")
# direction in which a metric is bad: +1 = high is bad, -1 = low is bad
IQM_BAD_DIRECTION = {"cjv": +1, "efc": +1, "qi_2": +1, "cnr": -1, "snr_total": -1, "fber": -1}


def euler_from_subject(subject_dir: Path) -> dict:
    """Euler numbers from recon-all.log (last occurrence) plus FreeSurfer run status."""
    scripts = Path(subject_dir) / "scripts"
    out = {"euler_lh": np.nan, "euler_rh": np.nan,
           "fs_done": (scripts / "recon-all.done").exists(),
           "fs_error": (scripts / "recon-all.error").exists()}
    log = scripts / "recon-all.log"
    if log.exists():
        hits = EULER_RE.findall(log.read_text(errors="ignore"))
        if hits:
            out["euler_lh"], out["euler_rh"] = map(int, hits[-1])
    out["euler_total"] = out["euler_lh"] + out["euler_rh"]
    return out


def collect_euler(fs_root: Path, subjects: dict[str, str]) -> pd.DataFrame:
    """subjects: participant_id -> subject directory name relative to fs_root."""
    rows = []
    for pid, rel in subjects.items():
        d = Path(fs_root) / rel
        rows.append({"participant_id": pid, "fs_path": str(d), **(euler_from_subject(d) if d.is_dir() else {})})
    return pd.DataFrame(rows)


def collect_mriqc(mriqc_dir: Path) -> pd.DataFrame:
    """One row per subject from MRIQC T1w IQM json files."""
    rows = []
    for f in sorted(Path(mriqc_dir).rglob("sub-*_T1w.json")):
        pid = f.name.split("_")[0]
        if f.name != f"{pid}_T1w.json":  # skip extra runs (e.g. sub-audicog16_run-2_T1w); primary T1w only
            continue
        js = json.loads(f.read_text())
        rows.append({"participant_id": pid, **{k: js.get(k) for k in IQM_BAD_DIRECTION}})
    df = pd.DataFrame(rows, columns=["participant_id", *IQM_BAD_DIRECTION])
    # MRIQC returns fber = -1 when the background is zeroed (defaced / pre-masked images): not a measurement
    df["fber"] = pd.to_numeric(df["fber"], errors="coerce").where(lambda v: v > 0)
    return df


def robust_z(x: pd.Series, groups: pd.Series) -> pd.Series:
    def _z(v: pd.Series) -> pd.Series:
        med = v.median()
        mad = 1.4826 * (v - med).abs().median()
        return (v - med) / mad if mad and np.isfinite(mad) else v * np.nan
    return x.groupby(groups).transform(_z)


def apply_rules(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Add robust z-scores, reasons and an automatic status (include / flag / exclude)."""
    df = df.copy()
    k = cfg["mad_k"]
    reasons: list[list[str]] = [[] for _ in range(len(df))]

    def add(mask: pd.Series, text: str) -> None:
        for i in np.flatnonzero(mask.fillna(False).to_numpy()):
            reasons[i].append(text)

    exclude = pd.Series(False, index=df.index)
    if "fs_error" in df:
        e = df["fs_error"].fillna(False).astype(bool)
        add(e, "recon-all error")
        exclude |= e
    if "euler_total" in df:
        floor = df["euler_total"] < cfg["euler_abs_floor"]
        add(floor, f"Euler<{cfg['euler_abs_floor']}")
        exclude |= floor.fillna(False)
        df["euler_z"] = robust_z(df["euler_total"], df["site"])
        add(df["euler_z"] < -k, f"Euler site-z<-{k}")

    n_iqm = pd.Series(0, index=df.index)
    for m, sign in IQM_BAD_DIRECTION.items():
        if m in df and df[m].notna().any():
            df[f"{m}_z"] = robust_z(pd.to_numeric(df[m], errors="coerce"), df["site"])
            bad = (sign * df[f"{m}_z"]) > k
            n_iqm += bad.fillna(False).astype(int)
            add(bad, f"{m} site-z {'>' if sign > 0 else '<-'}{k}")
    df["n_iqm_outliers"] = n_iqm

    flag = (df.get("euler_z", pd.Series(np.nan, index=df.index)) < -k).fillna(False)
    flag |= n_iqm >= cfg["min_iqm_outliers"]
    if "legacy_drop" in df:
        ld = df["legacy_drop"].fillna(False).astype(bool)
        add(ld, "excluded in v1 (reason unrecorded)")
        flag |= ld
    df["qc_reasons"] = ["; ".join(r) for r in reasons]
    df["qc_auto"] = np.where(exclude, "exclude", np.where(flag, "flag", "include"))
    return df


def merge_manual(new: pd.DataFrame, previous: pd.DataFrame | None) -> pd.DataFrame:
    """Keep manual review columns from a previous qc_exclusions.tsv; final = manual if set, else auto."""
    new = new.copy()
    for c in ("manual_decision", "manual_note", "reviewer"):
        new[c] = ""
    if previous is not None and len(previous):
        prev = previous.set_index("participant_id")
        for c in ("manual_decision", "manual_note", "reviewer"):
            if c in prev:
                new[c] = new["participant_id"].map(prev[c]).fillna("").astype(str)
    manual = new["manual_decision"].str.strip().str.lower()
    new["qc_status"] = np.where(manual.isin(["include", "exclude"]), manual, new["qc_auto"])
    return new
