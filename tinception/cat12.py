"""CAT12 helpers: which T1w each subject is processed from (same inputs as FreeSurfer v2), and readers for
the CAT12 xml outputs (report/cat_<sub>_T1w.xml, label/catROI_<sub>_T1w.xml) -> tidy feature rows.

Column naming: cat12_<measure> for global/quality values; <tissue>__<region> (tissue gm, wm, csf) for ROI
volumes in cm^3. Quality marks follow CAT12: 1 = A (best) ... 6 = F; rating % = 105 - 10 * mark.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from .features import clean_name

QUALITY = ["IQR", "SIQR", "NCR", "ICR", "res_RMS"]  # <qualityratings> marks (lower is better)
TISSUES = {"Vgm": "gm", "Vwm": "wm", "Vcsf": "csf"}   # <data> fields of catROI xml read as volumes


def t1_input(pid: str, site: str, t1cfg: dict, rawdata: Path, derivatives: Path) -> Path:
    """T1w file for one subject following analysis.yaml -> t1_inputs (N4 sites, acq overrides)."""
    if site in t1cfg.get("n4_sites", []):
        return derivatives / t1cfg["n4_dir"] / f"{pid}_T1w.nii.gz"
    acq = (t1cfg.get("acq") or {}).get(pid)
    return rawdata / pid / "anat" / f"{pid}_{f'acq-{acq}_' if acq else ''}T1w.nii.gz"


def _vec(text: str | None) -> np.ndarray:
    """CAT12 xml vectors look like '[1.2;3.4;5]' or '[1.2 3.4 5]'."""
    return np.array([float(v) for v in re.split(r"[;\s]+", (text or "").strip("[] \n")) if v])


def _num(parent: ET.Element | None, tag: str) -> float:
    text = parent.findtext(tag) if parent is not None else None
    return float(text) if text else np.nan


def read_report(xml: Path) -> dict:
    """Image quality marks, TIV and absolute tissue volumes (cm^3) from report/cat_*.xml. Missing fields
    (e.g. a failed subject with ignoreErrors=1 whose report holds only an error block) become NaN."""
    root = ET.parse(xml).getroot()
    q, s = root.find("qualityratings"), root.find("subjectmeasures")
    out = {f"cat12_{k}": _num(q, k) for k in QUALITY}
    out["cat12_IQR_pct"] = 105 - 10 * out["cat12_IQR"]
    out["cat12_TIV"] = _num(s, "vol_TIV")
    cgw = _vec(s.findtext("vol_abs_CGW") if s is not None else None)
    csf, gm, wm = cgw[:3] if len(cgw) >= 3 else (np.nan,) * 3
    out |= {"cat12_CSF": csf, "cat12_GM": gm, "cat12_WM": wm, "cat12_WMH": _num(s, "vol_abs_WMH")}
    return out


def read_rois(xml: Path) -> dict[str, dict]:
    """{atlas: {<tissue>__<region>: volume}} from label/catROI_*.xml (all atlases present in the file)."""
    out = {}
    for atlas in ET.parse(xml).getroot():
        names = [clean_name(n.text or "") for n in atlas.find("names")]
        if len(set(names)) != len(names):
            raise ValueError(f"{xml}: {atlas.tag} region names are not unique after cleaning")
        row = {}
        for t in atlas.find("data"):
            if t.tag not in TISSUES:
                continue
            vals = _vec(t.text)
            if len(vals) != len(names):
                raise ValueError(f"{xml}: {atlas.tag}/{t.tag} has {len(vals)} values for {len(names)} names")
            row |= {f"{TISSUES[t.tag]}__{n}": v for n, v in zip(names, vals, strict=True)}
        out[atlas.tag] = row
    return out


def collect(subjects: dict[str, Path]) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Global/quality table (output-status columns cat12_mwp1, cat12_roi) and one ROI table per atlas."""
    glob, rois = [], {}
    for pid, d in subjects.items():
        rep, roi = d / "report" / f"cat_{pid}_T1w.xml", d / "label" / f"catROI_{pid}_T1w.xml"
        row = {"participant_id": pid, "cat12_mwp1": (d / "mri" / f"mwp1{pid}_T1w.nii").is_file(),
               "cat12_roi": roi.is_file()}
        row |= read_report(rep) if rep.is_file() else {f"cat12_{k}": np.nan for k in QUALITY}
        glob.append(row)
        if roi.is_file():
            for atlas, vals in read_rois(roi).items():
                rois.setdefault(atlas, []).append({"participant_id": pid} | vals)
    return pd.DataFrame(glob), {a: pd.DataFrame(r) for a, r in rois.items()}
