"""FreeSurfer outputs -> tidy feature tables (one row per subject, one column per feature).

Column naming: <measure>__<hemi>_<region> for cortical parcellations (measure: thick, area, vol),
<structure> for aseg / subsegmentation volumes (mm^3). Global measures (eTIV, total GM, ...) are in
aseg_global. Readers only need a subject directory, so they work on the cluster and the SSD alike.

Auditory ROIs (exploratory, assembled from the tables above):
  Heschl's gyrus  = Destrieux G_temp_sup-G_T_transv;  planum temporale = Destrieux G_temp_sup-Plan_tempo
  HCP-MMP1 auditory: A1, LBelt, MBelt, PBelt, RI, A4, A5;  MGN = thalamic nuclei MGN (hemisphere-specific)
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

# cortical parcellations: stats file stem -> short atlas key
CORTICAL = {
    "aparc": "dk",
    "aparc.a2009s": "destrieux",
    "aparc.DKTatlas": "dkt",
    "Schaefer2018_400Parcels_7Networks_order": "schaefer400",
    "HCPMMP1": "hcpmmp",
}
MEASURES = {"ThickAvg": "thick", "SurfArea": "area", "GrayVol": "vol"}
# global measures are stored under the short name in the 2nd field of '# Measure <long>, <short>, ...'
ASEG_GLOBAL = {"eTIV", "BrainSegVol", "BrainSegVolNotVent", "TotalGrayVol", "CortexVol", "SubCortGrayVol",
               "CerebralWhiteMatterVol", "lhCortexVol", "rhCortexVol", "SupraTentorialVol",
               "VentricleChoroidVol", "MaskVol"}
SUBSEG = {
    "hippo": ("{h}.hippoSfVolumes-T1.v22.txt", True),
    "amyg": ("{h}.amygNucVolumes-T1.v22.txt", True),
    "thal": ("ThalamicNuclei.v13.T1.volumes.txt", False),
}
AUDITORY = {
    "HG": ("destrieux", "G_temp_sup-G_T_transv"),
    "PT": ("destrieux", "G_temp_sup-Plan_tempo"),
    "TTS": ("destrieux", "S_temporal_transverse"),
    "A1": ("hcpmmp", "A1"), "LBelt": ("hcpmmp", "LBelt"), "MBelt": ("hcpmmp", "MBelt"),
    "PBelt": ("hcpmmp", "PBelt"), "RI": ("hcpmmp", "RI"), "A4": ("hcpmmp", "A4"), "A5": ("hcpmmp", "A5"),
}


def _clean(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", name)


def read_aseg(subject_dir: Path) -> tuple[dict, dict]:
    """(structure volumes, global measures) from stats/aseg.stats."""
    f = Path(subject_dir) / "stats" / "aseg.stats"
    vols, glob = {}, {}
    if not f.exists():
        return vols, glob
    bv = Path(subject_dir) / "stats" / "brainvol.stats"  # FS 8 keeps TotalGray/Cortex/... here
    lines = f.read_text().splitlines() + (bv.read_text().splitlines() if bv.exists() else [])
    for line in lines:
        if line.startswith("# Measure"):
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 4 and parts[1] in ASEG_GLOBAL:
                glob[parts[1]] = float(parts[3])
        elif not line.startswith("#") and len(line.split()) >= 5:
            p = line.split()
            vols[_clean(p[4])] = float(p[3])  # columns: Index SegId NVoxels Volume_mm3 StructName ...
    return vols, glob


def read_cortical(subject_dir: Path, stem: str) -> dict:
    """measure__hemi_region from ?h.<stem>.stats (ThickAvg, SurfArea, GrayVol)."""
    out = {}
    for h in ("lh", "rh"):
        f = Path(subject_dir) / "stats" / f"{h}.{stem}.stats"
        if not f.exists():
            continue
        cols = None
        for line in f.read_text().splitlines():
            if line.startswith("# ColHeaders"):
                cols = line.split()[2:]
            elif line.startswith("# Measure Cortex, MeanThickness"):
                out[f"thick__{h}_MeanThickness"] = float(line.split(",")[3])
            elif line.startswith("# Measure Cortex, WhiteSurfArea"):
                out[f"area__{h}_WhiteSurfArea"] = float(line.split(",")[3])
            elif cols and not line.startswith("#") and line.strip():
                row = dict(zip(cols, line.split(), strict=False))
                name = _clean(row["StructName"])
                for src, short in MEASURES.items():
                    if src in row:
                        out[f"{short}__{h}_{name}"] = float(row[src])
    return out


def read_subseg(subject_dir: Path) -> dict:
    """hippocampal subfields, amygdala nuclei and thalamic nuclei volumes (legacy-tool outputs)."""
    out = {}
    mri = Path(subject_dir) / "mri"
    for key, (pattern, per_hemi) in SUBSEG.items():
        for h in ("lh", "rh") if per_hemi else ("",):
            f = mri / pattern.format(h=h)
            if not f.exists():
                continue
            for line in f.read_text().splitlines():
                p = line.split()
                if len(p) == 2:
                    out[f"{key}__{h + '_' if h else ''}{_clean(p[0])}"] = float(p[1])
    return out


def subject_features(subject_dir: Path) -> dict[str, dict]:
    """All tables for one subject: {table_name: {feature: value}}."""
    vols, glob = read_aseg(subject_dir)
    tables = {"aseg": vols, "aseg_global": glob, "subseg": read_subseg(subject_dir)}
    for stem, key in CORTICAL.items():
        tables[key] = read_cortical(subject_dir, stem)
    tables["auditory"] = auditory_rois(tables)
    return tables


def auditory_rois(tables: dict[str, dict]) -> dict:
    out = {}
    for roi, (atlas, label) in AUDITORY.items():
        t = tables.get(atlas, {})
        for h in ("lh", "rh"):
            hcp_label = f"{h[0].upper()}_{label}_ROI" if atlas == "hcpmmp" else label
            for meas in ("thick", "area", "vol"):
                v = t.get(f"{meas}__{h}_{hcp_label}")
                if v is not None:
                    out[f"{meas}__{h}_{roi}"] = v
    sub = tables.get("subseg", {})
    for side, h in (("Left", "lh"), ("Right", "rh")):
        mgn = sub.get(f"thal__{side}-MGN")
        if mgn is not None:
            out[f"vol__{h}_MGN"] = mgn
    return out


def fsaverage5_matrix(subject_dir: Path, measure: str, fwhm: int = 10) -> np.ndarray | None:
    """qcache map on fsaverage, restricted to the fsaverage5 vertices (icosahedral nesting: the first
    10242 vertices of fsaverage are exactly fsaverage5). Returns (20484,) lh then rh, or None."""
    import nibabel as nib

    parts = []
    for h in ("lh", "rh"):
        f = Path(subject_dir) / "surf" / f"{h}.{measure}.fwhm{fwhm}.fsaverage.mgh"
        if not f.exists():
            return None
        parts.append(np.asarray(nib.load(str(f)).get_fdata()).ravel()[:10242])
    return np.concatenate(parts).astype(np.float32)


def collect(subjects: dict[str, Path], vertex_measures=("thickness", "area", "volume"),
            fwhm: int = 10) -> tuple[dict[str, pd.DataFrame], dict[str, np.ndarray], list[str]]:
    """subjects: participant_id -> FreeSurfer subject dir. Returns (tables, vertex arrays, vertex ids)."""
    rows: dict[str, list] = {}
    vert: dict[str, list] = {m: [] for m in vertex_measures}
    vids: list[str] = []
    for pid, d in subjects.items():
        for name, feats in subject_features(d).items():
            rows.setdefault(name, []).append({"participant_id": pid, **feats})
        arrs = {m: fsaverage5_matrix(d, m, fwhm) for m in vertex_measures}
        if all(a is not None for a in arrs.values()):
            vids.append(pid)
            for m, a in arrs.items():
                vert[m].append(a)
    tables = {k: pd.DataFrame(v) for k, v in rows.items()}
    verts = {m: np.vstack(v) if v else np.empty((0, 20484), np.float32) for m, v in vert.items()}
    return tables, verts, vids
