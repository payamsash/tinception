"""Assemble the unified BIDS dataset (rawdata/) from all image sources.

One copy per image. Subject labels come from tinception.ids. The ID map records, for every
subject, where the image came from and which legacy FreeSurfer directory belongs to it.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from . import config, ids

# DROPS holds subjects excluded in v1 without a recorded reason; route them back to their site.
DROPS_PATTERNS = [
    (re.compile(r"^\d{4}$"), "tinmeg"),
    (re.compile(r"__?(TI|NT)-(HL|HA)$"), "cogtail"),
    (re.compile(r"^ICCAC_"), "iccac"),
    (re.compile(r"^(K\d+_controls|Triple_\d+_patients)$"), "triple"),
    (re.compile(r"^MP\d+_FIL$"), "london_2"),
    (re.compile(r"^indiTMS_"), "inditms"),
]


@dataclass
class ImageSource:
    site: str
    original_id: str
    suffix: str                      # e.g. T1w, T2w, acq-highreshippo_T2w, dwi, task-rest_run-1_bold
    datatype: str                    # anat | dwi | func
    path: Path
    sidecars: list[Path] = field(default_factory=list)
    legacy_drop: bool = False
    fs_dir: str | None = None

    @property
    def participant_id(self) -> str:
        return ids.bids_id(self.site, self.original_id)

    def target(self, rawdata: Path) -> Path:
        sub = self.participant_id
        return rawdata / sub / self.datatype / f"{sub}_{self.suffix}.nii.gz"


def drops_site(original_id: str) -> str | None:
    for pat, site in DROPS_PATTERNS:
        if pat.search(original_id):
            return site
    return None


def _stem(p: Path) -> str:
    return p.name.removesuffix(".nii.gz")


def legacy_sources(fs_root: Path | None = None) -> list[ImageSource]:
    """T1/T2 originals in vbm_per_site/<site>/ (+ DROPS, andes_t2); antinomics comes from its BIDS."""
    root = config.paths().sources["legacy_t1"]
    fs_root = fs_root or config.paths().freesurfer
    out: list[ImageSource] = []
    for site_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        name = site_dir.name
        if name in ("antinomics", "antinomics_t2"):
            continue
        for img in sorted(site_dir.glob("*.nii.gz")):
            if img.name.startswith("._"):
                continue
            oid = _stem(img)
            if name == "DROPS":
                site = drops_site(oid)
                if site is None:
                    raise ValueError(f"cannot route DROPS/{img.name} to a site")
                out.append(ImageSource(site, oid, "T1w", "anat", img, legacy_drop=True))
            elif name.endswith("_t2"):
                out.append(ImageSource(name.removesuffix("_t2"), oid, "T2w", "anat", img))
            else:
                out.append(ImageSource(name, ids.canonical_original_id(name, oid), "T1w", "anat", img))
    for s in out:
        if s.suffix == "T1w" and (fs_root / s.original_id).is_dir():
            s.fs_dir = s.original_id
    return out


def bids_sources(site: str, root: Path, suffix_map: dict[str, str | None] | None = None) -> list[ImageSource]:
    """Images from an existing BIDS dataset. The unified dataset has no sessions, so
    ``suffix_map`` renames session-specific suffixes (value None = skip the image)."""
    suffix_map = suffix_map or {}
    out = []
    for sub in sorted(p for p in root.glob("sub-*") if p.is_dir()):
        oid = sub.name.removeprefix("sub-")
        for img in sorted(sub.rglob("*.nii.gz")):
            if img.name.startswith("._"):
                continue
            datatype = img.parent.name
            suffix = _stem(img).removeprefix(f"{sub.name}_")
            if suffix in suffix_map:
                suffix = suffix_map[suffix]
                if suffix is None:
                    continue
            elif suffix.startswith("ses-"):
                raise ValueError(f"unmapped session file {img}")
            side = [img.with_name(_stem(img) + ext) for ext in (".json", ".bvec", ".bval")]
            out.append(ImageSource(site, oid, suffix, datatype, img, [s for s in side if s.exists()]))
    return out


AUDICOG_DIR = re.compile(r"_AUDICOG_Sujet(\d+)$")
AUDICOG_T1 = re.compile(r"^S(\d+)_T1w$")


def audicog_sources(roots: list[Path]) -> list[ImageSource]:
    """AUDICOG (Hobeika, Paris CENIR Prisma): raw per-series folders <date>_AUDICOG_Sujet<n>/S<k>_T1w/v_*.nii.
    Only T1w is used. If a subject has several T1w series, the first is `T1w`, later ones `run-<i>_T1w`."""
    out = []
    for root in roots:
        for subj in sorted(p for p in Path(root).glob("*_AUDICOG_Sujet*") if p.is_dir()):
            m = AUDICOG_DIR.search(subj.name)
            if not m:
                continue
            series = sorted((int(AUDICOG_T1.match(d.name)[1]), d) for d in subj.iterdir()
                            if d.is_dir() and AUDICOG_T1.match(d.name))
            for i, (_, d) in enumerate(series, start=1):
                nii = [f for f in d.glob("v_*.nii") if not f.name.startswith("._")]
                if len(nii) != 1:
                    raise ValueError(f"expected one v_*.nii in {d}, found {len(nii)}")
                side = nii[0].with_suffix(".json")
                out.append(ImageSource("audicog", m[1], "T1w" if i == 1 else f"run-{i}_T1w", "anat",
                                       nii[0], [side] if side.exists() else []))
    return out


def all_sources(include: set[str] | None = None) -> list[ImageSource]:
    src = config.paths().sources
    items = legacy_sources()
    fs_root = config.paths().freesurfer
    ant = bids_sources("antinomics", src["antinomics_bids"])
    for s in ant:
        if s.suffix == "T1w" and (fs_root / "antinomics_subjects" / s.original_id).is_dir():
            s.fs_dir = f"antinomics_subjects/{s.original_id}"
    items += ant
    # Husain datasets: T1w only at ses-A; T2w is a low-res 2D TSE overlay (ses-B repeat dropped)
    husain = {"ses-A_T1w": "T1w", "ses-A_T2w": "acq-lowres2d_T2w", "ses-B_T2w": None}
    items += bids_sources("uiuc", src["uiuc"], husain)
    items += bids_sources("whasc", src["whasc"], husain)
    items += audicog_sources([src[k] for k in ("audicog_first25", "audicog_last25") if k in src])
    if include:
        items = [s for s in items if s.datatype in include]
    seen: dict[Path, ImageSource] = {}
    for s in items:
        t = s.target(Path("/"))
        if t in seen:
            raise ValueError(f"target collision: {s.path} and {seen[t].path} -> {t}")
        seen[t] = s
    return items


def sidecar_from_registry(site: str) -> dict:
    meta = config.sites()[site]
    js = {"Manufacturer": meta.get("vendor"), "MagneticFieldStrength": meta.get("field_T"),
          "PulseSequenceType": meta.get("sequence"),
          "TinceptionNote": "Reconstructed from config/sites.yaml; no original DICOM sidecar available."}
    if meta.get("model"):
        js["ManufacturersModelName"] = meta["model"]
    return {k: v for k, v in js.items() if v is not None}


def _copy_data(src: Path, dst: Path) -> None:
    """Copy contents + mtime only: copying xattrs makes macOS write '._*' files on exFAT.
    An uncompressed .nii source written to a .nii.gz target is gzip-compressed."""
    if src.name.endswith(".nii") and dst.name.removesuffix(".part").endswith(".nii.gz"):
        with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
    else:
        shutil.copyfile(src, dst)
    st = src.stat()
    os.utime(dst, (st.st_atime, st.st_mtime))


def copy_source(s: ImageSource, rawdata: Path, overwrite: bool = False) -> Path:
    dst = s.target(rawdata)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if overwrite or not dst.exists():
        tmp = dst.with_name(dst.name + ".part")  # atomic: an interrupted copy never looks complete
        _copy_data(s.path, tmp)
        tmp.rename(dst)
    base = str(dst).removesuffix(".nii.gz")
    if s.sidecars:
        for sc in s.sidecars:
            ext = "".join(sc.suffixes[-1:])
            _copy_data(sc, Path(base + ext))
    elif not Path(base + ".json").exists():
        with open(base + ".json", "w") as f:
            json.dump(sidecar_from_registry(s.site), f, indent=2)
    return dst


def dataset_description() -> dict:
    return {
        "Name": "Tinception: multi-site structural MRI in chronic tinnitus",
        "BIDSVersion": "1.9.0",
        "DatasetType": "raw",
        "GeneratedBy": [{"Name": "tinception.bids", "Description": "Merged from site datasets; see code/id_map.tsv"}],
    }
