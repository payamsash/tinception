"""CAT12 helpers: which T1w each subject is processed from (same inputs as FreeSurfer v2)."""

from __future__ import annotations

from pathlib import Path


def t1_input(pid: str, site: str, t1cfg: dict, rawdata: Path, derivatives: Path) -> Path:
    """T1w file for one subject following analysis.yaml -> t1_inputs (N4 sites, acq overrides)."""
    if site in t1cfg.get("n4_sites", []):
        return derivatives / t1cfg["n4_dir"] / f"{pid}_T1w.nii.gz"
    acq = (t1cfg.get("acq") or {}).get(pid)
    return rawdata / pid / "anat" / f"{pid}_{f'acq-{acq}_' if acq else ''}T1w.nii.gz"
