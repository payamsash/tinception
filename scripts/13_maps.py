"""Vertex-wise and VBM group maps with PALM (analysis plan 6.6; user decisions D1-D6, 2026-10-31).

Each analysis (analysis.yaml maps_palm.analyses) is the H1 main model at every vertex / voxel: tinnitus + PTA4 +
PTA_HF beyond PTA4 + age spline + sex + Euler + site dummies [+ log eTIV for area, log CAT12 TIV for VBM]; the
hyperacusis map is the within-patient model (hyperacusis_z + PTA4 + covariates). PALM: TFCE (2D on the fsaverage5
surface, both hemispheres as two modalities with FWE across them; 3D for VBM), permutations within site
(exchangeability blocks, variance groups per site, z statistics), two-tailed FWE.

  prepare   -> <derivatives>/maps_palm/<analysis>/{design.csv, contrasts.csv, eb.csv, inputs, masks}, tasks.tsv
  fit --task i (1-based; SLURM array) -> <analysis>/palm_*
  collect   -> <results>/13_maps/{summary.tsv, summary.json}

Example (cluster, repo root):  .venv/bin/python scripts/13_maps.py prepare && sbatch --array=1-8 slurm/13_maps.sbatch
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import group as g

ID = "participant_id"
HEMIS = ("lh", "rh")
NV = 10242   # fsaverage5 vertices per hemisphere


def subjects_table(master: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Rows of one analysis sample with the design variables (PTA_HF_res is added in cmd_prepare, on the rows
    finally analysed)."""
    m = master.copy()
    if sample in ("primary", "clean"):
        m = m[m.in_case_control.astype(bool) & m.PTA_HF.notna()]
        if sample == "clean":
            m = m[m.group4_hf.isin(["CTRL", "TINN"])]
        m = m.assign(tinnitus=m.group.eq("TI").astype(int)).reset_index(drop=True)
    elif sample == "hyperacusis":
        m = m[m.in_patients.astype(bool) & m.hyperacusis_z.notna()].reset_index(drop=True)
    else:
        raise ValueError(sample)
    return m.assign(log_eTIV=np.log(m.eTIV), log_TIV=np.log(m.cat12_TIV))


def design(d: pd.DataFrame, terms: list[str], extra: list[str], age_df: int) -> pd.DataFrame:
    """Design matrix (intercept, terms, age spline, sex, Euler, extra, site dummies); continuous columns standardized."""
    import patsy

    d = d.copy()
    for c in [*terms, "euler_total", *extra]:
        if c != "tinnitus":
            d[c] = (d[c] - d[c].mean()) / d[c].std()
    rhs = " + ".join([*terms, f"cr(age, df={age_df}, constraints='center')", "C(sex)", "euler_total", *extra,
                      "C(site)"])
    return patsy.dmatrix(rhs, d, return_type="dataframe")


def _mgh(data: np.ndarray, path: Path) -> None:
    """vertices x subjects (or vertices) -> MGH (vertices, 1, 1, subjects)."""
    a = data.astype(np.float32)
    a = a.reshape(a.shape[0], 1, 1) if a.ndim == 1 else a.reshape(a.shape[0], 1, 1, -1)
    nib.save(nib.MGHImage(a, np.eye(4)), str(path))


def cmd_prepare(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    mcfg, age_df = cfg["maps_palm"], cfg["group"]["age_df"]
    root = p.derivatives / "maps_palm"
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    feat = p.derivatives / "features"
    vsub = (feat / "vertex_subjects.txt").read_text().split()
    vpos = {s: i for i, s in enumerate(vsub)}
    vbm_cache = {}   # fwhm -> (subject ids, 4D stack, affine): the clean sample reuses the main stack
    keep = np.zeros(2 * NV, bool)
    keep[np.load(p.derivatives / "normative" / "maps" / "data" / "vertex_thickness_index.npy")] = True  # 07b mask
    fs5 = p.freesurfer / "fsaverage5" / "surf"
    # PALM recognizes FreeSurfer surfaces by extension only (.white, .pial, ...): copy the average white surface
    srf = {h: root / "surf" / f"{h}.white_avg.white" for h in HEMIS}
    (root / "surf").mkdir(parents=True, exist_ok=True)
    for h in HEMIS:
        if not srf[h].exists():
            shutil.copy(fs5 / f"{h}.white_avg", srf[h])
    only = set(args.only.split(",")) if args.only else None
    old = None
    if only:   # partial rebuild: every other analysis must already be prepared; task numbers stay in config order
        if not (root / "tasks.tsv").exists():
            raise SystemExit("--only needs an existing tasks.tsv (run a full prepare first)")
        old = pd.read_csv(root / "tasks.tsv", sep="\t")
        missing = set(mcfg["analyses"]) - only - set(old.analysis)
        if missing:
            raise SystemExit(f"--only: analyses never prepared: {sorted(missing)}")
    tasks = []
    for name, a in mcfg["analyses"].items():
        if only and name not in only:   # keep the previous task row; inputs untouched
            tasks.append(old[old.analysis == name].iloc[0].to_dict() | {"task": len(tasks) + 1})
            continue
        out = root / name
        out.mkdir(parents=True, exist_ok=True)
        d = subjects_table(master, a["sample"])
        mod = a["modality"]
        if mod in ("thickness", "area"):
            d = d[d[ID].isin(vsub)]
        else:
            d = d[d.vbm_ok.astype(bool) & d.cat12_TIV.notna()]
        hyper = a["sample"] == "hyperacusis"
        terms = ["hyperacusis_z", "PTA4"] if hyper else ["tinnitus", "PTA4", "PTA_HF_res"]
        extra = {"thickness": [], "area": ["log_eTIV"]}.get(mod, ["log_TIV"])
        need = [t for t in terms if t != "PTA_HF_res"] + (["PTA_HF"] if not hyper else [])
        d = d.dropna(subset=[*need, "age", "sex", "euler_total", "site", *extra]).reset_index(drop=True)
        if not hyper:   # PTA_HF beyond PTA4, residualized on exactly the analysed rows (as 08_group)
            d = g.residualize(d, "PTA_HF", ["PTA4"], "PTA_HF_res")
        small = d.site.value_counts()
        if (small < mcfg["min_site_n"]).any():   # each site is its own variance group (-vg auto)
            raise SystemExit(f"{name}: sites below {mcfg['min_site_n']} subjects: {small[small < mcfg['min_site_n']].to_dict()}")
        for f in out.glob("palm_*"):   # outputs of an earlier design must not be collected with this one
            f.unlink()
        X = design(d, terms, extra, age_df)
        X.to_csv(out / "design.csv", header=False, index=False)
        (out / "design_columns.txt").write_text("\n".join(X.columns))
        C = np.zeros((len(a["contrasts"]), X.shape[1]))
        for i, c in enumerate(a["contrasts"]):
            C[i, list(X.columns).index(c)] = 1.0
        pd.DataFrame(C).to_csv(out / "contrasts.csv", header=False, index=False)
        pd.Series(pd.factorize(d.site)[0] + 1).to_csv(out / "eb.csv", header=False, index=False)
        d[[ID, "site"]].to_csv(out / "subjects.tsv", sep="\t", index=False)
        vbm = mod.startswith("vbm")
        n_perm, tail = (mcfg["vbm_n_perm"], mcfg["vbm_accel_tail"]) if vbm else (mcfg["n_perm"], mcfg["accel_tail"])
        common = ["-d", "design.csv", "-t", "contrasts.csv", "-eb", "eb.csv", "-within", "-vg", "auto", "-zstat",
                  "-T", "-n", str(n_perm), "-twotail", "-o", "palm"] + (["-accel", "tail"] if tail else [])
        if mod in ("thickness", "area"):
            X_all = np.load(feat / f"vertex_fsavg5_fwhm10_{mod}.npy", mmap_mode="r")
            x = np.asarray(X_all[[vpos[s] for s in d[ID]]], np.float64)
            x[:, ~keep] = 0.0
            if mod == "area":   # positive inside the mask (07b vertex_keep)
                x[:, keep] = np.log(x[:, keep])
            if not np.isfinite(x[:, keep]).all():
                raise SystemExit(f"{name}: non-finite values inside the vertex mask")
            argv = []
            for h, sl in zip(HEMIS, (slice(0, NV), slice(NV, 2 * NV)), strict=True):
                _mgh(x[:, sl].T, out / f"{h}.mgh")
                _mgh(keep[sl].astype(np.float32), out / f"{h}_mask.mgh")
                argv += ["-i", f"{h}.mgh", "-s", str(srf[h]), str(fs5 / f"{h}.white.avg.area.mgh"),
                         "-m", f"{h}_mask.mgh"]
            argv += ["-tfce2D", "-corrmod", *common]
            n_feat = int(keep.sum())
        else:
            fwhm = mod.replace("vbm", "")
            ids = list(d[ID])
            cached = vbm_cache.get(fwhm)
            if cached and set(ids) <= set(cached[0]):
                pos = {s: i for i, s in enumerate(cached[0])}
                vol, affine = cached[1][..., [pos[s] for s in ids]], cached[2]
            else:
                files = [p.derivatives / "cat12-12.9" / s / "mri" / f"s{fwhm}mwp1{s}_T1w.nii" for s in ids]
                absent = [str(f) for f in files if not f.exists()]
                if absent:
                    raise SystemExit(f"{name}: {len(absent)} CAT12 images missing, e.g. {absent[:3]}")
                imgs = [nib.load(str(f)) for f in files]
                affine = imgs[0].affine
                bad = [str(f) for f, im in zip(files, imgs, strict=True) if not np.allclose(im.affine, affine)]
                if bad:
                    raise SystemExit(f"{name}: images on a different grid: {bad[:3]}")
                vol = np.stack([np.asarray(im.dataobj, np.float32) for im in imgs], axis=-1)
                vbm_cache = {fwhm: (ids, vol, affine)}
            mask = (vol.mean(axis=-1) > mcfg["vbm_mask_threshold"]) & (vol.std(axis=-1) > 0)
            nib.save(nib.Nifti1Image(vol, affine), str(out / "gm.nii"))
            nib.save(nib.Nifti1Image(mask.astype(np.uint8), affine), str(out / "mask.nii"))
            argv = ["-i", "gm.nii", "-m", "mask.nii", *common]
            n_feat = int(mask.sum())
        (out / "palm_args.json").write_text(json.dumps(argv))
        tasks.append({"task": len(tasks) + 1, "analysis": name, "n": len(d), "n_sites": int(d.site.nunique()),
                      "features": n_feat, "contrasts": "+".join(a["contrasts"])})
        print(f"{name}: n={len(d)} ({d.site.nunique()} sites), {n_feat} features, design {X.shape}", flush=True)
    pd.DataFrame(tasks).to_csv(root / "tasks.tsv", sep="\t", index=False)


def cmd_fit(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    mcfg = cfg["maps_palm"]
    root = p.derivatives / "maps_palm"
    row = pd.read_csv(root / "tasks.tsv", sep="\t").iloc[args.task - 1]
    out = root / row.analysis
    argv = json.loads((out / "palm_args.json").read_text())
    if args.n_perm:  # pilot / timing runs
        i = argv.index("-n")
        argv[i + 1] = str(args.n_perm)
        argv[argv.index("-o") + 1] = f"pilot{args.n_perm}"
    call = "palm(" + ", ".join(f"'{a}'" for a in argv) + ")"
    cmd = f"addpath('{mcfg['palm']}'); cd('{out}'); {call}; exit"
    print(f"{row.analysis}: matlab -batch \"{cmd}\"", flush=True)
    subprocess.run(["bash", "-lc", f"module load {mcfg['matlab_module']} && matlab -batch \"{cmd}\""], check=True)


def cmd_collect(args: argparse.Namespace) -> None:
    """Per analysis x contrast (x hemisphere): minimum TFCE FWE p (vertex: corrected across both hemispheres,
    PALM 'mfwep'; VBM: 'fwep') and the number of significant vertices / voxels at p < .05, inside the mask."""
    import re

    p = config.paths()
    root = p.derivatives / "maps_palm"
    out = p.results / "13_maps"
    out.mkdir(parents=True, exist_ok=True)
    tasks = pd.read_csv(root / "tasks.tsv", sep="\t")
    rows = []
    for t in tasks.itertuples():
        d = root / t.analysis
        vertex = (d / "lh_mask.mgh").exists()
        # PALM omits _c<k> when there is a single contrast (and _m<k> with a single modality)
        pat = "palm_tfce_zvstat_mfwep_m*.mgz" if vertex else "palm_tfce_zvstat_fwep*.nii"
        files = sorted(d.glob(pat))
        if not files:
            rows.append({"analysis": t.analysis, "status": "no output"})
            continue
        for f in files:
            mm = re.search(r"fwep(?:_m(\d))?(?:_c(\d))?\.", f.name)
            hemi = HEMIS[int(mm.group(1)) - 1] if mm.group(1) else "volume"
            ci = int(mm.group(2)) - 1 if mm.group(2) else 0
            mask = np.asarray(nib.load(str(d / (f"{hemi}_mask.mgh" if vertex else "mask.nii"))).dataobj).ravel() > 0
            v = np.asarray(nib.load(str(f)).dataobj, float).ravel()[mask]
            rows.append({"analysis": t.analysis, "contrast": t.contrasts.split("+")[ci],
                         "hemisphere": hemi, "n": t.n, "min_fwep": float(v.min()), "n_sig_05": int((v < 0.05).sum()),
                         "n_tested": int(v.size), "file": f.name, "status": "ok"})
    s = pd.DataFrame(rows)
    s.to_csv(out / "summary.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps(g.json_clean(s.to_dict("records")), indent=2, allow_nan=False))
    print(s.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("prepare")
    pr.add_argument("--only", default=None, help="comma-separated analyses to (re)build; others keep their inputs")
    pr.set_defaults(func=cmd_prepare)
    f = sub.add_parser("fit")
    f.add_argument("--task", type=int, required=True, help="1-based row of tasks.tsv (SLURM_ARRAY_TASK_ID)")
    f.add_argument("--n-perm", type=int, default=None, help="override for pilot/timing runs (outputs pilot<n>_*)")
    f.set_defaults(func=cmd_fit)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
