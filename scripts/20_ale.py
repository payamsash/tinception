"""Coordinate-based meta-analysis of published tinnitus grey-matter morphometry (analysis plan 6.10 B2; details fixed
2026-11-04 before the paywalled coordinates were extracted).

Input: data/meta/coordinates_final.tsv (reconciled double extraction; one row per focus) with columns
  experiment_id (one per independent subject group, Turkeltaub 2012), pmid, study, modality, direction_tinnitus (+/-),
  space (MNI/Talairach), correction (corrected/uncorrected/svc), x, y, z, n_total
Conditions: primary (all foci, any direction); corrected only; VBM only; decreases; increases (each run only with
>= ale.min_experiments experiments, otherwise reported as underpowered). ALE (NiMARE; sample size = n_total;
Talairach -> MNI by NiMARE's Lancaster transform), cluster-forming p < .001, cluster-level FWE p < .05 (size), Monte
Carlo iterations from analysis.yaml. Per significant cluster: jackknife contribution of each experiment and the
fail-safe N (Acar et al. 2018: null experiments with random foci in the ALE mask added until the cluster is lost).
-> <results>/20_ale/{summary.json, clusters_<cond>.tsv, jackknife_<cond>.tsv, z_<cond>.nii.gz}

  run [--conditions ...] [--n-iters N]   (cluster: slurm/20_ale.sbatch)
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import stats as st

from tinception import config
from tinception.stats import group as g

TARGET = "z_desc-size_level-cluster_corr-FWE_method-montecarlo"
LOGP = "logp_desc-size_level-cluster_corr-FWE_method-montecarlo"
SPACE = {"MNI": "MNI", "TALAIRACH": "TAL", "TAL": "TAL"}
ALLOWED = {"modality": {"VBM", "surface"}, "direction_tinnitus": {"+", "-"},
           "correction": {"corrected", "uncorrected"}}
MIN_RUNNABLE = 2   # ALE needs >= 2 experiments even as a descriptive run


def load(path) -> pd.DataFrame:
    c = pd.read_csv(path, sep="\t")
    need = {"experiment_id", "pmid", "study", "modality", "direction_tinnitus", "space", "correction", "x", "y", "z",
            "n_total"}
    if need - set(c.columns):
        raise SystemExit(f"{path}: missing columns {sorted(need - set(c.columns))}")
    c["space"] = c.space.str.upper().map(SPACE)
    if c.space.isna().any() or c[["x", "y", "z", "n_total"]].isna().any().any():
        raise SystemExit("unknown space or missing coordinate / n")
    for col, ok in ALLOWED.items():   # the condition filters match these values exactly
        bad = sorted(set(c[col].astype(str)) - ok)
        if bad:
            raise SystemExit(f"{col}: values {bad} not in {sorted(ok)}")
    if (c.groupby("experiment_id").n_total.nunique() > 1).any():
        raise SystemExit("an experiment has more than one n_total")
    return c


def subset(c: pd.DataFrame, cond: str) -> pd.DataFrame:
    return {"primary": c, "corrected": c[c.correction.eq("corrected")], "vbm": c[c.modality.eq("VBM")],
            "decreases": c[c.direction_tinnitus.eq("-")], "increases": c[c.direction_tinnitus.eq("+")]}[cond]


def dataset(c: pd.DataFrame):
    from nimare.dataset import Dataset

    if (c.groupby("experiment_id").space.nunique() > 1).any():
        raise SystemExit("an experiment mixes MNI and Talairach foci (NiMARE takes one space per experiment)")
    src = {str(e): {"contrasts": {"1": {"coords": {"space": x.space.iloc[0], "x": x.x.tolist(), "y": x.y.tolist(),
                                                     "z": x.z.tolist()},
                                          "metadata": {"sample_sizes": [int(x.n_total.iloc[0])]}}}}
           for e, x in c.groupby("experiment_id")}
    return Dataset(src)


def run_ale(c: pd.DataFrame, acfg: dict, n_iters: int, n_cores: int, seed: int):
    from nimare.correct import FWECorrector
    from nimare.meta.cbma.ale import ALE

    np.random.seed(seed)   # NiMARE's Monte Carlo null draws from the global numpy RNG
    res = ALE().fit(dataset(c))
    return FWECorrector(method="montecarlo", voxel_thresh=acfg["voxel_thresh"], n_iters=n_iters,
                        n_cores=n_cores, vfwe_only=False).transform(res)


def significant(cres, acfg: dict) -> np.ndarray:
    """Voxels in clusters with Monte Carlo cluster-size FWE p < cluster_p, from the -log10 p map (p is discrete,
    k / n_iters: a cluster exactly at p = cluster_p is not significant)."""
    logp = np.asarray(cres.get_map(LOGP).dataobj)
    return logp > -np.log10(acfg["cluster_p"]) + 1e-9


def clusters(cres, acfg: dict) -> tuple[pd.DataFrame, np.ndarray]:
    """Significant clusters as a peak table and a label image array."""
    from scipy import ndimage

    img = cres.get_map(TARGET)
    lab, n = ndimage.label(significant(cres, acfg))
    logp = np.asarray(cres.get_map(LOGP).dataobj)
    zmap = np.asarray(cres.get_map("z").dataobj)
    rows = []
    for k in range(1, n + 1):
        m = lab == k
        ijk = np.argwhere(m & (zmap == zmap[m].max()))[0]
        xyz = nib.affines.apply_affine(img.affine, ijk)
        rows.append({"cluster": k, "n_voxels": int(m.sum()),
                     "volume_mm3": float(m.sum() * np.prod(img.header.get_zooms()[:3])),
                     "peak_x": float(xyz[0]), "peak_y": float(xyz[1]), "peak_z": float(xyz[2]),
                     "peak_ale_z": float(zmap[m].max()), "cluster_p_fwe": float(10 ** -logp[m].max())})
    return pd.DataFrame(rows), lab


def nimare_ids(tab: pd.DataFrame, lab: np.ndarray, nim: pd.DataFrame, affine) -> list:
    """Map NiMARE's cluster IDs (jackknife table columns) to our cluster numbers via the peak location."""
    inv = np.linalg.inv(affine)
    ids = {}
    for r in nim.itertuples():
        i, j, k = np.rint(nib.affines.apply_affine(inv, [r.X, r.Y, r.Z])).astype(int)
        if 0 <= i < lab.shape[0] and 0 <= j < lab.shape[1] and 0 <= k < lab.shape[2] and lab[i, j, k]:
            ids.setdefault(int(lab[i, j, k]), str(r[1]))
    return [ids.get(int(k)) for k in tab.cluster]


def failsafe(c: pd.DataFrame, lab: np.ndarray, ks: list[int], acfg: dict, n_cores: int, seed: int) -> dict:
    """Fail-safe N per cluster (Acar et al. 2018): smallest number m of null experiments (random foci in the ALE mask;
    foci count and n drawn from the included experiments) whose addition makes the cluster non-significant. Nulls
    are added cumulatively (the first m of one fixed set, so 'lost' is monotone in expectation); each probe m is
    fitted once for all clusters (cache); grid n, 2n, 4n, ... up to and including the cap max_factor * n, then
    bisection to a resolution of n/10."""
    from nimare.utils import get_template

    tmpl = get_template("mni152_2mm", mask="brain")   # the ALE analysis mask (NiMARE default)
    vox = np.argwhere(np.asarray(tmpl.dataobj) > 0)
    rng = np.random.default_rng(seed)
    per = c.groupby("experiment_id").agg(nf=("x", "size"), n=("n_total", "first"))
    n_exp, cap = len(per), int(acfg["failsafe_max_factor"] * len(per))
    rows = []
    for j in range(cap):
        p = per.iloc[rng.integers(n_exp)]
        xyz = nib.affines.apply_affine(tmpl.affine, vox[rng.integers(len(vox), size=int(p.nf))])
        rows += [{"experiment_id": f"null{j:04d}", "space": "MNI", "x": v[0], "y": v[1], "z": v[2],
                  "n_total": int(p.n)} for v in xyz]
    nulls = pd.DataFrame(rows)
    cache: dict[int, np.ndarray] = {}

    def sig(m: int) -> np.ndarray:
        if m not in cache:
            cc = pd.concat([c[["experiment_id", "space", "x", "y", "z", "n_total"]],
                            nulls[nulls.experiment_id < f"null{m:04d}"]], ignore_index=True)
            cache[m] = significant(run_ale(cc, acfg, acfg["failsafe_iters"], n_cores, seed), acfg)
        return cache[m]

    lost = lambda m, k: not sig(m)[lab == k].any()
    grid, m = [], n_exp
    while m < cap:
        grid.append(m)
        m *= 2
    grid.append(cap)
    out = {}
    for k in ks:
        hit = next((m for m in grid if lost(m, k)), None)
        if hit is None:
            out[k] = {"failsafe_n": None, "failsafe_ratio": None,
                      "failsafe_note": f"cluster survives {cap} null experiments (cap {acfg['failsafe_max_factor']}x)"}
            continue
        lo = ([0] + grid)[grid.index(hit)]
        hi = hit
        while hi - lo > max(1, n_exp // 10):
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if lost(mid, k) else (mid, hi)
        out[k] = {"failsafe_n": hi, "failsafe_ratio": hi / n_exp,
                  "failsafe_note": "Acar 2018: robust if clearly larger than the number of included experiments"}
    return out


def cmd_run(args: argparse.Namespace) -> None:
    from nimare.diagnostics import Jackknife

    cfg, p = config.analysis(), config.paths()
    acfg = cfg["ale"]
    out = args.out or p.results / "20_ale"
    out.mkdir(parents=True, exist_ok=True)
    n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    c = load(args.coords or Path(__file__).resolve().parents[1] / "data" / "meta" / "coordinates_final.tsv")
    n_iters = args.n_iters or acfg["n_iters"]
    summary = {"n_experiments": int(c.experiment_id.nunique()), "n_foci": len(c), "n_iters": n_iters,
               "conditions": {}}

    def save() -> None:   # after every condition: a later failure / time-out keeps what finished
        (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))

    for cond in args.conditions or acfg["conditions"]:
        s = subset(c, cond)
        info = {"n_experiments": int(s.experiment_id.nunique()), "n_foci": len(s),
                "n_subjects": int(s.groupby("experiment_id").n_total.first().sum())}
        if info["n_experiments"] < acfg["min_experiments"]:
            info["status"] = f"underpowered (< {acfg['min_experiments']} experiments)"
            if not args.descriptive or info["n_experiments"] < MIN_RUNNABLE:
                summary["conditions"][cond] = info
                save()
                print(cond, info)
                continue
            info["status"] += ": run as descriptive only (plan 6.10 B2)"
        else:
            info["status"] = "powered"
        cres = run_ale(s, acfg, n_iters, n_cores, cfg["seed"])
        cres.get_map("z").to_filename(out / f"z_{cond}.nii.gz")
        cres.get_map(TARGET).to_filename(out / f"zcluster_{cond}.nii.gz")
        tab, lab = clusters(cres, acfg)
        info["n_clusters"] = len(tab)
        if len(tab):
            # only FWE-significant clusters (the cluster z map is > 0 for every p < .001 cluster)
            jk = Jackknife(target_image=TARGET, voxel_thresh=None, n_cores=n_cores,
                           target_threshold=float(st.norm.isf(acfg["cluster_p"])) + 1e-9).transform(cres)
            nim = jk.tables[f"{TARGET}_tab-clust"]
            jk.tables[f"{TARGET}_diag-Jackknife_tab-counts_tail-positive"].to_csv(
                out / f"jackknife_{cond}.tsv", sep="\t", index=False)   # NiMARE 0.22 table names
            nim.to_csv(out / f"nimare_clusters_{cond}.tsv", sep="\t", index=False)
            tab["nimare_cluster_id"] = nimare_ids(tab, lab, nim, cres.get_map(TARGET).affine)
            if args.failsafe and cond == "primary":
                fs = failsafe(s, lab, tab.cluster.tolist(), acfg, n_cores, cfg["seed"])
                for k, v in fs.items():
                    for kk, vv in v.items():
                        tab.loc[tab.cluster == k, kk] = vv
        tab.to_csv(out / f"clusters_{cond}.tsv", sep="\t", index=False)
        info["clusters"] = tab.to_dict("records")
        summary["conditions"][cond] = info
        save()
        print(cond, json.dumps(g.json_clean(info), indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--coords", default=None)
    r.add_argument("--out", type=Path, default=None)
    r.add_argument("--conditions", nargs="*", default=None)
    r.add_argument("--n-iters", type=int, default=None, help="override (pilot)")
    r.add_argument("--failsafe", action="store_true", help="fail-safe N for significant primary clusters")
    r.add_argument("--descriptive", action="store_true",
                   help="also run conditions below ale.min_experiments, labelled descriptive (underpowered)")
    r.set_defaults(func=cmd_run)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
