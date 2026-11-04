"""Leave-site-out prediction (analysis plan 6.8, added 2026-11-03): can brain features (150 ROIs) predict tinnitus
status or distress in unseen sites, beyond the confounds? Positive controls: age and PTA4.

  run --task <name>  (one SLURM array task per prediction task) -> <results>/14_prediction/<task>.json,
                     predictions_<task>.tsv
  collect            -> <results>/14_prediction/{scores.tsv, summary.json} (Holm over the primary family)
  --features pyment-<model>  (plan 6.11 B): the frozen 64-d pyment embeddings (scripts/19b_pyment_embed.py) instead
                     of the 150 ROIs, same pipeline; tasks = prediction.embedding_tasks (adds the P3_sex positive
                     control); outputs in <results>/14_prediction_pyment-<model>/

Example (cluster, repo root):  sbatch --array=1-6 slurm/14_prediction.sbatch ; then 14_prediction.py collect
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import group as g
from tinception.stats import predict as pr

ID = "participant_id"


def task_data(name: str, master: pd.DataFrame, f: pd.DataFrame, cols: list[str], age_df: int):
    """(X, y, confounds, site, kind, analysed rows) for one task."""
    import patsy

    m = master.merge(f, on=ID)
    prim = m[m.in_case_control.astype(bool) & m.PTA_HF.notna()]
    spec = {"T1_tinnitus": (prim, "tinnitus", "classification", True, True),
            "T2_TNHL_vs_HLOS": (prim[prim.group4.isin(["TNHL", "HLOS"])], "tinnitus", "classification", True, True),
            "T3_TINN_vs_CTRL": (prim[prim.group4.isin(["TINN", "CTRL"])], "tinnitus", "classification", True, True),
            "T4_distress": (m[m.in_patients.astype(bool) & m.distress_z.notna()], "distress_z", "regression", True,
                            False),
            "P1_age": (prim, "age", "regression", False, False),
            "P2_PTA4": (prim, "PTA4", "regression", True, False),
            "P3_sex": (prim, "male", "classification", True, False)}   # sex beyond eTIV (plan 6.11 B)
    d, target, kind, use_age, hearing = spec[name]
    d = d.assign(tinnitus=d.group.eq("TI").astype(int), log_eTIV=np.log(d.eTIV),
                 male=d.sex.map({"M": 1.0, "F": 0.0}))
    need = [target, *cols, "group", "age", "sex", "euler_total", "log_eTIV", "site"] + (
        ["PTA4", "PTA_HF"] if hearing else [])
    d = d.dropna(subset=need).reset_index(drop=True)
    terms = (["C(sex)"] if target != "male" else []) + ["euler_total"]
    if hearing:
        d = g.residualize(d, "PTA_HF", ["PTA4"], "PTA_HF_res")
        terms += ["PTA4", "PTA_HF_res"]
    rhs = g.rhs(terms, age_df, True) if use_age else " + ".join([*terms, "log_eTIV"])   # same age term as 08
    C = np.asarray(patsy.dmatrix(rhs, d, return_type="dataframe").iloc[:, 1:], float)
    C = (C - C.mean(axis=0)) / np.where(C.std(axis=0) > 0, C.std(axis=0), 1)
    site = d.site.to_numpy()
    X = pr.site_standardize(d[cols].to_numpy(float), site)
    return X, d[target].to_numpy(float), C, site, kind, d


def setup(args: argparse.Namespace):
    """(config, paths, output dir, task list) for the chosen feature set."""
    cfg, p = config.analysis(), config.paths()
    pcfg = cfg["prediction"]
    if args.features == "roi":
        return cfg, p, p.results / "14_prediction", pcfg["tasks"]
    return cfg, p, p.results / f"14_prediction_{args.features}", pcfg["embedding_tasks"]


def features(args: argparse.Namespace, cfg: dict, p) -> tuple[pd.DataFrame, list[str]]:
    if args.features == "roi":
        return g.roi_features(p.derivatives / "features", cfg["normative"], cfg["group"])
    f = pd.read_csv(p.derivatives / "pyment" / f"embeddings_{args.features.removeprefix('pyment-')}.tsv", sep="\t")
    cols = [c for c in f.columns if c.startswith("emb") and f[c].std() > 0]   # drop dead (all-constant) units
    return f[[ID, *cols]], cols


def cmd_run(args: argparse.Namespace) -> None:
    cfg, p, out, tasks = setup(args)
    pcfg = cfg["prediction"]
    out.mkdir(parents=True, exist_ok=True)
    name = args.task or tasks[int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1]
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    f, cols = features(args, cfg, p)
    X, y, C, site, kind, d = task_data(name, master, f, cols, cfg["group"]["age_df"])
    if kind == "classification":
        y = y.astype(int)
    pred = pr.loso_predict(X, y, C, site, kind)
    sc = pr.score(y, pred, kind, site)
    key = "auc" if kind == "classification" else "r"
    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    n_perm = args.n_perm if args.n_perm is not None else pcfg["n_perm"]
    null = pr.permutation_null(X, y, C, site, kind, n_perm, cfg["seed"], n_jobs)
    fin = null[np.isfinite(null)]
    res = {"task": name, "features": args.features, "kind": kind, "n": len(y), "n_sites": int(pd.Series(site).nunique()),
           "n_class1": int(y.sum()) if kind == "classification" else None, "n_features": X.shape[1],
           "n_confounds": C.shape[1], **sc, "null_mean": float(fin.mean()) if len(fin) else None,
           "null_sd": float(fin.std()) if len(fin) else None,
           "null_95": float(np.quantile(fin, 0.95)) if len(fin) else None, "n_perm": n_perm,
           "n_perm_finite": len(fin), "p_perm": pr.perm_p(sc[key], null)}
    pd.DataFrame({ID: d[ID], "site": site, "y": y, "pred": pred}).to_csv(out / f"predictions_{name}.tsv", sep="\t",
                                                                         index=False)
    (out / f"{name}.json").write_text(json.dumps(g.json_clean(res), indent=2, allow_nan=False))
    print(json.dumps(g.json_clean(res), indent=1))


def cmd_collect(args: argparse.Namespace) -> None:
    import pingouin as pg

    cfg, _, out, tasks = setup(args)
    pcfg = cfg["prediction"]
    missing = [t for t in tasks if not (out / f"{t}.json").exists()]
    rows = [json.loads((out / f"{t}.json").read_text()) for t in tasks if t not in missing]
    s = pd.DataFrame(rows)
    pilot = s[s.n_perm != pcfg["n_perm"]].task.tolist() if len(s) else []
    if set(missing + pilot) & set(pcfg["primary_family"]):
        raise SystemExit(f"primary family incomplete (missing {missing}, pilot runs {pilot}): no Holm correction")
    prim = s.task.isin(pcfg["primary_family"])
    s["p_holm"] = np.nan
    s.loc[prim, "p_holm"] = pg.multicomp(s.loc[prim, "p_perm"].to_numpy(), method="holm")[1]
    s.to_csv(out / "scores.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps(g.json_clean({"scores": s.to_dict("records"), "missing": missing}),
                                                 indent=2, allow_nan=False))
    print(s.round(4).to_string(index=False), "\nmissing:", missing)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", default="roi", choices=["roi", "pyment-multi-2025", "pyment-reg-2025"])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--task", default=None, help="task name (default: from SLURM_ARRAY_TASK_ID)")
    r.add_argument("--n-perm", type=int, default=None, help="override (pilot)")
    r.set_defaults(func=cmd_run)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
