"""Leave-one-site-out (LOSO) transfer check of the ROI normative models (analysis plan 6.3, robustness).

For each site: fit the reference model on all OTHER sites' reference subjects, adapt it to the left-out site with
that site's controls only (PCNtoolkit transfer, as for a new scanner), and score the site's subjects. Compare with the
full-reference z-scores (scripts/07_normative.py): does a subject's deviation profile depend on its own site having
been part of the reference?

Subcommands
  tasks    reference x model x feature set x site -> <derivatives>/normative/loso/tasks.tsv
  fit      one task -> loso/<reference>/<model>/<set>/<site>.tsv (participant_id, role, fold, features)
  collect  per reference x model: merge sets and sites, compare with <reference>/<model>/z_all.tsv ->
           loso/<reference>/<model>/{z_loso.tsv, agreement.tsv, summary_site.tsv, summary.json}

Example (cluster, repo root)
  TINCEPTION_ENV=cluster uv run scripts/07c_normative_loso.py tasks
  T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/loso/tasks.tsv
  sbatch --array=1-$(($(wc -l < $T)-1)) slurm/07c_normative_loso.sbatch $T
  TINCEPTION_ENV=cluster uv run scripts/07c_normative_loso.py collect
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from tinception import config
from tinception.normative import core


def _dirs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    p = config.paths()
    norm = args.normative or p.derivatives / "normative"
    return args.master or p.root / "tinception_master.tsv", p.derivatives / "features", norm


def cmd_tasks(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, _, norm = _dirs(args)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    cc = master[master.in_case_control.astype(bool)]
    rows = []
    for r in cfg["loso"]["references"]:
        spec = cfg["references"][r]
        sites = sorted(cc.loc[cc[spec["column"]].eq(spec["value"]), "site"].unique())
        rows += [(r, m, s, site) for m in cfg["models"] for s in cfg["feature_sets"] for site in sites]
    (norm / "loso").mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame(rows, columns=["reference", "model", "feature_set", "site"])
    t.to_csv(norm / "loso" / "tasks.tsv", sep="\t", index=False)
    print(f"{len(t)} tasks -> {norm / 'loso' / 'tasks.tsv'}")


def cmd_fit(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, feat_dir, norm = _dirs(args)
    row = pd.read_csv(args.task_file, sep="\t").iloc[args.task - 1]
    spec = cfg["feature_sets"][row.feature_set]
    covs = list(cfg["models"][row.model]) + ([cfg["tiv_column"]] if spec["tiv"] else [])
    ref_spec = cfg["references"][row.reference]
    mcfg = core.reference_cfg(cfg, ref_spec)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    feats = core.load_feature_set(feat_dir, spec)
    features = [c for c in feats if c != core.ID]
    df, ref = core.design(master, feats, covs, cfg["batch_effects"], ref_spec)
    z = core.loso_transfer(df, ref, row.site, features, covs, cfg["batch_effects"], mcfg, config.analysis()["seed"],
                           cfg["loso"]["adapt_folds"])
    dest = norm / "loso" / row.reference / row.model / row.feature_set
    dest.mkdir(parents=True, exist_ok=True)
    z.to_csv(dest / f"{row.site}.tsv", sep="\t", index=False, na_rep="n/a")
    (dest / "settings.json").write_text(json.dumps(core.model_settings(mcfg, covs), indent=2))
    c = z[z.role == "control"][features].to_numpy()
    print(f"{row.reference}/{row.model}/{row.feature_set} LOSO {row.site}: {int((z.role == 'control').sum())} controls"
          f" (held-out z mean {np.nanmean(c):+.3f} sd {np.nanstd(c):.3f}), {int((z.role == 'scored').sum())} others; "
          f"rows without z {int(z[features].isna().all(axis=1).sum())}")


def cmd_collect(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, _, norm = _dirs(args)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)[[core.ID, "site"]]
    tasks = pd.read_csv(norm / "loso" / "tasks.tsv", sep="\t")
    for (ref_name, model), t in tasks.groupby(["reference", "model"], sort=False):
        d = norm / "loso" / ref_name / model
        missing = [f"{r.feature_set}/{r.site}" for r in t.itertuples() if not (d / r.feature_set / f"{r.site}.tsv").exists()]
        if missing:
            print(f"SKIP {ref_name}/{model}: {len(missing)} of {len(t)} tasks missing, e.g. {missing[:3]}")
            continue
        per_set = []
        for s in cfg["feature_sets"]:
            # a site can have no scored subjects for a model (e.g. triple has no PTA_HF -> no rows for Model B)
            parts = [pd.read_csv(d / s / f"{site}.tsv", sep="\t", na_values="n/a") for site in t.site.unique()]
            p = pd.concat([q for q in parts if len(q)], ignore_index=True)
            num = [c for c in p if c not in (core.ID, "role")]
            p[num] = p[num].apply(pd.to_numeric)
            per_set.append(p)
        roles = pd.concat([p[[core.ID, "role"]] for p in per_set]).drop_duplicates(core.ID)  # any set that scored it
        z = roles
        for p in per_set:
            z = z.merge(p.drop(columns=["role", "fold"]), on=core.ID, how="outer", validate="one_to_one")
        z.to_csv(d / "z_loso.tsv", sep="\t", index=False, na_rep="n/a")
        full = pd.read_csv(norm / ref_name / model / "z_all.tsv", sep="\t", na_values="n/a")
        features = [c for c in full if c not in (core.ID, "reference")]
        agree = core.profile_agreement(full, z, features, cfg["extreme_z"])
        agree = agree.merge(z[[core.ID, "role"]], on=core.ID).merge(master, on=core.ID, how="left")
        agree.to_csv(d / "agreement.tsv", sep="\t", index=False, na_rep="n/a")
        rows = []
        for site, a in agree.groupby("site"):
            ctrl = z[z.role.eq("control") & z[core.ID].isin(a[core.ID])][features].to_numpy()
            subj = np.nanmean(ctrl, axis=1)
            subj = subj[np.isfinite(subj)]
            se = subj.std(ddof=1) / np.sqrt(len(subj)) if len(subj) > 1 else np.nan
            rows.append({"site": site, "n_controls": int(a.role.eq("control").sum()),
                         "n_scored": int(a.role.eq("scored").sum()),
                         "profile_r_median": a.profile_r.median(), "profile_r_p10": a.profile_r.quantile(0.1),
                         "extreme_spearman": a[["n_extreme_full", "n_extreme_loso"]].corr("spearman").iloc[0, 1],
                         "loso_ctrl_mean": subj.mean() if len(subj) else np.nan, "loso_ctrl_mean_se": se,
                         "loso_ctrl_sd": np.nanstd(ctrl, ddof=1)})
        site_tab = pd.DataFrame(rows)
        lc = cfg["loso"]
        site_tab["ok"] = (site_tab.profile_r_median >= lc["min_profile_r"]) & \
            ((site_tab.loso_ctrl_mean.abs() - 1.96 * site_tab.loso_ctrl_mean_se) <= 0)
        site_tab.to_csv(d / "summary_site.tsv", sep="\t", index=False)
        summ = {"n_subjects": len(agree), "profile_r_median": float(agree.profile_r.median()),
                "extreme_spearman": float(agree[["n_extreme_full", "n_extreme_loso"]].corr("spearman").iloc[0, 1]),
                "sites_not_ok": site_tab.site[~site_tab.ok].tolist()}
        summ = {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in summ.items()}
        (d / "summary.json").write_text(json.dumps(summ, indent=2, allow_nan=False))
        print(f"{ref_name}/{model}: {len(agree)} subjects; median profile r {summ['profile_r_median']}, "
              f"extreme-count Spearman {summ['extreme_spearman']}; sites not ok: {summ['sites_not_ok']}")
        print(site_tab.round(3).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None, help="default: <root>/tinception_master.tsv")
    ap.add_argument("--normative", type=Path, default=None, help="default: <derivatives>/normative")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("tasks").set_defaults(func=cmd_tasks)
    f = sub.add_parser("fit")
    f.add_argument("--task-file", type=Path, required=True)
    f.add_argument("--task", type=int, required=True, help="1-based row of --task-file (SLURM_ARRAY_TASK_ID)")
    f.set_defaults(func=cmd_fit)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
