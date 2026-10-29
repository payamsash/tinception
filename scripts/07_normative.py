"""Normative models (PCNtoolkit 1.3 BLR) on FreeSurfer ROI features; settings in analysis.yaml -> normative.

Subcommands
  tasks    write the task list (reference x model x feature set) for slurm/07_normative.sbatch
  fit      run one task (or an ad-hoc pilot with --features): cross-fitted z-scores of the reference sample,
           full-reference z-scores for everyone else, fold statistics; the full model is saved for transfer
  collect  per reference x model: merge the feature sets (refuses sets fitted with different model settings),
           calibration of the held-out reference z-scores (pooled, per site, residual site/sex effects in mean and
           spread), per-feature calibration (failing features are excluded from H2), extreme-deviation counts

Outputs: <derivatives>/normative/<reference>/<model>/<feature set>/{z.tsv, stats.tsv, settings.json, model/}
         <derivatives>/normative/<reference>/<model>/{z_all.tsv, calibration_site.tsv, calibration_batch_tests.tsv,
         feature_calibration.tsv, calibration.json, extreme.tsv}

Example (cluster, repo root)
  TINCEPTION_ENV=cluster uv run scripts/07_normative.py tasks
  T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/tasks.tsv
  sbatch --array=1-$(($(wc -l < $T)-1)) slurm/07_normative.sbatch $T
  TINCEPTION_ENV=cluster uv run scripts/07_normative.py collect
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tinception import config
from tinception.normative import core


def _paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    p = config.paths()
    return (args.master or p.root / "tinception_master.tsv", args.features_dir or p.derivatives / "features",
            args.out or p.derivatives / "normative")


def cmd_tasks(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    _, _, out = _paths(args)
    out.mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame([(r, m, s) for r in cfg["references"] for m in cfg["models"] for s in cfg["feature_sets"]],
                     columns=["reference", "model", "feature_set"])
    t.to_csv(out / "tasks.tsv", sep="\t", index=False)
    print(f"{len(t)} tasks -> {out / 'tasks.tsv'}")


def cmd_fit(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, feat_dir, out = _paths(args)
    if args.task_file:
        row = pd.read_csv(args.task_file, sep="\t").iloc[args.task - 1]
        args.reference, args.model, args.feature_set = row.reference, row.model, row.feature_set
    spec = cfg["feature_sets"][args.feature_set]
    covs = list(cfg["models"][args.model]) + ([cfg["tiv_column"]] if spec["tiv"] else [])
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    feats = core.load_feature_set(feat_dir, spec)
    if args.features:  # pilot: a subset of this feature set
        missing = sorted(set(args.features) - set(feats.columns))
        if missing:
            raise SystemExit(f"not in {args.feature_set}: {missing}")
        feats = feats[[core.ID, *args.features]]
    features = [c for c in feats if c != core.ID]
    ref_spec = cfg["references"][args.reference]
    cfg = core.reference_cfg(cfg, ref_spec)
    df, ref = core.design(master, feats, covs, cfg["batch_effects"], ref_spec)
    dest = out / args.reference / args.model / (args.feature_set if not args.features else f"{args.feature_set}_pilot")
    dest.mkdir(parents=True, exist_ok=True)
    print(f"{args.reference}/{args.model}/{args.feature_set}: {len(features)} features, covariates age+{covs}, "
          f"{int(ref.sum())} reference / {len(df)} scored subjects")
    z, stats = core.cross_fit(df, ref, features, covs, cfg["batch_effects"], cfg, config.analysis()["seed"],
                              save_dir=dest / "model")
    z.insert(1, "reference", ref.to_numpy())
    z.to_csv(dest / "z.tsv", sep="\t", index=False, na_rep="n/a")
    stats.to_csv(dest / "stats.tsv", sep="\t", index=False, na_rep="n/a")
    (dest / "settings.json").write_text(json.dumps(core.model_settings(cfg, covs), indent=2))
    unscored = z[core.ID][z[features].isna().all(axis=1)].tolist()
    if unscored:
        print(f"WARNING: {len(unscored)} subjects not scored (site/sex level absent from training): {unscored[:10]}")
    bad = [f for f in features if z.loc[z[features].notna().any(axis=1), f].isna().any()]
    held = z[z.reference]
    print(f"-> {dest}  held-out reference z: mean {held[features].stack().mean():+.3f}, "
          f"sd {held[features].stack().std():.3f}; features with missing z: {bad or 'none'}")


def cmd_collect(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, _, out = _paths(args)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    master = master[[core.ID, *cfg["batch_effects"], *{r["column"] for r in cfg["references"].values()}]]
    for ref_name in cfg["references"]:
        for model in cfg["models"]:
            d = out / ref_name / model
            sets = [s for s in cfg["feature_sets"] if (d / s / "z.tsv").exists()]
            if not sets:
                continue
            if len(sets) < len(cfg["feature_sets"]):  # a partial merge would silently drop features
                print(f"SKIP {ref_name}/{model}: missing feature sets {sorted(set(cfg['feature_sets']) - set(sets))}")
                continue
            settings = {s: json.loads((d / s / "settings.json").read_text()) if (d / s / "settings.json").exists()
                        else None for s in sets}
            core_settings = {s: {k: v for k, v in (st or {}).items() if k != "covariates"} for s, st in settings.items()}
            if None in settings.values() or len({json.dumps(v, sort_keys=True) for v in core_settings.values()}) > 1:
                print(f"SKIP {ref_name}/{model}: feature sets fitted with different or unrecorded settings: {settings}")
                continue
            parts = [pd.read_csv(d / s / "z.tsv", sep="\t", na_values="n/a") for s in sets]
            # subjects can differ between sets (missing features) -> outer join; folds are per set (see <set>/z.tsv),
            # and every reference row of every set is out-of-fold, so reference membership comes from the master
            z = parts[0].drop(columns=["fold", "reference"])
            for p in parts[1:]:
                z = z.merge(p.drop(columns=["fold", "reference"]), on=core.ID, how="outer", validate="one_to_one")
            features = [c for c in z if c != core.ID]
            spec = cfg["references"][ref_name]
            z = z.merge(master, on=core.ID, how="left", validate="one_to_one")
            z.insert(1, "reference", z[spec["column"]].eq(spec["value"]))
            z = z.drop(columns=[c for c in master if c not in (core.ID, *cfg["batch_effects"])])
            z.drop(columns=cfg["batch_effects"]).to_csv(d / "z_all.tsv", sep="\t", index=False, na_rep="n/a")
            core.extreme_counts(z, features, cfg["extreme_z"]).to_csv(d / "extreme.tsv", sep="\t", index=False)
            held = z[z.reference]
            if held[features].notna().sum().sum() == 0:
                print(f"SKIP calibration {ref_name}/{model}: no held-out reference z-scores")
                continue
            pooled = core.write_calibration(d, held, held[cfg["batch_effects"]], features, cfg["calibration"], "roi",
                                            {"n_reference": len(held), "n_scored": len(z), "n_features": len(features),
                                             "sets": sets, "settings": settings[sets[0]] | {"covariates": None}})
            print(f"{ref_name}/{model}: {len(features)} features, {len(z)} subjects; held-out mean {pooled['mean']}, "
                  f"sd {pooled['sd']}; batch tests {pooled['batch_tests_run']}, p<.05 {pooled['batch_frac_p05']}, "
                  f"FDR hits {pooled['batch_fdr_hit_list'][:5]}; features failing {pooled['features_failing_list']} "
                  f"-> {'OK' if pooled['ok'] else 'CHECK'}; sites failing mean CI {pooled['sites_failing_mean']}, "
                  f"skew {pooled['sites_failing_skew']}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None, help="default: <root>/tinception_master.tsv")
    ap.add_argument("--features-dir", type=Path, default=None, help="default: <derivatives>/features")
    ap.add_argument("--out", type=Path, default=None, help="default: <derivatives>/normative")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("tasks").set_defaults(func=cmd_tasks)
    f = sub.add_parser("fit")
    f.add_argument("--task-file", type=Path, default=None)
    f.add_argument("--task", type=int, default=None, help="1-based row of --task-file (SLURM_ARRAY_TASK_ID)")
    f.add_argument("--reference", default="CO")
    f.add_argument("--model", default="A")
    f.add_argument("--feature-set", default="dk_thick")
    f.add_argument("--features", nargs="*", default=None, help="pilot: subset of the feature set")
    f.set_defaults(func=cmd_fit)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    args = ap.parse_args()
    if args.cmd == "fit" and args.task_file and not args.task:
        ap.error("--task-file needs --task")
    args.func(args)


if __name__ == "__main__":
    main()
