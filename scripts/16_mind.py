"""MIND networks and gradients (analysis plan 6.12; user decisions 2026-11-03). MIND matrices are computed with the
authors' implementation (isebenius/MIND, analysis.yaml mind.code), unchanged.

  subjects           -> rawdata/code/mind_subjects.tsv (every subject with FreeSurfer vertex features)
  compute --task i   one subject (SLURM array over rawdata/code/mind_subjects.tsv), every parcellation;
                     -> <derivatives>/mind/<parcellation>/<sub>.npy (float32) + regions.json
  collect            -> <derivatives>/mind/<parcellation>.npz (subjects x regions x regions) + completeness
  analyze            -> <results>/16_mind/ (regional degree, global MIND, Yeo-7 summaries, gradients; H1 models)

Example (cluster):  16_mind.py subjects; L=<rawdata>/code/mind_subjects.tsv;
                    sbatch --array=1-$(($(wc -l < $L)-1))%200 slurm/16_mind.sbatch
"""

from __future__ import annotations

import argparse
import json
import re
import sys

import numpy as np
import pandas as pd

from tinception import config

ID = "participant_id"
YEO = re.compile(r"7Networks_[LR]H_([A-Za-z]+)_")   # Schaefer names as returned by MIND: lh_7Networks_LH_Vis_1


def yeo_network(region: str) -> str:
    m = YEO.search(region)
    if not m:
        raise ValueError(f"not a Schaefer 7-network region name: {region}")
    return m.group(1)


def cmd_subjects(args: argparse.Namespace) -> None:
    p = config.paths()
    subs = (p.derivatives / "features" / "vertex_subjects.txt").read_text().split()
    out = p.rawdata / "code" / "mind_subjects.tsv"
    pd.DataFrame({ID: subs}).to_csv(out, sep="\t", index=False)
    print(f"{len(subs)} subjects -> {out}")


def cmd_compute(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    mcfg = cfg["mind"]
    sys.path.insert(1, mcfg["code"])
    from MIND import compute_MIND

    subs = pd.read_csv(p.rawdata / "code" / "mind_subjects.tsv", sep="\t")[ID].tolist()
    sub = subs[args.task - 1]
    root = p.derivatives / "mind"
    sdir = p.freesurfer / sub
    failed = []
    for key, parc in mcfg["parcellations"].items():   # a failing parcellation does not stop the others
        out = root / key / f"{sub}.npy"
        if out.exists():
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        missing = [h for h in ("lh", "rh") if not (sdir / "label" / f"{h}.{parc}.annot").exists()]
        if missing:
            failed.append(f"{key}: missing {missing} {parc}.annot")
            continue
        M = compute_MIND(str(sdir), mcfg["features"], parc)
        # region names are saved per subject (no shared file written concurrently); collect checks the order
        (root / key / f"{sub}.regions.json").write_text(json.dumps(list(M.index)))
        tmp = out.with_suffix(".tmp.npy")
        np.save(tmp, M.to_numpy(np.float32))
        tmp.rename(out)
        print(f"{sub} {key}: {M.shape}", flush=True)
    if failed:
        raise SystemExit(f"{sub}: " + "; ".join(failed))


def cmd_collect(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    root = p.derivatives / "mind"
    subs = pd.read_csv(p.rawdata / "code" / "mind_subjects.tsv", sep="\t")[ID].tolist()
    for key in cfg["mind"]["parcellations"]:
        have = [s for s in subs if (root / key / f"{s}.npy").exists()]
        missing = sorted(set(subs) - set(have))
        if not have:
            print(f"{key}: no subjects computed, skipped")
            continue
        names = {}
        for s in have:   # per-subject region files (compute) or the legacy shared regions.json
            f = root / key / f"{s}.regions.json"
            names[s] = tuple(json.loads(f.read_text())) if f.exists() else None
        legacy = root / key / "regions.json"
        ref = next((v for v in names.values() if v), None) or tuple(json.loads(legacy.read_text()))
        if any(v is not None and v != ref for v in names.values()):
            raise SystemExit(f"{key}: region order differs between subjects")
        X = np.stack([np.load(root / key / f"{s}.npy") for s in have])
        good = np.array([np.isfinite(m).all() for m in X])
        bad = [s for s, ok in zip(have, good, strict=True) if not ok]
        np.savez_compressed(root / f"{key}.npz", mind=X[good], subjects=np.array(have)[good], regions=np.array(ref))
        print(f"{key}: {int(good.sum())}/{len(subs)} subjects, {X.shape[1]} regions, missing {len(missing)} "
              f"{missing[:5]}, dropped non-finite {len(bad)} {bad[:5]}")


def _gradients(M: np.ndarray, ref: np.ndarray | None, gcfg: dict) -> np.ndarray:
    """Diffusion-map gradients (regions x n) of one MIND matrix; Procrustes-aligned to `ref` if given."""
    from brainspace.gradient import GradientMaps

    gm = GradientMaps(n_components=gcfg["n_components"], approach=gcfg["approach"], kernel=gcfg["kernel"],
                      random_state=0, alignment="procrustes" if ref is not None else None)
    gm.fit(M, sparsity=gcfg["sparsity"], reference=ref)
    return gm.aligned_ if ref is not None else gm.gradients_


def cmd_analyze(args: argparse.Namespace) -> None:
    from tinception.stats import evidence as ev
    from tinception.stats import group as g

    cfg, p = config.analysis(), config.paths()
    mcfg, gcfg = cfg["mind"], cfg["group"]
    root = p.derivatives / "mind"
    out = p.results / "16_mind"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    prim = master[master.in_case_control.astype(bool) & master.PTA_HF.notna()].copy()
    prim = prim.assign(tinnitus=prim.group.eq("TI").astype(int)).reset_index(drop=True)
    prim_ids = set(prim[ID])
    terms = ["tinnitus", "PTA4", "PTA_HF_res", *[f"C({c})" if c == "sex" else c for c in gcfg["covariates"]]]
    cont = ["PTA4", "PTA_HF_res", "euler_total"]
    con = {"tinnitus": {"tinnitus": 1.0}}
    rows, gsum = [], {}

    def family(df: pd.DataFrame, cols: list[str], name: str, key: str) -> None:
        d = g.prepare(df, cols, set())
        rhs = {y: g.rhs(terms, gcfg["age_df"], False) for y in d.attrs["outcomes"]}
        r = g.run_family(d, rhs, con, name, gcfg["tost_d"], cont)
        rows.append(r.assign(parcellation=key))

    for key in mcfg["parcellations"]:
        z = np.load(root / f"{key}.npz")
        X, subs, regions = z["mind"].astype(float), list(z["subjects"]), [str(r) for r in z["regions"]]
        keep = [i for i, s_ in enumerate(subs) if s_ in prim_ids]
        X, subs = X[keep], [subs[i] for i in keep]
        n_r = X.shape[1]
        off = ~np.eye(n_r, dtype=bool)
        deg = X.sum(axis=2) / (n_r - 1)                      # diagonal of MIND is 0
        cols = {ID: subs, **{f"deg__{r}": deg[:, j] for j, r in enumerate(regions)},
                "mind_global": X[:, off].mean(axis=1)}
        if key.startswith("schaefer"):                        # Yeo-7 summaries from the Schaefer labels
            net = np.array([yeo_network(r) for r in regions])
            if len(set(net)) != 7:
                raise SystemExit(f"{key}: expected 7 Yeo networks, got {sorted(set(net))}")
            for n_ in sorted(set(net)):
                i = net == n_
                cols[f"net_degree__{n_}"] = deg[:, i].mean(axis=1)
                w = X[:, i][:, :, i]
                cols[f"net_within__{n_}"] = w[:, ~np.eye(i.sum(), dtype=bool)].mean(axis=1)
                cols[f"net_between__{n_}"] = X[:, i][:, :, ~i].mean(axis=(1, 2))
        # gradients: template = mean MIND of the controls; every subject Procrustes-aligned to it
        is_co = np.isin(subs, prim.loc[prim.tinnitus == 0, ID])
        ref = _gradients(X[is_co].mean(axis=0), None, mcfg["gradients"])
        G = np.stack([_gradients(m, ref, mcfg["gradients"]) for m in X])           # subjects x regions x n
        cols["grad_dispersion"] = np.linalg.norm(G - G.mean(axis=1, keepdims=True), axis=2).mean(axis=1)
        cols["grad_g1_range"] = G[:, :, 0].max(axis=1) - G[:, :, 0].min(axis=1)
        if key == mcfg["primary"]:
            cols |= {f"g1__{r}": G[:, j, 0] for j, r in enumerate(regions)}
        feats = pd.DataFrame(cols)
        # PTA_HF beyond PTA4 residualized on exactly the analysed rows (as 08_group / 13_maps)
        df = g.residualize(prim.merge(feats, on=ID).reset_index(drop=True), "PTA_HF", ["PTA4"], "PTA_HF_res")
        deg_cols = [c for c in feats if c.startswith("deg__")]
        family(df, deg_cols, "degree", key)
        family(df, ["mind_global", "grad_dispersion", "grad_g1_range"], "global", key)
        if key.startswith("schaefer"):
            family(df, [c for c in feats if c.startswith("net_")], "yeo7", key)
        if key == mcfg["primary"]:
            family(df, [c for c in feats if c.startswith("g1__")], "g1_position", key)
        gsum[key] = {"n": len(df), "n_regions": n_r, "mean_mind": float(feats.mind_global.mean())}
        print(f"{key}: n={len(df)}, {n_r} regions", flush=True)
    m = pd.concat(rows, ignore_index=True)
    conv = m.converged.astype(bool)
    m["bf01_cauchy"] = [np.exp(-ev.log_bf10(b, s_)) if c else np.nan for b, s_, c in zip(m.beta, m.se, conv,
                                                                                         strict=True)]
    m.to_csv(out / "models.tsv", sep="\t", index=False, na_rep="n/a")
    q_thr = gcfg["fdr_q"]
    smry = m.groupby(["parcellation", "family"]).agg(
        n_tests=("p", "size"), hits=("q", lambda q: int((q < q_thr).sum())), min_q=("q", "min"),
        max_abs_d=("beta", lambda b: float(b.abs().max())), equivalent=("p_tost", lambda x: int((x < 0.05).sum())),
        bf01_median=("bf01_cauchy", "median"), nonconverged=("converged", lambda c: int((~c.astype(bool)).sum())))
    smry.reset_index().to_csv(out / "summary.tsv", sep="\t", index=False)
    summary = {"samples": gsum, "summary": smry.reset_index().to_dict("records"),
               "hits": m[m.q < q_thr][["parcellation", "family", "outcome", "beta", "p", "q"]].to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(smry.round(3).to_string())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    c.add_argument("--task", type=int, required=True, help="1-based subject index (SLURM_ARRAY_TASK_ID)")
    c.set_defaults(func=cmd_compute)
    sub.add_parser("subjects").set_defaults(func=cmd_subjects)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    sub.add_parser("analyze").set_defaults(func=cmd_analyze)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
