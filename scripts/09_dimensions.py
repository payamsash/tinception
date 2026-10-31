"""H3 (analysis plan 6.5): brain-clinical dimensions (PLS) and discrete subtypes in patients; settings in
analysis.yaml -> h3. Deviations: CO-reference normative z (Model A primary, Model B secondary; decision D1).

Outputs in <results>/09_h3/: pls.tsv (every model x clinical set x centering: mode-1 covariance, permutation p,
held-out leave-one-site-out within-site r and p; centering "site" = post-hoc sensitivity with brain and clinical
variables centred within site, added 2026-10-31), pls_weights.tsv (brain and clinical weights with bootstrap ratios), subtypes.json
(PCA on held-out controls, Gaussian mixtures, no-cluster null, stability, leave-one-site-out ARI, site association,
clinical validation), summary.json.
--part confound (analysis plan 6.9 A): confound-regressed (site fixed effects + age spline + sex), multi-mode PLS in
patients (all clinical sets), hearing-only PLS (PTA4, PTA_HF) in patients and in held-out controls, and the patients'
hearing mode projected onto the controls -> <results>/09_h3/confound/{modes.tsv, weights.tsv, projection.tsv,
summary.json}.

Example (cluster, repo root):  sbatch slurm/09_dimensions.sbatch [--part confound]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats as st

from tinception import config
from tinception.normative import core as ncore
from tinception.stats import group as g
from tinception.stats import h3

ID = "participant_id"


def deviations(norm: Path, model: str, ncfg: dict) -> tuple[pd.DataFrame, list[str]]:
    """z_all of one CO-reference model and its calibrated F1 + F3 features."""
    ok = pd.read_csv(norm / "CO" / model / "feature_calibration.tsv", sep="\t", index_col=0).ok
    feats = ncore.f1f3(ok, ncfg)
    return pd.read_csv(norm / "CO" / model / "z_all.tsv", sep="\t", na_values="n/a"), feats


def pls_block(d: pd.DataFrame, feats: list[str], clin: list[str], hcfg: dict, seed: int,
              centering: str = "none") -> tuple[dict, pd.DataFrame]:
    import gemmr

    x = d.dropna(subset=[*feats, *clin, "site"])
    X, Y, site = x[feats].to_numpy(float), x[clin].to_numpy(float), x.site.to_numpy()
    if centering == "site":  # only within-site brain-clinical covariation
        X, Y = h3.center_within(X, site), h3.center_within(Y, site)
    res = h3.pls_test(X, Y, site, hcfg["n_perm"], seed)
    rep = h3.pls_loso(X, Y, site, hcfg["n_perm"], seed + 1)
    bsr_u, bsr_v = h3.pls_bootstrap(X, Y, site, hcfg["n_boot"], seed + 2, res["U"], res["V"])
    try:
        n_req = int(np.ceil(np.nanmax(list(gemmr.pls_sample_size(X, Y, rs=(0.3,)).values()))))
    except (ValueError, RuntimeError, KeyError, IndexError) as err:  # gemmr's interpolation can fail outside its simulated range
        n_req = f"n/a ({str(err)[:60]})"
    row = {"centering": centering, "n": len(x), "n_sites": int(x.site.nunique()), "p_brain": len(feats), "clinical": "+".join(clin),
           "s1": res["s1"], "cov_explained": res["cov_explained"], "p_perm": res["p_perm"],
           "score_r": res["score_r"], "r_heldout": rep["r_heldout"], "p_heldout": rep["p_perm"],
           "gemmr_n_required_r0.3": n_req, "supported": bool(res["p_perm"] < 0.05 and rep["p_perm"] < 0.05)}
    w = pd.concat([pd.DataFrame({"side": "brain", "variable": feats, "weight": res["U"], "bsr": bsr_u}),
                   pd.DataFrame({"side": "clinical", "variable": clin, "weight": res["V"], "bsr": bsr_v})])
    return row, w.assign(clinical_set="+".join(clin), centering=centering)


def subtypes(z: pd.DataFrame, feats: list[str], master: pd.DataFrame, hcfg: dict, seed: int) -> dict:
    zz = z.merge(master[[ID, "site", "group", "hyperacusis_z", "affect_z", "log_duration", "distress_z"]], on=ID)
    zz = zz.dropna(subset=[*feats, "site"])
    ctrl, pat = zz[zz.reference.astype(bool)], zz[zz.group.eq("TI")]
    S, k_pc = h3.pca_from_controls(ctrl[feats].to_numpy(float), pat[feats].to_numpy(float), hcfg["pca_var"])
    ks = range(min(hcfg["ks"]), max(hcfg["ks"]) + 1)
    bic, models = h3.gmm_bic(S, ks, seed)
    best = min(bic, key=bic.get)
    obs, null = h3.cluster_null(S, bic, ks, hcfg["null_sims"], seed)
    p_null = float((np.sum(null >= obs) + 1) / (len(null) + 1))
    out = {"n_patients": len(pat), "n_controls_pca": len(ctrl), "pca_components": k_pc, "bic": bic, "best_k": best,
           "null_stat": obs, "null_p": p_null}
    k = max(2, best)  # criteria are evaluated for the best multi-cluster solution even if k = 1 wins
    labels = models[k].predict(S)
    site = pat.site.to_numpy()
    out |= {"k_evaluated": k, "cluster_sizes": np.bincount(labels, minlength=k).tolist(),
            "stability_ari": h3.stability_ari(S, labels, k, hcfg["stability_boot"], seed),
            "loso_ari": h3.loso_ari(S, site, labels, k, seed), "site_cramers_v": h3.cramers_v(labels, site)}
    val = {}  # clusters with < 5 patients with the variable -> not tested (reported via cluster_sizes)
    for v in ("hyperacusis_z", "affect_z", "log_duration", "distress_z"):
        x = pat[v].to_numpy(float)
        grp = [x[(labels == c) & np.isfinite(x)] for c in range(k)]
        if all(len(a) >= 5 for a in grp):
            val[v] = {"p_kruskal": float(st.kruskal(*grp).pvalue), "n": int(sum(len(a) for a in grp))}
    if val:  # Holm across the tested validation variables
        q = pg.multicomp([v["p_kruskal"] for v in val.values()], method="holm")[1]
        for v, qv in zip(val.values(), q, strict=True):
            v["p_holm"] = float(qv)
    out["clinical_validation"] = val
    crit = {"beats_null": p_null < 0.05 and best >= 2,
            "stable": out["stability_ari"] >= hcfg["min_stability_ari"],
            "replicates_loso": out["loso_ari"] >= hcfg["min_loso_ari"],
            "not_site": bool(out["site_cramers_v"] < hcfg["max_site_cramers_v"]),   # NaN (undefined) fails
            "clinically_distinct": any(v["p_holm"] < 0.05 for v in val.values())}
    out |= {"criteria": crit, "subtypes_claimed": bool(all(crit.values()))}
    return out


def confounds(d: pd.DataFrame, age_df: int) -> np.ndarray:
    """Age (centred natural cubic spline) + sex design, no intercept (site enters as fixed effect separately)."""
    import patsy

    return np.asarray(patsy.dmatrix(f"cr(age, df={age_df}, constraints='center') + C(sex)", d,
                                    return_type="dataframe").iloc[:, 1:], float)


def confound_block(x: pd.DataFrame, feats: list[str], clin: list[str], hcfg: dict, age_df: int,
                   seed: int) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Confound-regressed multi-mode PLS of one sample: per-mode test, LOSO replication, bootstrap ratios."""
    X, Y, site = x[feats].to_numpy(float), x[clin].to_numpy(float), x.site.to_numpy()
    C = confounds(x, age_df)
    Xr, Yr = h3.confound_residuals(X, C, site), h3.confound_residuals(Y, C, site)
    res = h3.pls_modes_test(Xr, Yr, site, hcfg["n_perm"], seed)
    rep = h3.pls_loso_confound(X, Y, C, site, res["U"], res["V"], hcfg["n_perm"], seed + 1)
    bu, bv = h3.pls_bootstrap_modes(Xr, Yr, site, hcfg["n_boot"], seed + 2, res["U"], res["V"])
    m = len(res["s"])
    modes = pd.DataFrame({"mode": np.arange(1, m + 1), "n": len(x), "n_sites": x.site.nunique(), "s": res["s"],
                          "cov_explained": res["cov_explained"], "p_perm": res["p_perm"], "score_r": res["score_r"],
                          "r_heldout": rep["r_heldout"], "p_heldout": rep["p_perm"]})
    modes["supported"] = (modes.p_perm < 0.05) & (modes.p_heldout < 0.05)
    w = pd.concat([pd.DataFrame({"mode": j + 1, "side": side, "variable": names, "weight": W[:, j], "bsr": B[:, j]})
                   for j in range(m) for side, names, W, B in (("brain", feats, res["U"], bu),
                                                                ("clinical", clin, res["V"], bv))])
    return modes, w, {"U": res["U"], "V": res["V"]}


def run_confound(args: argparse.Namespace, cfg: dict, p, master: pd.DataFrame, out: Path) -> None:
    hcfg, ncfg, seed, age_df = cfg["h3"], cfg["normative"], cfg["seed"], cfg["group"]["age_df"]
    out = out / "confound"
    out.mkdir(parents=True, exist_ok=True)
    norm = p.derivatives / "normative"
    hearing = ["PTA4", "PTA_HF"]
    modes, weights, proj = [], [], []
    cols = [ID, "site", "group", "age", "sex", *hcfg["clinical"], *hcfg["clinical_extra"]]
    for model in (hcfg["model_primary"], hcfg["model_secondary"]):
        z, feats = deviations(norm, model, ncfg)
        d = z.merge(master[cols], on=ID)
        pat, ctrl = d[d.group.eq("TI")], d[d.reference.astype(bool)]
        sets = [hcfg["clinical"]] + [[*hcfg["clinical"], e] for e in hcfg["clinical_extra"]] + [hearing]
        for i, clin in enumerate(sets):
            x = pat.dropna(subset=[*feats, *clin, "site", "age", "sex"])
            md, w, fit = confound_block(x, feats, clin, hcfg, age_df, seed + 10 * i)
            tag = {"model": model, "sample": "patients", "clinical": "+".join(clin)}
            modes.append(md.assign(**tag))
            weights.append(w.assign(**tag))
            print(f"{model} patients {tag['clinical']}:\n" + md.round(4).to_string(index=False))
            if clin == hearing:
                pat_fit = fit
        xc = ctrl.dropna(subset=[*feats, *hearing, "site", "age", "sex"])
        md, w, _ = confound_block(xc, feats, hearing, hcfg, age_df, seed + 100)
        tag = {"model": model, "sample": "controls", "clinical": "+".join(hearing)}
        modes.append(md.assign(**tag))
        weights.append(w.assign(**tag))
        print(f"{model} controls {tag['clinical']}:\n" + md.round(4).to_string(index=False))
        # patients' hearing mode 1 applied to the controls (fixed weights)
        Cc, sc = confounds(xc, age_df), xc.site.to_numpy()
        Xc = h3.confound_residuals(xc[feats].to_numpy(float), Cc, sc)
        Yc = h3.confound_residuals(xc[hearing].to_numpy(float), Cc, sc)
        pr = h3.project_r(Xc, Yc, sc, pat_fit["U"][:, 0], pat_fit["V"][:, 0], hcfg["n_perm"], seed + 200)
        proj.append({"model": model, "weights_from": "patients hearing-only mode 1", "applied_to": "controls", **pr})
        print(f"{model} patients' hearing mode -> controls: r {pr['r']:.3f} p {pr['p_perm']:.4f} (n {pr['n']})")
    modes, weights, proj = pd.concat(modes, ignore_index=True), pd.concat(weights), pd.DataFrame(proj)
    modes.to_csv(out / "modes.tsv", sep="\t", index=False, na_rep="n/a")
    weights.to_csv(out / "weights.tsv", sep="\t", index=False)
    proj.to_csv(out / "projection.tsv", sep="\t", index=False)
    core = "+".join(hcfg["clinical"])
    clin_w = weights[(weights.side == "clinical")]
    summary = {"supported_modes": modes[modes.supported][["model", "sample", "clinical", "mode", "p_perm",
                                                          "r_heldout", "p_heldout"]].to_dict("records"),
               "core_mode1": modes[(modes.clinical == core) & (modes["mode"] == 1)][
                   ["model", "p_perm", "r_heldout", "p_heldout", "supported"]].to_dict("records"),
               "clinical_weights_supported": clin_w.merge(modes[modes.supported][["model", "sample", "clinical", "mode"]],
                                                          on=["model", "sample", "clinical", "mode"])[
                   ["model", "sample", "clinical", "mode", "variable", "weight", "bsr"]].to_dict("records"),
               "projection_to_controls": proj.to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(json.dumps(g.json_clean(summary), indent=2, default=str))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default: <results>/09_h3")
    ap.add_argument("--part", choices=["main", "confound"], default="main")
    args = ap.parse_args()
    cfg, p = config.analysis(), config.paths()
    hcfg, ncfg, seed = cfg["h3"], cfg["normative"], cfg["seed"]
    if hcfg["features"] != "F1F3":
        raise SystemExit(f"h3.features={hcfg['features']!r}: only F1F3 is implemented")
    out = args.out or p.results / "09_h3"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(args.master or p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    master = master.assign(log_duration=np.log1p(master.duration_months))
    if args.part == "confound":
        run_confound(args, cfg, p, master, out)
        return
    norm = p.derivatives / "normative"
    rows, weights, subs = [], [], {}
    for model in (hcfg["model_primary"], hcfg["model_secondary"]):
        z, feats = deviations(norm, model, ncfg)
        d = z.merge(master[[ID, "site", "group", *hcfg["clinical"], *hcfg["clinical_extra"]]], on=ID)
        d = d[d.group.eq("TI")]
        sets = [hcfg["clinical"]] + [[*hcfg["clinical"], e] for e in hcfg["clinical_extra"]]
        for centering in ("none", "site"):
            for i, clin in enumerate(sets):
                r, w = pls_block(d, feats, clin, hcfg, seed + 10 * i, centering)
                rows.append({"model": model, "primary": model == hcfg["model_primary"] and i == 0
                             and centering == "none", **r})
                weights.append(w.assign(model=model))
                print(f"{model} {centering} {r['clinical']}: n={r['n']} cov {r['cov_explained']:.2f} "
                      f"p {r['p_perm']:.4f} | held-out r {r['r_heldout']:.3f} p {r['p_heldout']:.4f}")
        subs[model] = subtypes(z, feats, master, hcfg, seed)
    pls = pd.DataFrame(rows)
    pls.to_csv(out / "pls.tsv", sep="\t", index=False, na_rep="n/a")
    pd.concat(weights).to_csv(out / "pls_weights.tsv", sep="\t", index=False)
    (out / "subtypes.json").write_text(json.dumps(g.json_clean(subs), indent=2, allow_nan=False))
    prim = pls[pls.primary].iloc[0]
    summary = {"H3_dimension_supported": bool(prim.supported),
               "primary": prim.drop(labels=["primary"]).to_dict(),
               "secondary_and_exploratory": pls[~pls.primary][["model", "centering", "clinical", "n", "p_perm", "r_heldout",
                                                               "p_heldout", "supported"]].to_dict("records"),
               "subtypes_claimed": {m: v["subtypes_claimed"] for m, v in subs.items()},
               "subtypes": {m: {k: v[k] for k in ("best_k", "null_p", "stability_ari", "loso_ari", "site_cramers_v",
                                                  "criteria")} for m, v in subs.items()}}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(json.dumps(g.json_clean(summary), indent=2, default=str))


if __name__ == "__main__":
    main()
