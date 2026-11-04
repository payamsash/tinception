"""External normative reference (analysis plan 6.8, added 2026-11-03): pretrained braincharts models transferred to
our sites with symmetric 2-fold cross-fitting of the site adaptation (CO controls).

  prepare  -> <derivatives>/normative_braincharts/{subjects.tsv, responses_<set>.tsv, sets.json}
  (transfer: slurm/07d_braincharts.sbatch runs scripts/braincharts_transfer.py in the PCNtoolkit 0.35 env)
  analyze  -> <results>/07d_braincharts/{calibration.tsv, agreement.tsv, burden.tsv, effects.tsv, summary.json}
"""

from __future__ import annotations

import argparse
import json
import re

import numpy as np
import pandas as pd

from tinception import config
from tinception.normative import core as ncore
from tinception.stats import group as g
from tinception.stats import h2
from tinception.stats import robust as rb

ID = "participant_id"
SETS = {"dk_thick": ["lifespan_DK_46K_59sites", "site_ids_ct_dk_59sites.txt"],
        "dk_area": ["lifespan_sa_DK_46K_59sites", "site_ids_ct_dk_59sites.txt"],
        "aseg": ["lifespan_57K_82sites", "site_ids_ct_82sites.txt"]}


def to_bc(col: str) -> str:
    """Our feature name -> braincharts IDP name."""
    m = re.match(r"^(thick|area)__(lh|rh)_(.+)$", col)
    if m:
        return f"{'L' if m.group(2) == 'lh' else 'R'}_{m.group(3)}"
    return col.replace("-Thalamus", "-Thalamus-Proper")


def from_bc(idp: str, set_name: str) -> str:
    if set_name == "aseg":
        return idp.replace("-Thalamus-Proper", "-Thalamus")
    pre = "thick" if set_name == "dk_thick" else "area"
    return f"{pre}__{'lh' if idp.startswith('L_') else 'rh'}_{idp[2:]}"


def cmd_prepare(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    out = p.derivatives / "normative_braincharts"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    m = master[master.in_case_control.astype(bool)][[ID, "site", "age", "sex", "group"]].dropna().copy()
    m = m.assign(sex=m.sex.map({"F": 0, "M": 1}), control=m.group.eq("CO").astype(int),
                 sitenum=pd.factorize(m.site)[0] + 1000, site=("tinception_" + m.site))
    rng = np.random.default_rng(cfg["seed"])
    m["fold"] = 0
    for idx in m.groupby(["site", "group", "sex"]).groups.values():   # alternate within cells, random order/start
        order = rng.permutation(np.asarray(idx))
        m.loc[order, "fold"] = (np.arange(len(order)) + rng.integers(2)) % 2
    m[[ID, "age", "sex", "site", "sitenum", "fold", "control"]].to_csv(out / "subjects.tsv", sep="\t", index=False)
    fsets = cfg["normative"]["feature_sets"]
    feat = p.derivatives / "features"
    for name in SETS:
        t = ncore.load_feature_set(feat, fsets[name])   # raises on an empty match
        cols = [c for c in t.columns if c != ID]
        t.rename(columns={c: to_bc(c) for c in cols}).to_csv(out / f"responses_{name}.tsv", sep="\t", index=False,
                                                             na_rep="n/a")
        print(f"{name}: {len(cols)} features")
    (out / "sets.json").write_text(json.dumps(SETS, indent=2))
    cc = pd.crosstab(m[m.control == 1].site, m[m.control == 1].fold).reindex(sorted(m.site.unique()), fill_value=0)
    if (cc.to_numpy() == 0).any():
        raise SystemExit(f"a site x fold cell has no adaptation controls:\n{cc}")
    print(f"subjects {len(m)} ({m.control.sum()} controls), folds {m.fold.value_counts().to_dict()}, "
          f"min controls per site x fold {int(cc.to_numpy().min())}")


def cmd_analyze(args: argparse.Namespace) -> None:
    cfg, p = config.analysis(), config.paths()
    gcfg, ncfg = cfg["group"], cfg["normative"]
    d = p.derivatives / "normative_braincharts"
    out = p.results / "07d_braincharts"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    master = master[[ID, "site", "group", "PTA4", "PTA_HF", "euler_total", "eTIV"]]
    log = json.loads((d / "transfer_log.json").read_text())
    failed = {k: v["failures"] for k, v in log.items() if v["failures"]}
    if failed:
        raise SystemExit(f"transfer failures (see transfer_log.json): { {k: list(v)[:5] for k, v in failed.items()} }")
    z = None
    for name in SETS:
        zs = pd.read_csv(d / f"z_{name}.tsv", sep="\t", na_values="n/a")
        zs = zs.rename(columns={c: from_bc(c, name) for c in zs.columns if c != ID})
        z = zs if z is None else z.merge(zs, on=ID, how="outer")
    feats = [c for c in z.columns if c != ID]
    cover = z[feats].notna().mean()
    if (cover < 1).any():
        raise SystemExit(f"incomplete braincharts z: {cover[cover < 1].round(3).to_dict()}")
    z = z.merge(master, on=ID).assign(log_eTIV=lambda x: np.log(x.eTIV))
    ctrl, pat = z[z.group.eq("CO")], z[z.group.eq("TI")]
    fcal = ncore.feature_calibration(ctrl, feats, ncfg["calibration"]["feature"])
    fcal.to_csv(out / "calibration.tsv", sep="\t")
    held = ctrl[feats].to_numpy(float)
    pooled = {"mean": float(np.nanmean(held)), "sd": float(np.nanstd(held, ddof=1)),
              "n_features_failing": int((~fcal.ok).sum())}
    site_cal = ctrl.groupby("site")[feats].apply(lambda x: pd.Series({"mean": np.nanmean(x.to_numpy(float)),
                                                                      "sd": np.nanstd(x.to_numpy(float), ddof=1)}))
    # agreement with the in-house CO Model A deviations
    za = pd.read_csv(p.derivatives / "normative" / "CO" / "A" / "z_all.tsv", sep="\t", na_values="n/a")
    both = z[[ID, "group", *feats]].merge(za[[ID, *[f for f in feats if f in za]]], on=ID, suffixes=("", "_inhouse"))
    agree = []
    for f in feats:
        if f"{f}_inhouse" not in both:
            continue
        for grp, x in (("controls", both[both.group.eq("CO")]), ("patients", both[both.group.eq("TI")])):
            ok = x[[f, f"{f}_inhouse"]].dropna()
            fam = "dk_thick" if f.startswith("thick__") else "dk_area" if f.startswith("area__") else "aseg"
            agree.append({"feature": f, "family": fam, "group": grp, "n": len(ok),
                          "r": float(np.corrcoef(ok[f], ok[f"{f}_inhouse"])[0, 1])})
    agree = pd.DataFrame(agree)
    agree.to_csv(out / "agreement.tsv", sep="\t", index=False)
    # H2a burden on calibrated DK thickness + aseg
    # braincharts has no TIV term: volumes enter the burden with log eTIV as covariate; thickness alone as check
    f13 = ncore.f1f3(fcal.ok, ncfg)
    zz = z.assign(tinnitus=z.group.eq("TI").astype(int))
    nb = []
    for label, fl, covs in (("F1F3", f13, ("euler_total", "log_eTIV")),
                            ("F1_thickness", [f for f in f13 if f.startswith("thick__")], ("euler_total",))):
        b = h2.counts(zz, fl, ncfg["extreme_z"]).merge(zz[[ID, "tinnitus", "site", "euler_total", "log_eTIV"]], on=ID)
        rows = pd.DataFrame([h2.nb_burden(b, o, covariates=covs) for o in ("n_pos", "n_neg")])
        if rows.p_one_sided.isna().any():
            raise SystemExit(f"NB burden {label} did not converge: {rows[['outcome', 'error']].to_dict('records')}")
        rows["p_holm"] = _holm(rows.p_one_sided.to_numpy())   # Holm over the two directions (as 6.3)
        nb.append(rows.assign(features=label, n_features=len(fl)))
    nb = pd.concat(nb, ignore_index=True)
    nb.to_csv(out / "burden.tsv", sep="\t", index=False)
    # tinnitus effect on the braincharts z, primary sample (with PTA_HF), H1 main hearing terms
    prim = g.residualize(zz[zz.PTA_HF.notna()].reset_index(drop=True), "PTA_HF", ["PTA4"], "PTA_HF_res")
    safe = {f: f"f{i}" for i, f in enumerate(feats)}
    ps = prim.rename(columns=safe)
    for c in ("PTA4", "PTA_HF_res", "euler_total", "log_eTIV"):
        ps[c] = (ps[c] - ps[c].mean()) / ps[c].std()
    eff = []
    for f in feats:
        vol = not f.startswith("thick__")   # braincharts z has no TIV term: log eTIV for area / volume (as H1, D1)
        rhs = "tinnitus + PTA4 + PTA_HF_res + euler_total" + (" + log_eTIV" if vol else "") + " + C(site)"
        try:
            e = rb.wls_effect(ps, safe[f], rhs, "tinnitus")
        except (ValueError, KeyError, np.linalg.LinAlgError):
            e = {"beta": np.nan, "se": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan, "n": 0}
        fam = "dk_thick" if f.startswith("thick__") else "dk_area" if f.startswith("area__") else "aseg"
        eff.append({"family": fam, "outcome": f, **e})
    eff = pd.DataFrame(eff)
    eff["q"] = eff.groupby("family").p.transform(g.bh)
    apr = eff[eff.outcome.isin(gcfg["apriori"])].copy()
    apr_missing = [c for c in gcfg["apriori"] if c not in set(apr.outcome)]   # no braincharts model (MGN)
    apr["family"] = "apriori"
    apr["q"] = g.bh(apr.p)
    eff = pd.concat([eff, apr], ignore_index=True)
    eff.to_csv(out / "effects.tsv", sep="\t", index=False, na_rep="n/a")
    q_thr = gcfg["fdr_q"]
    summary = {"n_controls": len(ctrl), "n_patients": len(pat), "calibration_pooled": pooled,
               "site_mean_range": [float(site_cal["mean"].min()), float(site_cal["mean"].max())],
               "site_sd_range": [float(site_cal["sd"].min()), float(site_cal["sd"].max())],
               "agreement_median_r": {f"{f} | {g_}": float(v) for (f, g_), v in
                                      agree.groupby(["family", "group"]).r.median().items()},
               "agreement_min_r": {f"{f} | {g_}": float(v) for (f, g_), v in
                                   agree.groupby(["family", "group"]).r.min().items()},
               "apriori_family_size": len(apr), "apriori_missing": apr_missing,
               "burden": nb.to_dict("records"), "n_features_burden": len(f13),
               "effects_hits": eff[eff.q < q_thr][["family", "outcome", "beta", "p", "q"]].to_dict("records"),
               "effects_min_q": eff.groupby("family").q.min().to_dict(),
               "apriori": apr[["outcome", "beta", "ci_lo", "ci_hi", "p", "q"]].to_dict("records"),
               "apriori_max_abs_beta": float(apr.beta.abs().max())}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(json.dumps(g.json_clean({k: v for k, v in summary.items() if k != "apriori"}), indent=1, default=str))


def _holm(p: np.ndarray) -> np.ndarray:
    import pingouin as pg

    return pg.multicomp(p, method="holm")[1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prepare").set_defaults(func=cmd_prepare)
    sub.add_parser("analyze").set_defaults(func=cmd_analyze)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
