"""Step 08 group models (analysis plan 6.1); settings in analysis.yaml -> group.

Subcommands
  h1   H1a-H1d at ROI level -> <results>/08_group/h1/:
         models.tsv     every outcome x model x contrast (beta = adjusted d, CI, p, q within family x contrast,
                        partial r, TOST p, convergence)
         h1d.tsv        H1d: regions with a tinnitus effect in the PTA4-only model (q < .05, any family; one row per
                        region): attenuation after PTA_HF within the primary sample (and after PTA_EHF in its subset),
                        equivalence, controls-only PTA_HF effect, classification (point estimates)
         tiers.tsv      three-tier verdict per subregion (specific / consistent with a whole-structure effect)
         spin.json      H1c spatial correspondence of DK thickness maps (spin test)
         summary.json   verdict per hypothesis
  h2       H2 (plan 6.3/6.4) on the CO-reference normative z (calibrated features only) -> <results>/08_group/h2/:
           h2a.tsv (NB burden tests: models A/B4/B x F1+F3 / all ROI / vertex / voxel, with and without london_1),
           h2b_regions.tsv, summary.json (H2a, H2b network overlap vs per-patient spin null, ladder A -> B4 -> B)
  table1   TI vs CO descriptives per sample (pingouin: Welch t + Hedges g, chi-square + Cramer V) -> h1/table1.tsv
  h1d-tasks / h1d-boot / h1d-ci   bootstrap 95% CIs of the H1d attenuations (5,000 resamples within site; cluster)
Hearing terms: PTA4 (overall hearing) and PTA_HF_res = PTA_HF residualized on PTA4 within the sample (r = 0.86), with a
joint Wald test of both ("hearing_joint"); the tinnitus coefficient is the same as with raw PTA_HF.

Samples (all in_case_control): primary = with PTA_HF (main model, factorial); pta4only = all (PTA4-only model, H1a);
clean = primary with group4_hf in {CTRL, TINN}; controls = primary CO (PTA_HF map); ehf = primary with PTA_EHF.

Example (cluster, repo root):  TINCEPTION_ENV=cluster uv run scripts/08_group.py h1
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

from tinception import config
from tinception.normative import core as ncore
from tinception.stats import group as g
from tinception.stats import h2


def covariates() -> list[str]:
    """Standard covariates besides the age spline (analysis.yaml group.covariates; sex is categorical)."""
    return [f"C({c})" if c == "sex" else c for c in config.analysis()["group"]["covariates"]]


def load(args: argparse.Namespace) -> tuple[pd.DataFrame, dict, dict]:
    """Case-control subjects with every outcome; outcome -> volumetric flag; subregion tier info."""
    p, cfg = config.paths(), config.analysis()
    gcfg, fsets = cfg["group"], cfg["normative"]["feature_sets"]
    master = pd.read_csv(args.master or p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    feat_dir = args.features_dir or p.derivatives / "features"
    m = master[master.in_case_control.astype(bool)][[g.ID, "site", "sex", "age", "group", "group4", "group4_hf",
                                                     "PTA4", "PTA_HF", "PTA_EHF", "HL", "euler_total", "eTIV"]].copy()
    m["tinnitus"] = m.group.eq("TI").astype(int)
    m["HL"] = m.HL.astype(float)
    families, volumetric = {}, set()
    tables = {}
    for name in gcfg["families"]:
        spec = fsets[name]
        t = tables.setdefault(spec["table"], pd.read_csv(feat_dir / f"{spec['table']}.tsv", sep="\t", na_values="n/a"))
        cols = ncore.select_features(t, spec["include"], spec.get("exclude"))
        families[name] = cols
        if spec["tiv"]:
            volumetric |= set(cols)
    sub = tables.setdefault("subseg", pd.read_csv(feat_dir / "subseg.tsv", sep="\t", na_values="n/a"))
    tiers = {}
    for struct, spec in gcfg["subregions"].items():
        for c, info in g.subregion_frame(sub, spec["pattern"], spec["whole"], spec["parts_exclude"]).items():
            tiers[c] = info | {"structure": struct}
    wholes = sorted({v["whole"] for v in tiers.values()})
    volumetric |= set(tiers) | set(wholes)
    apriori = gcfg["apriori"]
    volumetric |= {c for c in apriori if not c.startswith("thick__")}
    cols = {g.ID: sub[g.ID]}
    for i, info in enumerate(tiers.values()):  # tier 3 covariate: log(whole - subregion), same segmentation
        info["rest_col"] = f"rest{i}"
        cols[f"rest{i}"] = np.log(info["rest"].where(info["rest"] > 0))
    rest = pd.DataFrame(cols)
    df = m
    for t in [*tables.values(), rest]:
        df = df.merge(t, on=g.ID, how="left", validate="one_to_one")
    df = df.assign(log_eTIV=np.log(df.eTIV)).copy()
    return df, {"families": families, "apriori": apriori, "tiers": tiers, "wholes": wholes,
                "volumetric": volumetric}, gcfg


def sample(df: pd.DataFrame, name: str) -> pd.DataFrame:
    s = {"pta4only": df, "primary": df[df.PTA_HF.notna()], "ehf": df[df.PTA_HF.notna() & df.PTA_EHF.notna()]}
    s["clean"] = s["primary"][s["primary"].group4_hf.isin(["CTRL", "TINN"])]
    s["controls"] = s["primary"][s["primary"].group.eq("CO")]
    out = s[name].reset_index(drop=True)
    # PTA_HF beyond PTA4 (r = 0.86): residualized within the analysis sample (plan 6.1, 2026-10-30)
    return g.residualize(out, "PTA_HF", ["PTA4"], "PTA_HF_res") if out.PTA_HF.notna().any() else out


CONTINUOUS = ("PTA4", "PTA_HF", "PTA_HF_res", "PTA_EHF", "euler_total", "log_eTIV")
HEARING = ["PTA4", "PTA_HF_res"]           # PTA4 = overall hearing; PTA_HF_res = high-frequency loss beyond PTA4


def run(df_s: pd.DataFrame, outcomes: list[str], terms: list[str], contrasts: dict, family: str, spec: dict,
        gcfg: dict, rest: bool = False, tost: bool = True, joints: dict | None = None,
        age_df: int | None = None) -> pd.DataFrame:
    """Log-transform volumetric outcomes and fit one family under one model (standardization inside fit)."""
    d = g.prepare(df_s, outcomes, spec["volumetric"])
    inv = {v: k for k, v in d.attrs["outcomes"].items()}
    rhs = {inv[o]: g.rhs(list(terms) + ([spec["tiers"][o]["rest_col"]] if rest else []), age_df or gcfg["age_df"],
                         o in spec["volumetric"]) for o in outcomes}
    cont = list(CONTINUOUS) + ([spec["tiers"][o]["rest_col"] for o in outcomes] if rest else [])
    return g.run_family(d, rhs, contrasts, family, gcfg["tost_d"] if tost else None, cont, joints)


def models_h1(df: pd.DataFrame, spec: dict, gcfg: dict) -> pd.DataFrame:
    cov = covariates()
    one = {"tinnitus": {"tinnitus": 1.0}}
    hear = {"tinnitus": {"tinnitus": 1.0}, "PTA4": {"PTA4": 1.0}, "PTA_HF_beyond_PTA4": {"PTA_HF_res": 1.0}}
    hjoint = {"hearing_joint": HEARING}
    factorial = {"TINN_vs_CTRL": {"tinnitus": 1.0}, "TNHL_vs_HLOS": {"tinnitus": 1.0, "tinnitus:HL": 1.0},
                 "HLOS_vs_CTRL": {"HL": 1.0}, "tinnitus_main": {"tinnitus": 1.0, "tinnitus:HL": 0.5},
                 "interaction": {"tinnitus:HL": 1.0}}
    prim, pta4, clean = (sample(df, s) for s in ("primary", "pta4only", "clean"))
    designs = {"main": (prim, ["tinnitus", *HEARING, *cov], hear, hjoint),
               "pta4only": (pta4, ["tinnitus", "PTA4", *cov], {"tinnitus": {"tinnitus": 1.0}, "PTA4": {"PTA4": 1.0}},
                            None),
               "factorial": (prim, ["tinnitus * HL", "PTA_HF_res", *cov], factorial, None),
               "clean_hearing": (clean, ["tinnitus", *HEARING, *cov], one, None)}
    fitted = set().union(*spec["families"].values())
    extra = [c for c in spec["apriori"] if c not in fitted]  # a priori measures not in any family (MGN)
    res = []
    for model, (d, terms, con, jt) in designs.items():
        for fam, cols in spec["families"].items():
            res.append(run(d, cols, terms, con, fam, spec, gcfg, joints=jt).assign(model=model))
        res.append(run(d, spec["wholes"] + extra, terms, con, "subregion_whole+", spec, gcfg, joints=jt)
                   .assign(model=model))
    out = pd.concat(res, ignore_index=True)
    # a priori family assembled from those fits, FDR within it
    apr = out[out.outcome.isin(spec["apriori"])].drop_duplicates(["model", "outcome", "contrast"]).copy()
    apr["family"] = "apriori"
    apr["q"] = apr.groupby("model", group_keys=False).apply(g.requery)
    keep = out.family != "subregion_whole+"
    whole = out[~keep & out.outcome.isin(spec["wholes"])].copy()
    whole["family"] = "subregion_whole"
    whole["q"] = whole.groupby("model", group_keys=False).apply(g.requery)
    structs = {}
    for c, info in spec["tiers"].items():
        structs.setdefault(f"{info['structure']}_{info['hemi']}", []).append(c)
    tiers = []
    for fam, cols in structs.items():  # three-tier rule: tier 2 absolute, tier 3 + rest of structure
        tiers.append(run(prim, cols, designs["main"][1], hear, f"tier2_{fam}", spec, gcfg).assign(model="main"))
        tiers.append(run(prim, cols, designs["main"][1], hear, f"tier3_{fam}", spec, gcfg, rest=True)
                     .assign(model="main_specificity"))
    # sensitivity (plan 6.8, added 2026-10-30): a priori family without Euler; H1a with linear age
    sens = [run(prim, spec["apriori"], ["tinnitus", *HEARING, *[c for c in cov if c != "euler_total"]], hear,
                "apriori", spec, gcfg).assign(model="main_no_euler"),
            run(pta4, gcfg["h1a"], ["tinnitus", "PTA4", *cov], {"PTA4": {"PTA4": 1.0}}, "h1a", spec, gcfg,
                tost=False, age_df=1).assign(model="pta4only_age_linear")]
    return pd.concat([out[keep], whole, apr, *tiers, *sens], ignore_index=True)


def sel(models: pd.DataFrame, model: str, fam: str, con: str) -> pd.DataFrame:
    return models[(models.model == model) & (models.family == fam) & (models.contrast == con)].set_index("outcome")


def h1d_table(df: pd.DataFrame, models: pd.DataFrame, spec: dict, gcfg: dict, ctrl_fit: pd.DataFrame) -> pd.DataFrame:
    """H1d rows, one per region with q < .05 for tinnitus in the PTA4-only model in any family."""
    cov, one = covariates(), {"tinnitus": {"tinnitus": 1.0}}
    sig = models[(models.model == "pta4only") & (models.contrast == "tinnitus") & (models.q < gcfg["fdr_q"])]
    sig = sig.sort_values("q").drop_duplicates("outcome")
    prim, ehf = sample(df, "primary"), sample(df, "ehf")
    full = ["tinnitus", "PTA4", "PTA_HF", *cov]
    rows = []
    for r in sig.itertuples():
        main = models[(models.model == "main") & (models.outcome == r.outcome) & (models.contrast == "tinnitus")
                      & ~models.family.isin(["apriori"])].iloc[0]
        base = run(prim, [r.outcome], ["tinnitus", "PTA4", *cov], one, "h1d", spec, gcfg, tost=False).iloc[0]
        e0 = run(ehf, [r.outcome], full, one, "ehf", spec, gcfg, tost=False).iloc[0]
        e1 = run(ehf, [r.outcome], [*full, "PTA_EHF"], one, "ehf", spec, gcfg, tost=False).iloc[0]
        c = ctrl_fit.loc[r.outcome]
        att = 1 - main.beta / base.beta
        indep = main.q < gcfg["fdr_q"] and att < 0.5
        hf = main.q >= gcfg["fdr_q"] and main.p_tost < 0.05 and c.p < 0.05 and np.sign(c.beta) == np.sign(base.beta)
        rows.append({"family": r.family, "outcome": r.outcome, "beta_pta4only_full_set": r.beta,
                     "q_pta4only_full_set": r.q, "beta_pta4only_primary": base.beta, "beta_main": main.beta,
                     "q_main": main.q, "p_tost_main": main.p_tost, "attenuation": att,
                     "controls_PTA_HF_beta": c.beta, "controls_PTA_HF_p": c.p,
                     "beta_ehf_without": e0.beta, "beta_ehf_with": e1.beta, "attenuation_ehf": 1 - e1.beta / e0.beta,
                     "n_ehf": int(e1.n),
                     "class": "hearing-independent" if indep else "explained by HF hearing" if hf else "ambiguous"})
    return pd.DataFrame(rows)


def tier_verdicts(models: pd.DataFrame, gcfg: dict) -> pd.DataFrame:
    t2 = models[models.family.str.startswith("tier2_") & (models.contrast == "tinnitus")].set_index("outcome")
    t3 = models[models.family.str.startswith("tier3_") & (models.contrast == "tinnitus")].set_index("outcome")
    v = pd.DataFrame({"structure": t2.family.str.replace("tier2_", ""), "beta_tier2": t2.beta, "q_tier2": t2.q,
                      "beta_tier3": t3.beta.reindex(t2.index), "q_tier3": t3.q.reindex(t2.index)})
    q = gcfg["fdr_q"]
    v["verdict"] = np.where(v.q_tier3 < q, "specific",
                            np.where(v.q_tier2 < q, "consistent with a whole-structure effect", "none"))
    return v.reset_index()


def cmd_h1(args: argparse.Namespace) -> None:
    df, spec, gcfg = load(args)
    cfg = config.analysis()
    out = (args.out or config.paths().results / "08_group") / "h1"
    out.mkdir(parents=True, exist_ok=True)
    cov, q = covariates(), gcfg["fdr_q"]
    ctrl = sample(df, "controls")
    models = models_h1(df, spec, gcfg)
    # controls-only PTA_HF effect, fitted once for every outcome (spin map + H1d)
    every = sorted(set().union(*spec["families"].values(), spec["wholes"], spec["apriori"]))
    ctrl_fit = run(ctrl, every, ["PTA_HF", *cov], {"PTA_HF": {"PTA_HF": 1.0}}, "controls", spec, gcfg,
                   tost=False).set_index("outcome")
    models = pd.concat([models, ctrl_fit.reset_index().assign(model="controls_PTA_HF")], ignore_index=True)
    models.to_csv(out / "models.tsv", sep="\t", index=False, na_rep="n/a")
    h1d = h1d_table(df, models, spec, gcfg, ctrl_fit)
    h1d.to_csv(out / "h1d.tsv", sep="\t", index=False, na_rep="n/a")
    tiers = tier_verdicts(models, gcfg)
    tiers.to_csv(out / "tiers.tsv", sep="\t", index=False, na_rep="n/a")

    # H1a: PTA4 effect on transverse temporal thickness, one-sided negative, PTA4-only model
    h1a = sel(models, "pta4only", "dk_thick", "PTA4").loc[gcfg["h1a"]]
    h1a_p = st.norm.cdf(h1a.z)
    # H1c spatial: DK thickness maps, tinnitus (main) vs PTA4 (main) and vs PTA_HF (controls only)
    coords, hemi, names = g.parcel_centroids(config.paths().freesurfer / "fsaverage")
    cols = [f"thick__{n}" for n in names]
    t_tin = sel(models, "main", "dk_thick", "tinnitus").loc[cols].z.to_numpy()
    maps = {"tinnitus_vs_PTA4": sel(models, "main", "dk_thick", "PTA4").loc[cols].z.to_numpy(),
            "tinnitus_vs_PTA_HF_controls": ctrl_fit.loc[cols].z.to_numpy()}
    spin = {k: g.spin_correlation(t_tin, v, coords, hemi, gcfg["spin_n"], cfg["seed"]) for k, v in maps.items()}
    for v in spin.values():  # not spatially equivalent: upper 95% bound (spin-null width) below spatial_r_max
        v["not_equivalent"] = bool(v["ci_hi"] < gcfg["spatial_r_max"])
    (out / "spin.json").write_text(json.dumps(g.json_clean(spin), indent=2, allow_nan=False))

    apr = sel(models, "main", "apriori", "tinnitus")
    fac = models[(models.model == "factorial") & (models.family == "apriori")]
    cln = sel(models, "clean_hearing", "apriori", "tinnitus")
    mgn = tiers.set_index("outcome")
    summary = {
        "samples": {s: len(sample(df, s)) for s in ("primary", "pta4only", "clean", "controls", "ehf")},
        "nonconverged_fits": models.loc[~models.converged.astype(bool), ["model", "family", "outcome"]]
        .drop_duplicates().to_dict("records"),
        "H1a": {"outcomes": gcfg["h1a"], "beta": h1a.beta.tolist(), "p_one_sided": h1a_p.tolist(),
                "supported": bool((h1a_p < 0.05).all() and (h1a.beta < 0).all()),
                "sensitivity_age_linear": sel(models, "pta4only_age_linear", "h1a", "PTA4")
                .pipe(lambda t: {"beta": t.beta.tolist(), "p_one_sided": st.norm.cdf(t.z).tolist()}),
                "hearing_joint_main": models[(models.model == "main") & (models.family == "dk_thick")
                                             & (models.contrast == "hearing_joint")].set_index("outcome")
                .loc[gcfg["h1a"], ["chi2", "p"]].to_dict("index")},
        "H1b": {"significant": apr.index[apr.q < q].tolist(),
                "equivalent_to_null": apr.index[(apr.q >= q) & (apr.p_tost < 0.05)].tolist(),
                "mgn_specific": {o: mgn.loc[o, "verdict"] for o in apr.index[apr.q < q] if o in mgn.index},
                "supported": bool((apr.q < q).any())},
        "H1c": {"factorial_tinnitus_main": fac[(fac.contrast == "tinnitus_main") & (fac.q < q)].outcome.tolist(),
                "factorial_interaction": fac[(fac.contrast == "interaction") & (fac.q < q)].outcome.tolist(),
                "clean_hearing": cln.index[cln.q < q].tolist(),
                "clean_hearing_equivalent_to_null": cln.index[(cln.q >= q) & (cln.p_tost < 0.05)].tolist(),
                "spatial": spin},
        "H1d": {"n_tested": len(h1d), "classes": h1d["class"].value_counts().to_dict() if len(h1d) else {},
                "bootstrap_ci": "pending: h1d-tasks / h1d-boot / h1d-ci"},
        "subregions": tiers.verdict.value_counts().to_dict(),
        "sensitivity_no_euler_H1b": sel(models, "main_no_euler", "apriori", "tinnitus")
        .pipe(lambda t: {"significant": t.index[t.q < q].tolist(), "max_abs_d": float(t.beta.abs().max())}),
        "hearing_joint_q05": {f: int(((models.model == "main") & (models.family == f)
                                      & (models.contrast == "hearing_joint") & (models.q < q)).sum())
                              for f in [*spec["families"], "subregion_whole"]},
        "tinnitus_hits_main": {f: sel(models, "main", f, "tinnitus").pipe(lambda t: t.index[t.q < q].tolist())
                               for f in [*spec["families"], "subregion_whole"]},
    }
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(json.dumps(g.json_clean(summary), indent=2))


def _held_out(master: pd.DataFrame, z: pd.DataFrame) -> pd.DataFrame:
    """Scored subjects with tinnitus flag: patients (TI) and held-out reference controls (CO)."""
    d = z.merge(master, on=g.ID, how="left", validate="one_to_one")
    d = d[d.reference.astype(bool) | d.group.eq("TI")].copy()
    d["tinnitus"] = d.group.eq("TI").astype(int)
    return d


def _wilson(k: int, n: int) -> tuple[float, float]:
    from statsmodels.stats.proportion import proportion_confint

    lo, hi = proportion_confint(k, n, method="wilson")
    return float(lo * 100), float(hi * 100)


def cmd_h2(args: argparse.Namespace) -> None:
    """H2 (plan 6.3/6.4) on the CO-reference normative z-scores; calibrated features only."""
    import pingouin as pg
    from neuromaps.nulls.spins import gen_spinsamples

    cfg = config.analysis()
    gcfg, ncfg = cfg["group"], cfg["normative"]
    thr, seed = ncfg["extreme_z"], cfg["seed"]
    p = config.paths()
    norm = p.derivatives / "normative"
    out = (args.out or p.results / "08_group") / "h2"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(args.master or p.root / "tinception_master.tsv", sep="\t", na_values="n/a",
                         low_memory=False)[[g.ID, "site", "group", "euler_total"]]
    ok = {m: pd.read_csv(norm / "CO" / m / "feature_calibration.tsv", sep="\t", index_col=0).ok for m in ncfg["models"]}
    zs = {m: _held_out(master, pd.read_csv(norm / "CO" / m / "z_all.tsv", sep="\t", na_values="n/a"))
          for m in ncfg["models"]}
    aseg = ncore.select_features(pd.DataFrame(columns=list(ok["B"].index)), ncfg["feature_sets"]["aseg"]["include"])

    def feats(model: str, which: str) -> list[str]:
        good = ok[model][ok[model]].index.tolist()  # features failing calibration are excluded (plan 6.3)
        if which == "all":
            return good
        return [f for f in good if f.startswith("thick__") or f in aseg]   # F1 + F3

    # H2a: burden, patients vs held-out controls (Model B primary; positive and negative separately, Holm over 2)
    rows = []
    common_ids = set(zs["B"][g.ID]) & set(zs["B4"][g.ID]) & set(zs["A"][g.ID])  # model rows comparable (Model B set)
    for model in ncfg["models"]:
        for which in ("F1F3", "all"):
            fl = feats(model, which)
            c = h2.counts(zs[model], fl, thr).merge(zs[model][[g.ID, "tinnitus", "site", "euler_total"]], on=g.ID)
            c = c[c[g.ID].isin(common_ids)]
            for subset, d in (("all_sites", c), ("without_london_1", c[c.site != "london_1"])):
                for o in ("n_pos", "n_neg", "n_any"):
                    rows.append({"model": model, "features": which, "n_features_set": len(fl), "subset": subset,
                                 **h2.nb_burden(d, o)})
    maps = p.derivatives / "normative" / "maps" / "CO" / "B"
    for st_ in ncfg["maps"]["sets"]:  # secondary: map level, Model B
        if not (maps / st_ / "z.npy").exists():
            continue
        z = np.load(maps / st_ / "z.npy")
        subj = pd.read_csv(maps / st_ / "subjects.tsv", sep="\t")
        good = pd.read_csv(maps / st_ / "feature_calibration.tsv", sep="\t").ok.to_numpy()
        zz = pd.DataFrame(z[:, good]).assign(**{g.ID: subj[g.ID].to_numpy(), "reference": subj.reference.to_numpy()})
        d = _held_out(master, zz)
        c = h2.counts(d, list(range(int(good.sum()))), thr).merge(d[[g.ID, "tinnitus", "site", "euler_total"]], on=g.ID)
        for o in ("n_pos", "n_neg", "n_any"):
            rows.append({"model": "B", "features": st_, "n_features_set": int(good.sum()), "subset": "all_sites",
                         **h2.nb_burden(c, o)})
    h2a = pd.DataFrame(rows)
    prim = (h2a.model == "B") & (h2a.features == "F1F3") & (h2a.subset == "all_sites") & h2a.outcome.isin(["n_pos", "n_neg"])
    h2a["p_holm_primary"] = np.nan
    pp = h2a.loc[prim, "p_one_sided"]
    if len(pp) != 2 or pp.isna().any():
        raise SystemExit(f"H2a primary needs two converged tests (n_pos, n_neg); got {h2a.loc[prim].to_dict('records')}")
    h2a.loc[prim, "p_holm_primary"] = pg.multicomp(pp.to_numpy(), method="holm")[1]
    h2a.to_csv(out / "h2a.tsv", sep="\t", index=False, na_rep="n/a")

    # H2b: region vs network overlap among patients (Model B, DK thickness), spin null per patient
    coords, hemi, names = g.parcel_centroids(p.freesurfer / "fsaverage")
    keep = [i for i, n in enumerate(names) if ok["B"].get(f"thick__{n}", False)]
    coords, hemi, names = coords[keep], hemi[keep], [names[i] for i in keep]
    net_map = h2.dk_to_yeo(p.freesurfer / "fsaverage", names)
    nets = sorted({v for v in net_map.values() if v != "none"})
    net_of = np.array([nets.index(net_map[n]) if net_map[n] in nets else -1 for n in names])
    spins = gen_spinsamples(coords, hemi, n_rotate=gcfg["spin_n"], method="vasa", seed=seed)
    cols = [f"thick__{n}" for n in names]
    zB = zs["B"].dropna(subset=cols)  # complete DK thickness maps only (NaN must not count as "not extreme")
    h2b, regions = {}, []
    for direction, sign in (("positive", 1), ("negative", -1)):
        for grp, d in (("patients", zB[zB.tinnitus == 1]), ("controls", zB[zB.tinnitus == 0])):
            flags = (sign * d[cols].to_numpy(float)) > thr
            reg, net = h2.overlap(flags, net_of, nets)
            k, n = int(flags[:, reg.argmax()].sum()), len(flags)
            res = {"n": n, "max_region": names[int(reg.argmax())], "max_region_pct": float(reg.max()),
                   "max_region_ci": _wilson(k, n), "network_pct": dict(zip(nets, net.tolist(), strict=True))}
            if grp == "patients":
                null = h2.spin_null_overlap(flags, net_of, nets, spins, gcfg["spin_n"], seed)
                pnet = (np.sum(null >= net[None, :], axis=0) + 1) / (null.shape[0] + 1)
                res |= {"network_null_mean_pct": dict(zip(nets, null.mean(axis=0).tolist(), strict=True)),
                        "network_p_spin": dict(zip(nets, pnet.tolist(), strict=True)),
                        "network_q_spin": dict(zip(nets, g.bh(pd.Series(pnet)).tolist(), strict=True)),
                        "region_low": bool(_wilson(k, n)[1] < 10)}
            h2b[f"{direction}_{grp}"] = res
            regions.append(pd.DataFrame({"region": names, "network": [net_map[x] for x in names], "pct": reg,
                                         "direction": direction, "group": grp}))
    pd.concat(regions).to_csv(out / "h2b_regions.tsv", sep="\t", index=False)

    # Ladder A -> B4 -> B (plan 6.4): patients' burden over the same subjects and features (calibrated in all three)
    common = sorted(set(feats("A", "F1F3")) & set(feats("B4", "F1F3")) & set(feats("B", "F1F3")))
    ids = pd.Series(sorted(common_ids))  # subjects scored by all three models
    burden = {m: h2.counts(zs[m][zs[m][g.ID].isin(ids)], common, thr).set_index(g.ID).n_any.reindex(ids)
              for m in ncfg["models"]}
    if any(b.isna().any() for b in burden.values()):
        raise SystemExit("ladder: a subject lacks a burden under one of the models")
    site = zs["B"].set_index(g.ID).site.reindex(ids).to_numpy()
    pat = zs["B"].set_index(g.ID).tinnitus.reindex(ids).to_numpy() == 1  # ids come from Model B, so no NaN
    ladder = {}
    for grp, mask in (("patients", pat), ("controls", ~pat)):
        a, b4, b = (burden[m].to_numpy(float)[mask] for m in ("A", "B4", "B"))
        ladder[grp] = {"mean_A": a.mean(), "mean_B4": b4.mean(), "mean_B": b.mean(),
                       "A_minus_B4": (a - b4).mean(),
                       "A_minus_B4_ci": h2.bootstrap_mean_diff(a, b4, site[mask], gcfg["bootstrap_n"], seed),
                       "B4_minus_B": (b4 - b).mean(),
                       "B4_minus_B_ci": h2.bootstrap_mean_diff(b4, b, site[mask], gcfg["bootstrap_n"], seed + 1)}
    sub = h2a[prim]
    summary = {"threshold": thr, "H2a_primary": sub[["outcome", "irr", "ci_lo", "ci_hi", "p_one_sided", "p_holm_primary",
                                                    "mean_patients", "mean_controls", "n"]].to_dict("records"),
               "H2a_supported": bool((sub.p_holm_primary < 0.05).any()),
               "H2a_nonconverged": h2a.loc[~h2a.converged.astype(bool), ["model", "features", "subset", "outcome"]]
               .to_dict("records"),
               "H2a_without_london_1": h2a[(h2a.model == "B") & (h2a.features == "F1F3") & (h2a.subset == "without_london_1")
                                           & h2a.outcome.isin(["n_pos", "n_neg"])][["outcome", "irr", "p_one_sided"]]
               .to_dict("records"),
               "H2b": h2b, "ladder_F1F3_any": ladder, "n_features_ladder": len(common)}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(json.dumps(g.json_clean(summary), indent=2))


def cmd_table1(args: argparse.Namespace) -> None:
    """Descriptive TI vs CO comparison per analysis sample (pingouin: Welch t-test with Hedges g; chi-square)."""
    import pingouin as pg

    df, _, _ = load(args)
    out = (args.out or config.paths().results / "08_group") / "h1"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in ("primary", "pta4only", "clean"):
        d = sample(df, s)
        ti, co = d[d.tinnitus == 1], d[d.tinnitus == 0]
        for v in ("age", "PTA4", "PTA_HF", "PTA_EHF", "euler_total", "eTIV"):
            a, b = ti[v].dropna(), co[v].dropna()
            if len(a) < 2 or len(b) < 2:
                continue
            t = pg.ttest(a, b, correction=True).iloc[0]
            rows.append({"sample": s, "variable": v, "TI_n": len(a), "TI_mean": a.mean(), "TI_sd": a.std(),
                         "CO_n": len(b), "CO_mean": b.mean(), "CO_sd": b.std(), "test": "Welch t", "stat": t["T"],
                         "p": t["p_val"] if "p_val" in t else t["p-val"],
                         "effect": pg.compute_effsize(a, b, eftype="hedges"), "effect_type": "Hedges g"})
        for v in ("sex", "site"):
            _, _, stats = pg.chi2_independence(d, x="tinnitus", y=v)
            r = stats[stats.test == "pearson"].iloc[0]
            rows.append({"sample": s, "variable": v, "TI_n": len(ti), "CO_n": len(co), "test": "chi-square",
                         "stat": r["chi2"], "p": r["pval"], "effect": r["cramer"], "effect_type": "Cramer V",
                         "TI_mean": ti[v].eq("F").mean() if v == "sex" else np.nan,
                         "CO_mean": co[v].eq("F").mean() if v == "sex" else np.nan})
    t1 = pd.DataFrame(rows)
    t1.to_csv(out / "table1.tsv", sep="\t", index=False, na_rep="n/a")
    print(t1.round(3).to_string(index=False))


def cmd_h1d_tasks(args: argparse.Namespace) -> None:
    out = (args.out or config.paths().results / "08_group") / "h1"
    h1d = pd.read_csv(out / "h1d.tsv", sep="\t", na_values="n/a")
    t = pd.DataFrame([(o, k) for o in h1d.outcome for k in range(args.chunks)], columns=["outcome", "chunk"])
    t.to_csv(out / "h1d_tasks.tsv", sep="\t", index=False)
    print(f"{len(t)} tasks ({len(h1d)} regions x {args.chunks} chunks) -> {out / 'h1d_tasks.tsv'}")


def cmd_h1d_boot(args: argparse.Namespace) -> None:
    df, spec, gcfg = load(args)
    cfg = config.analysis()
    out = (args.out or config.paths().results / "08_group") / "h1"
    tasks = pd.read_csv(out / "h1d_tasks.tsv", sep="\t")
    row = tasks.iloc[args.task - 1]
    n_chunks = int(tasks.chunk.max()) + 1
    n = int(np.ceil(gcfg["bootstrap_n"] / n_chunks))
    cov = covariates()
    full = ["tinnitus", "PTA4", "PTA_HF", *cov]
    vol = row.outcome in spec["volumetric"]
    draws = {}
    for name, (s, base, top) in {"attenuation": ("primary", ["tinnitus", "PTA4", *cov], full),
                                 "attenuation_ehf": ("ehf", full, [*full, "PTA_EHF"])}.items():
        d = g.prepare(sample(df, s), [row.outcome], spec["volumetric"])
        draws[name] = g.bootstrap_attenuation(d, "y0", g.rhs(base, gcfg["age_df"], vol), g.rhs(top, gcfg["age_df"], vol),
                                              list(CONTINUOUS), n, cfg["seed"] + 1000 * int(row.chunk))
    (out / "boot").mkdir(exist_ok=True)
    np.savez(out / "boot" / f"{row.outcome}_{int(row.chunk):03d}.npz", **draws)
    print(f"{row.outcome} chunk {row.chunk}: {n} resamples, NaN {int(np.isnan(draws['attenuation']).sum())}")


def cmd_h1d_ci(args: argparse.Namespace) -> None:
    out = (args.out or config.paths().results / "08_group") / "h1"
    h1d = pd.read_csv(out / "h1d.tsv", sep="\t", na_values="n/a")
    for col in ("attenuation", "attenuation_ehf"):
        lo, hi, nb = [], [], []
        for o in h1d.outcome:
            files = sorted((out / "boot").glob(f"{o}_*.npz"))
            v = np.concatenate([np.load(f)[col] for f in files]) if files else np.array([])
            v = v[np.isfinite(v)]
            lo.append(np.percentile(v, 2.5) if len(v) else np.nan)
            hi.append(np.percentile(v, 97.5) if len(v) else np.nan)
            nb.append(len(v))
        h1d[f"{col}_ci_lo"], h1d[f"{col}_ci_hi"], h1d[f"{col}_n_boot"] = lo, hi, nb
    h1d.to_csv(out / "h1d.tsv", sep="\t", index=False, na_rep="n/a")
    print(h1d[["outcome", "attenuation", "attenuation_ci_lo", "attenuation_ci_hi", "attenuation_n_boot", "class"]]
          .to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None)
    ap.add_argument("--features-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default: <results>/08_group")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("h1").set_defaults(func=cmd_h1)
    sub.add_parser("table1").set_defaults(func=cmd_table1)
    sub.add_parser("h2").set_defaults(func=cmd_h2)
    t = sub.add_parser("h1d-tasks")
    t.add_argument("--chunks", type=int, default=10, help="bootstrap chunks per region")
    t.set_defaults(func=cmd_h1d_tasks)
    b = sub.add_parser("h1d-boot")
    b.add_argument("--task", type=int, required=True, help="1-based row of h1d_tasks.tsv (SLURM_ARRAY_TASK_ID)")
    b.set_defaults(func=cmd_h1d_boot)
    sub.add_parser("h1d-ci").set_defaults(func=cmd_h1d_ci)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
