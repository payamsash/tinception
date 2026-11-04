"""Pretrained brain-age gap (analysis plan 6.11 A; fixed 2026-11-04 before running).

Inputs: <derivatives>/pyment/predictions_{multi-2025,reg-2025}.csv (slurm/19_brainage_*.sbatch: the author's
containers, FastSurfer preprocessing, both published brain-age models).
Primary model = the one whose raw gap (predicted - age) is less correlated with age in the controls (author
recommendation: least biased). Bias correction: gap regressed on age + age^2 in the controls of the other
site-stratified fold (2-fold cross-fitting) and the fit removed from everyone. Test: corrected gap with the H1 main
model (tinnitus + PTA4 + PTA_HF beyond PTA4 + age spline + sex + Euler, site random intercept; outcome standardized,
so beta = adjusted d), TOST |d| < .2 and BF01; also the PTA4-only sample. Quality: MAE and r(predicted, age) in
controls, overall and per site. -> <results>/19_brainage/{quality.tsv, models.tsv, summary.json}
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import evidence as ev
from tinception.stats import group as g

ID = "participant_id"
MODELS = ("multi-2025", "reg-2025")


def load_predictions(d) -> pd.DataFrame:
    out = None
    for m in MODELS:
        f = d / f"predictions_{m}.csv"
        p = pd.read_csv(f)
        p = pd.DataFrame({ID: "sub-" + p.subject.astype(str), f"pred_{m}": p.age.astype(float)})
        if p[ID].duplicated().any():
            raise SystemExit(f"{f}: duplicated subjects")
        out = p if out is None else out.merge(p, on=ID, how="outer", validate="one_to_one")
    return out


def cross_fit_correct(df: pd.DataFrame, pred: str, seed: int) -> pd.Series:
    """Bias-corrected gap: (pred - age) minus the fit of gap ~ age + age^2 estimated in the controls of the other
    fold (folds alternate within site x group in random order)."""
    rng = np.random.default_rng(seed)
    fold = pd.Series(0, index=df.index)
    for idx in df.groupby(["site", "group"]).groups.values():
        o = rng.permutation(np.asarray(idx))
        fold.loc[o] = (np.arange(len(o)) + rng.integers(2)) % 2
    gap = df[pred] - df.age
    X = np.column_stack([np.ones(len(df)), df.age, df.age ** 2])
    out = pd.Series(np.nan, index=df.index)
    for k in (0, 1):
        tr = (fold != k) & df.group.eq("CO") & gap.notna()
        b = np.linalg.lstsq(X[tr.to_numpy()], gap[tr].to_numpy(), rcond=None)[0]
        te = (fold == k).to_numpy()
        out[te] = gap[te] - X[te] @ b
    return out


def main() -> None:
    cfg, p = config.analysis(), config.paths()
    gcfg = cfg["group"]
    out = p.results / "19_brainage"
    out.mkdir(parents=True, exist_ok=True)
    pr = load_predictions(p.derivatives / "pyment")
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    d = master[master.in_case_control.astype(bool)].merge(pr, on=ID, how="left")
    if d.age.isna().any():
        raise SystemExit(f"{int(d.age.isna().sum())} case-control subjects without age")
    missing = int(d[[f"pred_{m}" for m in MODELS]].isna().any(axis=1).sum())
    d = d.dropna(subset=[f"pred_{m}" for m in MODELS]).reset_index(drop=True)

    # quality and model choice (controls)
    co = d[d.group.eq("CO")]
    qual, bias = [], {}
    for m in MODELS:
        gap = co[f"pred_{m}"] - co.age
        bias[m] = float(np.corrcoef(gap, co.age)[0, 1])
        qual.append({"model": m, "site": "all", "n": len(co), "mae": float(np.abs(gap).mean()),
                     "r_pred_age": float(np.corrcoef(co[f"pred_{m}"], co.age)[0, 1]), "mean_gap": float(gap.mean()),
                     "r_gap_age": bias[m]})
        for s_, x in co.groupby("site"):
            gx = x[f"pred_{m}"] - x.age
            qual.append({"model": m, "site": s_, "n": len(x), "mae": float(np.abs(gx).mean()),
                         "r_pred_age": float(np.corrcoef(x[f"pred_{m}"], x.age)[0, 1]) if len(x) > 2 else np.nan,
                         "mean_gap": float(gx.mean()), "r_gap_age": np.nan})
    qual = pd.DataFrame(qual)
    qual.to_csv(out / "quality.tsv", sep="\t", index=False)
    primary = min(bias, key=lambda k: abs(bias[k]))

    # corrected gaps and the H1-style test
    d = d.assign(tinnitus=d.group.eq("TI").astype(int),
                 **{f"gap_{m}": cross_fit_correct(d, f"pred_{m}", cfg["seed"]) for m in MODELS})
    cov = [f"C({c})" if c == "sex" else c for c in gcfg["covariates"]]
    rows = []
    for sample in ("primary", "pta4only"):
        s_ = d[d.PTA_HF.notna()] if sample == "primary" else d
        s_ = s_.reset_index(drop=True)
        terms = ["tinnitus", "PTA4", *cov]
        cont = ["PTA4", "euler_total"]
        if sample == "primary":
            s_ = g.residualize(s_, "PTA_HF", ["PTA4"], "PTA_HF_res")
            terms.insert(2, "PTA_HF_res")
            cont.append("PTA_HF_res")
        cols = [f"gap_{m}" for m in MODELS]
        pp = g.prepare(s_, cols, set())
        rhs = {y: g.rhs(terms, gcfg["age_df"], False) for y in pp.attrs["outcomes"]}
        r = g.run_family(pp, rhs, {"tinnitus": {"tinnitus": 1.0}}, sample, gcfg["tost_d"], cont)
        inv = {v: k for k, v in pp.attrs["outcomes"].items()}
        for _, x in r.iterrows():
            m = x.outcome.replace("gap_", "")
            sd = float(g.fit_frame(pp, inv[x.outcome], rhs[inv[x.outcome]])[x.outcome].std())   # analysed rows
            rows.append({"sample": sample, "model": m, "primary_model": m == primary, "n": int(x.n), "d": x.beta,
                         "ci_lo": x.ci_lo, "ci_hi": x.ci_hi, "p": x.p, "p_tost": x.p_tost,
                         "bf01": float(np.exp(-ev.log_bf10(x.beta, x.se))), "years_per_d": sd,
                         "diff_years": x.beta * sd, "converged": bool(x.converged)})
    res = pd.DataFrame(rows)
    res.to_csv(out / "models.tsv", sep="\t", index=False)
    summary = {"n_case_control_with_predictions": len(d), "missing_predictions": missing, "primary_model": primary,
               "r_gap_age_controls": bias, "quality_all": qual[qual.site == "all"].to_dict("records"),
               "tests": res.to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(qual[qual.site == "all"].round(3).to_string(index=False))
    print("primary:", primary)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
