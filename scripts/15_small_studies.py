"""Small-study resampling (analysis plan 6.10 A): how often do typical small tinnitus studies drawn from our data
report "significant" regional differences, which regions, with how much effect-size inflation, and how many vanish
when hearing is controlled?

Outputs in <results>/15_small_studies/: conditions.tsv (one row per n x sampling mode), regions.tsv (per condition
x region: hit rate, sign consistency, |d| of hits, full-sample d), literature.tsv, summary.json.

Example (cluster, repo root):  TINCEPTION_ENV=cluster .venv/bin/python scripts/15_small_studies.py
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import group as g
from tinception.stats import resample as rs

ID = "participant_id"


def covariates(d: pd.DataFrame) -> dict[str, np.ndarray]:
    """Per-subject numpy columns used by every draw (built once)."""
    return {"tin": d.tinnitus.to_numpy(float), "age": d.age.to_numpy(float), "sex": d.sex.eq("M").to_numpy(float),
            "tiv": d.log_eTIV.to_numpy(float), "pta4": d.PTA4.to_numpy(float), "ptahf": d.PTA_HF.to_numpy(float),
            "site": pd.factorize(d.site)[0]}


def design(cv: dict, rows: np.ndarray, hearing: bool, site_dummies: bool, tiv: bool) -> np.ndarray:
    """[1, tinnitus, age, sex, (log eTIV), (PTA4, PTA_HF beyond PTA4 within the draw), (site dummies)]; tinnitus
    is column 1. Site dummies only if the draw spans > 1 site."""
    n = len(rows)
    cols = [np.ones(n), cv["tin"][rows], cv["age"][rows], cv["sex"][rows]]
    if tiv:
        cols.append(cv["tiv"][rows])
    if hearing:
        pta4, hf = cv["pta4"][rows], cv["ptahf"][rows]
        A = np.column_stack([np.ones(n), pta4])
        cols += [pta4, hf - A @ np.linalg.lstsq(A, hf, rcond=None)[0]]   # as g.residualize, within the draw
    X = np.column_stack(cols)
    if site_dummies:
        codes = np.unique(cv["site"][rows], return_inverse=True)[1]
        if codes.max() > 0:
            X = np.column_stack([X, np.eye(codes.max() + 1)[codes][:, 1:]])
    return X


def tests(cv: dict, rows: np.ndarray, Y: dict, hearing: bool, site_dummies: bool) -> tuple[np.ndarray, np.ndarray]:
    """t and p of tinnitus for every ROI: thickness without, area / volume with log eTIV (as H1)."""
    t, p = [], []
    for tiv, y in ((False, Y["thick"][rows]), (True, Y["vol"][rows])):
        ti, pi = rs.ols_t(design(cv, rows, hearing, site_dummies, tiv), y)
        t.append(ti)
        p.append(pi)
    return np.concatenate(t), np.concatenate(p)


def main() -> None:
    cfg, p = config.analysis(), config.paths()
    scfg, seed = cfg["small_studies"], cfg["seed"]
    out = p.results / "15_small_studies"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    f, cols = g.roi_features(p.derivatives / "features", cfg["normative"], cfg["group"])
    thick = [c for c in cols if c.startswith("thick__")]
    vol = [c for c in cols if c not in thick]
    order = thick + vol
    unknown = {k: [r for r in v if r not in order] for k, v in scfg["literature_regions"].items()}
    if any(unknown.values()):
        raise SystemExit(f"literature regions not among the ROIs: {unknown}")
    lit = {k: [order.index(r) for r in v] for k, v in scfg["literature_regions"].items()}
    d = master[master.in_case_control.astype(bool) & master.PTA_HF.notna()].merge(f, on=ID)
    d = d.assign(tinnitus=d.group.eq("TI").astype(int), log_eTIV=np.log(d.eTIV))
    d = d.dropna(subset=[*cols, "group", "age", "sex", "eTIV", "PTA4", "PTA_HF", "site"]).reset_index(drop=True)
    cv = covariates(d)
    Yall = {"thick": d[thick].to_numpy(float), "vol": d[vol].to_numpy(float)}
    everyone = np.arange(len(d))
    full = {}
    for hearing in (False, True):   # full-sample reference, same OLS models with site dummies
        t, _ = tests(cv, everyone, Yall, hearing, True)
        full[hearing] = rs.d_from_t(t, int(d.tinnitus.sum()), int((1 - d.tinnitus).sum()))
    full_abs = np.abs(full[False])
    idx_ti = {s: np.flatnonzero((d.site == s) & (d.tinnitus == 1)) for s in d.site.unique()}
    idx_co = {s: np.flatnonzero((d.site == s) & (d.tinnitus == 0)) for s in d.site.unique()}
    all_ti, all_co = np.flatnonzero(d.tinnitus == 1), np.flatnonzero(d.tinnitus == 0)
    cnt = d.groupby("site").tinnitus.agg(["sum", "count"])
    m, N = len(order), scfg["n_draws"]
    crow, rrow, lrow = [], [], []
    for ci, (n, mode) in enumerate([(n, mo) for n in scfg["n_per_group"] for mo in scfg["modes"]]):
        rng = np.random.default_rng([seed, n, scfg["modes"].index(mode)])   # independent stream per condition
        within = mode == "within_site"
        ok = cnt.index[(cnt["sum"] >= n) & (cnt["count"] - cnt["sum"] >= n)].tolist() if within else None
        if within and not ok:
            print(f"n={n} within_site: no site with {n} per group, skipped", flush=True)
            continue
        hit_p, hit_q, pos, dsum = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)
        any_p = any_q = n_p = n_q = vanish = n_nan = 0
        d_hits, d_full_hits = [], []
        for _ in range(N):
            if within:
                s_ = ok[rng.integers(len(ok))]
                rows = np.concatenate([rng.choice(idx_ti[s_], n, replace=False), rng.choice(idx_co[s_], n, replace=False)])
            else:
                rows = np.concatenate([rng.choice(all_ti, n, replace=False), rng.choice(all_co, n, replace=False)])
            t1, p1 = tests(cv, rows, Yall, False, not within)
            _, p2 = tests(cv, rows, Yall, True, not within)
            n_nan += int(np.isnan(p1).any())
            sig = p1 < 0.05
            sig_q = rs.bh(p1) < 0.05
            dv = np.abs(rs.d_from_t(t1, n, n))
            hit_p += sig
            hit_q += sig_q
            pos += sig & (t1 > 0)
            dsum += np.where(sig, dv, 0)
            any_p += sig.any()
            any_q += sig_q.any()
            n_p += sig.sum()
            n_q += sig_q.sum()
            vanish += (sig & ~(p2 < 0.05)).sum()
            d_hits.extend(dv[sig].tolist())
            d_full_hits.extend(full_abs[sig].tolist())   # the same regions' full-sample |d| (paired)
        dh, dfh = np.array(d_hits), np.array(d_full_hits)
        crow.append({"n_per_group": n, "mode": mode, "n_draws": N, "eligible_sites": len(ok) if within else
                     d.site.nunique(), "draws_with_nan_p": n_nan,
                     "share_any_p05": any_p / N, "share_any_q05": any_q / N, "mean_n_p05": n_p / N,
                     "mean_n_q05": n_q / N,
                     "median_abs_d_of_hits": float(np.median(dh)) if len(dh) else np.nan,
                     "median_full_abs_d_same_regions": float(np.median(dfh)) if len(dfh) else np.nan,
                     "median_inflation_ratio": float(np.median(dh / np.maximum(dfh, 1e-3))) if len(dh) else np.nan,
                     "median_full_abs_d_all_rois": float(np.median(full_abs)),
                     "share_hits_vanish_with_hearing": vanish / n_p if n_p else np.nan})
        for i, r in enumerate(order):
            rrow.append({"n_per_group": n, "mode": mode, "region": r, "hit_rate_p05": hit_p[i] / N,
                         "hit_rate_q05": hit_q[i] / N, "share_positive_of_hits": pos[i] / hit_p[i] if hit_p[i]
                         else np.nan, "mean_abs_d_of_hits": dsum[i] / hit_p[i] if hit_p[i] else np.nan,
                         "full_d": float(full[False][i]), "full_d_hearing": float(full[True][i])})
        for k, ii in lit.items():
            lrow.append({"n_per_group": n, "mode": mode, "source": k, "n_regions": len(ii),
                         "mean_hit_rate_p05": float(np.mean(hit_p[ii] / N)), "max_hit_rate_p05": float(np.max(hit_p[ii] / N))})
        print(json.dumps(g.json_clean(crow[-1])), flush=True)
    cond, reg, lt = pd.DataFrame(crow), pd.DataFrame(rrow), pd.DataFrame(lrow)
    cond.to_csv(out / "conditions.tsv", sep="\t", index=False)
    reg.to_csv(out / "regions.tsv", sep="\t", index=False)
    lt.to_csv(out / "literature.tsv", sep="\t", index=False)
    top = reg.sort_values("hit_rate_p05", ascending=False).groupby(["n_per_group", "mode"]).head(5)
    summary = {"n_primary": len(d), "n_rois": len(order), "conditions": cond.to_dict("records"),
               "top_regions": top[["n_per_group", "mode", "region", "hit_rate_p05", "share_positive_of_hits",
                                   "mean_abs_d_of_hits", "full_d"]].to_dict("records"),
               "literature": lt.to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(cond.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
