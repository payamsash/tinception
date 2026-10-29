"""Robustness tools for step 08c (analysis plan 6.8): within-site matching, overlap weights, ComBat-GAM, weighted /
fixed-site regression with site-clustered SEs (t(G - 1) reference), and a DerSimonian-Laird meta-analysis of per-site
effects (pooled d, tau^2, I^2).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as st


def match_within_site(df: pd.DataFrame, group: str = "tinnitus", exact: tuple[str, ...] = ("site", "sex"),
                      on: tuple[str, ...] = ("age", "PTA4"), caliper: float = 0.25, seed: int = 0) -> pd.DataFrame:
    """1:1 greedy nearest-neighbour matching of group == 1 to group == 0 within cells of `exact`, on the
    standardized `on` variables, without replacement; a pair is accepted if their Euclidean distance in SD units is
    at most `caliper`.
    Treated subjects are visited in random order. Returns the matched rows with a `pair` id."""
    rng = np.random.default_rng(seed)
    d = df.dropna(subset=[group, *exact, *on]).copy()
    z = (d[list(on)] - d[list(on)].mean()) / d[list(on)].std(ddof=1)
    pairs = []
    for _, cell in d.groupby(list(exact)):
        t = cell.index[cell[group] == 1].to_numpy()
        c = list(cell.index[cell[group] == 0])
        for i in rng.permutation(t):
            if not c:
                break
            dist = np.sqrt(((z.loc[c].to_numpy() - z.loc[i].to_numpy()) ** 2).sum(axis=1))
            j = int(np.argmin(dist))
            if dist[j] <= caliper:
                pairs.append((i, c.pop(j)))
    rows = [i for p in pairs for i in p]
    out = d.loc[rows].copy()
    out["pair"] = np.repeat(np.arange(len(pairs)), 2)
    return out


def overlap_weights(df: pd.DataFrame, group: str, covariates: list[str]) -> pd.Series:
    """Overlap weights from a logistic propensity model (group ~ covariates): treated 1 - e, controls e."""
    import statsmodels.formula.api as smf

    one_sided = df.groupby("site")[group].nunique() if "site" in df else pd.Series(dtype=int)
    if (one_sided < 2).any():
        raise ValueError(f"sites with only one group cannot be weighted: {one_sided.index[one_sided < 2].tolist()}")
    ps = smf.logit(f"{group} ~ {' + '.join(covariates)}", df).fit(disp=0, maxiter=200).predict(df)
    if not np.isfinite(ps).all() or ((ps < 1e-6) | (ps > 1 - 1e-6)).any():
        raise ValueError("propensity scores at 0/1 (separation)")
    return pd.Series(np.where(df[group] == 1, 1 - ps, ps), index=df.index)


def combat(df: pd.DataFrame, outcomes: list[str], site: str, linear: list[str], smooth: list[str]) -> pd.DataFrame:
    """ComBat-GAM (neuroHarmonize): remove site effects from `outcomes` while preserving the covariates (linear
    terms and smooth terms such as age). Rows with missing values in any column are dropped."""
    from neuroHarmonize import harmonizationLearn

    d = df.dropna(subset=[*outcomes, site, *linear, *smooth]).copy()
    cov = d[[site, *linear, *smooth]].rename(columns={site: "SITE"}).reset_index(drop=True)
    _, adj = harmonizationLearn(d[outcomes].to_numpy(float), cov, smooth_terms=list(smooth))
    d[outcomes] = adj
    return d


def wls_effect(df: pd.DataFrame, y: str, rhs: str, term: str, weights: pd.Series | None = None,
               site: str | None = "site") -> dict:
    """(Weighted) least squares y ~ rhs. With `site`: SEs clustered by site and a t(G - 1) reference (few
    clusters); with site=None (e.g. within one site): HC3 SEs and a t(residual df) reference.
    Returns beta, se, CI, p and n."""
    import statsmodels.formula.api as smf

    d = df.copy()
    w = np.ones(len(d)) if weights is None else weights.reindex(d.index).to_numpy(float)
    model = smf.wls(f"{y} ~ {rhs}", d, weights=w, missing="drop")
    if site is None:
        res = model.fit(cov_type="HC3")
        df_t = max(int(res.df_resid), 1)
    else:
        groups = pd.factorize(d.loc[model.data.row_labels, site])[0]
        res = model.fit(cov_type="cluster", cov_kwds={"groups": groups})
        df_t = max(len(set(groups)) - 1, 1)
    b, se = float(res.params[term]), float(res.bse[term])
    tcrit = st.t.ppf(0.975, df_t)
    return {"beta": b, "se": se, "ci_lo": b - tcrit * se, "ci_hi": b + tcrit * se,
            "p": float(2 * st.t.sf(abs(b / se), df_t)), "n": int(res.nobs)}


def dl_meta(effects: np.ndarray, ses: np.ndarray) -> dict:
    """DerSimonian-Laird random-effects meta-analysis: pooled effect (95% CI), tau^2, I^2 (%), Q and its p."""
    e, v = np.asarray(effects, float), np.asarray(ses, float) ** 2
    ok = np.isfinite(e) & np.isfinite(v) & (v > 0)
    e, v = e[ok], v[ok]
    k = len(e)
    if k < 2:
        return {"k": k, "pooled": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "tau2": np.nan, "i2": np.nan,
                "q": np.nan, "p_q": np.nan}
    w = 1 / v
    fixed = (w * e).sum() / w.sum()
    q = float((w * (e - fixed) ** 2).sum())
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c)
    wr = 1 / (v + tau2)
    pooled = (wr * e).sum() / wr.sum()
    se = np.sqrt(1 / wr.sum())
    return {"k": k, "pooled": float(pooled), "ci_lo": float(pooled - 1.96 * se), "ci_hi": float(pooled + 1.96 * se),
            "tau2": float(tau2), "i2": float(max(0.0, (q - (k - 1)) / q) * 100) if q > 0 else 0.0, "q": q,
            "p_q": float(st.chi2.sf(q, k - 1))}
