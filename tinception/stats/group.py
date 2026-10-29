"""Group models for H1 (analysis plan 6.1): linear mixed models per outcome, effect sizes, FDR, contrasts,
equivalence tests, the three-tier subregion rule and spin tests.

Conventions (decisions D1/D2, 2026-10-29): volumes and areas are log-transformed with log eTIV as covariate,
thickness is not; every outcome is standardized within its analysis set, so a binary group coefficient is an
adjusted Cohen's d; continuous covariates are standardized; age enters as a natural cubic spline; site is a random
intercept (statsmodels MixedLM, REML).
"""

from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats as st

ID = "participant_id"


def zscore(x: pd.Series) -> pd.Series:
    return (x - x.mean()) / x.std(ddof=1)


def prepare(df: pd.DataFrame, outcomes: list[str], volumetric: set[str]) -> pd.DataFrame:
    """Analysis frame: log-transform volumetric outcomes (D1); columns renamed y0, y1, ... for formulas (mapping in
    .attrs['outcomes']). Standardization (D2) happens in fit(), on the complete cases actually analysed."""
    out = df.copy()
    names = {}
    for i, c in enumerate(outcomes):
        v = out[c].astype(float)
        if c in volumetric:
            if (v <= 0).any():
                raise ValueError(f"{c}: non-positive values cannot be log-transformed")
            v = np.log(v)
        names[f"y{i}"] = c
        out[f"y{i}"] = v
    out.attrs["outcomes"] = names
    return out


def rhs(terms: list[str], age_df: int, log_etiv: bool) -> str:
    """Right-hand side: terms + age (natural cubic spline with age_df df; age_df = 1: linear) + log eTIV for
    volumetric outcomes."""
    age = "age" if age_df == 1 else f"cr(age, df={age_df}, constraints=\"center\")"
    return " + ".join([*terms, age, *(["log_eTIV"] if log_etiv else [])])


def residualize(df: pd.DataFrame, y: str, on: list[str], name: str) -> pd.DataFrame:
    """Add `name` = residual of y regressed (OLS, with intercept) on `on`, within df (NaN where incomplete).
    Used for PTA_HF beyond PTA4: the span of the hearing terms, hence every other coefficient, is unchanged, but the
    hearing coefficients become interpretable (PTA4 = overall hearing, residual = extra high-frequency loss)."""
    out = df.copy()
    ok = out[[y, *on]].notna().all(axis=1)
    X = np.column_stack([np.ones(ok.sum()), out.loc[ok, on].to_numpy(float)])
    beta, *_ = np.linalg.lstsq(X, out.loc[ok, y].to_numpy(float), rcond=None)
    out[name] = np.nan
    out.loc[ok, name] = out.loc[ok, y].to_numpy(float) - X @ beta
    return out


OPTIMIZERS = ("nm", "cg", "lbfgs")


def fit_frame(df: pd.DataFrame, y: str, rhs_: str, continuous=(), groups: str = "site") -> pd.DataFrame:
    """Complete cases of every formula variable (and `groups`), with y and the `continuous` covariates standardized
    on exactly those rows (D2). Shared by the mixed models and the robustness regressions."""
    used = [c for c in dict.fromkeys(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", f"{y} + {rhs_}")) if c in df]
    d = df.dropna(subset=[*used, groups]).copy()
    for c in [y, *[c for c in continuous if c in used]]:
        d[c] = zscore(d[c].astype(float))
    return d


def fit(df: pd.DataFrame, y: str, rhs_: str, continuous: list[str] = (), groups: str = "site"):
    """Mixed model y ~ rhs with a random intercept per `groups`, on the complete cases of every formula variable;
    y and the `continuous` covariates are standardized on exactly those rows (D2: a binary group coefficient is an
    adjusted Cohen's d). Returns (result, converged)."""
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import ConvergenceWarning

    d = fit_frame(df, y, rhs_, continuous, groups)
    res = None
    with warnings.catch_warnings():
        # boundary / singular-covariance warnings are expected when a site variance is truly ~0; validity is judged
        # below (converged, finite log-likelihood and finite fixed-effect SEs)
        warnings.simplefilter("ignore", ConvergenceWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        for method in OPTIMIZERS:  # L-BFGS can stall at site variance 0 with llf = inf (2026-10-29), so it is last
            try:
                r = smf.mixedlm(f"{y} ~ {rhs_}", d, groups=d[groups]).fit(reml=True, method=method)
            except (np.linalg.LinAlgError, ValueError):
                continue
            if r.converged and np.isfinite(r.llf) and np.isfinite(r.bse_fe).all():
                res = r
                break
            res = res or r
    if res is None:
        raise np.linalg.LinAlgError(f"{y}: no optimizer could fit the mixed model")
    converged = bool(res.converged and np.isfinite(res.llf) and np.isfinite(res.bse_fe).all())
    res.fixed_site = smf.ols(f"{y} ~ {rhs_} + C({groups})", d).fit()  # cross-check: site as fixed effect
    return res, converged


def contrast(res, weights: dict[str, float]) -> dict:
    """Wald estimate of a linear combination of fixed effects, with 95% CI, two-sided p and partial r."""
    names = list(res.fe_params.index)
    L = np.array([weights.get(n, 0.0) for n in names])
    missing = set(weights) - set(names)
    if missing:
        raise KeyError(f"terms not in model: {sorted(missing)}")
    beta = float(L @ res.fe_params.to_numpy())
    se = float(np.sqrt(L @ res.cov_params().loc[names, names].to_numpy() @ L))
    z = beta / se
    df_resid = res.nobs - len(names)
    out = {"beta": beta, "se": se, "z": z, "p": 2 * st.norm.sf(abs(z)), "ci_lo": beta - 1.96 * se,
           "ci_hi": beta + 1.96 * se, "partial_r": z / np.sqrt(z ** 2 + df_resid), "n": int(res.nobs),
           "site_var": float(res.cov_re.iloc[0, 0])}
    fx = getattr(res, "fixed_site", None)
    if fx is not None:  # same contrast with site as a fixed effect (cross-check of the mixed model)
        Lf = np.array([weights.get(n, 0.0) for n in fx.params.index])
        out["beta_fixed_site"] = float(Lf @ fx.params.to_numpy())
        out["se_fixed_site"] = float(np.sqrt(Lf @ fx.cov_params().to_numpy() @ Lf))
    return out


def joint(res, terms: list[str]) -> dict:
    """Joint Wald chi-square test that all `terms` are 0 (e.g. every hearing term)."""
    names = list(res.fe_params.index)
    idx = [names.index(t) for t in terms]
    b = res.fe_params.to_numpy()[idx]
    V = res.cov_params().loc[names, names].to_numpy()[np.ix_(idx, idx)]
    chi2 = float(b @ np.linalg.solve(V, b))
    return {"chi2": chi2, "df": len(terms), "p": float(st.chi2.sf(chi2, len(terms))), "n": int(res.nobs)}


def tost(beta: float, se: float, bound: float) -> float:
    """Two one-sided tests for |effect| < bound (equivalence); returns the TOST p (max of the two)."""
    return float(max(st.norm.sf((beta + bound) / se), st.norm.sf((bound - beta) / se)))


def bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values via pingouin (NaN p stays NaN)."""
    import pingouin as pg

    q = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.any():
        q[ok] = pg.multicomp(p[ok].to_numpy(), method="fdr_bh")[1]
    return q


def run_family(df: pd.DataFrame, rhs_by_outcome: dict[str, str], contrasts: dict[str, dict[str, float]],
               family: str, tost_bound: float | None = None, continuous: list[str] = (),
               joints: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Fit every outcome (prepared column y<i>) and evaluate each named contrast (and joint Wald tests); BH q within
    family x contrast, over converged fits only (a non-converged fit keeps its estimate but gets no q)."""
    rows = []
    for y, r in rhs_by_outcome.items():
        res, conv = fit(df, y, r, continuous)
        base = {"family": family, "outcome": df.attrs["outcomes"][y], "converged": conv}
        for name, w in contrasts.items():
            e = contrast(res, w)
            if tost_bound is not None:
                e["p_tost"] = tost(e["beta"], e["se"], tost_bound)
            rows.append(base | {"contrast": name, **e})
        for name, terms in (joints or {}).items():
            rows.append(base | {"contrast": name, **joint(res, terms)})
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = out.p.where(out.converged).groupby(out.contrast, group_keys=False).apply(bh)
    return out


def requery(rows: pd.DataFrame) -> pd.Series:
    """BH q recomputed within a new family (e.g. the a priori family assembled from other families' fits)."""
    return rows.p.where(rows.converged).groupby(rows.contrast, group_keys=False).apply(bh)


def subregion_frame(table: pd.DataFrame, pattern: str, whole: str, parts_exclude: str) -> dict[str, dict]:
    """{subregion column: {"hemi", "whole", "rest" series}} for the three-tier rule: rest = whole - subregion,
    both from the same segmentation."""
    out = {}
    for c in table:
        m = re.match(pattern, c)
        if not m or re.search(parts_exclude, c):
            continue
        w = whole.format(h=m.group(1))
        out[c] = {"hemi": m.group(1), "whole": w, "rest": table[w] - table[c]}
    return out


def parcel_centroids(fsaverage, annot: str = "aparc") -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Centroids (on the sphere) of every annotation parcel of both hemispheres, hemisphere ids (0 lh, 1 rh) and
    names '<hemi>_<parcel>' (medial wall / unknown / corpus callosum dropped), for spin tests."""
    import nibabel as nib

    coords, hemi, names = [], [], []
    for k, h in enumerate(("lh", "rh")):
        xyz, _ = nib.freesurfer.read_geometry(str(fsaverage / "surf" / f"{h}.sphere"))
        labels, _, pnames = nib.freesurfer.read_annot(str(fsaverage / "label" / f"{h}.{annot}.annot"))
        for i, n in enumerate(pnames):
            n = n.decode() if isinstance(n, bytes) else n
            if n in ("unknown", "corpuscallosum", "Medial_wall") or not (labels == i).any():
                continue
            c = xyz[labels == i].mean(axis=0)
            coords.append(c / np.linalg.norm(c) * 100)
            hemi.append(k)
            names.append(f"{h}_{n}")
    return np.array(coords), np.array(hemi), names


def spin_correlation(x: np.ndarray, y: np.ndarray, coords: np.ndarray, hemi: np.ndarray, n: int,
                     seed: int) -> dict:
    """Spearman r of two parcel maps with a spin-permutation p (Vasa method; rotations of x, two-sided)."""
    from neuromaps.nulls.spins import gen_spinsamples

    r = st.spearmanr(x, y).statistic
    spins = gen_spinsamples(coords, hemi, n_rotate=n, method="vasa", seed=seed)
    null = np.array([st.spearmanr(x[spins[:, k]], y).statistic for k in range(n)])
    # 95% CI of r whose width comes from the spin null (spatial autocorrelation): sd of Fisher-z under the null
    sd_z = np.std(np.arctanh(np.clip(null, -0.999999, 0.999999)), ddof=1)
    zr = np.arctanh(np.clip(r, -0.999999, 0.999999))
    return {"r": float(r), "p_spin": float((np.sum(np.abs(null) >= abs(r)) + 1) / (n + 1)), "n_spin": n,
            "ci_lo": float(np.tanh(zr - 1.96 * sd_z)), "ci_hi": float(np.tanh(zr + 1.96 * sd_z)),
            "n_eff": float(1 / sd_z ** 2 + 3)}


def json_clean(o):
    """Recursively replace float NaN/inf by None so json.dumps(..., allow_nan=False) writes strict JSON."""
    if isinstance(o, dict):
        return {k: json_clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def bootstrap_attenuation(df: pd.DataFrame, y: str, base_rhs: str, full_rhs: str, continuous: list[str], n: int,
                          seed: int, groups: str = "site") -> np.ndarray:
    """Attenuation 1 - beta_tinnitus(full) / beta_tinnitus(base) in n bootstrap resamples drawn within each site
    (stratified), both models refitted on the same resample (analysis plan 6.1, H1d)."""
    rng = np.random.default_rng(seed)
    idx = [np.flatnonzero(df[groups].to_numpy() == s) for s in pd.unique(df[groups])]
    out = np.full(n, np.nan)
    for b in range(n):
        d = df.iloc[np.concatenate([rng.choice(i, len(i)) for i in idx])].reset_index(drop=True)
        try:
            b0 = fit(d, y, base_rhs, continuous)[0].fe_params["tinnitus"]
            b1 = fit(d, y, full_rhs, continuous)[0].fe_params["tinnitus"]
            out[b] = 1 - b1 / b0
        except (np.linalg.LinAlgError, ValueError):
            continue
    return out
