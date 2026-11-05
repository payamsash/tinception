"""Small-study resampling helpers: OLS t equals statsmodels, BH equals pingouin, d from t."""

import numpy as np
import pingouin as pg
import statsmodels.api as sm

from tinception.stats import resample as rs


def test_ols_t_and_bh():
    rng = np.random.default_rng(0)
    X = np.column_stack([np.ones(80), rng.integers(0, 2, 80), rng.normal(size=80)])
    Y = rng.normal(size=(80, 5)) + X[:, [1]] * 0.5
    t, p = rs.ols_t(X, Y)
    ref = sm.OLS(Y[:, 2], X).fit()
    assert np.isclose(t[2], ref.tvalues[1]) and np.isclose(p[2], ref.pvalues[1])
    pv = rng.uniform(size=20)
    assert np.allclose(rs.bh(pv), pg.multicomp(pv, method="fdr_bh")[1])
    assert np.isclose(rs.d_from_t(np.array([2.0]), 30, 30)[0], 2 * np.sqrt(2 / 30))


def test_rank_deficient_and_nan_bh():
    rng = np.random.default_rng(1)
    X = np.column_stack([np.ones(40), rng.integers(0, 2, 40), np.zeros(40)])     # constant (all-zero) column
    t, _ = rs.ols_t(X, rng.normal(size=(40, 3)))
    ref = sm.OLS(rng.normal(size=40), X[:, :2]).fit()
    assert np.isfinite(t).all() and ref.df_resid == 38
    Xc = np.column_stack([np.ones(40), np.r_[np.ones(20), np.zeros(20)], np.r_[np.ones(20), np.zeros(20)]])
    assert np.isnan(rs.ols_t(Xc, rng.normal(size=(40, 2)))[0]).all()            # tinnitus = a site dummy
    q = rs.bh(np.array([0.01, np.nan, 0.04, 0.5]))
    assert np.isnan(q[1]) and np.allclose(q[[0, 2, 3]], pg.multicomp(np.array([0.01, 0.04, 0.5]), method="fdr_bh")[1])
