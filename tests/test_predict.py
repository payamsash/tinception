"""Leave-site-out prediction tools on synthetic data: signal is recovered, a confound-only signal is removed by the
in-fold deconfounding, and permutation p is small only with signal."""

import numpy as np

from tinception.stats import predict as pr


def _data(n=240, p=20, effect=1.0, confound=0.0, seed=0):
    rng = np.random.default_rng(seed)
    site = np.repeat(np.arange(4), n // 4)
    c = rng.normal(size=n)
    y = (rng.normal(size=n) + confound * c > 0).astype(int)
    X = rng.normal(size=(n, p)) + site[:, None] * 2.0       # site offsets (removed by within-site standardization)
    X[:, :3] += effect * y[:, None] + confound * c[:, None] * 2
    return pr.site_standardize(X, site), y, c[:, None], site


def test_signal_recovered_and_permutation():
    X, y, C, site = _data(effect=1.0)
    s = pr.score(y, pr.loso_predict(X, y, C, site, "classification"), "classification", site)
    assert s["auc"] > 0.8
    null = pr.permutation_null(X, y, C, site, "classification", 19, 0)
    assert s["auc"] > null.max()


def test_confound_only_signal_removed():
    X, y, C, site = _data(effect=0.0, confound=1.5, seed=1)
    raw = pr.score(y, pr.loso_predict(X, y, np.empty((len(y), 0)), site, "classification"), "classification",
                   site)["auc"]
    adj = pr.score(y, pr.loso_predict(X, y, C, site, "classification"), "classification", site)["auc"]
    # confound-driven signal is gone (cross-validated confound regression can push AUC below .5 when the confound is
    # strongly tied to y: Snoek et al. 2019; this makes the one-sided permutation test conservative, not liberal)
    assert raw > 0.7 and adj < 0.6


def test_regression():
    rng = np.random.default_rng(2)
    site = np.repeat(np.arange(4), 60)
    X = rng.normal(size=(240, 15))
    y = X[:, 0] * 0.8 + rng.normal(0, 0.6, 240)
    s = pr.score(y, pr.loso_predict(X, y, None, site, "regression"), "regression", site)
    assert s["r"] > 0.6 and s["r2"] > 0.3


def test_site_level_target_differences_do_not_bias_within_site_r():
    # target differs strongly between sites, features carry only within-site signal: within-site r is positive
    rng = np.random.default_rng(3)
    site = np.repeat(np.arange(6), 50)
    y = rng.normal(size=300) + site * 3.0
    X = rng.normal(size=(300, 10))
    X[:, 0] += (y - site * 3.0)
    X = pr.site_standardize(X, site)
    s = pr.score(y, pr.loso_predict(X, y, None, site, "regression"), "regression", site)
    assert s["r"] > 0.5


def test_single_class_sites_excluded_and_nan_p():
    y = np.array([0, 1, 0, 1, 1, 1])
    site = np.array([0, 0, 0, 0, 1, 1])
    s = pr.score(y, np.array([.1, .9, .2, .8, .5, .6]), "classification", site)
    assert s["n_auc"] == 4 and s["auc"] == 1.0
    assert np.isnan(pr.perm_p(np.nan, np.zeros(5))) and pr.perm_p(1.0, np.array([0.0, np.nan])) == 0.5
