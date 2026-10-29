import numpy as np
import pandas as pd

from tinception.stats import robust as rb


def _df(n=600, seed=0, d=0.0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({"site": rng.choice(list("abcdef"), n), "sex": rng.choice(["F", "M"], n),
                       "age": rng.uniform(20, 80, n), "tinnitus": rng.integers(0, 2, n)})
    df["PTA4"] = 0.4 * df.age + rng.normal(0, 8, n) + 3 * df.tinnitus
    off = dict(zip("abcdef", rng.normal(0, 0.8, 6), strict=True))
    df["y"] = d * df.tinnitus - 0.02 * df.age + df.site.map(off) + rng.normal(0, 1, n)
    return df


def test_matching_is_within_site_and_sex_and_balances():
    m = rb.match_within_site(_df(), caliper=0.25)
    pairs = m.groupby("pair")
    assert (pairs.site.nunique() == 1).all() and (pairs.sex.nunique() == 1).all()
    assert (pairs.tinnitus.sum() == 1).all()
    smd = (m[m.tinnitus == 1].PTA4.mean() - m[m.tinnitus == 0].PTA4.mean()) / m.PTA4.std()
    assert abs(smd) < 0.1


def test_overlap_weights_range():
    df = _df()
    w = rb.overlap_weights(df, "tinnitus", ["age", "PTA4", "C(site)"])
    assert ((w > 0) & (w < 1)).all()


def test_wls_effect_cluster_and_hc3():
    df = _df(d=0.5, seed=2)
    e = rb.wls_effect(df, "y", "tinnitus + age + C(site)", "tinnitus")
    assert 0.3 < e["beta"] < 0.7 and e["ci_lo"] < e["beta"] < e["ci_hi"] and e["p"] < 0.01
    one = df[df.site == "a"]
    h = rb.wls_effect(one, "y", "tinnitus + age", "tinnitus", site=None)
    assert np.isfinite(h["se"]) and h["n"] == len(one)


def test_combat_removes_site_means_and_keeps_effect():
    df = _df(d=0.5, seed=3)
    h = rb.combat(df.assign(sex_M=df.sex.eq("M").astype(int)), ["y"], "site", ["tinnitus", "sex_M"], ["age"])
    resid = h.y - 0.5 * h.tinnitus + 0.02 * h.age
    assert resid.groupby(h.site).mean().std() < 0.1
    e = rb.wls_effect(h, "y", "tinnitus + age", "tinnitus")
    assert 0.3 < e["beta"] < 0.7


def test_dl_meta_homogeneous_and_heterogeneous():
    hom = rb.dl_meta(np.array([0.2, 0.21, 0.19, 0.2]), np.array([0.1, 0.1, 0.1, 0.1]))
    assert np.isclose(hom["pooled"], 0.2, atol=0.01) and hom["tau2"] == 0 and hom["i2"] == 0
    het = rb.dl_meta(np.array([-0.8, 0.9, -0.7, 1.0]), np.array([0.1, 0.1, 0.1, 0.1]))
    assert het["tau2"] > 0.3 and het["i2"] > 90
    assert np.isnan(rb.dl_meta(np.array([0.1]), np.array([0.1]))["pooled"])


def test_caliper_is_euclidean_in_sd_units():
    df = pd.DataFrame({"site": "a", "sex": "F", "tinnitus": [1, 0, 1, 0], "age": [50.0, 50.0, 60.0, 80.0],
                       "PTA4": [10.0, 10.5, 20.0, 20.0]})
    m = rb.match_within_site(df, caliper=0.25, seed=0)
    assert len(m) == 2 and set(m.age) == {50.0}                    # the 60 vs 80 pair is far beyond 0.25 SD


def test_overlap_weights_refuses_one_group_sites():
    df = _df()
    df.loc[df.site == "a", "tinnitus"] = 1
    import pytest

    with pytest.raises(ValueError, match="only one group"):
        rb.overlap_weights(df, "tinnitus", ["age", "C(site)"])
