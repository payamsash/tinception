import numpy as np
import pandas as pd
import pytest

from tinception.stats import group as g


def _data(n=600, seed=0, d_tin=0.5, interaction=0.0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({"participant_id": [f"s{i}" for i in range(n)], "site": rng.choice(list("abcdef"), n),
                       "sex": rng.choice(["F", "M"], n), "age": rng.uniform(20, 80, n),
                       "tinnitus": rng.integers(0, 2, n), "HL": rng.integers(0, 2, n).astype(float),
                       "euler_total": rng.normal(-50, 10, n), "eTIV": rng.normal(1.5e6, 1e5, n)})
    site_off = dict(zip("abcdef", rng.normal(0, 0.5, 6), strict=True))
    df["thick"] = (d_tin * df.tinnitus + interaction * df.tinnitus * df.HL - 0.01 * df.age
                   + df.site.map(site_off) + rng.normal(0, 1, n))
    df["vol"] = np.exp(8 + 0.1 * df.thick + 0.5 * np.log(df.eTIV / 1.5e6))
    df["log_eTIV"] = np.log(df.eTIV)
    return df


def test_prepare_logs_volumes_and_rejects_nonpositive():
    df = _data()
    d = g.prepare(df, ["thick", "vol"], {"vol"})
    assert d.attrs["outcomes"] == {"y0": "thick", "y1": "vol"}
    assert np.allclose(d.y0, df.thick) and np.allclose(d.y1, np.log(df.vol))
    with pytest.raises(ValueError, match="non-positive"):
        g.prepare(df.assign(vol=0.0), ["vol"], {"vol"})


def test_fit_standardizes_on_the_analysed_rows():
    d = g.prepare(_data(), ["thick"], set())
    d.loc[:99, "euler_total"] = np.nan                     # 100 incomplete rows
    res, _ = g.fit(d, "y0", g.rhs(["tinnitus", "euler_total"], 3, False), ["euler_total"])
    used = res.model.data.frame
    assert len(used) == len(d) - 100
    assert np.isclose(used.y0.std(), 1) and np.isclose(used.euler_total.std(), 1)


def test_effect_recovered_and_null_is_equivalent():
    d = g.prepare(_data(d_tin=0.6), ["thick"], set())
    res, conv = g.fit(d, "y0", g.rhs(["tinnitus", "C(sex)", "euler_total"], 3, False), ["euler_total"])
    e = g.contrast(res, {"tinnitus": 1.0})
    assert conv and e["beta"] > 0.3 and e["p"] < 1e-6 and e["ci_lo"] < e["beta"] < e["ci_hi"]
    d0 = g.prepare(_data(d_tin=0.0, n=2000, seed=3), ["thick"], set())
    e0 = g.contrast(g.fit(d0, "y0", g.rhs(["tinnitus", "C(sex)", "euler_total"], 3, False), ["euler_total"])[0],
                    {"tinnitus": 1.0})
    assert e0["p"] > 0.01 and g.tost(e0["beta"], e0["se"], 0.2) < 0.05


def test_factorial_contrasts_combine_terms():
    d = g.prepare(_data(d_tin=0.2, interaction=0.6, n=2000, seed=1), ["thick"], set())
    res, _ = g.fit(d, "y0", g.rhs(["tinnitus * HL", "C(sex)"], 3, False))
    b = res.fe_params
    tnhl = g.contrast(res, {"tinnitus": 1.0, "tinnitus:HL": 1.0})
    assert np.isclose(tnhl["beta"], b["tinnitus"] + b["tinnitus:HL"])
    assert tnhl["beta"] > g.contrast(res, {"tinnitus": 1.0})["beta"]
    with pytest.raises(KeyError):
        g.contrast(res, {"nope": 1.0})


def test_bh_and_tost():
    p = pd.Series([0.01, 0.04, np.nan, 0.03])
    q = g.bh(p)
    assert np.isnan(q[2]) and np.allclose(q.dropna(), [0.03, 0.04, 0.04])
    assert g.tost(0.0, 0.05, 0.2) < 0.001 and g.tost(0.25, 0.05, 0.2) > 0.5


def test_subregion_frame_rest_is_whole_minus_part():
    t = pd.DataFrame({"hippo__lh_CA1": [2.0, 3.0], "hippo__lh_tail": [1.0, 1.0], "hippo__lh_fissure": [0.1, 0.1],
                      "hippo__lh_Whole_hippocampus": [10.0, 12.0]})
    f = g.subregion_frame(t, r"^hippo__(lh|rh)_", "hippo__{h}_Whole_hippocampus", "Whole_|fissure")
    assert set(f) == {"hippo__lh_CA1", "hippo__lh_tail"}
    assert f["hippo__lh_CA1"]["rest"].tolist() == [8.0, 9.0] and f["hippo__lh_CA1"]["hemi"] == "lh"


def test_spin_correlation():
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(80, 3))
    coords = xyz / np.linalg.norm(xyz, axis=1, keepdims=True) * 100
    hemi = np.repeat([0, 1], 40)
    x = coords[:, 2] + rng.normal(0, 5, 80)                 # smooth map
    same = g.spin_correlation(x, x + rng.normal(0, 1, 80), coords, hemi, 200, seed=1)
    rand = g.spin_correlation(x, rng.normal(size=80), coords, hemi, 200, seed=1)
    assert same["r"] > 0.9 and same["p_spin"] < 0.05 and rand["p_spin"] > 0.05
    assert same["ci_lo"] < same["r"] < same["ci_hi"] and rand["ci_lo"] < 0 < rand["ci_hi"]


def test_requery_and_nonconverged_excluded_from_fdr():
    rows = pd.DataFrame({"contrast": ["t"] * 3, "p": [0.01, 0.02, 0.001], "converged": [True, True, False]})
    q = g.requery(rows)
    assert np.isnan(q[2]) and np.allclose(q[:2], [0.02, 0.02])


def test_json_clean():
    import json

    o = g.json_clean({"a": float("nan"), "b": [np.float64(1.5), np.inf], "c": np.int64(3), "d": np.bool_(True)})
    assert json.dumps(o, allow_nan=False) == '{"a": null, "b": [1.5, null], "c": 3, "d": true}'


def test_bootstrap_attenuation_runs_within_sites():
    d = g.prepare(_data(n=300).assign(PTA4=lambda x: x.tinnitus * 0.5 + np.random.default_rng(0).normal(size=300)),
                  ["thick"], set())
    a = g.bootstrap_attenuation(d, "y0", g.rhs(["tinnitus"], 3, False), g.rhs(["tinnitus", "PTA4"], 3, False),
                                ["PTA4"], 5, seed=0)
    assert a.shape == (5,) and np.isfinite(a).all()


def test_mixed_model_recovers_site_variance_and_matches_fixed_site():
    d = g.prepare(_data(n=900, seed=5), ["thick"], set())
    res, conv = g.fit(d, "y0", g.rhs(["tinnitus", "C(sex)"], 3, False))
    e = g.contrast(res, {"tinnitus": 1.0})
    assert conv and np.isfinite(res.llf) and e["site_var"] > 0.05      # simulated site SD 0.5 on raw scale
    assert abs(e["beta"] - e["beta_fixed_site"]) < 0.02 and abs(e["se"] - e["se_fixed_site"]) < 0.01


def test_residualized_pta_hf_leaves_tinnitus_coefficient_unchanged():
    rng = np.random.default_rng(9)
    df = _data(n=500, seed=9)
    df["PTA4"] = rng.normal(20, 10, 500) + 3 * df.tinnitus
    df["PTA_HF"] = 0.9 * df.PTA4 + rng.normal(0, 4, 500)
    df = g.residualize(df, "PTA_HF", ["PTA4"], "PTA_HF_res")
    assert abs(np.corrcoef(df.PTA_HF_res, df.PTA4)[0, 1]) < 1e-10
    d = g.prepare(df, ["thick"], set())
    cont = ["PTA4", "PTA_HF", "PTA_HF_res"]
    raw = g.contrast(g.fit(d, "y0", g.rhs(["tinnitus", "PTA4", "PTA_HF"], 3, False), cont)[0], {"tinnitus": 1.0})
    res_fit = g.fit(d, "y0", g.rhs(["tinnitus", "PTA4", "PTA_HF_res"], 3, False), cont)[0]
    resid = g.contrast(res_fit, {"tinnitus": 1.0})
    assert np.isclose(raw["beta"], resid["beta"], atol=1e-4) and np.isclose(raw["se"], resid["se"], atol=1e-4)
    j = g.joint(res_fit, ["PTA4", "PTA_HF_res"])
    assert j["df"] == 2 and 0 <= j["p"] <= 1


def test_bh_matches_statsmodels_and_linear_age():
    from statsmodels.stats.multitest import multipletests

    p = pd.Series(np.random.default_rng(1).uniform(size=30))
    assert np.allclose(g.bh(p), multipletests(p, method="fdr_bh")[1])
    assert g.rhs(["tinnitus"], 1, True) == "tinnitus + age + log_eTIV"


def test_spin_correlation_spins_matches_manual():
    from scipy import stats as st

    from tinception.stats import group as g

    rng = np.random.default_rng(0)
    x, y = rng.normal(size=30), rng.normal(size=30)
    spins = np.array([rng.permutation(30) for _ in range(50)]).T
    r = g.spin_correlation_spins(x, y, spins)
    null = np.array([st.spearmanr(x[spins[:, k]], y).statistic for k in range(50)])
    assert np.isclose(r["r"], st.spearmanr(x, y).statistic)
    assert np.isclose(r["p_spin"], (np.sum(np.abs(null) >= abs(r["r"])) + 1) / 51)
