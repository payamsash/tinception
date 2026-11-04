"""Network context tools on synthetic data: vectorized Spearman with spin nulls equals scipy; a map built from one
seed's connectivity profile makes that seed the top epicentre; subject epicentre flags; ENIGMA name mapping."""

import numpy as np
import pandas as pd
from scipy import stats as st

from tinception.stats import network as net


def _toy(n=40, seed=0):
    rng = np.random.default_rng(seed)
    names = [f"lh_r{i}" for i in range(n // 2)] + [f"rh_r{i}" for i in range(n // 2)]
    C = np.abs(rng.normal(size=(n, n)))
    C = (C + C.T) / 2
    np.fill_diagonal(C, 0)
    spins = np.array([rng.permutation(n) for _ in range(500)]).T   # regions x n (random permutations as "spins")
    conn = {"fc_ctx": pd.DataFrame(C, index=names, columns=names),
            "fc_sctx": pd.DataFrame(np.abs(rng.normal(size=(3, n))), index=["Left-Thalamus", "Left-Putamen",
                                                                            "Left-Caudate"], columns=names)}
    return names, conn, spins, rng


def test_spearman_null_matches_scipy():
    rng = np.random.default_rng(1)
    x, Y = rng.normal(size=30), rng.normal(size=(4, 30))
    spins = np.array([rng.permutation(30) for _ in range(20)]).T
    r, null = net.spearman_null(x, Y, spins)
    assert np.allclose(r, [st.spearmanr(x, y).statistic for y in Y])
    assert np.isclose(null[3, 2], st.spearmanr(x[spins[:, 3]], Y[2]).statistic)


def test_epicentre_recovers_seed():
    names, conn, spins, rng = _toy()
    prof = net.seed_profiles(names, conn, "fc")
    x = pd.Series(prof.loc["lh_r3"].fillna(prof.loc["lh_r3"].mean()).to_numpy() + rng.normal(0, 0.05, len(names)),
                  index=names)
    e = net.epicentres(x, prof, spins)
    top = e.sort_values("r", ascending=False).iloc[0]
    assert top.seed == "lh_r3" and top.p_spin < 0.01
    assert np.isnan(prof.loc["lh_r3", "lh_r3"]) and len(prof) == len(names) + 3


def test_subject_epicentres_flags_sign():
    names, conn, spins, rng = _toy(seed=2)
    prof = net.seed_profiles(names, conn, "fc")
    base = prof.loc["rh_r5"].fillna(prof.loc["rh_r5"].mean()).to_numpy()
    Z = np.vstack([base + rng.normal(0, 0.05, len(names)), -base + rng.normal(0, 0.05, len(names)),
                   rng.normal(size=len(names))])
    f = net.subject_epicentres(Z, {"fc": prof}, spins)["fc"]
    s = list(prof.index).index("rh_r5")
    assert f[0, s] == 1 and f[1, s] == -1 and f.shape == (3, len(prof))


def test_enigma_names_cover_dk():
    conn = net.connectomes()
    assert conn["fc_ctx"].shape == (68, 68) and "lh_transversetemporal" in conn["fc_ctx"].index
    assert set(conn["sc_sctx"].index) == set(net.SUBCORTEX.values())


def test_nan_maps_not_significant_and_subset_subcortex():
    names, conn, spins, _ = _toy(seed=4)
    prof = net.seed_profiles(names, conn, "fc")
    flat = np.zeros((1, len(names)))
    assert (net.subject_epicentres(flat, {"fc": prof}, spins)["fc"] == 0).all()
    assert np.isnan(net.perm_p(np.array([np.nan]), np.zeros((10, 1)))[0])
    sub = pd.Series([1.0, 2.0, 3.0], index=["Left-Thalamus", "Left-Putamen", "Left-Caudate"])
    conn["sc_sctx"] = conn["fc_sctx"]
    h = net.subcortical_hub_test(sub, conn, 50, 0)
    assert len(h) == 2 and h[0]["n_regions"] == 3


def test_spin_correlation_spins_nan_map():
    from tinception.stats import group as g

    rng = np.random.default_rng(0)
    spins = np.array([rng.permutation(20) for _ in range(30)]).T
    r = g.spin_correlation_spins(np.ones(20), rng.normal(size=20), spins)
    assert np.isnan(r["r"]) and np.isnan(r["p_spin"])
