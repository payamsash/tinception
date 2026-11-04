"""Co-deviation tools on synthetic data: centring per cell, corr = numpy, global test detects a coupling difference
and is null without one, edge test finds the changed pair, matrix spin correlation recovers an embedded connectome."""

import numpy as np

from tinception.stats import codeviation as cd


def _groups(n=200, p=12, couple=0.0, seed=0):
    rng = np.random.default_rng(seed)
    site = np.tile(np.arange(4), n // 4)
    is_a = np.repeat([True, False], n // 2)
    Z = rng.normal(size=(n, p)) + site[:, None] * 0.5          # site offsets must not create coupling
    shared = rng.normal(size=n)
    Z[is_a, 0] += couple * shared[is_a]
    Z[is_a, 1] += couple * shared[is_a]
    return Z, is_a, site


def test_center_and_corr():
    Z, is_a, site = _groups()
    Zc = cd.center(Z, site, is_a)
    for s in range(4):
        for a in (True, False):
            assert np.allclose(Zc[(site == s) & (is_a == a)].mean(axis=0), 0)
    assert np.allclose(cd.corr(Z), np.corrcoef(Z, rowvar=False))


def test_global_test_power_and_null():
    Z, is_a, site = _groups(couple=1.5)
    res = cd.global_test(Z, is_a, site, 199, 0)
    assert res["p_frobenius"] < 0.05 and res["mean_r_diff"] > 0
    Z0, a0, s0 = _groups(couple=0.0, seed=3)
    assert cd.global_test(Z0, a0, s0, 199, 0)["p_frobenius"] > 0.05


def test_edge_tests_find_pair():
    Z, is_a, site = _groups(couple=1.5)
    Zc = cd.center(Z, site, is_a)
    e = cd.edge_tests(Zc[is_a], Zc[~is_a], [f"r{i}" for i in range(Z.shape[1])])
    top = e.sort_values("p").iloc[0]
    assert {top.region_a, top.region_b} == {"r0", "r1"} and top.q < 0.05


def test_matrix_spin_corr_recovers_connectome():
    rng = np.random.default_rng(0)
    p = 20
    C = np.abs(rng.normal(size=(p, p)))
    C = (C + C.T) / 2
    R = C + rng.normal(0, 0.1, (p, p))
    R = (R + R.T) / 2
    spins = np.array([rng.permutation(p) for _ in range(200)]).T
    res = cd.matrix_spin_corr(R, C, spins)
    assert res["r"] > 0.8 and res["p_spin"] < 0.01


def test_nan_statistics_give_nan_p():
    Z, is_a, site = _groups(couple=1.5)
    Z[:, 3] = 1.0                                   # constant region
    res = cd.global_test(Z, is_a, site, 19, 0)
    assert np.isnan(res["frobenius"]) and np.isnan(res["p_frobenius"]) and np.isnan(res["p_mean_r_diff"])
    p = 10
    rng = np.random.default_rng(0)
    spins = np.array([rng.permutation(p) for _ in range(20)]).T
    out = cd.matrix_spin_corr(np.ones((p, p)), rng.normal(size=(p, p)), spins)
    assert np.isnan(out["r"]) and np.isnan(out["p_spin"])


def test_center_cells_matches_pandas_and_edge_df():
    Z, is_a, site = _groups()
    assert np.allclose(cd.center(Z, site, is_a), cd.center_cells(Z, cd.cells(site, is_a)))
    Zc = cd.center(Z, site, is_a)
    e1 = cd.edge_tests(Zc[is_a], Zc[~is_a], [f"r{i}" for i in range(Z.shape[1])])
    e4 = cd.edge_tests(Zc[is_a], Zc[~is_a], [f"r{i}" for i in range(Z.shape[1])], 4, 4)
    assert (np.abs(e4.z) <= np.abs(e1.z) + 1e-12).all()   # fewer df -> larger SE -> smaller |z|
