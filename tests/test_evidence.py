"""Bayes factor from estimate + SE and minimum detectable effect."""

import numpy as np
from scipy import integrate
from scipy import stats as st

from tinception.stats import evidence as ev


def test_normal_prior_closed_form_matches_integral():
    beta, se, tau = 0.05, 0.07, 0.2
    num = integrate.quad(lambda d: st.norm.pdf(beta, d, se) * st.norm.pdf(d, 0, tau), -2, 2)[0]
    assert np.isclose(ev.bf10(beta, se, "normal", tau), num / st.norm.pdf(beta, 0, se), rtol=1e-6)


def test_bf_direction_and_mde():
    assert ev.bf10(0.0, 0.07) < 1 / 3          # precise null estimate favours H0
    assert ev.bf10(0.5, 0.07) > 100            # large effect favours H1
    assert np.isclose(ev.mde(0.1), 0.2802, atol=1e-3)
    assert np.isnan(ev.bf10(np.nan, 0.1))
    assert ev.category(12) == "strong H0" and ev.category(0.2) == "moderate H1"


def test_cauchy_matches_pingouin_jzs_and_edge_cases():
    import pingouin as pg
    import pytest

    n1, n2, d = 394, 352, 0.12
    se = np.sqrt(1 / n1 + 1 / n2)
    assert np.isclose(ev.bf10(d, se), float(pg.bayesfactor_ttest(d / se, n1, n2, r=0.707)), rtol=2e-3)
    assert np.isfinite(ev.log_bf10(5.0, 0.1)) and ev.log_bf10(5.0, 0.1) > 100       # no underflow
    with pytest.raises(ValueError):
        ev.bf10(np.nan, 0.1, prior="Normal")
    assert ev.category(np.inf) == "extreme H0" and ev.category(0.0) == "extreme H1"
    assert ev.mde(0.1, 0.05 / 16) > ev.mde(0.1)
