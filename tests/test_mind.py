"""MIND helpers: Schaefer network parsing from MIND region names."""

import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("mind16", Path(__file__).parents[1] / "scripts" / "16_mind.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_yeo_network_parsing():
    assert m.yeo_network("lh_7Networks_LH_Vis_1") == "Vis"
    assert m.yeo_network("rh_7Networks_RH_Default_pCunPCC_2") == "Default"
    assert m.yeo_network("lh_7Networks_LH_SalVentAttn_FrOperIns_3") == "SalVentAttn"
    with pytest.raises(ValueError):
        m.yeo_network("lh_bankssts")
