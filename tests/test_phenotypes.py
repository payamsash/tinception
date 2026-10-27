import numpy as np
import pandas as pd
import pytest

from tinception import bids, ids
from tinception.io import merge_on_id
from tinception.phenotypes import audiometry, harmonize

CFG = {"pta4_freqs": [500, 1000, 2000, 4000], "pta_hf_freqs": [4000, 8000],
       "pta_hf468_freqs": [4000, 6000, 8000], "pta_ehf_freqs": [10000, 12500, 14000, 16000],
       "hl_threshold_db": 25}


def test_bids_ids_are_site_prefixed_and_alphanumeric():
    assert ids.bids_id("uiuc", "sub-108") == "sub-uiuc108"
    assert ids.bids_id("whasc", "108") == "sub-whasc108"
    assert ids.bids_id("cogtail", "0dexlx_TI-HL") == "sub-cogtail0dexlxTIHL"
    assert ids.bids_id("london_1", "AG_e1") == "sub-london1AGe1"


def test_tinmeg_ids_are_zero_padded():
    assert ids.canonical_original_id("tinmeg", 539) == "0539"
    assert ids.canonical_original_id("tinspect", "MR3976_ACFK") == "MR3976_ACFK"


def test_pta_binaural_and_better_ear():
    df = pd.DataFrame({
        "thr_R_500": [10], "thr_R_1000": [10], "thr_R_2000": [20], "thr_R_4000": [40], "thr_R_8000": [60],
        "thr_L_500": [30], "thr_L_1000": [30], "thr_L_2000": [40], "thr_L_4000": [60], "thr_L_8000": [80],
    })
    out = audiometry.summarise(df, CFG)
    assert out.loc[0, "PTA4"] == pytest.approx((20 + 20 + 30 + 50) / 4)
    assert out.loc[0, "PTA4_better_ear"] == pytest.approx(20.0)
    assert out.loc[0, "PTA4_asymmetry"] == pytest.approx(20.0)
    assert out.loc[0, "PTA_HF"] == pytest.approx((50 + 70) / 2)
    assert np.isnan(out.loc[0, "PTA_HF468"])  # 6 kHz never measured -> missing, not partial
    assert bool(out.loc[0, "HL"]) is True


def test_reported_binaural_threshold_is_not_overwritten():
    df = pd.DataFrame({"thr_B_500": [5.0], "thr_R_500": [50.0], "thr_L_500": [50.0]})
    assert audiometry.fill_binaural(df).loc[0, "thr_B_500"] == 5.0


def test_distress_z_within_instrument_patients_only():
    df = pd.DataFrame({"group": ["TI", "TI", "TI", "TI", "CO"],
                       "THI": [10, 30, np.nan, np.nan, 50],
                       "TFI": [np.nan, np.nan, 20, 40, np.nan]})
    d = harmonize.distress_score(df)
    assert list(d.distress_instrument[:4]) == ["THI", "THI", "TFI", "TFI"]
    assert np.isnan(d.distress_z[4])
    assert d.distress_z[0] == pytest.approx(-d.distress_z[1])
    assert d.distress_z[2] == pytest.approx(d.distress_z[0])


def test_merge_on_id_is_strict():
    a = pd.DataFrame({"participant_id": ["s1", "s2"], "x": [1, 2]})
    b = pd.DataFrame({"participant_id": ["s1"], "y": [3]})
    with pytest.raises(KeyError):
        merge_on_id(a, b, require="left")
    with pytest.raises(ValueError):
        merge_on_id(a, pd.DataFrame({"participant_id": ["s1", "s1"], "y": [1, 2]}))
    with pytest.raises(ValueError):
        merge_on_id(a, pd.DataFrame({"participant_id": ["s1"], "x": [9]}))


@pytest.mark.parametrize("oid,site", [
    ("0838", "tinmeg"), ("2yt56w__TI-HA", "cogtail"), ("ICCAC_794kGN", "iccac"),
    ("K20_controls", "triple"), ("Triple_06_patients", "triple"), ("MP03906_FIL", "london_2"),
    ("indiTMS_1045i", "inditms"),
])
def test_drops_routing(oid, site):
    assert bids.drops_site(oid) == site


def test_antinomics_override_is_per_ear():
    # 'a': both ears in .mat; 'b': only L in .mat; 'c': no .mat -> legacy ears minus the known-bad 250/2000 Hz,
    # legacy 12500 relabelled 12000; other sites untouched
    df = pd.DataFrame({"site": ["antinomics"] * 3 + ["uiuc"], "original_id": ["a", "b", "c", "a"],
                       "thr_R_250": [1.0, 2.0, 2.0, 3.0], "thr_R_2000": [1.0, 2.0, 2.0, 3.0],
                       "thr_R_4000": [9.0, 7.0, 7.0, 3.0], "thr_R_12500": [5.0, 6.0, 6.0, np.nan],
                       "thr_L_4000": [9.0, 9.0, 9.0, 3.0]})
    aud = pd.DataFrame({"site": "antinomics", "original_id": ["a", "b"],
                        "thr_R_4000": [10.0, np.nan], "thr_R_12000": [11.0, np.nan], "thr_L_4000": [20.0, 21.0]})
    out = harmonize.apply_audiogram_overrides(df, aud)
    a, b, c, u = (out.iloc[i] for i in range(4))
    assert (a.thr_R_4000, a.thr_R_12000, a.thr_L_4000) == (10.0, 11.0, 20.0)
    assert np.isnan(a.thr_R_250) and np.isnan(a.thr_R_12500)
    assert (b.thr_R_4000, b.thr_R_12000, b.thr_L_4000) == (7.0, 6.0, 21.0)
    for r in (b, c):
        assert np.isnan(r.thr_R_250) and np.isnan(r.thr_R_2000) and np.isnan(r.thr_R_12500)
    assert (c.thr_R_4000, c.thr_R_12000, c.thr_L_4000) == (7.0, 6.0, 9.0)
    assert [r.audiogram_source for r in (a, b, c, u)] == [
        "antinomics_mat", "antinomics_mat+legacy", "site_master", "site_master"]
    assert u.thr_R_250 == 3.0
