from pathlib import Path

from tinception import cat12

T1CFG = {"n4_sites": ["london_2"], "n4_dir": "n4_london2", "acq": {"sub-london1AGe1": "echomean"}}


def test_t1_input_rules():
    raw, der = Path("/raw"), Path("/der")
    assert cat12.t1_input("sub-london2MP1FIL", "london_2", T1CFG, raw, der) == der / "n4_london2/sub-london2MP1FIL_T1w.nii.gz"
    assert cat12.t1_input("sub-london1AGe1", "london_1", T1CFG, raw, der) == \
        raw / "sub-london1AGe1/anat/sub-london1AGe1_acq-echomean_T1w.nii.gz"
    assert cat12.t1_input("sub-uiuc1", "uiuc", T1CFG, raw, der) == raw / "sub-uiuc1/anat/sub-uiuc1_T1w.nii.gz"
