from pathlib import Path

import numpy as np
import pytest

from tinception import cat12

T1CFG = {"n4_sites": ["london_2"], "n4_dir": "n4_london2", "acq": {"sub-london1AGe1": "echomean"}}


def test_t1_input_rules():
    raw, der = Path("/raw"), Path("/der")
    assert cat12.t1_input("sub-london2MP1FIL", "london_2", T1CFG, raw, der) == der / "n4_london2/sub-london2MP1FIL_T1w.nii.gz"
    assert cat12.t1_input("sub-london1AGe1", "london_1", T1CFG, raw, der) == \
        raw / "sub-london1AGe1/anat/sub-london1AGe1_acq-echomean_T1w.nii.gz"
    assert cat12.t1_input("sub-uiuc1", "uiuc", T1CFG, raw, der) == raw / "sub-uiuc1/anat/sub-uiuc1_T1w.nii.gz"


REPORT = """<S><qualityratings><IQR>1.84</IQR><NCR>1.6</NCR><res_RMS>2</res_RMS></qualityratings>
<subjectmeasures><vol_abs_CGW>[395.0 670.5 610.1 0 0]</vol_abs_CGW><vol_abs_WMH>2.7</vol_abs_WMH>
<vol_TIV>1674.8</vol_TIV></subjectmeasures></S>"""
ROI = """<S><suit><ids>[1;2]</ids><names><item>Left I IV</item><item>Vermis VI</item></names>
<data><Vgm>[2.9;6.6]</Vgm><Vwm>[0.8;0.1]</Vwm></data></suit></S>"""


def test_read_cat12_xml(tmp_path):
    d = tmp_path / "sub-x"
    (d / "report").mkdir(parents=True)
    (d / "label").mkdir()
    (d / "mri").mkdir()
    (d / "report" / "cat_sub-x_T1w.xml").write_text(REPORT)
    (d / "label" / "catROI_sub-x_T1w.xml").write_text(ROI)
    (d / "mri" / "mwp1sub-x_T1w.nii").write_text("")
    glob, rois = cat12.collect({"sub-x": d, "sub-y": tmp_path / "sub-y"})
    x = glob.set_index("participant_id").loc["sub-x"]
    assert (x.cat12_IQR, x.cat12_TIV, x.cat12_CSF, x.cat12_GM, x.cat12_WM, x.cat12_WMH) == (
        1.84, 1674.8, 395.0, 670.5, 610.1, 2.7)
    assert abs(x.cat12_IQR_pct - 86.6) < 1e-9 and x.cat12_mwp1
    assert not glob.set_index("participant_id").loc["sub-y", "cat12_mwp1"]
    assert rois["suit"].iloc[0].to_dict() == {"participant_id": "sub-x", "gm__Left_I_IV": 2.9,
                                              "gm__Vermis_VI": 6.6, "wm__Left_I_IV": 0.8, "wm__Vermis_VI": 0.1}


def test_cat12_partial_outputs(tmp_path):
    # failed subject: report with only an error block, no ROI xml -> NaN row, not a crash
    d = tmp_path / "sub-z"
    (d / "report").mkdir(parents=True)
    (d / "report" / "cat_sub-z_T1w.xml").write_text("<S><error>segmentation failed</error></S>")
    glob, rois = cat12.collect({"sub-z": d})
    z = glob.iloc[0]
    assert np.isnan(z.cat12_IQR) and np.isnan(z.cat12_TIV) and not z.cat12_roi and not rois


def test_cat12_roi_reads_only_tissues_and_rejects_name_clashes(tmp_path):
    f = tmp_path / "roi.xml"
    f.write_text("<S><a><names><item>R 1</item></names><data><Vgm>[1]</Vgm><Vcsf>[2]</Vcsf>"
                 "<Vother>[9]</Vother></data></a></S>")
    assert cat12.read_rois(f) == {"a": {"gm__R_1": 1.0, "csf__R_1": 2.0}}
    f.write_text("<S><a><names><item>R-1</item><item>R 1</item><item>R/1</item></names>"
                 "<data><Vgm>[1;2;3]</Vgm></data></a></S>")
    with pytest.raises(ValueError, match="not unique"):
        cat12.read_rois(f)
