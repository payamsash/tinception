from tinception import features

ASEG = """# Measure EstimatedTotalIntraCranialVol, eTIV, Estimated Total Intracranial Volume, 1500000.0, mm^3
# Measure BrainSegNotVent, BrainSegVolNotVent, Brain Segmentation Volume Without Ventricles, 1000000.0, mm^3
# ColHeaders  Index SegId NVoxels Volume_mm3 StructName normMean normStdDev normMin normMax normRange
  1   4     6621     6620.9  Left-Lateral-Ventricle   30.0 10.0 5 90 85
"""
APARC = """# Measure Cortex, MeanThickness, Mean Thickness, 2.5, mm
# ColHeaders StructName NumVert SurfArea GrayVol ThickAvg ThickStd MeanCurv GausCurv FoldInd CurvInd
transversetemporal 900 300 1000 2.6 0.5 0.1 0.0 1 0.1
"""


def test_parsers(tmp_path):
    (tmp_path / "stats").mkdir()
    (tmp_path / "mri").mkdir()
    (tmp_path / "stats" / "aseg.stats").write_text(ASEG)
    (tmp_path / "stats" / "lh.aparc.stats").write_text(APARC)
    (tmp_path / "mri" / "ThalamicNuclei.v13.T1.volumes.txt").write_text("Left-MGN 105.3\nRight-MGN 113.9\n")
    t = features.subject_features(tmp_path)
    assert t["aseg_global"]["eTIV"] == 1500000.0
    assert t["aseg"]["Left-Lateral-Ventricle"] == 6620.9
    assert t["dk"]["thick__lh_transversetemporal"] == 2.6 and t["dk"]["area__lh_transversetemporal"] == 300
    assert t["dk"]["thick__lh_MeanThickness"] == 2.5
    assert t["auditory"]["vol__lh_MGN"] == 105.3
