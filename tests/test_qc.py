import numpy as np
import pandas as pd

from tinception import qc

CFG = {"mad_k": 3.5, "min_iqm_outliers": 2, "euler_abs_floor": -434}


def test_euler_parsed_from_last_log_occurrence(tmp_path):
    s = tmp_path / "sub-x" / "scripts"
    s.mkdir(parents=True)
    (s / "recon-all.log").write_text(
        "orig.nofix lheno =  -80, rheno = -60\n...\norig.nofix lheno =  -40, rheno = -32\n")
    (s / "recon-all.done").write_text("")
    out = qc.euler_from_subject(tmp_path / "sub-x")
    assert (out["euler_lh"], out["euler_rh"], out["euler_total"]) == (-40, -32, -72)
    assert out["fs_done"] and not out["fs_error"]


def _cohort(n=20):
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "participant_id": [f"sub-{i}" for i in range(n)], "site": ["a"] * n,
        "euler_total": rng.normal(-60, 10, n).round(), "fs_error": False, "legacy_drop": False,
        "cjv": rng.normal(0.4, 0.02, n), "cnr": rng.normal(3, 0.1, n),
        "snr_total": rng.normal(12, 0.5, n), "efc": rng.normal(0.5, 0.01, n),
        "fber": rng.normal(2000, 50, n), "qi_2": rng.normal(0.005, 0.0005, n),
    })


def test_rules_exclude_flag_include():
    df = _cohort()
    df.loc[0, "euler_total"] = -600                 # below absolute floor -> exclude
    df.loc[1, "euler_total"] = -150                 # site outlier -> flag
    df.loc[2, ["cjv", "cnr"]] = [0.6, 2.0]           # two IQM outliers -> flag
    df.loc[3, "cjv"] = 0.6                           # one IQM outlier only -> include
    df.loc[4, "legacy_drop"] = True                  # v1 drop -> flag for review
    out = qc.apply_rules(df, CFG)
    assert list(out.qc_auto[:5]) == ["exclude", "flag", "flag", "include", "flag"]
    assert (out.qc_auto[5:] == "include").all()


def test_manual_decisions_survive_rerun():
    auto = pd.DataFrame({"participant_id": ["a", "b"], "qc_auto": ["flag", "include"]})
    prev = pd.DataFrame({"participant_id": ["a"], "manual_decision": ["include"],
                         "manual_note": ["ok on visual"], "reviewer": ["PS"]})
    out = qc.merge_manual(auto, prev)
    assert list(out.qc_status) == ["include", "include"]
    assert out.loc[0, "manual_note"] == "ok on visual"


def test_fber_sentinel_is_missing(tmp_path):
    import json
    d = tmp_path / "sub-a" / "anat"
    d.mkdir(parents=True)
    (d / "sub-a_T1w.json").write_text(json.dumps({"cjv": 0.5, "cnr": 2, "snr_total": 9, "efc": 0.5,
                                                  "fber": -1.0, "qi_2": 0.003}))
    out = qc.collect_mriqc(tmp_path)
    assert out.loc[0, "participant_id"] == "sub-a" and np.isnan(out.loc[0, "fber"])


def test_extra_runs_ignored(tmp_path):
    import json
    d = tmp_path / "sub-b" / "anat"
    d.mkdir(parents=True)
    for name in ("sub-b_T1w.json", "sub-b_run-2_T1w.json"):
        (d / name).write_text(json.dumps({"cjv": 0.5, "cnr": 2, "snr_total": 9, "efc": 0.5, "fber": 10, "qi_2": 0.003}))
    assert len(qc.collect_mriqc(tmp_path)) == 1


def test_strict_rule_column():
    df = _cohort()
    df.loc[0, "euler_total"] = -350          # < -300 strict floor, above -434 primary floor
    cfg = {**CFG, "strict": {"euler_abs_floor": -300, "euler_site_z": -4, "iqm_k": 4, "min_iqm_outliers": 2}}
    out = qc.apply_rules(df, cfg)
    assert out.loc[0, "qc_auto"] != "exclude" and bool(out.loc[0, "qc_strict_exclude"])
    assert not out.qc_strict_exclude[5:].any()
