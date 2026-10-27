"""Build the analysis master table: one row per subject with design, demographics, hearing,
clinical scores, TIV, QC metrics and analysis-inclusion flags.

Inputs : rawdata/participants.tsv, QC dir (qc_metrics.tsv, qc_exclusions.tsv), a features dir with
         aseg_global.tsv (FreeSurfer eTIV) and cat12_global.tsv (CAT12 TIV, tissue volumes, IQR); on the
         cluster that is <derivatives>/features, locally a copy of those two tables
Outputs: <root>/tinception_master.tsv (+ .xlsx with a data dictionary sheet)

Inclusion flags (recomputed every run; QC manual decisions flow in from qc_exclusions.tsv):
  core_complete      age, sex, group and PTA4 present
  qc_ok              qc_status != exclude  (flags pending visual review count as ok until decided)
  site_design        case-control (site recruited both groups) | tinnitus-only
  in_case_control    core_complete & qc_ok & case-control site  -> group comparisons, normative models
  in_patients        core_complete & qc_ok & group == TI        -> within-patient analyses (all sites)
  vbm_ok             qc_ok & CAT12 IQR not excluded (analysis.yaml cat12.iqr_max_mark) -> VBM / CAT12 ROIs
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tinception import config

GROUPS = {
    "id": ["participant_id", "site", "original_id", "legacy_id"],
    "scanner": ["scanner_vendor", "scanner_model", "field_T", "sequence"],
    "design": ["site_design", "group", "tinnitus", "group4", "group4_source", "group4_hf", "group4_recruited",
               "hearing_aid"],
    "demographics": ["age", "sex", "education_years"],
    "hearing": ["PTA4", "pta_source", "PTA_HF", "PTA_HF468", "PTA_EHF", "PTA4_better_ear",
                "PTA4_asymmetry", "HL", "HL_HF", "PTA_reported"],
    "tinnitus": ["THI", "TFI", "TQ_goebel", "distress_instrument", "distress_raw", "distress_z",
                 "duration_months", "duration_note", "laterality", "laterality_asym", "laterality_raw",
                 "loudness_vas100", "loudness_nrs10", "TPFQ"],
    "hyperacusis": ["hyperacusis_instrument", "hyperacusis_raw", "hyperacusis_z", "hyperacusis_khalfa", "IHS"],
    "affect": ["HADS_A", "HADS_D", "HADS", "BDI", "BAI", "MDI", "STAI", "PSQ_index", "PHQ15", "ISI", "affect_z"],
    "other": ["handedness", "education_level"],
    "imaging": ["has_T1w", "has_T2w", "fs_done", "eTIV", "cat12_TIV", "cat12_GM", "cat12_WM", "cat12_CSF",
                "cat12_WMH"],
    "qc": ["euler_lh", "euler_rh", "euler_total", "euler_z", "cjv", "cnr", "snr_total", "efc", "fber",
           "qi_2", "n_iqm_outliers", "qc_auto", "qc_status", "qc_reasons", "qc_strict_exclude", "fs_input_note",
           "cat12_IQR", "cat12_IQR_pct", "cat12_IQR_z", "cat12_iqr_flag", "cat12_iqr_exclude",
           "manual_decision", "manual_note",
           "legacy_drop", "legacy_drop_reason"],
    "flags": ["core_complete", "qc_ok", "in_case_control", "in_patients", "vbm_ok", "phenotype_status",
              "corrections"],
}


CAT12_COLS = ["cat12_TIV", "cat12_GM", "cat12_WM", "cat12_CSF", "cat12_WMH", "cat12_IQR", "cat12_IQR_pct",
              "cat12_IQR_z", "cat12_iqr_flag", "cat12_iqr_exclude"]


def read_features(features: Path) -> pd.DataFrame:
    """FreeSurfer eTIV (aseg_global.tsv) + CAT12 global/quality columns (cat12_global.tsv), one row per subject."""
    fs = pd.read_csv(features / "aseg_global.tsv", sep="\t", na_values="n/a", usecols=["participant_id", "eTIV"])
    fs["eTIV"] = fs["eTIV"].round(0)
    cat = pd.read_csv(features / "cat12_global.tsv", sep="\t", na_values="n/a", usecols=["participant_id", *CAT12_COLS])
    return fs.merge(cat, on="participant_id", how="outer", validate="one_to_one")


def legacy_drop_reason(df: pd.DataFrame) -> pd.Series:
    reasons = []
    for r in df.itertuples():
        if not r.legacy_drop:
            reasons.append("")
            continue
        miss = [v for v, x in (("age", r.age), ("sex", r.sex), ("group", r.group), ("PTA", r.PTA4))
                if pd.isna(x)]
        if r.phenotype_status == "no_phenotype":
            reasons.append("no phenotype row")
        elif miss:
            reasons.append("missing " + "/".join(miss))
        else:
            reasons.append("no missing core data (likely v1 ID mismatch, e.g. double underscore)")
    return pd.Series(reasons, index=df.index)


def build(participants: Path, qc_dir: Path, features: Path) -> pd.DataFrame:
    p = pd.read_csv(participants, sep="\t", na_values="n/a")
    p = p[p.has_T1w.astype(bool)].copy()
    sites = config.sites()
    p["scanner_vendor"] = p.site.map(lambda s: sites[s].get("vendor"))
    p["scanner_model"] = p.site.map(lambda s: sites[s].get("model"))
    p["field_T"] = p.site.map(lambda s: sites[s].get("field_T"))
    p["sequence"] = p.site.map(lambda s: sites[s].get("sequence"))
    both = p.groupby("site").group.agg(lambda g: {"CO", "TI"} <= set(g.dropna()))
    p["site_design"] = p.site.map(lambda s: "case-control" if both[s] else "tinnitus-only")

    q = pd.read_csv(qc_dir / "qc_metrics.tsv", sep="\t", na_values="n/a")
    qcols = [c for c in GROUPS["qc"] + ["fs_done"] if c in q and c not in (
        "legacy_drop", "legacy_drop_reason", "fs_input_note", "qc_status", "manual_decision", "manual_note")]
    p = p.merge(q[["participant_id", *qcols]], on="participant_id", how="left", validate="one_to_one")
    e = pd.read_csv(qc_dir / "qc_exclusions.tsv", sep="\t", keep_default_na=False)
    p = p.merge(e[["participant_id", "qc_status", "manual_decision", "manual_note"]],
                on="participant_id", how="left", validate="one_to_one")
    p = p.merge(read_features(features), on="participant_id", how="left", validate="one_to_one")

    p["fs_input_note"] = p.participant_id.map(config.analysis()["qc"].get("fs_input_notes", {}))
    p["legacy_drop_reason"] = legacy_drop_reason(p)
    p["core_complete"] = p[["age", "sex", "group", "PTA4"]].notna().all(axis=1)
    p["qc_ok"] = p.qc_status.fillna("") != "exclude"
    if "qc_strict_exclude" not in p:
        raise SystemExit("qc_metrics.tsv has no qc_strict_exclude: rerun `03_qc.py aggregate` with analysis.yaml qc.strict")
    p["qc_strict_exclude"] = p["qc_strict_exclude"].fillna(True).astype(bool) | ~p["qc_ok"]
    p["in_case_control"] = p.core_complete & p.qc_ok & (p.site_design == "case-control")
    p["in_patients"] = p.core_complete & p.qc_ok & (p.group == "TI")
    p["vbm_ok"] = p.qc_ok & ~p["cat12_iqr_exclude"].fillna(True).astype(bool)  # no CAT12 output -> not ok
    for col, ok in (("eTIV", "qc_ok"), ("cat12_TIV", "vbm_ok")):
        miss = p.participant_id[p[ok] & p[col].isna()]
        if len(miss):
            print(f"WARNING: {col} missing for {len(miss)} {ok} subjects (stale features copy?): {miss.tolist()[:10]}")

    # single-frequency thresholds right after the hearing summaries: thr_<ear>_<Hz>, ears B, R, L
    thr = sorted((c for c in p if c.startswith("thr_")),
                 key=lambda c: ("BRL".index(c.split("_")[1]), int(c.split("_")[2])))
    cols = []
    for g, cs in GROUPS.items():
        cols += [c for c in cs if c in p]
        if g == "hearing":
            cols += thr
    return p[cols].sort_values(["site", "participant_id"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--participants", type=Path, default=None)
    ap.add_argument("--qc-dir", type=Path, default=Path("qc/final"))
    ap.add_argument("--features", type=Path, default=Path("qc/features"),
                    help="dir with aseg_global.tsv + cat12_global.tsv (default: local copy in qc/features)")
    ap.add_argument("--out", type=Path, default=None, help="output directory (default: data root)")
    args = ap.parse_args()

    out = args.out or config.paths().root
    m = build(args.participants or config.paths().rawdata / "participants.tsv", args.qc_dir, args.features)
    m.to_csv(out / "tinception_master.tsv", sep="\t", index=False, na_rep="n/a")
    dictionary = pd.DataFrame([(g, c) for g, cs in GROUPS.items() for c in cs if c in m]
                              + [("audiogram", c) for c in m if c.startswith("thr_")],
                              columns=["group", "column"])
    with pd.ExcelWriter(out / "tinception_master.xlsx") as xw:
        m.to_excel(xw, sheet_name="master", index=False)
        dictionary.to_excel(xw, sheet_name="columns", index=False)

    print(f"{len(m)} subjects -> {out}/tinception_master.tsv/.xlsx")
    print(pd.crosstab(m.site, [m.site_design, m.in_case_control], margins=True).to_string())
    print("\ncore_complete:", int(m.core_complete.sum()), " qc_ok:", int(m.qc_ok.sum()),
          " in_case_control:", int(m.in_case_control.sum()), " in_patients:", int(m.in_patients.sum()),
          " vbm_ok:", int(m.vbm_ok.sum()))
    cc = m[m.in_case_control]
    print("\ncase-control set by group4:", cc.group4.value_counts().to_dict())


if __name__ == "__main__":
    main()
