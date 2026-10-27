"""Build the analysis master table: one row per subject with design, demographics, hearing,
clinical scores, TIV, QC metrics and analysis-inclusion flags.

Inputs : rawdata/participants.tsv, QC dir (qc_metrics.tsv, qc_exclusions.tsv), eTIV tables
Outputs: <root>/tinception_master.tsv (+ .xlsx with a data dictionary sheet)

Inclusion flags (recomputed every run; QC manual decisions flow in from qc_exclusions.tsv):
  core_complete      age, sex, group and PTA4 present
  qc_ok              qc_status != exclude  (flags pending visual review count as ok until decided)
  site_design        case-control (site recruited both groups) | tinnitus-only
  in_case_control    core_complete & qc_ok & case-control site  -> group comparisons, normative models
  in_patients        core_complete & qc_ok & group == TI        -> within-patient analyses (all sites)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
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
    "imaging": ["has_T1w", "has_T2w", "fs_location", "fs_dir", "eTIV"],
    "qc": ["euler_lh", "euler_rh", "euler_total", "euler_z", "cjv", "cnr", "snr_total", "efc", "fber",
           "qi_2", "n_iqm_outliers", "qc_auto", "qc_status", "qc_reasons", "qc_strict_exclude", "fs_input_note",
           "manual_decision", "manual_note",
           "legacy_drop", "legacy_drop_reason"],
    "flags": ["core_complete", "qc_ok", "in_case_control", "in_patients", "phenotype_status", "corrections"],
}


def read_tiv(files: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(f, sep="\t", header=None, usecols=[0, 1], names=["participant_id", "eTIV"])
              for f in files if f.exists()]
    t = pd.concat(frames).dropna(subset=["eTIV"]).drop_duplicates("participant_id", keep="first")
    t["eTIV"] = pd.to_numeric(t["eTIV"], errors="coerce").round(0)
    return t


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


def build(participants: Path, qc_dir: Path, tiv_files: list[Path]) -> pd.DataFrame:
    p = pd.read_csv(participants, sep="\t", na_values="n/a")
    p = p[p.has_T1w.astype(bool)].copy()
    sites = config.sites()
    p["scanner_vendor"] = p.site.map(lambda s: sites[s].get("vendor"))
    p["scanner_model"] = p.site.map(lambda s: sites[s].get("model"))
    p["field_T"] = p.site.map(lambda s: sites[s].get("field_T"))
    p["sequence"] = p.site.map(lambda s: sites[s].get("sequence"))
    both = p.groupby("site").group.agg(lambda g: {"CO", "TI"} <= set(g.dropna()))
    p["site_design"] = p.site.map(lambda s: "case-control" if both[s] else "tinnitus-only")
    p["fs_location"] = np.where(p.fs_dir.notna() & ~p.participant_id.isin(_cluster_ids(tiv_files)),
                                "ssd_legacy_name", "cluster_bids_name")

    q = pd.read_csv(qc_dir / "qc_metrics.tsv", sep="\t", na_values="n/a")
    qcols = [c for c in GROUPS["qc"] if c in q and c not in ("legacy_drop", "legacy_drop_reason", "fs_input_note",
                                                                "qc_status", "manual_decision", "manual_note")]
    p = p.merge(q[["participant_id", *qcols]], on="participant_id", how="left", validate="one_to_one")
    e = pd.read_csv(qc_dir / "qc_exclusions.tsv", sep="\t", keep_default_na=False)
    p = p.merge(e[["participant_id", "qc_status", "manual_decision", "manual_note"]],
                on="participant_id", how="left", validate="one_to_one")
    p = p.merge(read_tiv(tiv_files), on="participant_id", how="left", validate="one_to_one")

    p["fs_input_note"] = p.participant_id.map(config.analysis()["qc"].get("fs_input_notes", {}))
    p["legacy_drop_reason"] = legacy_drop_reason(p)
    p["core_complete"] = p[["age", "sex", "group", "PTA4"]].notna().all(axis=1)
    p["qc_ok"] = p.qc_status.fillna("") != "exclude"
    if "qc_strict_exclude" not in p:
        raise SystemExit("qc_metrics.tsv has no qc_strict_exclude: rerun `03_qc.py aggregate` with analysis.yaml qc.strict")
    p["qc_strict_exclude"] = p["qc_strict_exclude"].fillna(True).astype(bool) | ~p["qc_ok"]
    p["in_case_control"] = p.core_complete & p.qc_ok & (p.site_design == "case-control")
    p["in_patients"] = p.core_complete & p.qc_ok & (p.group == "TI")

    # single-frequency thresholds right after the hearing summaries: thr_<ear>_<Hz>, ears B, R, L
    thr = sorted((c for c in p if c.startswith("thr_")),
                 key=lambda c: ("BRL".index(c.split("_")[1]), int(c.split("_")[2])))
    cols = []
    for g, cs in GROUPS.items():
        cols += [c for c in cs if c in p]
        if g == "hearing":
            cols += thr
    return p[cols].sort_values(["site", "participant_id"]).reset_index(drop=True)


def _cluster_ids(tiv_files: list[Path]) -> set[str]:
    f = [x for x in tiv_files if "cluster" in x.name and x.exists()]
    return set(pd.read_csv(f[0], sep="\t", header=None)[0]) if f else set()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--participants", type=Path, default=None)
    ap.add_argument("--qc-dir", type=Path, default=Path("qc/final"))
    ap.add_argument("--tiv", type=Path, nargs="*", default=[Path("qc/tiv_cluster.tsv"), Path("qc/tiv_ssd.tsv")])
    ap.add_argument("--out", type=Path, default=None, help="output directory (default: data root)")
    args = ap.parse_args()

    out = args.out or config.paths().root
    m = build(args.participants or config.paths().rawdata / "participants.tsv", args.qc_dir, args.tiv)
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
          " in_case_control:", int(m.in_case_control.sum()), " in_patients:", int(m.in_patients.sum()))
    cc = m[m.in_case_control]
    print("\ncase-control set by group4:", cc.group4.value_counts().to_dict())


if __name__ == "__main__":
    main()
