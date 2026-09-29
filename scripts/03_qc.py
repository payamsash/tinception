"""Quality control: collect Euler numbers + MRIQC IQMs and write qc_metrics / qc_exclusions.

Subcommands
  euler      collect FreeSurfer Euler numbers + run status from one FreeSurfer root
             (run on the cluster for subjects there, locally on the SSD for the rest)
  aggregate  merge euler_*.tsv + MRIQC jsons with participants.tsv, apply site-wise rules

Examples
  # cluster: 570 subjects in BIDS-named FreeSurfer dirs
  python scripts/03_qc.py euler --fs-root $FS --naming bids --out qc/euler_cluster.tsv
  # local SSD: the 304 complete legacy subjects (legacy dir names from participants.fs_dir)
  python scripts/03_qc.py euler --fs-root $SSD_FS --naming legacy --only-missing qc/euler_cluster.tsv --out qc/euler_ssd.tsv
  python scripts/03_qc.py aggregate --euler qc/euler_*.tsv --mriqc $MRIQC_DIR --out qc/

Outputs (aggregate)
  qc_metrics.tsv       all metrics + robust site-wise z-scores
  qc_exclusions.tsv    participant_id, site, qc_auto, qc_reasons, manual_decision/manual_note/reviewer,
                       qc_status (= manual decision if given, else automatic). Edit manual_* columns by
                       hand after visual review; re-running aggregate keeps them.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import pandas as pd

from tinception import config, qc


def cmd_euler(args: argparse.Namespace) -> None:
    part = pd.read_csv(args.participants or config.paths().rawdata / "participants.tsv",
                       sep="\t", na_values="n/a")
    part = part[part.has_T1w.astype(bool)]
    if args.naming == "bids":
        subjects = {p: p for p in part.participant_id}
    else:
        subjects = dict(zip(part.participant_id, part.fs_dir, strict=True))
        subjects = {p: d for p, d in subjects.items() if isinstance(d, str)}
    if args.only_missing:
        have = pd.read_csv(args.only_missing, sep="\t")
        done = set(have.loc[have.euler_total.notna(), "participant_id"])
        subjects = {p: d for p, d in subjects.items() if p not in done}
    root = Path(args.fs_root)
    subjects = {p: d for p, d in subjects.items() if (root / d).is_dir()}
    eul = qc.collect_euler(root, subjects)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    eul.to_csv(args.out, sep="\t", index=False, na_rep="n/a")
    print(f"{len(eul)} subjects, Euler found for {eul.euler_total.notna().sum()} -> {args.out}")


def cmd_aggregate(args: argparse.Namespace) -> None:
    cfg = config.analysis()["qc"]
    part = pd.read_csv(args.participants or config.paths().rawdata / "participants.tsv",
                       sep="\t", na_values="n/a")
    df = part.loc[part.has_T1w.astype(bool), ["participant_id", "site", "group", "legacy_drop"]]

    if args.euler:
        eul = pd.concat([pd.read_csv(f, sep="\t", na_values="n/a") for f in args.euler])
        eul = eul.sort_values("euler_total").drop_duplicates("participant_id", keep="last")
        df = df.merge(eul, on="participant_id", how="left", validate="one_to_one")
    if args.mriqc:
        iqm = qc.collect_mriqc(args.mriqc)
        df = df.merge(iqm, on="participant_id", how="left", validate="one_to_one")

    df = qc.apply_rules(df, cfg)
    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "qc_metrics.tsv", sep="\t", index=False, na_rep="n/a")

    excl_file = args.out / "qc_exclusions.tsv"
    prev = pd.read_csv(excl_file, sep="\t", keep_default_na=False) if excl_file.exists() else None
    cols = ["participant_id", "site", "group", "qc_auto", "qc_reasons", "euler_total", "n_iqm_outliers"]
    excl = qc.merge_manual(df[cols], prev)
    excl.to_csv(excl_file, sep="\t", index=False, na_rep="n/a")
    with open(args.out / "qc_provenance.json", "w") as f:
        json.dump({"date": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"), "rules": cfg,
                   "euler_files": [str(e) for e in args.euler or []],
                   "mriqc_dir": str(args.mriqc) if args.mriqc else None}, f, indent=2)

    print(pd.crosstab(excl.site, excl.qc_status, margins=True).to_string())
    missing = df.euler_total.isna().sum() if "euler_total" in df else len(df)
    print(f"\nEuler missing: {missing}; MRIQC missing: {df['cjv'].isna().sum() if 'cjv' in df else len(df)}")
    print(f"-> {args.out}/qc_metrics.tsv, qc_exclusions.tsv (edit manual_decision after visual review)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--participants", type=Path, default=None)
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("euler")
    e.add_argument("--fs-root", type=Path, required=True)
    e.add_argument("--naming", choices=["bids", "legacy"], default="bids")
    e.add_argument("--only-missing", type=Path, default=None,
                   help="skip subjects that already have Euler in this tsv")
    e.add_argument("--out", type=Path, required=True)
    e.set_defaults(func=cmd_euler)

    a = sub.add_parser("aggregate")
    a.add_argument("--euler", type=Path, nargs="*", default=[])
    a.add_argument("--mriqc", type=Path, default=None)
    a.add_argument("--out", type=Path, default=Path("qc"))
    a.set_defaults(func=cmd_aggregate)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
