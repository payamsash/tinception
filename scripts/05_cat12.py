"""CAT12.9 volume pipeline (secondary VBM + cerebellar ROIs).

Subcommands
  list     write the CAT12 input list (participant_id, t1): every subject with a T1w, using the same
           inputs as FreeSurfer v2 (analysis.yaml -> t1_inputs: N4 for london_2, echo-mean for two
           london_1 subjects, skipped subjects). QC exclusions are applied at analysis time, not here.
  collect  CAT12 xml outputs -> <derivatives>/features/: cat12_global.tsv (quality marks, TIV, tissue volumes,
           site-wise robust z of the IQR mark and the analysis.yaml cat12 QC rule) and cat12_<atlas>.tsv
           (GM/WM/CSF ROI volumes, cm^3).

Example (cluster, repo root)
  TINCEPTION_ENV=cluster uv run scripts/05_cat12.py list
  L=/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/cat12_list.tsv
  sbatch --array=1-$(($(wc -l < $L)-1))%100 slurm/05_cat12.sbatch $L
  TINCEPTION_ENV=cluster uv run scripts/05_cat12.py collect
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tinception import cat12, config, qc


def cmd_list(args: argparse.Namespace) -> None:
    p = config.paths()
    t1cfg = config.analysis()["t1_inputs"]
    part = pd.read_csv(p.rawdata / "participants.tsv", sep="\t", na_values="n/a")
    part = part[part.has_T1w.astype(bool) & ~part.participant_id.isin(t1cfg.get("skip", []))]
    out = pd.DataFrame({
        "participant_id": part.participant_id,
        "t1": [str(cat12.t1_input(pid, site, t1cfg, p.rawdata, p.derivatives))
               for pid, site in zip(part.participant_id, part.site, strict=True)],
    })
    missing = out[~out.t1.map(lambda s: Path(s).is_file())]
    if len(missing) and not args.allow_missing:
        raise SystemExit(f"{len(missing)} input files missing, e.g.\n{missing.head().to_string(index=False)}")
    dest = args.out or p.rawdata / "code" / "cat12_list.tsv"
    out.to_csv(dest, sep="\t", index=False)
    print(f"{len(out)} subjects ({len(missing)} missing inputs) -> {dest}")


def cmd_collect(args: argparse.Namespace) -> None:
    p = config.paths()
    cfg = config.analysis()
    root = args.root or p.derivatives / f"cat12-{cfg['cat12']['version']}"
    out = args.out or p.derivatives / "features"
    out.mkdir(parents=True, exist_ok=True)
    subjects = {d.name: d for d in sorted(root.glob("sub-*")) if d.is_dir()}
    if not subjects:
        raise SystemExit(f"no CAT12 subject directories in {root}")
    glob, rois = cat12.collect(subjects)
    part = pd.read_csv(p.rawdata / "participants.tsv", sep="\t", na_values="n/a", usecols=["participant_id", "site"])
    glob = part.merge(glob, on="participant_id", how="right", validate="one_to_one")
    if glob.site.isna().any():
        print(f"WARNING: not in participants.tsv (no site, no IQR z): {glob.participant_id[glob.site.isna()].tolist()}")
    # higher mark = worse image: flag beyond the site-wise robust z used for MRIQC (unknown where z cannot be
    # computed, e.g. MAD = 0 or no site); exclude worse than grade C or without a report
    glob["cat12_IQR_z"] = qc.robust_z(glob["cat12_IQR"], glob["site"])
    glob["cat12_iqr_flag"] = (glob["cat12_IQR_z"] > cfg["qc"]["mad_k"]).astype("boolean").where(glob["cat12_IQR_z"].notna())
    glob["cat12_iqr_exclude"] = (glob["cat12_IQR"] > cfg["cat12"]["iqr_max_mark"]) | glob["cat12_IQR"].isna()
    glob.to_csv(out / "cat12_global.tsv", sep="\t", index=False, na_rep="n/a")
    for atlas, t in rois.items():
        t.to_csv(out / f"cat12_{atlas}.tsv", sep="\t", index=False, na_rep="n/a")
        print(f"cat12_{atlas:20s} {t.shape[0]:4d} subjects x {t.shape[1] - 1:4d} features")
    print(f"{len(glob)} subjects in {root}; report for {glob.cat12_IQR.notna().sum()}, "
          f"mwp1 for {glob.cat12_mwp1.sum()}, ROI xml for {glob.cat12_roi.sum()}; IQR exclude "
          f"{glob.cat12_iqr_exclude.sum()}, flag {glob.cat12_iqr_flag.sum()} (unknown {glob.cat12_iqr_flag.isna().sum()}) -> {out}")
    for atlas, t in rois.items():
        if len(t) != glob.cat12_roi.sum():
            print(f"WARNING: cat12_{atlas} has {len(t)} subjects, {glob.cat12_roi.sum()} have an ROI xml")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("list", help="write the CAT12 input list")
    a.add_argument("--out", type=Path, default=None, help="default: rawdata/code/cat12_list.tsv")
    a.add_argument("--allow-missing", action="store_true", help="write the list even if inputs are missing")
    a.set_defaults(func=cmd_list)
    c = sub.add_parser("collect", help="CAT12 xml outputs -> feature tables")
    c.add_argument("--root", type=Path, default=None, help="default: <derivatives>/cat12-<version>")
    c.add_argument("--out", type=Path, default=None, help="default: <derivatives>/features")
    c.set_defaults(func=cmd_collect)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
