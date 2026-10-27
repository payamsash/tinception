"""CAT12.9 volume pipeline (secondary VBM + cerebellar ROIs).

Subcommands
  list   write the CAT12 input list (participant_id, t1): every subject with a T1w, using the same
         inputs as FreeSurfer v2 (analysis.yaml -> t1_inputs: N4 for london_2, echo-mean for two
         london_1 subjects, skipped subjects). QC exclusions are applied at analysis time, not here.

Example (cluster, repo root)
  TINCEPTION_ENV=cluster uv run scripts/05_cat12.py list
  L=/shares/tinnitus.orl.med.uzh/tinception/rawdata/code/cat12_list.tsv
  sbatch --array=1-$(($(wc -l < $L)-1))%100 slurm/05_cat12.sbatch $L
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tinception import cat12, config


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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("list", help="write the CAT12 input list")
    a.add_argument("--out", type=Path, default=None, help="default: rawdata/code/cat12_list.tsv")
    a.add_argument("--allow-missing", action="store_true", help="write the list even if inputs are missing")
    a.set_defaults(func=cmd_list)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
