"""Extract FreeSurfer feature tables for all subjects (run on the cluster after slurm/06_annot_stats).

Outputs in <derivatives>/features/:
  aseg.tsv, aseg_global.tsv, subseg.tsv, dk.tsv, dkt.tsv, destrieux.tsv, schaefer400.tsv, hcpmmp.tsv,
  auditory.tsv                 one row per subject (participant_id + features)
  vertex_fsavg5_fwhm10_{thickness,area,volume}.npy + vertex_subjects.txt   (n_subjects x 20484, lh then rh)
  completeness.tsv             per subject: number of features per table (spot missing outputs)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tinception import config, features


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fs-root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fwhm", type=int, default=10)
    args = ap.parse_args()

    fs = args.fs_root or config.paths().freesurfer
    out = args.out or config.paths().derivatives / "features"
    out.mkdir(parents=True, exist_ok=True)
    subjects = {d.name: d for d in sorted(fs.glob("sub-*")) if d.is_dir() and "_failed" not in d.name}
    print(f"{len(subjects)} subjects in {fs}")

    tables, verts, vids = features.collect(subjects, fwhm=args.fwhm)
    comp = pd.DataFrame({"participant_id": list(subjects)})
    for name, t in tables.items():
        t.to_csv(out / f"{name}.tsv", sep="\t", index=False, na_rep="n/a")
        comp = comp.merge(t.set_index("participant_id").notna().sum(axis=1).rename(name).reset_index(),
                          on="participant_id", how="left")
        print(f"{name:12s} {t.shape[0]:4d} subjects x {t.shape[1] - 1:5d} features")
    comp.to_csv(out / "completeness.tsv", sep="\t", index=False)
    for m, a in verts.items():
        np.save(out / f"vertex_fsavg5_fwhm{args.fwhm}_{m}.npy", a)
    (out / "vertex_subjects.txt").write_text("\n".join(vids) + "\n")
    print(f"vertex matrices: {len(vids)} subjects x 20484 vertices -> {out}")


if __name__ == "__main__":
    main()
