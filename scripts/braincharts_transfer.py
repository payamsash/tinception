"""Transfer of pretrained braincharts models (Rutherford et al. 2022) to new sites with adaptation data: runs in the
separate PCNtoolkit 0.35 environment (braincharts models are 0.x pickles), standalone (no tinception imports).
Called by slurm/07d_braincharts.sbatch after `07d_braincharts.py prepare`.

Inputs (--data): subjects.tsv (participant_id, age, sex 0/1, site, sitenum, fold, control), responses_<set>.tsv
(participant_id + one column per braincharts IDP name), sets.json ({set: [model dir, site-id file]}).
For every set, IDP and fold k: adaptation = controls of the other fold, test = everyone in fold k.
Output: z_<set>.tsv (participant_id x IDP), transfer_log.json.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd


@contextmanager
def chdir(path):
    """contextlib.chdir is Python >= 3.11; this env is 3.10."""
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def main() -> None:
    from pcntoolkit.normative import predict
    from pcntoolkit.util.utils import create_design_matrix

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--braincharts", type=Path, required=True, help="braincharts repo (models/, docs/)")
    ap.add_argument("--xmin", type=float, default=-5)
    ap.add_argument("--xmax", type=float, default=110)
    args = ap.parse_args()
    args.data, args.braincharts = args.data.resolve(), args.braincharts.resolve()
    subj = pd.read_csv(args.data / "subjects.tsv", sep="\t")
    sets = json.loads((args.data / "sets.json").read_text())
    log = {}
    for name, (model, site_file) in sets.items():
        site_ids_tr = (args.braincharts / "docs" / site_file).read_text().splitlines()
        resp = pd.read_csv(args.data / f"responses_{name}.tsv", sep="\t").set_index("participant_id")
        d = subj.set_index("participant_id").join(resp, how="inner")
        # design matrices once per fold (age / sex / site do not depend on the IDP); rows subset per IDP below
        design = {}
        for k in sorted(d.fold.unique()):
            for tag, x in (("te", d[d.fold == k]), ("ad", d[(d.fold != k) & d.control.astype(bool)])):
                X = create_design_matrix(x[["age", "sex"]], site_ids=x.site, all_sites=site_ids_tr, basis="bspline",
                                         xmin=args.xmin, xmax=args.xmax)
                design[(k, tag)] = pd.DataFrame(np.asarray(X), index=x.index)
        Z = pd.DataFrame(np.nan, index=d.index, columns=resp.columns)
        fails = {}
        for idp in resp.columns:
            mdir = args.braincharts / "models" / model / idp / "Models"
            if not mdir.exists():
                fails[idp] = "no model"
                continue
            for k in sorted(d.fold.unique()):
                te = d[(d.fold == k) & d[idp].notna()]
                ad = d[(d.fold != k) & d.control.astype(bool) & d[idp].notna()]
                with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as w, chdir(w):
                    files = {}
                    for tag, x in (("te", te), ("ad", ad)):
                        files[f"cov_{tag}"] = f"cov_{tag}.txt"
                        np.savetxt(files[f"cov_{tag}"], design[(k, tag)].loc[x.index].to_numpy())
                        files[f"resp_{tag}"] = f"resp_{tag}.txt"
                        np.savetxt(files[f"resp_{tag}"], x[idp].to_numpy(float))
                        files[f"site_{tag}"] = f"site_{tag}.txt"
                        np.savetxt(files[f"site_{tag}"], x.sitenum.to_numpy(int))
                    try:
                        out = predict(files["cov_te"], alg="blr", respfile=files["resp_te"], model_path=str(mdir),
                                      adaptrespfile=files["resp_ad"], adaptcovfile=files["cov_ad"],
                                      adaptvargroupfile=files["site_ad"], testvargroupfile=files["site_te"],
                                      inputsuffix="estimate")   # braincharts files: NM_0_0_estimate.pkl
                        Z.loc[te.index, idp] = np.asarray(out[2], float).ravel()
                    except Exception as err:  # noqa: BLE001 - one failing IDP must not stop the others; logged
                        fails[f"{idp}|fold{k}"] = repr(err)[:200]
            print(f"{name} {idp}: finite {int(Z[idp].notna().sum())}/{len(Z)}", flush=True)
        Z.reset_index().to_csv(args.data / f"z_{name}.tsv", sep="\t", index=False, na_rep="n/a")
        log[name] = {"model": model, "site_file": site_file, "n": len(Z), "idps": len(resp.columns), "failures": fails}
    (args.data / "transfer_log.json").write_text(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
