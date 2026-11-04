"""Coordinated co-deviation (analysis plan 6.9 D; Hettwer et al. 2022): do patients' normative deviations co-vary
across regions differently from held-out controls', and does their co-deviation follow the normative connectome?

Per CO-reference model (analysis.yaml codeviation.models): F1 + F3 deviations, centred within site x group.
Outputs in <results>/12_codeviation/: global.tsv, edges_<model>.tsv, connectome.tsv, summary.json.

Example (cluster, repo root):  sbatch slurm/12_codeviation.sbatch
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tinception import config
from tinception.normative import core as ncore
from tinception.stats import codeviation as cd
from tinception.stats import group as g
from tinception.stats import network as net

ID = "participant_id"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default: <results>/12_codeviation")
    args = ap.parse_args()
    cfg, p = config.analysis(), config.paths()
    ccfg, ncfg, seed = cfg["codeviation"], cfg["normative"], cfg["seed"]
    q_thr = cfg["group"]["fdr_q"]
    out = args.out or p.results / "12_codeviation"
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(args.master or p.root / "tinception_master.tsv", sep="\t", na_values="n/a",
                         low_memory=False)[[ID, "site", "group"]]
    conn = net.connectomes()
    coords, hemi, names = g.parcel_centroids(p.freesurfer / "fsaverage")
    spins_for = net.spin_cache(coords, hemi, names, seed)
    glob, cmat, edge_hits = [], [], {}
    for model in ccfg["models"]:
        d = p.derivatives / "normative" / "CO" / model
        ok = pd.read_csv(d / "feature_calibration.tsv", sep="\t", index_col=0).ok
        feats = ncore.f1f3(ok, ncfg)
        z = pd.read_csv(d / "z_all.tsv", sep="\t", na_values="n/a").merge(master, on=ID)
        z = z[z.group.eq("TI") | z.reference.astype(bool)].dropna(subset=[*feats, "site"])
        is_pat, site = z.group.eq("TI").to_numpy(), z.site.to_numpy()
        Z = z[feats].to_numpy(float)
        res = cd.global_test(Z, is_pat, site, ccfg["n_perm"], seed)
        glob.append({"model": model, "n_patients": int(is_pat.sum()), "n_controls": int((~is_pat).sum()),
                     "n_regions": len(feats), **res})
        print(f"{model} global: {glob[-1]}", flush=True)
        Zc = cd.center(Z, site, is_pat)
        e = cd.edge_tests(Zc[is_pat], Zc[~is_pat], feats, int(pd.Series(site[is_pat]).nunique()),
                          int(pd.Series(site[~is_pat]).nunique()))
        e.to_csv(out / f"edges_{model}.tsv", sep="\t", index=False)
        hits = e[e.q < q_thr].sort_values("p")
        edge_hits[model] = {"n": len(hits), "top": hits.head(10).to_dict("records")}
        # connectome alignment of the cortical part (thickness)
        th = [f for f in feats if f.startswith("thick__")]
        cort = [f.replace("thick__", "") for f in th]
        ci = [feats.index(f) for f in th]
        Rp, Rc = cd.corr(Zc[is_pat][:, ci]), cd.corr(Zc[~is_pat][:, ci])
        diff = cd.fisher_z(Rp) - cd.fisher_z(Rc)
        sp = spins_for(cort, ccfg["spin_n"])
        for kind in ("fc", "sc"):
            C = conn[f"{kind}_ctx"].loc[cort, cort].to_numpy(float)
            for label, M in (("patients", Rp), ("controls", Rc), ("patients_minus_controls", diff)):
                cmat.append({"model": model, "connectivity": kind, "matrix": label, "n_regions": len(cort),
                             **cd.matrix_spin_corr(M, C, sp)})
        print(f"{model}: edges q<{q_thr}: {len(hits)}; connectome done", flush=True)
    glob, cmat = pd.DataFrame(glob), pd.DataFrame(cmat)
    glob.to_csv(out / "global.tsv", sep="\t", index=False)
    cmat.to_csv(out / "connectome.tsv", sep="\t", index=False)
    summary = {"q_threshold": q_thr, "global": glob.to_dict("records"), "connectome": cmat.to_dict("records"),
               "edges_significant": edge_hits}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(glob.round(4).to_string(index=False))
    print(cmat.round(3).to_string(index=False))
    print({m: v["n"] for m, v in edge_hits.items()})


if __name__ == "__main__":
    main()
