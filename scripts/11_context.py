"""Gradient contextualization (analysis plan 6.9 C; BrainSpace / Park et al. 2022): spatial correlation of the group
effect maps of 6.9 B with cortical hierarchy axes at DK level, with spin nulls.

Axes (data/context/dk_gradients.tsv, built by `prepare`):
  fc_gradient1_margulies2016  principal functional gradient (Margulies 2016; neuromaps, fsLR 32k -> fsaverage)
  sa_axis_sydnor2021          sensorimotor-association axis (Sydnor 2021; neuromaps)
  bigbrain_gradient_paquola   BigBrain microstructural profile covariance gradient (ENIGMA Toolbox, DK)
  fc_gradient1_hcp_dk         principal gradient of the ENIGMA HCP functional connectome at DK level (BrainSpace;
                              sign aligned to Margulies) - robustness of the first axis, outside the BH family
BH over the three primary axes per map.
Outputs: <results>/11_context/{gradients.tsv, summary.json}.

  prepare (local, needs internet):  python scripts/11_context.py prepare --fsaverage <dir with label/?h.aparc.annot>
  run (cluster):                    sbatch slurm/11_context.sbatch
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import group as g
from tinception.stats import network as net

REPO = Path(__file__).resolve().parents[1]
AXES_FILE = REPO / "data" / "context" / "dk_gradients.tsv"
ROBUSTNESS_AXES = ("fc_gradient1_hcp_dk",)
NEUROMAPS = {"fc_gradient1_margulies2016": ("margulies2016", "fcgradient01"),
             "sa_axis_sydnor2021": ("sydnor2021", "SAaxis")}


def _parcellate(lh: np.ndarray, rh: np.ndarray, fsaverage: Path) -> pd.Series:
    """Mean of a fsaverage (164k) vertex map within each DK parcel; exact zeros (medial-wall fill) and NaN ignored."""
    import nibabel as nib

    out = {}
    for h, data in (("lh", lh), ("rh", rh)):
        labels, _, names = nib.freesurfer.read_annot(str(fsaverage / "label" / f"{h}.aparc.annot"))
        for i, n in enumerate(names):
            n = n.decode() if isinstance(n, bytes) else n
            if n in ("unknown", "corpuscallosum"):
                continue
            v = data[labels == i]
            v = v[np.isfinite(v) & (v != 0)]
            if len(v):
                out[f"{h}_{n}"] = float(v.mean())
    return pd.Series(out)


def cmd_prepare(args: argparse.Namespace) -> None:
    import nibabel as nib
    from brainspace.gradient import GradientMaps
    from neuromaps.datasets import fetch_annotation, fetch_atlas
    from scipy.spatial import cKDTree

    conn = net.connectomes()
    names = list(conn["fc_ctx"].index)
    # fsLR 32k -> fsaverage 164k by nearest neighbour on neuromaps' registration spheres (fsLR sphere deformed to
    # fsaverage, as used by neuromaps.transforms.fslr_to_fsaverage; no Connectome Workbench needed); parcel means
    # barely depend on the interpolation. Check: the fsLR medial wall must land on the fsaverage medial wall.
    fslr, fsa = fetch_atlas("fsLR", "32k"), fetch_atlas("fsaverage", "164k")
    reg = [Path(fslr["sphere"][h]).parent / f"tpl-fsLR_space-fsaverage_den-32k_hemi-{H}_sphere.surf.gii"
           for h, H in ((0, "L"), (1, "R"))]
    nn = [cKDTree(nib.load(reg[h]).agg_data()[0]).query(nib.load(fsa["sphere"][h]).agg_data()[0])[1] for h in (0, 1)]
    for h in (0, 1):
        m1 = np.asarray(nib.load(fslr["medial"][h]).agg_data())[nn[h]] == 0
        m2 = np.asarray(nib.load(fsa["medial"][h]).agg_data()) == 0
        dice = 2 * (m1 & m2).sum() / (m1.sum() + m2.sum())
        if dice < 0.8:
            raise SystemExit(f"fsLR -> fsaverage misaligned (medial-wall Dice {dice:.2f}, hemisphere {h})")
    axes, prov = {}, {}
    for key, (src, desc) in NEUROMAPS.items():
        f = fetch_annotation(source=src, desc=desc, space="fsLR", den="32k")
        lh, rh = (np.asarray(nib.load(f[h]).agg_data(), float)[nn[h]] for h in (0, 1))
        axes[key] = _parcellate(lh, rh, args.fsaverage).reindex(names)
        prov[key] = {"neuromaps": [src, desc, "fsLR", "32k"],
                     "transform": "nearest neighbour, neuromaps tpl-fsLR_space-fsaverage_den-32k sphere -> "
                                  "fsaverage 164k sphere (medial-wall Dice >= 0.8 checked)",
                     "parcellation": "fsaverage aparc (DK), vertex mean, zeros/NaN ignored"}
    import enigmatoolbox

    bb = np.loadtxt(Path(enigmatoolbox.__file__).parent / "histology" / "bb_gradient_aparc.csv", delimiter=",")
    axes["bigbrain_gradient_paquola"] = pd.Series(bb, index=names)   # ENIGMA DK order = connectome order
    prov["bigbrain_gradient_paquola"] = {"enigmatoolbox": "histology/bb_gradient_aparc.csv (Paquola et al. 2019)"}
    gm = GradientMaps(n_components=2, approach="dm", kernel="normalized_angle", random_state=0)
    g1 = pd.Series(gm.fit(conn["fc_ctx"].to_numpy()).gradients_[:, 0], index=names)
    if np.corrcoef(g1, axes["fc_gradient1_margulies2016"])[0, 1] < 0:
        g1 = -g1
    axes["fc_gradient1_hcp_dk"] = g1
    prov["fc_gradient1_hcp_dk"] = {"brainspace": "GradientMaps dm normalized_angle, ENIGMA HCP fc_ctx (DK)",
                                   "sign": "aligned to Margulies 2016 G1"}
    df = pd.DataFrame(axes)
    if df.isna().any().any():
        raise SystemExit(f"missing parcels: {df.index[df.isna().any(axis=1)].tolist()}")
    AXES_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(AXES_FILE, sep="\t", index_label="region")
    AXES_FILE.with_suffix(".json").write_text(json.dumps(prov, indent=2))
    print(df.corr(method="spearman").round(2).to_string())
    print(f"-> {AXES_FILE}")


def cmd_run(args: argparse.Namespace) -> None:
    from statsmodels.stats.multitest import multipletests

    cfg, p = config.analysis(), config.paths()
    ncfg, seed = cfg["network"], cfg["seed"]
    out = args.out or p.results / "11_context"
    out.mkdir(parents=True, exist_ok=True)
    if not AXES_FILE.exists():
        raise SystemExit(f"{AXES_FILE} missing: run `11_context.py prepare` locally and sync data/ to the cluster")
    axes = pd.read_csv(AXES_FILE, sep="\t", index_col="region")
    coords, hemi, names = g.parcel_centroids(p.freesurfer / "fsaverage")
    spins_for = net.spin_cache(coords, hemi, names, seed)
    rows = []
    for name, (ctx, _) in net.effect_maps(p.results, cfg["h3"], ncfg).items():
        x = ctx.reindex(names).dropna()
        spins = spins_for(x.index, ncfg["spin_n"])
        for ax in axes:
            r = g.spin_correlation_spins(x.to_numpy(float), axes.loc[x.index, ax].to_numpy(float), spins)
            rows.append({"map": name, "axis": ax, "n_regions": len(x), **r})
        print(f"{name}: done", flush=True)
    res = pd.DataFrame(rows)
    res["primary_axis"] = ~res.axis.isin(ROBUSTNESS_AXES)
    prim = res[res.primary_axis]
    res["q"] = np.nan   # BH over the primary axes per map; robustness axes keep only p
    res.loc[prim.index, "q"] = prim.groupby("map").p_spin.transform(lambda v: multipletests(v, method="fdr_bh")[1])
    res.to_csv(out / "gradients.tsv", sep="\t", index=False)
    q_thr = cfg["group"]["fdr_q"]
    summary = {"q_threshold": q_thr, "hits": res[res.q < q_thr][["map", "axis", "r", "ci_lo", "ci_hi", "p_spin", "q"]]
               .to_dict("records"), "all": res[["map", "axis", "r", "p_spin", "q"]].to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(res.pivot(index="map", columns="axis", values="r").round(2).to_string())
    print(json.dumps(g.json_clean(summary["hits"]), indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("prepare")
    pr.add_argument("--fsaverage", type=Path, required=True, help="fsaverage dir with label/?h.aparc.annot")
    pr.set_defaults(func=cmd_prepare)
    ru = sub.add_parser("run")
    ru.add_argument("--out", type=Path, default=None, help="default <results>/11_context")
    ru.set_defaults(func=cmd_run)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
