"""Molecular, network and null-model context (analysis plan 6.13; fixed 2026-11-04 before running).

Runs in a separate NiSpace environment (NiSpace 0.1.0b2, commit 3b59468; brainsmash 0.11.0) on small tables exported
from the cluster (H1 coefficients, within-patient coefficients, Model A deviations, master columns, HCP degree):
NiSpace provides the DK-parcellated reference maps (PET, cortex features, ENIGMA thickness) and the null maps
(spins / Moran / brainsmash); correlations and permutations are computed here with the repo's helpers.

  PYTHONPATH=. python scripts/18_molecular_context.py --inputs <dir> --out <dir> [--n-null 5000] [--n-perm 5000]

Outputs: colocalization.tsv (target x reference x null), null_robustness.tsv, individual.tsv, pet_maps.tsv,
summary.json.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from tinception.stats.network import _zrank
from tinception.stats.resample import bh

NULLS = ("alexander_bloch", "moran", "burt2020")


def zrank(a: np.ndarray, axis: int = -1) -> np.ndarray:
    """Standardized ranks along `axis` (network._zrank with a default axis for 1-D maps)."""
    return _zrank(a, axis=axis)
GRADIENT_PRIMARY = ("fc_gradient1_margulies2016", "sa_axis_sydnor2021", "bigbrain_gradient_paquola")   # as 11_context
CORTICAL_TRACER = {"D23": "flb457"}   # Hansen 2022: FLB457 for cortical D2/3 (raclopride / fallypride are striatal)


def to_nsp(name: str) -> str:
    """thick__lh_bankssts / lh_bankssts -> hemi-L_lab-bankssts (NiSpace DesikanKilliany labels); fails loudly."""
    m = re.match(r"^(?:thick__)?(lh|rh)_(.+)$", name)
    if not m:
        raise ValueError(f"not a DK cortical name: {name}")
    return f"hemi-{'L' if m.group(1) == 'lh' else 'R'}_lab-{m.group(2)}"


def clean(o):
    """NaN / inf -> None and numpy scalars -> Python, for strict JSON."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


def null_maps(x: pd.Series, parc, n: int, seed: int) -> dict[str, np.ndarray]:
    """method -> (n x parcels) null maps of one target map (NiSpace generators; spins / Moran / brainsmash)."""
    from nispace.nulls import generate_null_maps

    X = pd.DataFrame([x.to_numpy(float)], columns=x.index, index=["m"])
    out = {}
    for m in NULLS:
        nm = generate_null_maps(m, X, parc, n_nulls=n, seed=seed, verbose=False)
        out[m] = np.asarray((nm[0] if isinstance(nm, tuple) else nm).to_array(), float).reshape(-1, len(x))
    return out


def corr_null(x: np.ndarray, Y: np.ndarray, nulls: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Spearman r of map x with every row of Y; two-sided p from the null maps of x. Rows of Y with missing parcels
    use their finite parcels (nulls rotated on all parcels, then masked); null maps with non-finite values inside
    the mask are dropped instead of being counted as non-exceedances. Returns r, p, number of valid nulls.
    The ranked full-mask nulls are computed once and reused for every complete row."""
    r, p, nv = np.full(len(Y), np.nan), np.full(len(Y), np.nan), np.zeros(len(Y), int)
    base = np.isfinite(x)
    ok_b = np.isfinite(nulls[:, base]).all(axis=1)
    zn_b = zrank(nulls[ok_b][:, base], axis=1)
    for i, y in enumerate(Y):
        m = np.isfinite(y) & base
        if (m == base).all():
            zn = zn_b
        else:
            okm = np.isfinite(nulls[:, m]).all(axis=1)
            zn = zrank(nulls[okm][:, m], axis=1)
        zy = zrank(y[m])
        r[i] = zy @ zrank(x[m]) / m.sum()
        rn = zn @ zy / m.sum()
        nv[i] = len(rn)
        if np.isfinite(r[i]) and len(rn):
            p[i] = (np.sum(np.abs(rn) >= abs(r[i])) + 1) / (len(rn) + 1)
    return r, p, nv


def pet_maps(X: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One map per PET target (Hansen 2022 procedure): same-tracer maps z-scored and averaged weighted by their
    sample size; per target the tracer with the largest total n, except FLB457 for D2/3 (cortical binding)."""
    meta = pd.DataFrame({"name": X.index,
                         "target": [re.search(r"target-([^_]+)", s).group(1) for s in X.index],
                         "tracer": [re.search(r"tracer-([^_]+)", s).group(1) for s in X.index],
                         "n": [int(re.search(r"_n-(\d+)", s).group(1)) for s in X.index]})
    rows, info = {}, []
    for target, g_ in meta.groupby("target"):
        tracer = CORTICAL_TRACER.get(target, g_.groupby("tracer").n.sum().idxmax())
        sel = g_[g_.tracer == tracer]
        V = X.loc[sel.name]
        Zs = V.sub(V.mean(axis=1), axis=0).div(V.std(axis=1), axis=0)
        w = sel.n.to_numpy(float)
        avg = Zs.mul(w, axis=0).sum(axis=0, min_count=1) / w.sum()
        rows[f"{target}_{tracer}"] = avg.where(Zs.notna().all())   # a parcel missing in any source stays missing
        info.append({"map": f"{target}_{tracer}", "target": target, "tracer": tracer, "n_maps": len(sel),
                     "n_subjects": int(w.sum()), "sources": ";".join(sel.name)})
    return pd.DataFrame(rows).T, pd.DataFrame(info)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-null", type=int, default=5000)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    from nispace.datasets import fetch_parcellation, fetch_reference

    args.out.mkdir(parents=True, exist_ok=True)
    parc = fetch_parcellation("DesikanKilliany", verbose=False)
    labels = [str(lab) for lab in parc.get_labels()]
    kw = {"parcellation": "DesikanKilliany", "print_references": False, "verbose": False}
    pet, pet_info = pet_maps(fetch_reference("pet", **kw)[labels])
    pet_info.to_csv(args.out / "pet_maps.tsv", sep="\t", index=False)
    refs = {"pet": pet, "cortexfeatures": fetch_reference("cortexfeatures", **kw)[labels],
            "enigmathick": fetch_reference("enigmathick", **kw)[labels]}

    # ---------------- target maps (DK-68 thickness coefficients) ----------------
    m = pd.read_csv(args.inputs / "h1_models.tsv", sep="\t", na_values="n/a")
    main = m[(m.model == "main") & (m.family == "dk_thick")]
    targets, zmaps = {}, {}
    for con in ("tinnitus", "PTA4", "PTA_HF_beyond_PTA4"):
        x = main[main.contrast == con].set_index("outcome")
        targets[con] = x.beta.rename(to_nsp)[labels]
        zmaps[con] = x.z.rename(to_nsp)[labels]
    ctrl = m[(m.model == "controls_PTA_HF") & m.outcome.str.startswith("thick__")].set_index("outcome")
    zmaps["PTA_HF_controls"] = ctrl.z.rename(to_nsp)[labels]
    a = pd.read_csv(args.inputs / "patients_assoc.tsv", sep="\t", na_values="n/a")
    hy = a[(a.predictor == "hyperacusis_z") & (a.hearing == "PTA4") & (a.family == "dk_thick")].set_index("outcome")
    targets["hyperacusis"] = hy.beta.rename(to_nsp)[labels]
    nulls = {k: null_maps(v, parc, args.n_null, args.seed) for k, v in targets.items()}
    znull = null_maps(zmaps["tinnitus"], parc, args.n_null, args.seed)   # H1c uses the z maps

    # ---------------- group-level colocalization ----------------
    rows = []
    for tname, x in targets.items():
        for rname, R in refs.items():
            for meth in NULLS:
                r, p, nv = corr_null(x.to_numpy(float), R.to_numpy(float), nulls[tname][meth])
                q = bh(p)
                rows += [{"target": tname, "reference_set": rname, "reference": ref, "null": meth, "r": r[i],
                          "p": p[i], "q": q[i], "n_nulls_valid": int(nv[i])} for i, ref in enumerate(R.index)]
    col = pd.DataFrame(rows)
    wide = col.pivot_table(index=["target", "reference_set", "reference"], columns="null", values="q")
    rr = col.groupby(["target", "reference_set", "reference"]).r.first()
    consensus = wide[(wide[list(NULLS)] < 0.05).all(axis=1)].join(rr).reset_index()
    col.to_csv(args.out / "colocalization.tsv", sep="\t", index=False)

    # ---------------- null-model robustness of earlier spatial results ----------------
    deg = pd.read_csv(args.inputs / "degree.tsv", sep="\t", index_col=0).rename(index=to_nsp).loc[labels]
    grad = pd.read_csv(REPO / "data" / "context" / "dk_gradients.tsv", sep="\t", index_col="region")
    grad = grad.rename(index=to_nsp).loc[labels, list(GRADIENT_PRIMARY)]
    rob = []
    for meth in NULLS:   # H1c: tinnitus z map vs PTA4 z map and vs the controls-only PTA_HF z map
        Y = np.vstack([zmaps["PTA4"].to_numpy(float), zmaps["PTA_HF_controls"].to_numpy(float)])
        r, p, _ = corr_null(zmaps["tinnitus"].to_numpy(float), Y, znull[meth])
        rob += [{"family": "H1c", "target": "tinnitus (z)", "test": t, "null": meth, "r": r[i], "p": p[i]}
                for i, t in enumerate(("vs PTA4 (z)", "vs PTA_HF controls (z)"))]
    for tname in ("tinnitus", "PTA_HF_beyond_PTA4", "PTA4"):
        for meth in NULLS:
            for fam, M in (("hub", deg), ("gradient", grad)):
                r, p, _ = corr_null(targets[tname].to_numpy(float), M.to_numpy(float).T, nulls[tname][meth])
                rob += [{"family": fam, "target": tname, "test": c, "null": meth, "r": r[i], "p": p[i]}
                        for i, c in enumerate(M.columns)]
    rob = pd.DataFrame(rob)
    rob["q"] = rob.groupby(["family", "target", "null"]).p.transform(bh)   # BH as in the original analyses
    rob.to_csv(args.out / "null_robustness.tsv", sep="\t", index=False)

    # ---------------- individual-level PET colocalization (JuSpace-style, Freedman-Lane) ----------------
    z = pd.read_csv(args.inputs / "zA.tsv", sep="\t", na_values="n/a")
    ms = pd.read_csv(args.inputs / "master.tsv", sep="\t", na_values="n/a")
    th = [c for c in z.columns if c.startswith("thick__")]
    d = z.merge(ms, on="participant_id")
    d = d[d.in_case_control.astype(bool) & (d.group.eq("TI") | d.reference.astype(bool))]
    d = d.dropna(subset=[*th, "PTA4", "PTA_HF", "euler_total", "site"]).reset_index(drop=True)
    A = np.column_stack([np.ones(len(d)), d.PTA4])
    hf_res = d.PTA_HF.to_numpy() - A @ np.linalg.lstsq(A, d.PTA_HF.to_numpy(), rcond=None)[0]
    C = np.column_stack([d.PTA4, hf_res, d.euler_total]).astype(float)
    Z = d[th].rename(columns=to_nsp)[labels].to_numpy(float)
    g = d.group.eq("TI").to_numpy()
    site = d.site.to_numpy()
    Xfull = np.column_stack([np.ones(len(d)), g.astype(float), C])
    Xred = np.column_stack([np.ones(len(d)), C])
    P = pet.to_numpy(float)
    masks = [np.isfinite(v) for v in P]
    PETz = [zrank(v[mk]) for v, mk in zip(P, masks, strict=True)]
    full = [i for i, mk in enumerate(masks) if mk.all()]
    part = [i for i, mk in enumerate(masks) if not mk.all()]
    PETfull = np.vstack([PETz[i] for i in full])

    def statistic(Y: np.ndarray) -> np.ndarray:
        """Covariate effects removed with the full model (group kept); every patient z-scored against the controls;
        per-patient Spearman r with every PET map -> patients x PET."""
        b = np.linalg.lstsq(Xfull, Y, rcond=None)[0]
        Ya = Y - C @ b[2:]
        ref = Ya[~g]
        zp = (Ya[g] - ref.mean(axis=0)) / ref.std(axis=0, ddof=1)
        out = np.empty((int(g.sum()), len(masks)))
        out[:, full] = zrank(zp, axis=1) @ PETfull.T / Z.shape[1]
        for i in part:
            out[:, i] = zrank(zp[:, masks[i]], axis=1) @ PETz[i] / masks[i].sum()
        return out

    obs = statistic(Z)
    mean_r = obs.mean(axis=0)
    gam = np.linalg.lstsq(Xred, Z, rcond=None)[0]
    F, Rres = Xred @ gam, Z - Xred @ gam                      # reduced model (covariates only)
    rng = np.random.default_rng(args.seed)
    idx = [np.flatnonzero(site == s_) for s_ in np.unique(site)]
    null = np.empty((args.n_perm, len(masks)))
    for k in range(args.n_perm):                              # Freedman-Lane: permute reduced-model residuals within site
        perm = np.arange(len(d))
        for i in idx:
            perm[i] = i[rng.permutation(len(i))]
        null[k] = statistic(F + Rres[perm]).mean(axis=0)
    p_ind = np.where(np.isfinite(mean_r), (np.sum(np.abs(null) >= np.abs(mean_r), axis=0) + 1) / (args.n_perm + 1),
                     np.nan)
    ind = pd.DataFrame({"pet": pet.index, "mean_r_patients": mean_r, "null_sd": null.std(axis=0), "p_perm": p_ind,
                        "q": bh(p_ind)})
    # within patients: individual colocalization vs distress / hyperacusis (Spearman, permutation within site)
    pat = d[g].reset_index(drop=True)
    for v in ("distress_z", "hyperacusis_z"):
        ok = pat[v].notna().to_numpy()
        yv, sv, cv = pat.loc[ok, v].to_numpy(float), pat.loc[ok, "site"].to_numpy(), obs[ok]
        zc = zrank(cv, axis=0)
        r_obs = zrank(yv) @ zc / len(yv)
        pidx = [np.flatnonzero(sv == s_) for s_ in np.unique(sv)]
        nn = np.empty((args.n_perm, cv.shape[1]))
        for k in range(args.n_perm):
            yp = yv.copy()
            for i in pidx:
                yp[i] = yv[rng.permutation(i)]
            nn[k] = zrank(yp) @ zc / len(yp)
        pv = np.where(np.isfinite(r_obs), (np.sum(np.abs(nn) >= np.abs(r_obs), axis=0) + 1) / (args.n_perm + 1), np.nan)
        ind[f"r_{v}"], ind[f"p_{v}"], ind[f"q_{v}"], ind[f"n_{v}"] = r_obs, pv, bh(pv), int(ok.sum())
    ind.to_csv(args.out / "individual.tsv", sep="\t", index=False)

    summary = {"n_null": args.n_null, "n_perm": args.n_perm, "pet_maps": pet_info.to_dict("records"),
               "n_reference_maps": {k: len(v) for k, v in refs.items()},
               "consensus_colocalizations_q05_all_nulls": consensus.to_dict("records"),
               "min_q_by_target_set_null": {f"{a_} | {b_} | {c_}": v for (a_, b_, c_), v in
                                            col.groupby(["target", "reference_set", "null"]).q.min().items()},
               "null_robustness": rob.to_dict("records"), "individual": ind.to_dict("records"),
               "n_patients": int(g.sum()), "n_controls": int((~g).sum())}
    (args.out / "summary.json").write_text(json.dumps(clean(summary), indent=2, allow_nan=False))
    print(consensus.to_string(index=False) if len(consensus) else "no consensus colocalization")
    print(rob.pivot_table(index=["family", "target", "test"], columns="null", values="q").round(4).to_string())
    print(ind.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
