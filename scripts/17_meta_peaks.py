"""Replication of published meta-analytic tinnitus VBM peaks in our data (analysis plan 6.10 B1, fixed 2026-11-04).

Per peak (analysis.yaml meta_peaks): mean CAT12 modulated GM (mwp1) in a sphere around the MNI coordinate, GM voxels
only (mean > maps_palm.vbm_mask_threshold); coordinates as published (e.g. the "R PCC" peak lies at x = 0); mixed
model (site random intercept) of the matching contrast — all / hearing-matched peaks: tinnitus vs controls
(primary sample); normal-hearing peak: TINN vs CTRL; hearing-loss peaks: TNHL vs HLOS (both by group4_hf) — with
age spline + sex + Euler + log CAT12 TIV, without and with PTA4 + PTA_HF beyond PTA4; TOST |d| < .2, BF01, sign
agreement with the published direction. -> <results>/17_meta_peaks/{peaks.tsv, summary.json}
"""

from __future__ import annotations

import json

import nibabel as nib
import numpy as np
import pandas as pd

from tinception import config
from tinception.stats import evidence as ev
from tinception.stats import group as g

ID = "participant_id"


def sphere_means(files: list, peaks: dict, radius: float, gm_threshold: float) -> tuple[np.ndarray, list[int]]:
    """subjects x peaks: mean voxel value within `radius` mm of each MNI coordinate, restricted to grey-matter voxels
    (mean value across subjects > gm_threshold, as the PALM VBM mask) so that spheres at the brain edge are not
    diluted by background / CSF. Only the bounding box of all spheres is read from each image.
    Returns the means and the number of GM voxels per sphere."""
    ref = nib.load(str(files[0]))
    ijk = np.indices(ref.shape[:3]).reshape(3, -1).T
    xyz = nib.affines.apply_affine(ref.affine, ijk)
    spheres = []
    for name, (x, y, z, *_rest) in peaks.items():
        m = np.linalg.norm(xyz - np.array([x, y, z]), axis=1) <= radius
        if m.sum() == 0:
            raise SystemExit(f"{name}: sphere outside the image grid")
        spheres.append(m.reshape(ref.shape[:3]))
    union = np.any(spheres, axis=0)
    lo, hi = np.argwhere(union).min(axis=0), np.argwhere(union).max(axis=0) + 1
    box = tuple(slice(a_, b_) for a_, b_ in zip(lo, hi, strict=True))
    sub = [sp[box] for sp in spheres]
    vols = []
    for f in files:
        img = nib.load(str(f))
        if img.shape[:3] != ref.shape[:3] or not np.allclose(img.affine, ref.affine):
            raise SystemExit(f"{f}: different grid")
        vols.append(np.asarray(img.dataobj[box], np.float32))
    V = np.stack(vols)                                       # subjects x box
    gm = V.mean(axis=0) > gm_threshold
    masks = [sp & gm for sp in sub]
    n_vox = [int(m.sum()) for m in masks]
    if min(n_vox) == 0:
        raise SystemExit(f"a sphere has no GM voxels: {dict(zip(peaks, n_vox, strict=True))}")
    return np.stack([V[:, m].mean(axis=1) for m in masks], axis=1), n_vox


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peak-set", default="peaks",
                    help="key under analysis.yaml meta_peaks (default: the published B1 peaks; ale_peaks = 6.10 B2 ALE "
                         "clusters -> <results>/17_meta_peaks_ale_peaks/)")
    ap.add_argument("--exclude-sites", nargs="*", default=[],
                    help="sites left out (independence check, e.g. sites whose subjects overlap a source study)")
    args = ap.parse_args()
    cfg, p = config.analysis(), config.paths()
    mcfg, gcfg = cfg["meta_peaks"], cfg["group"]
    out = p.results / ("17_meta_peaks" if args.peak_set == "peaks" else f"17_meta_peaks_{args.peak_set}")
    if args.exclude_sites:
        out = out.with_name(out.name + "_without_" + "_".join(sorted(args.exclude_sites)))
    out.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(p.root / "tinception_master.tsv", sep="\t", na_values="n/a", low_memory=False)
    d = master[master.in_case_control.astype(bool) & master.PTA_HF.notna() & master.vbm_ok.astype(bool)
               & master.cat12_TIV.notna() & ~master.site.isin(args.exclude_sites)].reset_index(drop=True)
    files = [p.derivatives / "cat12-12.9" / s_ / "mri" / f"mwp1{s_}_T1w.nii" for s_ in d[ID]]
    absent = [str(f) for f in files if not f.exists()]
    if absent:
        raise SystemExit(f"{len(absent)} mwp1 images missing, e.g. {absent[:2]}")
    peaks = mcfg[args.peak_set]
    vals, n_vox = sphere_means(files, peaks, mcfg["radius_mm"], cfg["maps_palm"]["vbm_mask_threshold"])
    d = pd.concat([d, pd.DataFrame(vals, columns=list(peaks))], axis=1)
    d = d.assign(log_TIV=np.log(d.cat12_TIV))
    # complete cases first, so counts, residualization and fits use the same rows; both hearing subgroups use the
    # PTA_HF-based groups (group4_hf), as the H1 clean-hearing contrast (normal hearing = TINN vs CTRL,
    # hearing loss = TNHL vs HLOS; complementary)
    d = d.dropna(subset=["age", "sex", "euler_total", "PTA4", "PTA_HF", "log_TIV", "site", *peaks]).reset_index(drop=True)
    samples = {"all": d.assign(tinnitus=d.group.eq("TI").astype(int)),
               "normal_hearing": d[d.group4_hf.isin(["CTRL", "TINN"])].assign(
                   tinnitus=lambda x: x.group4_hf.eq("TINN").astype(int)),
               "hearing_loss": d[d.group4_hf.isin(["TNHL", "HLOS"])].assign(
                   tinnitus=lambda x: x.group4_hf.eq("TNHL").astype(int))}
    samples = {k: g.residualize(v.reset_index(drop=True), "PTA_HF", ["PTA4"], "PTA_HF_res")
               for k, v in samples.items()}
    cov = [f"C({c})" if c == "sex" else c for c in gcfg["covariates"]] + ["log_TIV"]
    cont = ["PTA4", "PTA_HF_res", "euler_total", "log_TIV"]
    rows = []
    for name, (_x, _y, _z, direction, contrast) in peaks.items():
        s_ = samples[contrast]
        pr = g.prepare(s_, [name], set())
        for label, terms in (("no_hearing", ["tinnitus", *cov]), ("hearing", ["tinnitus", "PTA4", "PTA_HF_res", *cov])):
            rhs = {y: g.rhs(terms, gcfg["age_df"], False) for y in pr.attrs["outcomes"]}
            con = {"tinnitus": {"tinnitus": 1.0}} | ({"PTA4": {"PTA4": 1.0}, "PTA_HF_res": {"PTA_HF_res": 1.0}}
                                                      if label == "hearing" else {})
            rr = g.run_family(pr, rhs, con, contrast, gcfg["tost_d"], cont).set_index("contrast")
            r = rr.loc["tinnitus"]
            hear = {f"{h}_{k}": float(rr.loc[h, c]) for h in ("PTA4", "PTA_HF_res") if h in rr.index
                    for k, c in (("beta", "beta"), ("p", "p"))}   # hearing terms at the same sphere (descriptive)
            rows.append(hear | {"peak": name, "contrast": contrast, "model": label, "published_direction": int(direction),
                         "gm_voxels": n_vox[list(peaks).index(name)], "n": int(r.n), "n_tinnitus": int(s_.tinnitus.sum()), "d": float(r.beta), "ci_lo": float(r.ci_lo),
                         "ci_hi": float(r.ci_hi), "p": float(r.p), "p_tost": float(r.p_tost),
                         "bf01": float(np.exp(-ev.log_bf10(r.beta, r.se))), "same_sign": bool(np.sign(r.beta) == direction) if direction else None,   # 0 = mixed / no direction
                         "converged": bool(r.converged)})
    t = pd.DataFrame(rows)
    t["q"] = t.p.where(t.converged).groupby(t.model).transform(g.bh)   # converged fits only (as run_family)
    t.to_csv(out / "peaks.tsv", sep="\t", index=False)
    summary = {"radius_mm": mcfg["radius_mm"], "n_primary": len(d),
               "by_model": {m: {"n_peaks": len(x), "n_p05": int((x.p < .05).sum()), "n_q05": int((x.q < .05).sum()),
                                "n_equivalent": int((x.p_tost < .05).sum()), "n_same_sign": int(x.same_sign.eq(True).sum()),
                                "median_bf01": float(x.bf01.median()), "max_abs_d": float(x.d.abs().max())}
                            for m, x in t.groupby("model")},
               "peaks": t.to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False))
    print(t[["peak", "contrast", "model", "n", "d", "ci_lo", "ci_hi", "p", "q", "p_tost", "bf01", "same_sign"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
