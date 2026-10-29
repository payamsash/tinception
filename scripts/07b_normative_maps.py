"""Normative models at vertex and voxel level (phase 2): the ROI model of scripts/07_normative.py, batched.

Settings: analysis.yaml -> normative (model, covariate ladder, references) and normative.maps (sets, chunking).

Subcommands
  prepare  build the subject x feature matrices in <derivatives>/normative/maps/data/:
             vertex_<measure>: fsaverage5 maps (features/vertex_fsavg5_fwhm10_<measure>.npy), medial wall removed
             voxel_gm: CAT12 s6mwp1 resampled to <voxel_mm> mm, voxels with mean modulated GM > threshold
           -> <set>.npy (float32), <set>_subjects.txt, <set>_index.npy (original vertex / flat voxel index),
              voxel_gm_mask.nii.gz
  tasks    task list (reference x model x set x chunk) -> maps/tasks.tsv
  fit      one task: cross-fitted z for the reference sample, full-reference z for everyone else (no model saved)
           -> maps/<reference>/<model>/<set>/chunk_<k>.npz (z, ids, reference flag, folds, feature ids, settings)
              + chunk_<k>_stats.tsv
  collect  check that the chunks cover every feature once with identical settings, then assemble ->
           maps/<reference>/<model>/<set>/{z.npy, subjects.tsv, features.txt, settings.json, calibration*.tsv/json,
           feature_calibration.tsv, extreme.tsv}; map-level calibration criteria (normative.calibration.maps)
Feature ids are v<original index>: fsaverage5 vertex (0-20483, lh then rh) or flat voxel index in the mask grid.

Example (cluster, repo root)
  sbatch slurm/07b_prepare.sbatch                 # prepare (CAT12 resampling takes a while)
  TINCEPTION_ENV=cluster uv run scripts/07b_normative_maps.py tasks
  T=/shares/tinnitus.orl.med.uzh/tinception/derivatives/normative/maps/tasks.tsv
  sbatch --array=1-$(($(wc -l < $T)-1))%100 slurm/07b_normative_maps.sbatch $T
  sbatch slurm/07b_collect.sbatch
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tinception import config
from tinception.normative import core


def _dirs(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    p = config.paths()
    maps = args.out or p.derivatives / "normative" / "maps"
    return args.master or p.root / "tinception_master.tsv", maps, maps / "data"


def _load(data: Path, name: str) -> tuple[np.ndarray, list[str]]:
    x = np.load(data / f"{name}.npy", mmap_mode="r")
    subjects = (data / f"{name}_subjects.txt").read_text().split()
    if x.shape[0] != len(subjects):
        raise ValueError(f"{name}: {x.shape[0]} rows for {len(subjects)} subjects")
    return x, subjects


FALLBACK_RANK = {"refit": 1, "gaussian": 2, "failed": 3}
VERTEX_MEASURES = ("thickness", "area")  # share one vertex set (positivity checked on both)


def vertex_keep(feats: Path, fsaverage: Path) -> np.ndarray:
    """Cortex vertices (no medial wall) that are finite and > 0 for every scanned subject in thickness AND area: a
    zero marks a vertex outside that subject's cortex (medial-wall boundary), and zero-inflated vertices break the
    warped BLR fit. One mask for both measures keeps the vertex sets identical. Deliberately all scanned subjects
    (stricter than the analysed ones; 2026-10-28: 18,363 vs 18,376 vertices, i.e. 13 vertices lost)."""
    keep = core.fsavg5_cortex_mask(fsaverage)
    for m in VERTEX_MEASURES:
        x = np.load(feats / f"vertex_fsavg5_fwhm10_{m}.npy", mmap_mode="r")
        keep &= np.isfinite(x).all(axis=0) & (np.asarray(x) > 0).all(axis=0) & (x.std(axis=0) > 0)
    return keep


def prepare_vertex(spec: dict, data: Path) -> None:
    p = config.paths()
    feats = p.derivatives / "features"
    x = np.load(feats / f"vertex_fsavg5_fwhm10_{spec['measure']}.npy")
    subjects = (feats / "vertex_subjects.txt").read_text().split()
    cortex = core.fsavg5_cortex_mask(p.freesurfer / "fsaverage")
    keep = vertex_keep(feats, p.freesurfer / "fsaverage")
    name = f"vertex_{spec['measure']}"
    np.save(data / f"{name}.npy", x[:, keep].astype(np.float32))
    np.save(data / f"{name}_index.npy", np.flatnonzero(keep))
    (data / f"{name}_subjects.txt").write_text("\n".join(subjects) + "\n")
    print(f"{name}: {len(subjects)} subjects x {keep.sum()} vertices (cortex {cortex.sum()} of {cortex.size}; "
          f"{cortex.sum() - keep.sum()} cortex vertices dropped for values <= 0 in any subject)")


def prepare_voxel(spec: dict, data: Path, cfg: dict, master: pd.DataFrame) -> None:
    import nibabel as nib
    from nilearn.image import resample_img

    root = config.paths().derivatives / f"cat12-{config.analysis()['cat12']['version']}"
    eligible = master[master.in_case_control.astype(bool) & master.vbm_ok.astype(bool)].participant_id
    files = {pid: root / pid / "mri" / f"{spec['input']}{pid}_T1w.nii" for pid in eligible}  # mask from analysed
    files = {pid: f for pid, f in files.items() if f.is_file()}                              # subjects only
    if not files:
        raise SystemExit(f"no {spec['input']} files for eligible subjects under {root}")
    ref = nib.load(str(next(iter(files.values()))))
    vox = cfg["voxel_mm"]
    zoom = np.abs(np.diag(ref.affine)[:3])
    affine = ref.affine.copy()
    affine[:3, :3] = np.diag(np.sign(np.diag(ref.affine)[:3]) * vox)
    shape = tuple(int(np.ceil(s * z / vox)) for s, z in zip(ref.shape[:3], zoom, strict=True))
    imgs = []
    for f in files.values():
        img = nib.load(str(f))
        if not np.allclose(img.affine, ref.affine):
            raise ValueError(f"{f}: grid differs from {next(iter(files.values()))}")
        r = resample_img(img, target_affine=affine, target_shape=shape, interpolation="continuous",
                         force_resample=True, copy_header=True)
        imgs.append(np.asarray(r.dataobj, dtype=np.float32).ravel())
    x = np.vstack(imgs)
    keep = (x.mean(axis=0) > cfg["voxel_mask_threshold"]) & (x.std(axis=0) > 0)
    nib.save(nib.Nifti1Image(keep.reshape(shape).astype(np.uint8), affine), str(data / "voxel_gm_mask.nii.gz"))
    np.save(data / "voxel_gm.npy", x[:, keep])
    np.save(data / "voxel_gm_index.npy", np.flatnonzero(keep))
    (data / "voxel_gm_subjects.txt").write_text("\n".join(files) + "\n")
    print(f"voxel_gm: {len(files)} subjects x {keep.sum()} voxels ({vox} mm, mean GM > {cfg['voxel_mask_threshold']})")


def cmd_prepare(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]["maps"]
    master_f, _, data = _dirs(args)
    data.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    for name, spec in cfg["sets"].items():
        if args.sets and name not in args.sets:
            continue
        if spec["source"] == "vertex":
            prepare_vertex(spec, data)
        else:
            prepare_voxel(spec, data, cfg, master)


def _set_file(name: str, spec: dict) -> str:
    return f"vertex_{spec['measure']}" if spec["source"] == "vertex" else name


def cmd_tasks(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    _, maps, data = _dirs(args)
    rows = []
    for name, spec in cfg["maps"]["sets"].items():
        if args.sets and name not in args.sets:
            continue
        n = np.load(data / f"{_set_file(name, spec)}.npy", mmap_mode="r").shape[1]
        for k, _ in enumerate(core.chunks(n, cfg["maps"]["chunk"])):
            rows += [(r, m, name, k) for r in cfg["maps"]["references"] for m in cfg["models"]]
    t = pd.DataFrame(rows, columns=["reference", "model", "feature_set", "chunk"])
    dest = maps / (args.name or "tasks.tsv")  # never rewrite a task file a running array still reads by row
    t.to_csv(dest, sep="\t", index=False)
    print(f"{len(t)} tasks -> {dest}")


def cmd_fit(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, maps, data = _dirs(args)
    row = pd.read_csv(args.task_file, sep="\t").iloc[args.task - 1]
    spec = cfg["maps"]["sets"][row.feature_set]
    x, subjects = _load(data, _set_file(row.feature_set, spec))
    a, b = core.chunks(x.shape[1], cfg["maps"]["chunk"])[row.chunk]
    index = np.load(data / f"{_set_file(row.feature_set, spec)}_index.npy")
    features = [f"v{j}" for j in index[a:b]]
    feats = pd.DataFrame(np.asarray(x[:, a:b], dtype=np.float64), columns=features).assign(
        **{core.ID: subjects})[[core.ID, *features]]
    tiv = [spec.get("tiv_column", cfg["tiv_column"])] if spec["tiv"] else []
    covs = list(cfg["models"][row.model]) + tiv
    ref_spec = cfg["references"][row.reference]
    mcfg = core.reference_cfg(cfg, ref_spec)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)
    if spec["source"] == "voxel":  # CAT12 image-quality rule on top of the FreeSurfer QC in in_case_control
        master = master[master.vbm_ok.astype(bool)]
    df, ref = core.design(master, feats, covs, cfg["batch_effects"], ref_spec)
    z, stats = core.cross_fit(df, ref, features, covs, cfg["batch_effects"], mcfg, config.analysis()["seed"])
    dest = maps / row.reference / row.model / row.feature_set
    dest.mkdir(parents=True, exist_ok=True)
    np.savez(dest / f"chunk_{row.chunk:04d}.npz", z=z[features].to_numpy(np.float32), ids=z[core.ID].to_numpy(str),
             reference=ref.to_numpy(), fold=z["fold"].to_numpy(), columns=np.arange(a, b), n_features=x.shape[1],
             settings=json.dumps(core.model_settings(mcfg, covs), sort_keys=True))
    stats.to_csv(dest / f"chunk_{row.chunk:04d}_stats.tsv", sep="\t", index=False, na_rep="n/a")
    held = z.loc[ref.to_numpy(), features].to_numpy()
    print(f"{row.reference}/{row.model}/{row.feature_set} chunk {row.chunk} [{a}:{b}): {int(ref.sum())} reference / "
          f"{len(df)} scored; held-out z mean {np.nanmean(held):+.3f} sd {np.nanstd(held):.3f}; "
          f"columns with missing z {int(np.isnan(z[features].to_numpy()).any(axis=0).sum())}")


def cmd_collect(args: argparse.Namespace) -> None:
    cfg = config.analysis()["normative"]
    master_f, maps, data = _dirs(args)
    master = pd.read_csv(master_f, sep="\t", na_values="n/a", low_memory=False)[[core.ID, *cfg["batch_effects"]]]
    for ref_name in cfg["maps"]["references"]:
        for model in cfg["models"]:
            for name, spec in cfg["maps"]["sets"].items():
                d = maps / ref_name / model / name
                n = np.load(data / f"{_set_file(name, spec)}.npy", mmap_mode="r").shape[1]
                index = np.load(data / f"{_set_file(name, spec)}_index.npy")
                files = [d / f"chunk_{k:04d}.npz" for k in range(len(core.chunks(n, cfg["maps"]["chunk"])))]
                missing = [f.name for f in files if not f.exists()]
                if missing:
                    print(f"SKIP {ref_name}/{model}/{name}: {len(missing)} of {len(files)} chunks missing")
                    continue
                parts = [np.load(f) for f in files]
                ids = parts[0]["ids"]
                problems = []
                if not np.array_equal(np.concatenate([p["columns"] for p in parts]), np.arange(n)):
                    problems.append("chunks do not cover the current features exactly once (stale chunks?)")
                if any(int(p["n_features"]) != n for p in parts):
                    problems.append("chunks fitted on a different prepared matrix")
                if len({str(p["settings"]) for p in parts}) > 1:
                    problems.append("chunks fitted with different settings")
                if any(not np.array_equal(p["ids"], ids) for p in parts):
                    problems.append("chunks scored different subjects")
                if problems:
                    print(f"SKIP {ref_name}/{model}/{name}: {'; '.join(problems)}")
                    continue
                z = np.hstack([p["z"] for p in parts])
                reference = parts[0]["reference"]
                cols = [f"v{j}" for j in index]
                np.save(d / "z.npy", z)
                (d / "features.txt").write_text("\n".join(cols) + "\n")
                (d / "settings.json").write_text(json.dumps(json.loads(str(parts[0]["settings"])), indent=2))
                pd.DataFrame({core.ID: ids, "reference": reference}).to_csv(d / "subjects.tsv", sep="\t", index=False)
                full = pd.DataFrame(z, columns=cols).assign(**{core.ID: ids})
                core.extreme_counts(full, cols, cfg["extreme_z"]).to_csv(d / "extreme.tsv", sep="\t", index=False)
                held = full[reference].merge(master, on=core.ID, how="left", validate="one_to_one")
                st = pd.concat([pd.read_csv(f.with_name(f"{f.stem}_stats.tsv"), sep="\t", na_values="n/a")
                                for f in files])
                fb = st.dropna(subset=["fallback"]) if "fallback" in st else st.iloc[:0]
                # one count per feature: its worst outcome in any of the fits
                per_feature = fb.groupby("feature").fallback.agg(lambda v: max(v, key=FALLBACK_RANK.get)) \
                    if len(fb) else pd.Series(dtype=object)
                fallback = {k: int(v) for k, v in per_feature.value_counts().items()}
                pooled = core.write_calibration(d, held, held[cfg["batch_effects"]], cols, cfg["calibration"], "map",
                                                {"n_reference": int(reference.sum()), "n_scored": len(ids),
                                                 "n_features": n, "fallback_features": fallback})
                print(f"{ref_name}/{model}/{name}: {n} features, {len(ids)} subjects; held-out mean {pooled['mean']}, "
                      f"sd {pooled['sd']}; fallback {pooled['fallback_features']}; features failing "
                      f"{pooled['features_failing']}, batch FDR hits "
                      f"{pooled['batch_fdr_hits']} -> {'OK' if pooled['ok'] else 'CHECK'}; sites failing mean CI "
                      f"{pooled['sites_failing_mean']}, skew {pooled['sites_failing_skew']}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None, help="default: <root>/tinception_master.tsv")
    ap.add_argument("--out", type=Path, default=None, help="default: <derivatives>/normative/maps")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--sets", nargs="*", default=None, help="subset of normative.maps.sets")
    p.set_defaults(func=cmd_prepare)
    t = sub.add_parser("tasks")
    t.add_argument("--sets", nargs="*", default=None, help="subset of normative.maps.sets")
    t.add_argument("--name", default=None, help="task file name in maps/ (default tasks.tsv)")
    t.set_defaults(func=cmd_tasks)
    f = sub.add_parser("fit")
    f.add_argument("--task-file", type=Path, required=True)
    f.add_argument("--task", type=int, required=True, help="1-based row of --task-file (SLURM_ARRAY_TASK_ID)")
    f.set_defaults(func=cmd_fit)
    sub.add_parser("collect").set_defaults(func=cmd_collect)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
