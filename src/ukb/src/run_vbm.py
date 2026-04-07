import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from tqdm.auto import tqdm


def run(cmd, check=True, capture_output=False, text=True):
    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=text,
    )


def ensure_dx_folder(remote_folder):
    run(["dx", "mkdir", "-p", remote_folder])


def find_gm_pve_file(subject_dir: Path) -> Path | None:
    c = subject_dir / "T1_fast" / "T1_brain_pve_1.nii.gz"
    if c.exists():
        return c
    return None


def smooth_nifti_mm(in_file: Path, out_file: Path, fwhm_mm: float = 7.0):
    img = nib.load(str(in_file))
    data = img.get_fdata(dtype=np.float32)

    zooms = img.header.get_zooms()[:3]
    sigma_mm = fwhm_mm / 2.354820045
    sigma_vox = [sigma_mm / z for z in zooms]

    smoothed = gaussian_filter(data, sigma=sigma_vox)

    out_img = nib.Nifti1Image(smoothed.astype(np.float32), img.affine, img.header)
    nib.save(out_img, str(out_file))


def modulate_gm(gm_mni_file: Path, jac_file: Path, out_file: Path):
    gm_img = nib.load(str(gm_mni_file))
    jac_img = nib.load(str(jac_file))

    gm_data = gm_img.get_fdata(dtype=np.float32)
    jac_data = jac_img.get_fdata(dtype=np.float32)

    gm_mod = gm_data * jac_data
    gm_mod = np.nan_to_num(gm_mod, nan=0.0, posinf=0.0, neginf=0.0)

    out_img = nib.Nifti1Image(gm_mod.astype(np.float32), gm_img.affine, gm_img.header)
    nib.save(out_img, str(out_file))


def extract_subject_zip(local_zip: Path, extract_dir: Path):
    with zipfile.ZipFile(local_zip, "r") as zf:
        zf.extractall(extract_dir)


def find_subject_root(extract_dir: Path) -> Path | None:
    for cand in extract_dir.rglob("*"):
        if not cand.is_dir():
            continue
        has_transforms = (cand / "transforms").exists()
        has_t1 = (cand / "T1.nii.gz").exists() or (cand / "T1_fast").exists()
        if has_transforms and has_t1:
            return cand
    return None


def get_processed_subjects(out_project_dir: str) -> set[str]:
    try:
        result = run(["dx", "ls", out_project_dir], capture_output=True)
        return {
            line.rstrip("/").strip()
            for line in result.stdout.splitlines()
            if line.strip().endswith("/")
        }
    except subprocess.CalledProcessError:
        return set()


def upload_final_map(subject_id: str, final_map: Path, out_project_dir: str):
    remote_folder = f"{out_project_dir}/{subject_id}"
    ensure_dx_folder(remote_folder)
    run(["dx", "upload", str(final_map), "--path", f"{remote_folder}/"])


def process_subject(subject_id: str, subject_dir: Path, work_out_dir: Path, mni_ref: Path, fwhm_mm: float = 6.0) -> Path:
    gm_pve = find_gm_pve_file(subject_dir)
    if gm_pve is None:
        raise FileNotFoundError(f"GM PVE not found for {subject_id}")

    warp = subject_dir / "transforms" / "T1_to_MNI_warp_coef.nii.gz"
    premat = subject_dir / "transforms" / "T1_to_MNI_linear.mat"

    if not warp.exists():
        raise FileNotFoundError(f"Missing warp file: {warp}")
    if not premat.exists():
        raise FileNotFoundError(f"Missing affine file: {premat}")

    work_out_dir.mkdir(parents=True, exist_ok=True)

    gm_mni = work_out_dir / "gm_mni.nii.gz"
    jac = work_out_dir / "jacobian.nii.gz"
    gm_mod = work_out_dir / "gm_mod_mni.nii.gz"
    gm_mod_s = work_out_dir / "gm_mod_mni_s3.nii.gz"

    run([
        "applywarp",
        f"--ref={mni_ref}",
        f"--in={gm_pve}",
        f"--warp={warp}",
        f"--premat={premat}",
        f"--out={gm_mni}",
        "--interp=spline",
    ])

    run([
        "fnirtfileutils",
        f"--in={warp}",
        f"--ref={mni_ref}",
        f"--jac={jac}",
    ])

    modulate_gm(gm_mni, jac, gm_mod)
    smooth_nifti_mm(gm_mod, gm_mod_s, fwhm_mm=fwhm_mm)

    return gm_mod_s


def main(init_subject: int, end_subject: int):
    # ---------------- config ----------------
    df_fname = "./results/ukb_vol.csv"
    t1_code = "20252"          # change if needed
    extension = "2_0.zip"      # change if needed
    out_project_dir = "/VBM"
    fwhm_mm = 7.0

    workdir = Path.cwd()
    raw_dir = workdir / "raw_zips"
    extract_base = workdir / "extract"
    subjects_dir = workdir / "subjects"
    work_out_base = workdir / "work_out"
    logs_dir = workdir / "logs"

    for d in [raw_dir, extract_base, subjects_dir, work_out_base, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    fsl_dir = Path(os.environ.get("FSLDIR", "/usr/local/fsl"))
    mni_ref = fsl_dir / "data" / "standard" / "MNI152_T1_2mm_brain.nii.gz"
    if not mni_ref.exists():
        raise FileNotFoundError(f"MNI reference not found: {mni_ref}")

    # ---------------- subjects table ----------------
    run(["dx", "download", df_fname, "-o", "ukb_vol.csv"])
    df = pd.read_csv("ukb_vol.csv")
    df["eid"] = df["eid"].astype(str)

    all_subjects = df["eid"].tolist()

    processed_subjects = get_processed_subjects(out_project_dir)
    missing_subjects = [s for s in all_subjects if s not in processed_subjects]
    subjects = sorted(missing_subjects)[init_subject:end_subject]

    print(f"Found {len(processed_subjects)} processed subjects")
    print(f"Found {len(missing_subjects)} missing subjects")
    print(f"Running {len(subjects)} subjects")

    errors = []

    for subject in tqdm(subjects):
        folder_id = subject[:2]
        zip_name = f"{subject}_{t1_code}_{extension}"
        dx_path = f"/Bulk/Brain MRI/T1/{folder_id}/{zip_name}"

        local_zip = raw_dir / zip_name
        extract_dir = extract_base / f"extract_{subject}"
        local_subject_dir = subjects_dir / subject
        work_out_dir = work_out_base / subject

        for p in [local_zip]:
            if p.exists():
                p.unlink()
        for p in [extract_dir, local_subject_dir, work_out_dir]:
            if p.exists():
                shutil.rmtree(p, ignore_errors=True)

        extract_dir.mkdir(parents=True, exist_ok=True)

        try:
            run(["dx", "download", dx_path, "-o", str(local_zip)])
            extract_subject_zip(local_zip, extract_dir)

            found_subject_root = find_subject_root(extract_dir)
            if found_subject_root is None:
                raise FileNotFoundError(f"Could not locate extracted subject folder for {subject}")

            shutil.move(str(found_subject_root), str(local_subject_dir))

            final_map = process_subject(
                subject_id=subject,
                subject_dir=local_subject_dir,
                work_out_dir=work_out_dir,
                mni_ref=mni_ref,
                fwhm_mm=fwhm_mm,
            )

            upload_final_map(
                subject_id=subject,
                final_map=final_map,
                out_project_dir=out_project_dir,
            )

        except Exception as e:
            print(f"Error for {subject}: {e}")
            errors.append({"subject_id": subject, "error": str(e)})

        finally:
            if local_zip.exists():
                local_zip.unlink()
            shutil.rmtree(extract_dir, ignore_errors=True)
            shutil.rmtree(local_subject_dir, ignore_errors=True)
            shutil.rmtree(work_out_dir, ignore_errors=True)

    errors_df = pd.DataFrame(errors)
    errors_csv = logs_dir / f"vbm_errors_{init_subject}_{end_subject}.csv"
    errors_df.to_csv(errors_csv, index=False)

    ensure_dx_folder(f"{out_project_dir}/_logs")
    run(["dx", "upload", str(errors_csv), "--path", f"{out_project_dir}/_logs/"])

    print("done")


if __name__ == "__main__":
    main(init_subject=0, end_subject=1500)