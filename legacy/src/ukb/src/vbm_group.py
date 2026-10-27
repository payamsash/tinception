from pathlib import Path
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
import nibabel as nib

from nilearn.glm.second_level import non_parametric_inference, SecondLevelModel
from nilearn.image import resample_to_img
from nilearn.plotting import plot_stat_map


# --------------------------------------------------
# Config
# --------------------------------------------------
REMOTE_BASE = "/VBM"
OUTDIR = Path("VBM")
TMPDIR = Path("tmp_vbm")
OUTDIR.mkdir(exist_ok=True, parents=True)
TMPDIR.mkdir(exist_ok=True, parents=True)

# --------------------------------------------------
# 1) Load and prepare dataframe
# --------------------------------------------------
df_fname = "./GWAS/main/master_gwas.tsv"
subprocess.run(["dx", "download", df_fname, "-o", "master.tsv"], check=True)

df = pd.read_csv("master.tsv", sep="\t")
df = df[["subjects", "age", "sex", "srt", "eTIV", "tin_status"]].copy()
df["subjects"] = df["subjects"].astype(str)

# map group labels for FSGD
df["tin_status"] = df["tin_status"].map({0: "CO", 1: "TI"})

# optional but recommended: center continuous covariates
for col in ["age", "srt", "eTIV"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
    df[col] = df[col] - df[col].mean()


# --------------------------------------------------
# 2) Find common subjects
# --------------------------------------------------
result = subprocess.run(
    ["dx", "ls", REMOTE_BASE],
    capture_output=True,
    text=True,
    check=True
)

processed_subjects = {
    line.rstrip("/").strip()
    for line in result.stdout.splitlines()
    if line.strip().endswith("/")
}

common_subjects = [s for s in df["subjects"].tolist() if s in processed_subjects]
print(f"Found {len(common_subjects)} common subjects")

df = df[df["subjects"].isin(common_subjects)].copy()
df["subjects"] = pd.Categorical(df["subjects"], categories=common_subjects, ordered=True)
df = df.sort_values("subjects").reset_index(drop=True)
df["subjects"] = df["subjects"].astype(str)

# drop rows with missing values in required columns
df = df.dropna(subset=["subjects", "age", "sex", "srt", "eTIV", "tin_status"]).copy()


# --------------------------------------------------
# 3) Parallel download
# --------------------------------------------------
def download_one(subj: str):
    gm_remote = f"{REMOTE_BASE}/{subj}/gm_mod_mni_s3.nii.gz"


    gm_local = TMPDIR / f"{subj}.nii.gz"
    if not gm_local.exists():
        try:
            subprocess.run(
                ["dx", "download", gm_remote, "-o", str(gm_local)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return subj, True, None
        except subprocess.CalledProcessError as e:
            return subj, False, str(e)
        
    return subj, True, None

max_workers = min(16, max(4, (os.cpu_count() or 8)))
downloaded = []
failed = []

with ThreadPoolExecutor(max_workers=max_workers) as ex:
    futures = {ex.submit(download_one, subj): subj for subj in df["subjects"]}
    for fut in as_completed(futures):
        subj, ok, err = fut.result()
        if ok:
            downloaded.append(subj)
        else:
            failed.append((subj, err))

print(f"Downloaded {len(downloaded)} subjects")
if failed:
    print(f"Failed {len(failed)} subjects")
    print("Examples:", failed[:5])
    
downloaded_set = set(downloaded)
final_subjects = [s for s in df["subjects"].tolist() if s in downloaded_set]
df = df[df["subjects"].isin(final_subjects)].copy()
    
## add site to covariates
match_fname = f"./results/ukb_vol_harmonized.csv"
subprocess.run(["dx", "download", match_fname], check=True)
df_matched = pd.read_csv("ukb_vol_harmonized.csv")
df_matched["subject_id"] = df_matched["subject_id"].astype(str)
df_matched.rename(columns={"subject_id": "subjects"}, inplace=True)

df = df.merge(
            df_matched[["subjects", "SITE"]],
            on="subjects",
            how="inner"
            )


img_dir = Path("tmp_vbm")
out_dir = Path("VBM")
out_dir.mkdir(parents=True, exist_ok=True)
df = df.copy()
df["img_path"] = df["subjects"].astype(str).apply(lambda s: img_dir / f"{s}.nii.gz")
df = df[df["img_path"].apply(lambda p: p.exists())].copy()
df = df.sort_values("subjects").reset_index(drop=True)


df["tin_status_num"] = df["tin_status"].map({"CO": 0, "TI": 1})
df["sex_num"] = df["sex"]

for col in ["age", "srt", "eTIV"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
    df[col] = df[col] - df[col].mean()

site_dummies = pd.get_dummies(df["SITE"], prefix="SITE", drop_first=True)

design_matrix = pd.concat(
    [
        pd.Series(1.0, index=df.index, name="intercept"),
        df[["tin_status_num", "age", "sex_num", "srt", "eTIV"]],
        site_dummies,
    ],
    axis=1,
)
design_matrix = design_matrix.astype(float)
valid = design_matrix.notnull().all(axis=1)
df = df.loc[valid].reset_index(drop=True)
design_matrix = design_matrix.loc[valid].reset_index(drop=True)
second_level_input = df["img_path"].tolist()

mask_img = None  # or your explicit GM mask path

# --------------------------------------------------
# 3) Permutation TFCE helper
# --------------------------------------------------
def run_perm(name, contrast):
    outputs = non_parametric_inference(
        second_level_input=second_level_input,
        design_matrix=design_matrix,
        second_level_contrast=contrast,
        mask=mask_img,
        model_intercept=False,   # because intercept is already in design matrix
        n_perm=5000,
        two_sided_test=False,    # directional test
        tfce=True,
        smoothing_fwhm=None,
        n_jobs=-1,
        random_state=42,
        verbose=1,
    )
    for key, img in outputs.items():
        nib.save(img, out_dir / f"{name}_{key}.nii.gz")
    return outputs

# --------------------------------------------------
# 4) Run both directions
# --------------------------------------------------
perm_ti_gt_co = run_perm("TI_gt_CO", "tin_status_num")
perm_co_gt_ti = run_perm("CO_gt_TI", "-tin_status_num")

# --------------------------------------------------
# 5) Standard GLM for effect maps
# --------------------------------------------------
print("working on second level ...")
slm = SecondLevelModel(mask_img=mask_img, smoothing_fwhm=None)
slm = slm.fit(second_level_input, design_matrix=design_matrix)

# Effect maps (beta / contrast estimate)
effect_ti_gt_co = slm.compute_contrast(
    second_level_contrast="tin_status_num",
    output_type="effect_size",
)
effect_co_gt_ti = slm.compute_contrast(
    second_level_contrast="-tin_status_num",
    output_type="effect_size",
)

# Optional: uncorrected stat maps from the parametric GLM
stat_ti_gt_co = slm.compute_contrast(
    second_level_contrast="tin_status_num",
    output_type="stat",
)
stat_co_gt_ti = slm.compute_contrast(
    second_level_contrast="-tin_status_num",
    output_type="stat",
)

nib.save(effect_ti_gt_co, out_dir / "TI_gt_CO_effect_size.nii.gz")
nib.save(effect_co_gt_ti, out_dir / "CO_gt_TI_effect_size.nii.gz")
nib.save(stat_ti_gt_co, out_dir / "TI_gt_CO_parametric_stat.nii.gz")
nib.save(stat_co_gt_ti, out_dir / "CO_gt_TI_parametric_stat.nii.gz")
subprocess.run(["dx", "upload", "-r", "VBM/", "--dest", "/VBM_results/"], check=True)
