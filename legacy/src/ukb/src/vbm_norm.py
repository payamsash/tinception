import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"

from pathlib import Path
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel
)

import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import nibabel as nib
from multiprocessing import Pool

#os.environ["OMP_NUM_THREADS"] = "8"
#os.environ["MKL_NUM_THREADS"] = "8"
#os.environ["OPENBLAS_NUM_THREADS"] = "8"

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

# if sex is strings, uncomment and adapt:
# df["sex"] = df["sex"].map({"M": 0, "F": 1})

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
df.rename(columns={"subjects": "subject_id"}, inplace=True)

### bringing site information
match_fname = f"./results/ukb_vol_harmonized.csv"
subprocess.run(["dx", "download", match_fname], check=True)
df_matched = pd.read_csv("ukb_vol_harmonized.csv")
df_matched["subject_id"] = df_matched["subject_id"].astype(str)

## merging 
df = df.merge(df_matched[["subject_id", "SITE"]], on="subject_id", how="inner") 

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
    futures = {ex.submit(download_one, subj): subj for subj in df["subject_id"]}
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
final_subjects = [s for s in df["subject_id"].tolist() if s in downloaded_set]
df = df[df["subject_id"].isin(final_subjects)].copy()
# df = df.query('subject_id != "5223474"')

# --------------------------------------------------
# 4) Preparing for norm modeling
# --------------------------------------------------
subprocess.run(["dx", "download", "subcortical_mask_thr80.nii.gz"], check=True)
mask = nib.load('subcortical_mask_thr80.nii.gz').get_fdata() > 0
rows = []
for subj in tqdm(df["subject_id"]):
    gm_file = TMPDIR / f"{subj}.nii.gz"
    img = nib.load(str(gm_file))
    rows.append(img.get_fdata()[mask])

y = np.vstack(rows)
print(y.shape)
voxel_cols = [f"v{i}" for i in range(1, y.shape[1] + 1)]
y_df = pd.DataFrame(y, columns=voxel_cols)
df = pd.concat([df.reset_index(drop=True), y_df.reset_index(drop=True)], axis=1)
df.dropna(inplace=True)
print(df.head(5))
print(df.shape)

# IMPORTANT: avoid BLAS oversubscription
vbm_norm_dir = Path("./vbm_norm_ukb")

base_kwargs = {
    "covariates": ['age', 'sex', 'srt', 'eTIV'],
    "batch_effects": ["SITE"],
    "subject_ids": "subject_id"
}

# split response variables into chunks
all_response_vars = list(df.columns[7:])
n_jobs = 8
chunks = np.array_split(all_response_vars, n_jobs)

def run_chunk(args):
    i, response_subset = args

    print(f"Starting chunk {i} with {len(response_subset)} variables")

    kwargs = base_kwargs.copy()
    kwargs["response_vars"] = list(response_subset)

    norm_train = NormData.from_dataframe(
        name=f"train_chunk_{i}",
        dataframe=df[df['tin_status'] == "CO"],
        **kwargs
    )

    norm_test = NormData.from_dataframe(
        name=f"test_chunk_{i}",
        dataframe=df,
        **kwargs
    )

    template_blr = BLR(
        name=f"payam_blr_{i}",
        basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
        fixed_effect=True,
        heteroskedastic=True,
        warp_name="warpsinharcsinh"
    )

    model = NormativeModel(
        template_regression_model=template_blr,
        savemodel=False,
        evaluate_model=True,
        saveresults=True,
        saveplots=False,
        save_dir=str(vbm_norm_dir / f"norm_model_chunk_{i}"),
        inscaler="standardize",
        outscaler="none",
    )

    model.fit_predict(norm_train, norm_test)
    print(f"Finished chunk {i}")


if __name__ == "__main__":
    with Pool(n_jobs) as p:
        p.map(run_chunk, list(enumerate(chunks)))
    subprocess.run(["dx", "upload", "-r", "vbm_norm_ukb/"], check=True)