from pathlib import Path
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn.glm.second_level import non_parametric_inference, SecondLevelModel

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_dir = tinception_dir / "VBM" / "struc"
out_dir = tinception_dir / "vbm_results" / "no_TFCE"
fname_covars = tinception_dir / "VBM_design" / "covars.csv"
df = pd.read_csv(fname_covars)
df = df[["tinception_id", "age", "sex", "PTA", "site", "TIV", "group"]].dropna().copy()
df["img_path"] = df["tinception_id"].astype(str).apply(lambda s: vbm_dir / f"{s}_struc_GM_to_template_GM_mod.nii.gz")
df = df[df["img_path"].apply(lambda p: p.exists())].copy()
df = df.sort_values("tinception_id").reset_index(drop=True)


df["group_num"] = df["group"].map({"CO": 0, "TI": 1})
df["sex_num"] = df["sex"]

for col in ["age", "PTA", "TIV"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
    df[col] = df[col] - df[col].mean()

site_dummies = pd.get_dummies(df["site"], prefix="site", drop_first=True)

design_matrix = pd.concat(
    [
        pd.Series(1.0, index=df.index, name="intercept"),
        df[["group_num", "age", "sex_num", "PTA", "TIV"]],
        site_dummies,
    ],
    axis=1,
)
design_matrix = design_matrix.astype(float)
valid = design_matrix.notnull().all(axis=1)
df = df.loc[valid].reset_index(drop=True)
design_matrix = design_matrix.loc[valid].reset_index(drop=True)
second_level_input = df["img_path"].tolist()
mask_img = tinception_dir / "vbm_norm" / "GM_mask.nii.gz"

## Permutation helper without TFCE
def run_perm(name, contrast):
    logp_img = non_parametric_inference(
        second_level_input=second_level_input,
        design_matrix=design_matrix,
        second_level_contrast=contrast,
        mask=mask_img,
        model_intercept=False,   
        n_perm=5000,
        two_sided_test=False,    
        tfce=False,
        threshold=None,          
        smoothing_fwhm=None,
        n_jobs=1,
        random_state=42,
        verbose=1,
    )

    # returned image = -log10(FWER-corrected voxelwise p-values)
    nib.save(logp_img, out_dir / f"{name}_logp_max_t.nii.gz")
    return logp_img

perm_ti_gt_co_logp = run_perm("TI_gt_CO", "group_num")
perm_co_gt_ti_logp = run_perm("CO_gt_TI", "-group_num")

## Standard GLM for t maps and effect maps
print("working on second level ...")
slm = SecondLevelModel(
    mask_img=mask_img,
    smoothing_fwhm=None,
    )
slm = slm.fit(second_level_input, design_matrix=design_matrix)

# Effect maps
effect_ti_gt_co = slm.compute_contrast("group_num", output_type="effect_size")
effect_co_gt_ti = slm.compute_contrast("-group_num", output_type="effect_size")

# T maps
t_ti_gt_co = slm.compute_contrast("group_num", output_type="stat")
t_co_gt_ti = slm.compute_contrast("-group_num", output_type="stat")

nib.save(effect_ti_gt_co, out_dir / "TI_gt_CO_effect_size.nii.gz")
nib.save(effect_co_gt_ti, out_dir / "CO_gt_TI_effect_size.nii.gz")
nib.save(t_ti_gt_co, out_dir / "TI_gt_CO_t.nii.gz")
nib.save(t_co_gt_ti, out_dir / "CO_gt_TI_t.nii.gz")