from pathlib import Path
import pandas as pd
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel
)
import nibabel as nib

"""
This script performs normative modeling of cortical thickness (SBM) data.

- Loads FreeSurfer surface thickness maps for LH and RH hemispheres per subject.
- Merges thickness data with covariates (age, sex, PTA, TIV, site, group).
- Constructs PCNToolkit BLR-based normative models using control subjects (CO) as training.
- Predicts deviations for tinnitus subjects (TI) across the cortical surface.
- Saves model results to a specified SBM normative model directory.
"""

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
subjects_fs_dir = tinception_dir / "subjects_fs_dir"
vbm_norm_dir = tinception_dir / "vbm_norm"
sbm_norm_dir = tinception_dir / "sbm_norm"

df = pd.read_csv(vbm_norm_dir / "covars.csv")
df.sort_values(by="tinception_id", inplace=True)
df.drop(columns=["THI", "distance", "weights", "subclass"], inplace=True)
train_indices = df[df['group'] == "CO"].index.tolist()
test_indices = df[df['group'] == "TI"].index.tolist()
df['site'] = pd.factorize(df['site'])[0]

cols_lh = [f'v{i}_lh' for i in range(1, 10242 + 1)]
cols_rh = [f'v{i}_rh' for i in range(1, 10242 + 1)]
y_dfss = []

for subject in df['subject_ID']:
    surf_dir = subjects_fs_dir / subject / "surf"
    
    y_dfs = []
    for hemi, cols in zip(["lh", "rh"], [cols_lh, cols_rh]):
        fname = surf_dir / f"{hemi}.thickness.fwhm10.fsaverage5.mgh"
        img = nib.load(fname)
        data = img.get_fdata().squeeze().reshape(1, -1)
        y_df = pd.DataFrame(data, columns=cols)
        y_dfs.append(y_df)
    
    y_dfs = pd.concat(y_dfs, axis=1)
    
    y_dfs["subject_ID"] = subject
    y_dfss.append(y_dfs)

y_df = pd.concat(y_dfss, axis=0)
y_df = y_df.merge(
                df,
                on="subject_ID",
                how="inner"
)

cols = y_df.columns.tolist()
last8 = cols[-8:]
rest = cols[:-8]
new_order = last8 + rest
df = y_df[new_order]

kwargs = {
        "covariates": ['age', 'sex', 'PTA', 'TIV'],
        "batch_effects": ["site"],
        "response_vars": list(df.columns[8:]), 
        "subject_ids": "tinception_id"
        }
norm_train_all = NormData.from_dataframe(
                                        name="train",
                                        dataframe=df[df['group'] == "CO"],
                                        **kwargs
                                        )
norm_test_tinnitus = NormData.from_dataframe(
                                        name="test",
                                        dataframe=df,
                                        **kwargs
                                        )

template_blr = BLR(
            name="payam_blr",
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
            save_dir=str(sbm_norm_dir / "norm_model"),
            inscaler="standardize",
            outscaler="none",
        )

model.fit_predict(norm_train_all, norm_test_tinnitus)


