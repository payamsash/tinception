from pathlib import Path
import pandas as pd
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel
)
import nibabel as nib

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_norm_dir = tinception_dir / "vbm_norm"
df = pd.read_csv(vbm_norm_dir / "covars.csv")
df.sort_values(by="tinception_id", inplace=True)
df.drop(columns=["subject_ID", "THI", "distance", "weights", "subclass"], inplace=True)

train_indices = df[df['group'] == "CO"].index.tolist()
test_indices = df[df['group'] == "TI"].index.tolist()

covariates = ['age', 'sex', 'site', 'TIV']
df['site'] = pd.factorize(df['site'])[0]

# Load the 4D brain data and the mask
data = nib.load(vbm_norm_dir / 'GM_mod_merg_s3.nii.gz').get_fdata()
mask1 = nib.load(vbm_norm_dir / 'GM_mask.nii.gz').get_fdata() > 0
mask2 = nib.load(vbm_norm_dir / 'subcortical_mask_thr80.nii.gz').get_fdata() > 0
mask = mask1 & mask2
y = data[mask].T

col_names = [f'v{i}' for i in range(1, y.shape[1] + 1)]
y_df = pd.DataFrame(y, columns=col_names)
df = pd.concat([df, y_df], axis=1)

kwargs = {
        "covariates": ['age', 'sex', 'TIV'],
        "batch_effects": ["site"],
        "response_vars": col_names, 
        "subject_ids": "tinception_id"
        }
norm_train_all = NormData.from_dataframe(
                                        name="train",
                                        dataframe=df[df['group'] == "CO"],
                                        **kwargs
                                        )
norm_test_tinnitus = NormData.from_dataframe(
                                        name="test",
                                        dataframe=df[df['group'] == "TI"],
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
            savemodel=True,
            evaluate_model=True,
            saveresults=True,
            saveplots=False,
            save_dir=str(vbm_norm_dir / "norm_model_subcortical"),
            inscaler="standardize",
            outscaler="none",
        )

model.fit_predict(norm_train_all, norm_test_tinnitus)