from pathlib import Path
import pandas as pd
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel
)
import nibabel as nib
import re



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



###################

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
df = pd.read_csv(tinception_dir / "subcortical_roi" / "master.csv")


atlas_name = "Amygdala"
df = df.query(f'structure == "{atlas_name}"')
if not atlas_name == "Thalamic-nuclei":
    df["region"] = df["region"] + "-" + df["hemi"]
df.drop(columns="hemi", inplace=True)
index_cols = list(df.columns[2:]) # except region and volume

## long to wide
df_wide = df.pivot(index=index_cols, columns='region', values='volume')
df_wide = df_wide.reset_index()
df_wide.columns.name = None

## define covariates
demographic_cols = ["age", "sex", "PTA", "TIV"]
start_idx = 6 + len(demographic_cols)
response_cols = df_wide.columns[start_idx:].tolist()

## run norm model
kwargs = {
            "covariates": demographic_cols,
            "batch_effects": ["site"],
            "response_vars": response_cols, 
            "subject_ids": "tinception_id"
            }
saving_dir = tinception_dir / "subcortical_roi" / "norm_models" / atlas_name
saving_dir.mkdir(exist_ok=True)

## create norm objects
norm_train_all = NormData.from_dataframe(
                                        name="train",
                                        dataframe=df_wide.query('group == "CO"'),
                                        **kwargs
                                        )
norm_test_all = NormData.from_dataframe(
                                        name="test",
                                        dataframe=df_wide,
                                        **kwargs
                                        )

## define models
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
                save_dir=str(saving_dir),
                inscaler="standardize",
                outscaler="none",
                )
model.fit_predict(norm_train_all, norm_test_all)



####################
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
df_main = pd.read_csv(tinception_dir / "subcortical_roi" / "master.csv")

atlases_ukb = ["amygdalar_nuclei", "hippo_subfields", "thalamic_nuclei"]
atlases_main = ["Amygdala", "Hippo", "Thalamic-nuclei"]

idx = 2
if idx == 0:
    end_idx = -2
elif idx in [1, 2]:
    end_idx = None

atlas_main = atlases_main[idx]
atlas_ukb = atlases_ukb[idx]

ukb_norm_model_fname = tinception_dir / "subcortical_roi" / "norm_models_ukb" / atlas_ukb
model = NormativeModel.load(str(ukb_norm_model_fname))


## prepare df_main
mask = df_main["hemi"].notna()
df_main.loc[mask, "region"] = df_main["region"] + "-" + df_main["hemi"]
df_main.drop(columns=["hemi", "THI", "subject_ID"], inplace=True)
df_main = df_main.query(f'structure == "{atlas_main}"')

df_wide = df_main.pivot(
    index=['tinception_id', 'age', 'sex', 'site', 'PTA', 'TIV', 'group'], 
    columns='region', 
    values='volume'
).reset_index()
df_wide.columns.name = None


df_wide.columns = [col.replace("-", "_") for col in df_wide.columns]
def clean_region_name(col):
    metadata = ['tinception_id', 'age', 'sex', 'site', 'PTA', 'TIV', 'group']
    if col in metadata:
        return col
    
    new_name = col
    if col.endswith("_lh"):
        new_name = col.replace("_lh", "_left_hemisphere")
    elif col.endswith("_rh"):
        new_name = col.replace("_rh", "_right_hemisphere")
        
    return f"Volume_of_{new_name}"

df_wide.columns = [clean_region_name(c) for c in df_wide.columns]
df_wide.rename(columns={"site": "SITE", "PTA": "srt", "TIV": "eTIV"}, inplace=True)

## just a place holder (will delete later)
if idx == 0:
    df_wide["Volume_of_Whole_amygdala_left_hemisphere"] = df_wide['Volume_of_Accessory_Basal_nucleus_left_hemisphere']
    df_wide["Volume_of_Whole_amygdala_right_hemisphere"] = df_wide['Volume_of_Accessory_Basal_nucleus_left_hemisphere']
if idx == 1:
    df_wide.drop(columns=["Volume_of_Whole_hippocampus_left_hemisphere", 
                            "Volume_of_Whole_hippocampus_right_hemisphere"], inplace=True)

if idx == 2:
    def rename_thalamus_cols(col):
        pattern = r"^Volume_of_(Left|Right)_(.+)$"
        match = re.match(pattern, col)
        
        if match:
            side = match.group(1).lower()  # left / right
            region = match.group(2)
            
            # Remove parentheses like MV(Re) -> MVRe
            region = region.replace("(", "").replace(")", "")
            
            return f"Volume_of_{region}_{side}_hemisphere"
        
        return col  # keep other columns unchanged

    df_wide = df_wide.rename(columns=rename_thalamus_cols)
    df_wide["Volume_of_Whole_thalamus"] = df_wide["Volume_of_Whole_thalamus_left_hemisphere"]
    df_wide.drop(columns=["Volume_of_Whole_thalamus_left_hemisphere", 
                            "Volume_of_Whole_thalamus_right_hemisphere"], inplace=True)
    
## create norm class for main dataset
kwargs = {
            "covariates": ["age", "sex", "srt", "eTIV"],
            # "batch_effects": ["SITE"],
            "response_vars": list(df_wide.columns[7:end_idx]), 
            "subject_ids": "tinception_id"
            }
norm_data_main_co = NormData.from_dataframe(
                                        name="main",
                                        dataframe=df_wide.query('group == "CO"'),
                                        **kwargs
                                        )
norm_data_main_all = NormData.from_dataframe(
                                        name="main",
                                        dataframe=df_wide,
                                        **kwargs
                                        )
model.save_dir = str(tinception_dir / "subcortical_roi" / "new_norm_models" / atlas_main)
extended_model = model.extend_predict(
                                    extend_data=norm_data_main_co,
                                    predict_data=norm_data_main_all,
                                    save_dir=str(tinception_dir / "subcortical_roi" / "new_norm_models" / atlas_main)
                                    )