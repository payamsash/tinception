from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
from brainspace.gradient import GradientMaps

## paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
ssa_dir = tinception_dir / "ssa_data"
ssa_results_dir = tinception_dir / "ssa_results"
n_parcels = 1000
n_network = 7
approach = "dm"
n_components=10
df_master = pd.read_csv(tinception_dir / "vbm_norm" / "covars.csv")

## keep common subjects
subjects_with_ssa = []
for fname in sorted(ssa_dir.iterdir()):
    name = fname.name
    if (not name.startswith(".")) and (name.endswith(".csv")):
        subject = name.split("_Schaefer")[0]
        subjects_with_ssa.append(subject)

df_master = df_master.query('subject_ID == @subjects_with_ssa').copy()

df_master.reset_index(drop=True, inplace=True)
subject_ids = df_master['subject_ID'].values
co_idxs = df_master[df_master["group"] == "CO"].index.to_list()
ti_idxs = df_master[df_master["group"] == "TI"].index.to_list()



ref_subject = "0539"
ref_fname = ssa_dir / f"{ref_subject}_Schaefer2018_{n_parcels}Parcels_{n_network}Networks_order.csv"
df_ref = pd.read_csv(ref_fname)
ref_parcel_names = df_ref.iloc[:, 0].tolist()

connectomes = []
for subject in tqdm(df_master['subject_ID']):

    fname = ssa_dir / f"{subject}_Schaefer2018_{n_parcels}Parcels_{n_network}Networks_order.csv"
    df = pd.read_csv(fname, index_col=0)
    
    # Check if any parcels are missing before fixing
    missing = set(ref_parcel_names) - set(df.index)
    if missing:
        print(f"Subject {subject} missing {len(missing)} parcels. Patching...")

    # reindex to match rows and columns
    df_fixed = df.reindex(index=ref_parcel_names, columns=ref_parcel_names, fill_value=0)
    connectome = df_fixed.to_numpy()
    np.fill_diagonal(connectome, 1)
    connectomes.append(connectome)

reference = np.array(connectomes)[co_idxs].mean(axis=0)
## compute and align gradients
print("Computing gradients will take a while, be patient ...")
gm = GradientMaps(
                    n_components=n_components,
                    approach=approach,
                    kernel=None,
                    alignment='procrustes',
                    random_state=0
                    )
gm.fit(reference)
gm.fit(connectomes, reference=gm.aligned_)

## compute mean connectome per group
avg_connectome_co = np.array(connectomes)[co_idxs].mean(axis=0)
avg_connectome_ti = np.array(connectomes)[ti_idxs].mean(axis=0)

## save the avg conns
np.save(ssa_results_dir / f"avg_connectome_co.npy", avg_connectome_co)
np.save(ssa_results_dir / f"avg_connectome_ti.npy", avg_connectome_ti)

## prepare dfs grad and lambda
grads = gm.aligned_
lambdas = gm.lambdas_

data_swapped = np.transpose(grads, (0, 2, 1))
data_reshaped = data_swapped.reshape(-1, 1000)

components = np.arange(1, 11)
index = pd.MultiIndex.from_product([subject_ids, components], 
                                names=['subject_ID', 'component'])

df_grad = pd.DataFrame(data_reshaped, index=index, columns=ref_parcel_names)
df_grad.reset_index(inplace=True)

metadata_cols = ["tinception_id", "group", "sex", "age", "site", "PTA", "THI", "TIV"]
df_grad = df_grad.merge(
                        df_master[["subject_ID"] + metadata_cols],
                        on="subject_ID",
                        how="inner"
                        )
df_grad.drop(columns="subject_ID", inplace=True)
parcel_cols = [col for col in df_grad.columns if col not in metadata_cols + ["component"]]
df_grad = df_grad[metadata_cols + ["component"] + parcel_cols]
df_grad.rename(columns={"tinception_id": "subjects"}, inplace=True)

flat_values = lambdas.flatten()
repeated_subjects = np.repeat(subject_ids, 10)
tiled_components = np.tile(np.arange(1, 11), len(subject_ids))
df_lambda = pd.DataFrame({
                    'subject_ID': repeated_subjects,
                    'component': tiled_components,
                    'value': flat_values
                })
df_lambda = df_lambda.merge(df_master[["subject_ID", "group", "tinception_id"]],
                        on="subject_ID",
                        how="inner"
                        )
df_lambda.drop(columns="subject_ID", inplace=True)
df_lambda.rename(columns={"tinception_id": "subjects"}, inplace=True)

df_grad.to_csv(ssa_results_dir / "grads.csv", index=False)
df_lambda.to_csv(ssa_results_dir / "lambdas.csv", index=False)
