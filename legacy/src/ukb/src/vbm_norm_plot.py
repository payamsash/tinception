import pandas as pd
from functools import reduce
import numpy as np
import subprocess
from pathlib import Path
from tqdm import tqdm
import nibabel as nib
from nilearn.plotting import plot_stat_map, plot_glass_brain
import matplotlib.pyplot as plt

from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def add_grp_col(df_test, df_train):
    df_test['group'] = np.where(df_test['subject_ids'].isin(df_train['subject_ids']), 'CO', 'TI')
    cols = df_test.columns.tolist()
    cols.insert(1, cols.pop())
    df_test = df_test[cols]
    return df_test

for i in tqdm(range(72)):
    if i == 0:
        fname = f"vbm_norm_ukb/norm_model_chunk_{i}/results/Z_train_chunk_{i}.csv"
        subprocess.run(["dx", "download", fname], check=True)

    fname = f"vbm_norm_ukb/norm_model_chunk_{i}/results/Z_test_chunk_{i}.csv"
    subprocess.run(["dx", "download", fname], check=True)
    
    fname = f"vbm_norm_ukb/norm_model_chunk_{i}/results/statistics_train_chunk_{i}.csv"
    subprocess.run(["dx", "download", fname], check=True)

fname_train = f"Z_train_chunk_0.csv"
df_train = pd.read_csv(fname_train)

dfs = []
for i in tqdm(range(72)):
    
    fname_test = f"Z_test_chunk_{i}.csv"
    df_test = pd.read_csv(fname_test)
    df_test.drop(columns=["observations"], inplace=True)
    dfs.append(df_test)

merged_df = reduce(lambda left, right: pd.merge(left, right, on="subject_ids"), dfs)
df_test = add_grp_col(merged_df, df_train)

## sort
fixed_cols = df_test.columns[:2]
sorted_cols = sorted(
    df_test.columns[2:],
    key=lambda x: int(x[1:])  # remove 'v' and convert to int
)
df_test = df_test[list(fixed_cols) + sorted_cols]


v_cols = [c for c in df_test.columns if c.startswith('v')]

hits_mask = (df_test[v_cols] > 1.96).astype(int)
hits_mask['group'] = df_test['group']

counts = hits_mask.groupby('group').sum() 
totals = df_test.groupby('group')[v_cols].count()

co_hits_arr = counts.loc['CO'].values
ti_hits_arr = counts.loc['TI'].values
co_total = totals.loc['CO'].values
ti_total = totals.loc['TI'].values

results_list = []

for i, col in enumerate(v_cols):
    co_h = co_hits_arr[i]
    ti_h = ti_hits_arr[i]
    co_f = co_total[i] - co_h
    ti_f = ti_total[i] - ti_h
    
    # Contingency Table: [[Hits_CO, Fails_CO], [Hits_TI, Fails_TI]]
    table = [[co_h, co_f], [ti_h, ti_f]]
    
    # chi2_contingency returns: chi2, p-value, dof, expected
    # Using correction=True (Yates' correction) is safer for small counts
    _, p, _, _ = chi2_contingency(table, correction=True)
    
    results_list.append({
        'voxel': col,
        'CO_hits': co_h,
        'TI_hits': ti_h,
        'CO_rate': co_h / co_total[i],
        'TI_rate': ti_h / ti_total[i],
        'p_uncorrected': p
    })

analysis_df = pd.DataFrame(results_list)
reject, pvals_corrected, _, _ = multipletests(
    analysis_df['p_uncorrected'], 
    alpha=0.05, 
    method='fdr_bh'
)

analysis_df['p_fdr'] = pvals_corrected
analysis_df['is_significant'] = reject
analysis_df = analysis_df.sort_values('p_fdr')
# significant_voxels = analysis_df[analysis_df['is_significant'] == True]

significant_voxels = analysis_df.query('p_uncorrected < 0.05')
print(f"Found {len(significant_voxels)} significant voxels after FDR correction.")


############# plot it on the img
img_bg = "MNI152_T1_0.5mm.nii.gz"
#subprocess.run(["dx", "download", img_bg], check=True)

mask_path = "HO_sub_thr80_mask_2mm.nii.gz"
#subprocess.run(["dx", "download", mask_path], check=True)
mask_img = nib.load(mask_path)
mask = nib.load(mask_path).get_fdata() > 0
affine = mask_img.affine

mask_data = mask.astype(float)
sig_map_data = np.zeros(mask_data.shape)
significant_indices = [int(v.replace('v', '')) - 1 for v in significant_voxels['voxel']]

mask_1d_indices = np.where(mask.ravel())[0]
sig_global_indices = mask_1d_indices[significant_indices]
sig_map_data.ravel()[sig_global_indices] = 1
sig_nifti = nib.Nifti1Image(sig_map_data, affine)

kwargs = {
            "colorbar": False,
            "cbar_tick_format": "%.2g",
            "annotate": False,
            "draw_cross": False,
            "radiological": False,
            "cmap": 'Reds',
            "threshold": 0,
            "symmetric_cbar": False,
            "vmin": None,
            "vmax": None,
            "dim": -0.3,
            "black_bg": True,
            "cut_coords": (32, -28, -26)
        }
fig = plot_stat_map(
            stat_map_img=sig_nifti,
            bg_img=img_bg,
            display_mode="ortho",
            **kwargs
            )
fig.savefig(f"vbm_norm_ukb_0.05.pdf", dpi=600, bbox_inches='tight')






###############
dfs = []
for i in tqdm(range(72)):
    
    fname_train = f"statistics_train_chunk_{i}.csv"
    df_train = pd.read_csv(fname_train)
    dfs.append(df_train)

df = reduce(lambda left, right: pd.merge(left, right, on="statistic"), dfs)

df_t = df.T
df_t.columns = df_t.iloc[0]
df_t = df_t.drop(df_t.index[0])
df_t = df_t.reset_index().rename(columns={'index': 'voxel_id'})


## plot
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), facecolor='white')

metrics = ["Rho", "MSLL", "EXPV"]
cmaps = ['magma', 'RdBu_r', 'YlOrRd'] 

for i, metric in enumerate(metrics):
    # Data preparation
    data_col = pd.to_numeric(df_t[metric], errors='coerce').fillna(0)
    
    # Range and processing logic
    if metric == "Rho":
        v_values = data_col.values
        v_min, v_max, thr = 0.2, 0.7, 0.1
    elif metric == "MSLL":
        v_values = data_col.values
        # Symmetric range around 0 makes MSLL "Good vs Bad" fit clear
        v_min, v_max, thr = -0.5, 0.5, 0.01 
    elif metric == "EXPV":
        v_values = data_col.clip(lower=0).values
        v_min, v_max, thr = 0.05, 0.5, 0.01

    # Map to 3D volume
    stat_map_data = np.zeros(mask.shape)
    stat_map_data[mask] = v_values
    img = nib.Nifti1Image(stat_map_data, mask_img.affine)
    
    # Plotting
    display = plot_glass_brain(
        stat_map_img=img,
        display_mode="lyrz",
        colorbar=True,
        annotate=True,
        threshold=thr,         
        cmap=cmaps[i],
        vmin=v_min,
        vmax=v_max,
        plot_abs=False,        
        axes=axes[i],
        black_bg=False         
    )
    
    # Custom Labels instead of Titles (Paper Style)
    axes[i].text(-0.05, 0.5, metric, transform=axes[i].transAxes, 
                fontsize=16, fontweight='bold', va='center', ha='right')

# Final layout adjustments
plt.subplots_adjust(hspace=0.05, left=0.15)

# Save as Vector PDF for high quality
# output_file = tinception_dir / "plots" / "VBM_Normative_Metrics_Final.pdf"
plt.savefig("vbm_norm_ukb_quality.pdf", dpi=600, bbox_inches='tight')
plt.show()

