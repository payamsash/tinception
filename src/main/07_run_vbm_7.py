from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

import nibabel as nib
from nilearn.plotting import plot_stat_map

"""
This script identifies and visualizes statistically significant deviations in normative VBM data.

- Compares control (CO) and tinnitus (TI) groups for extreme voxel-wise deviations (>1.96 SD).
- Performs chi-square tests per voxel and corrects p-values using FDR (Benjamini-Hochberg).
- Generates a binary significance map and saves it as a NIfTI file.
- Plots significant voxels over an MNI152 template for visualization.
"""


## read the files
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_norm_dir = tinception_dir / "vbm_norm"
plots_dir = tinception_dir / "plots"
fsl_dir = Path("/Users/payamsadeghishabestari/fsl")

df_test = pd.read_csv(vbm_norm_dir / "results" / "Z_test.csv")
df_train = pd.read_csv(vbm_norm_dir / "results" / "Z_train.csv")

v_cols = [c for c in df_test.columns if c.startswith('v')]
other_cols = [c for c in df_test.columns if not c.startswith('v')]
v_cols_sorted = sorted(v_cols, key=lambda x: int(x[1:]))

# Re-index the dataframe
df_test = df_test[other_cols + v_cols_sorted]
df_train = df_train[other_cols + v_cols_sorted]

df_test['group'] = np.where(df_test['subject_ids'].isin(df_train['subject_ids']), 'CO', 'TI')
## bring grop to beginning
cols = df_test.columns.tolist()
cols.insert(2, cols.pop())
df_test = df_test[cols]
df_test.drop(columns=["observations", "subject_ids"], inplace=True)

## find significant extremes
v_cols = [c for c in df_test.columns if c.startswith('v')]

hits_mask = (df_test[v_cols].abs() > 1.96).astype(int)
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
significant_voxels = analysis_df[analysis_df['is_significant'] == True]
print(f"Found {len(significant_voxels)} significant voxels after FDR correction.")


### plotting
img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"
mask_img = nib.load(vbm_norm_dir / 'GM_mask.nii.gz')
mask1 = nib.load(vbm_norm_dir / 'GM_mask.nii.gz').get_fdata() > 0
mask2 = nib.load(vbm_norm_dir / 'subcortical_mask_thr80.nii.gz').get_fdata() > 0
mask = mask1 & mask2

mask_data = mask.astype(float) 
affine = mask_img.affine
sig_map_data = np.zeros(mask_data.shape)
significant_indices = [int(v.replace('v', '')) - 1 for v in significant_voxels['voxel']]

mask_1d_indices = np.where(mask.ravel())[0]
sig_global_indices = mask_1d_indices[significant_indices]
sig_map_data.ravel()[sig_global_indices] = 1
sig_nifti = nib.Nifti1Image(sig_map_data, affine)
mask_file = vbm_norm_dir / f"individuals_sig_mask.nii.gz"
nib.save(sig_nifti, mask_file)

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
fig.savefig(plots_dir / "vbm" / f"vbm_norm_brain.pdf", dpi=600, bbox_inches='tight')