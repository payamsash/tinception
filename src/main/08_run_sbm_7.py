from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from nilearn.plotting import plot_surf_stat_map
from nilearn import datasets

"""
SBM Normative Vertex-Wise Analysis and Visualization

This script performs vertex-wise analysis of surface-based morphometry (SBM) data 
using normative modeling. It identifies vertices with extreme deviations in 
tinnitus (TI) versus control (CO) participants, computes per-vertex chi-square 
statistics, applies FDR correction, and visualizes significant vertices on 
the fsaverage5 surface.

Steps:
1. Load normative Z-scores for training (CO) and test (TI + CO) subjects.
2. Identify vertices with extreme Z-scores (|Z| > 1.96) and compute hit rates.
3. Perform chi-square tests per vertex to compare TI vs CO extremes.
4. Apply FDR correction for multiple comparisons.
5. Map significant vertices to fsaverage5 surface meshes (LH and RH).
6. Visualize lateral surface maps with TI hit rates using a hot colormap.
7. Save the figure as sbm_norm_brain.pdf.

Inputs:
- Z_train.csv, Z_test.csv : vertex-wise normative Z-scores
- fsaverage5 surface meshes from nilearn.datasets

Outputs:
- Analysis dataframe with per-vertex statistics (optional)
- PDF figure of significant vertices on LH/RH surfaces
"""

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
sbm_norm_dir = tinception_dir / "sbm_norm"
plots_dir = tinception_dir / "plots"

df_test = pd.read_csv(sbm_norm_dir / "results" / "Z_test.csv")
df_train = pd.read_csv(sbm_norm_dir / "results" / "Z_train.csv")

v_lh_cols = [c for c in df_test.columns if c.startswith('v') and c.endswith("_lh")]
v_rh_cols = [c for c in df_test.columns if c.startswith('v') and c.endswith("_rh")]
other_cols = [c for c in df_test.columns if not c.startswith('v')]
v_lh_cols_sorted = sorted(v_lh_cols, key=lambda x: int(x.replace('v', '').replace('_lh', '')))
v_rh_cols_sorted = sorted(v_rh_cols, key=lambda x: int(x.replace('v', '').replace('_rh', '')))

# Re-index the dataframe
df_test = df_test[other_cols + v_lh_cols_sorted + v_rh_cols_sorted]
df_train = df_train[other_cols + v_lh_cols_sorted + v_rh_cols_sorted]

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
    
    if co_h == 0 and ti_h == 0:
        p = 1.0
    else:
        co_f = co_total[i] - co_h
        ti_f = ti_total[i] - ti_h
        table = [[co_h, co_f], [ti_h, ti_f]]
        try:
            _, p, _, _ = chi2_contingency(table, correction=True)
        except ValueError:
            p = 1.0
    
    results_list.append({
        'vertex': col,
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

raw_001 = analysis_df[analysis_df['p_uncorrected'] < 0.001]
print(f"Voxels with p < 0.001 (uncorrected): {len(raw_001)}")


## plotting
fsaverage = datasets.fetch_surf_fsaverage('fsaverage5')
n_vertices = 10242 
lh_stats = np.zeros(n_vertices)
rh_stats = np.zeros(n_vertices)

for _, row in raw_001.iterrows():
    name = row['vertex']
    val = row['TI_rate'] 
    idx = int(name.replace('v', '').split('_')[0]) - 1
    
    if '_lh' in name and idx < n_vertices:
        lh_stats[idx] = val
    elif '_rh' in name and idx < n_vertices:
        rh_stats[idx] = val


hemi_params = [
    (fsaverage.infl_left,  lh_stats, fsaverage.sulc_left,  'left'),
    (fsaverage.infl_right, rh_stats, fsaverage.sulc_right, 'right')
]

fig, axes = plt.subplots(1, 2, figsize=(10, 5), subplot_kw={'projection': '3d'}, layout="tight")

for (mesh, stats, sulc, side), ax in zip(hemi_params, axes):
    display = plot_surf_stat_map(
            mesh,
            stats,
            sulc,
            hemi=side,
            view='lateral',
            colorbar=False,
            threshold=0.001,
            cmap='hot',
            vmin=0, vmax=None,
            axes=ax
            )
fig.savefig(plots_dir / "sbm" / f"sbm_norm_brain.pdf", dpi=600, bbox_inches='tight')