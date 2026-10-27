### SBM quality
import pandas as pd
from functools import reduce
import numpy as np
import subprocess
from pathlib import Path
from tqdm import tqdm
import nibabel as nib
from nilearn.plotting import plot_stat_map

from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from nilearn import plotting, datasets
import matplotlib.pyplot as plt


for i in tqdm(range(8)):
    fname = f"sbm_norm_ukb/norm_model_chunk_{i}/results/statistics_train_chunk_{i}.csv"
    subprocess.run(["dx", "download", fname], check=True)
    

dfs = []
for i in tqdm(range(8)):
    fname_train = f"statistics_train_chunk_{i}.csv"
    df_train = pd.read_csv(fname_train)
    dfs.append(df_train)

df = reduce(lambda left, right: pd.merge(left, right, on="statistic"), dfs)

df_t = df.T
df_t.columns = df_t.iloc[0]
df_t = df_t.drop(df_t.index[0])
df_t = df_t.reset_index().rename(columns={'index': 'voxel_id'})


fsaverage = datasets.fetch_surf_fsaverage('fsaverage5')
n_vertices = 10242 # Standard for fsaverage5

# 2. Setup Figure: 3 rows (metrics) x 2 columns (hemispheres)
fig = plt.figure(figsize=(8, 8), facecolor='white')
metrics = ["Rho", "MSLL", "EXPV"]
cmaps = ['magma', 'RdBu_r', 'YlOrRd']

for row_idx, metric in enumerate(metrics):
    # Initialize blank stats for this metric
    lh_stats = np.zeros(n_vertices)
    rh_stats = np.zeros(n_vertices)
    
    # Process metric data
    data_col = pd.to_numeric(df_t[metric], errors='coerce').fillna(0)
    
    # Range & Threshold Logic
    if metric == "Rho":
        v_min, v_max, thr = -0.1, 0.4, 0
    elif metric == "MSLL":
        v_min, v_max, thr = -0.2, 0.2, 0
    elif metric == "EXPV":
        data_col = data_col.clip(lower=0)
        v_min, v_max, thr = 0, 0.1, 0

    # Map the voxel_id (v123_lh) to the correct vertex and hemi
    for val, name in zip(data_col, df_t['voxel_id']):
        try:
            # Parse 'v123_lh' -> index 122
            idx = int(name.replace('v', '').split('_')[0]) - 1
            if idx >= n_vertices: continue
            
            if '_lh' in name:
                lh_stats[idx] = val
            elif '_rh' in name:
                rh_stats[idx] = val
        except ValueError:
            continue

    # Plotting Hemispheres
    hemi_params = [
        (fsaverage.infl_left,  lh_stats, fsaverage.sulc_left,  'left',  1 + (row_idx*2)),
        (fsaverage.infl_right, rh_stats, fsaverage.sulc_right, 'right', 2 + (row_idx*2))
    ]

    for mesh, stats, sulc, side, plot_pos in hemi_params:
        ax = fig.add_subplot(3, 2, plot_pos, projection='3d')
        plotting.plot_surf_stat_map(
            mesh, stats, bg_map=sulc,
            hemi=side, view='lateral',
            colorbar=False,
            #colorbar=(side == 'right'), # Only show colorbar on the right side
            threshold=thr,
            cmap=cmaps[row_idx],
            vmin=v_min, vmax=v_max,
            axes=ax,
            alpha=0.7 # Makes the sulcal folds visible but not distracting
        )
        if side == 'left':
            ax.text2D(-0.1, 0.5, metric, transform=ax.transAxes, 
                    fontsize=16, fontweight='bold', va='center', rotation=90)
        
        ax.dist = 6

plt.subplots_adjust(
    left=0.05,
    right=0.95,
    bottom=0.05,
    top=0.95,
    wspace=-0.35,   # controls lh–rh distance
    hspace=-0.15
)
plt.savefig("sbm_norm_ukb_quality.pdf", dpi=600, bbox_inches='tight')
plt.show()