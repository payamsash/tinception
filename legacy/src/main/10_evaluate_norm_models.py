from pathlib import Path
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn import plotting, datasets
import matplotlib.pyplot as plt


## assign paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
fsl_dir = Path("/Users/payamsadeghishabestari/fsl")
vbm_norm_dir = tinception_dir / "vbm_norm"
sbm_norm_dir = tinception_dir / "sbm_norm"

############# VBM
## read mask images
img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"
mask_img = nib.load(vbm_norm_dir / 'GM_mask.nii.gz')
mask1 = nib.load(vbm_norm_dir / 'GM_mask.nii.gz').get_fdata() > 0
mask2 = nib.load(vbm_norm_dir / 'subcortical_mask_thr80.nii.gz').get_fdata() > 0
mask = mask1 & mask2

## read stat norm
df = pd.read_csv(vbm_norm_dir / "results" / "statistics_train.csv")
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
    display = plotting.plot_glass_brain(
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
output_file = tinception_dir / "plots" / "VBM_Normative_Metrics_Final.pdf"
plt.savefig(output_file, dpi=600, bbox_inches='tight')
plt.show()

############## SBM

df = pd.read_csv(sbm_norm_dir / "results" / "statistics_train.csv")
df_t = df.T
df_t.columns = df_t.iloc[0]
df_t = df_t.drop(df_t.index[0])
df_t = df_t.reset_index().rename(columns={'index': 'voxel_id'})

# 1. Setup Surface Data
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
        v_min, v_max, thr = 0.2, 0.7, 0.1
    elif metric == "MSLL":
        v_min, v_max, thr = -0.5, 0.5, 0.01
    elif metric == "EXPV":
        data_col = data_col.clip(lower=0)
        v_min, v_max, thr = 0.05, 0.4, 0.01

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
            # colorbar=(side == 'right'), # Only show colorbar on the right side
            threshold=thr,
            cmap=cmaps[row_idx],
            vmin=v_min, vmax=v_max,
            axes=ax,
            darkness=0.5 # Makes the sulcal folds visible but not distracting
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
plt.savefig(tinception_dir / "plots" / "SBM_Normative_Metrics.pdf", dpi=600, bbox_inches='tight')
plt.show()