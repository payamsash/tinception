from pathlib import Path
import numpy as np
from scipy.ndimage import label
import pandas as pd
import seaborn as sns
import nibabel as nib
from nilearn import image
from nilearn.plotting import plot_stat_map
import matplotlib.pyplot as plt

############ VBM stat ############
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
fsl_dir = Path("/Users/payamsadeghishabestari/fsl")
vbm_dir = tinception_dir / "vbm_results"
vbm_norm_dir = tinception_dir / "vbm_norm"
plots_dir = tinception_dir / "plots" / "vbm"
plots_dir.mkdir(exist_ok=True)

img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"
img_fname =  vbm_dir / "fslvbm_s3_tfce_corrp_tstat2.nii.gz" # TI > CO

p_thr = 0.005
thr = 1 - p_thr
'''
kwargs = {
            "colorbar": True,
            "cbar_tick_format": "%.2g",
            "annotate": False,
            "draw_cross": False,
            "threshold": thr,
            "radiological": False,
            "cmap": 'autumn',
            "symmetric_cbar": False,
            "vmin": thr,
            "vmax": 1,
            "dim": -0.3,
            "black_bg": True,
            "cut_coords": 3
        }
for disp_mode in ["x", "y", "z"]:

    fig = plot_stat_map(
                    stat_map_img=img_fname,
                    bg_img=img_bg,
                    display_mode=disp_mode,
                    **kwargs
                    )
    fig.savefig(plots_dir / f"Figure_{disp_mode}.pdf", dpi=600, bbox_inches='tight')

############ VBM norm ############
df = pd.read_csv(vbm_norm_dir / "norm_model_subcortical" / "results" / "Z_test.csv")
img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"

mask_path = vbm_norm_dir / 'GM_mask.nii.gz'
subcortical_mask_path = vbm_norm_dir / 'subcortical_mask_thr80.nii.gz'

mask_img = nib.load(mask_path)
mask_affine = mask_img.affine
mask_shape = mask_img.shape

mask1 = mask_img.get_fdata() > 0
mask2 = nib.load(subcortical_mask_path).get_fdata() > 0
final_mask = mask1 & mask2  

## read deviation scores
subject_ids = df["subject_ids"].values
df.drop(columns=["observations", "subject_ids"], inplace=True)
thr = 1.96
average_dev = df.mean(axis=0)
large_dev = (df > thr).sum(axis=0)
small_dev = (df < -thr).sum(axis=0)
total_dev = large_dev + small_dev
df_extreme = pd.DataFrame({
                "Average_Dev": average_dev,
                "Large_Dev": large_dev,
                "Small_Dev": small_dev,
                "Total_Dev": total_dev
                }).T

## reconstruct
data = df_extreme.iloc[3].values
reconstructed_vol = np.zeros(mask_shape)
reconstructed_vol[final_mask] = data
output_nii = nib.Nifti1Image(reconstructed_vol, mask_affine)
smoothed_nii = image.smooth_img(output_nii, fwhm=7)

## plot
kwargs = {
            "colorbar": True,
            "cbar_tick_format": "%.2g",
            "annotate": False,
            "draw_cross": False,
            "radiological": False,
            "cmap": 'magma',
            "threshold": 12,
            "symmetric_cbar": False,
            "vmin": None,
            "vmax": None,
            "dim": -0.3,
            "black_bg": True,
            "cut_coords": 3
        }

for disp_mode in ["x", "y", "z"]:
    fig = plot_stat_map(
                    stat_map_img=smoothed_nii,
                    bg_img=img_bg,
                    display_mode=disp_mode,
                    **kwargs
                    )
    fig.savefig(plots_dir / f"Figure_norm_{disp_mode}.pdf", dpi=600, bbox_inches='tight')
'''
############ VBM boxplots ############
## getting co and ti idxs
df = pd.read_csv(vbm_norm_dir / "covars.csv")
co_idxs = df[df["group"] == "CO"].index
ti_idxs = df[df["group"] == "TI"].index

print("Creating 3d binary mask ...")
corrp_file = nib.load(vbm_dir / "fslvbm_s3_tfce_corrp_tstat2.nii.gz")
corrp_data = corrp_file.get_fdata()
affine = corrp_file.affine
sig_mask_3d = (corrp_data > 1 -  p_thr).astype(int)

print("Separating clusters ...")
labeled_mask, num_clusters = label(sig_mask_3d)
print(f"Found {num_clusters} separate significant clusters.")

print("Extract subject data from large GM file")
gm_file = nib.load(vbm_norm_dir / "GM_mod_merg_s3.nii.gz")
gm_data = gm_file.get_fdata()

print("Computing t-val and cohen distance ...")
tval_file = nib.load(vbm_dir / "fslvbm_s3_tstat2.nii.gz")
tval_data = tval_file.get_fdata()

peak_ps, peak_ts, cohen_ds = [], [], []
df_plots = []
for cluster_id in range(1, num_clusters + 1):
    print(f"working on cluster {cluster_id} ... ")
    cluster_indices = np.where(labeled_mask == cluster_id)
    num_voxels = len(cluster_indices[0])

    ## transform voxel index to MNI world coordinates
    avg_i = np.mean(cluster_indices[0])
    avg_j = np.mean(cluster_indices[1])
    avg_k = np.mean(cluster_indices[2])
    avg_voxel = np.array([avg_i, avg_j, avg_k, 1])
    mni_coords = (affine @ avg_voxel)[:3]
    print(f"Cluster {cluster_id}: Size = {num_voxels} voxels, Avg MNI = {mni_coords}")


    ## compute tval and cohen d 
    peak_p = np.max(np.abs(corrp_data[cluster_indices]))
    peak_t = np.max(np.abs(tval_data[cluster_indices]))
    cohen_d = peak_t * np.sqrt((1/len(co_idxs)) + (1/len(ti_idxs)))

    cluster_values = gm_data[cluster_indices].mean(axis=0)
    df_plots.append(
        pd.DataFrame({
                    'gm_volume': cluster_values,
                    'group': ["CO" if i in co_idxs else "TI" for i in range(len(cluster_values))],
                    "cluster_id": [cluster_id] * len(cluster_values)
                        }
        ))
    peak_ps.append(peak_p)
    peak_ts.append(peak_t)
    cohen_ds.append(cohen_d)
    
df_plot = pd.concat(df_plots, axis=0)
label_map = {
    1: "Uvula-RH",
    2: "Subcallosal-RH",
    3: "Subcallosal-LH"
}
df_plot['region'] = df_plot['cluster_id'].map(label_map)
df_plot["group"] = df_plot["group"].map({"CO": "Tinnitus", "TI": "Control"}) # (I mistakenly swap 0 and 1 in vbm design matrix)
df_plot.drop(columns="cluster_id", inplace=True)

## plotting 
pal = ['#1f77b4', '#d62728']
order = ["Control", "Tinnitus"]

g = sns.FacetGrid(
                data=df_plot,
                row="region",
                row_order=label_map.values(),
                xlim=[0, 1],
                height=1.8,
                aspect=3
                )
g.map_dataframe(
            sns.stripplot,
            x="gm_volume",
            hue="group",
            palette=pal,
            linewidth=0,
            size=5.5,
            dodge=True,
            edgecolor=None,
            jitter=0.19,
            alpha=0.15,
            hue_order=order
)
g.map_dataframe(
            sns.boxplot,
            x="gm_volume",
            hue="group",
            palette=pal,
            width=0.8,
            dodge=True,
            linewidth=1.8,
            gap=0.3,
            fill=False,
            hue_order=order,
            showfliers=False
)

g.despine(left=True)
g.set(yticks=[])
g.add_legend(loc="upper right", bbox_to_anchor=(0.9, 1))
g.tight_layout()

for ax, pval, tval, cohen in zip(g.axes.flat, peak_ps, peak_ts, cohen_ds):
    full_title = ax.get_title()
    region_name = full_title.split('=')[-1].strip()
    new_title = (f"{region_name}\n t = {tval:.2f} | p = {1 - pval:.4f} | d = {cohen:.2f}")
        
    ax.set_title(new_title, fontstyle='italic', fontsize=8)

g.savefig(
        plots_dir / "boxplots.pdf",
        format="pdf",
        dpi=300,
        bbox_inches="tight"
        )