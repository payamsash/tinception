from pathlib import Path

import pandas as pd
import numpy as np
from scipy import ndimage
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from statsmodels.stats.multitest import multipletests
import pingouin as pg

import nibabel as nib

## path to masks
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
masks_dir = tinception_dir / "vbm_results" / "cluster_masks"
putamen_lh_mask = masks_dir / "cluster_putamen-lh_mask.nii.gz"
putamen_rh_mask = masks_dir / "cluster_putamen-rh_mask.nii.gz"
subcollasal_mask = masks_dir / "cluster_subcollasal_mask.nii.gz"

# Load image
img = nib.load(tinception_dir / "vbm_results" / "HF_TFCE" / "CO_gt_TI_main_HF_logp_max_tfce.nii.gz")
data = img.get_fdata()

thr = 1.3
binary = data > thr
labeled, n_clusters = ndimage.label(binary)
print(f"Found {n_clusters} clusters")

cluster_sizes = ndimage.sum(binary, labeled, index=range(1, n_clusters + 1))
largest_label = np.argmax(cluster_sizes) + 1 
largest_cluster = labeled == 8 # subcollasal

out_img = nib.Nifti1Image(largest_cluster.astype(np.uint8), img.affine, img.header)
nib.save(out_img, subcollasal_mask)

## load images and dfs
df_main = pd.read_csv(tinception_dir / "VBM_design" / "covars.csv")
df_regions = pd.read_csv(tinception_dir / "VBM_design" / "gm.csv")

subset_idx = df_main.index.to_list()
gms_fname = tinception_dir / "vbm_norm" / "GM_mod_merg_s3.nii.gz"
img = nib.load(gms_fname)
data = img.get_fdata()
n_subs = data.shape[3]


voxel_vol = np.prod(img.header.get_zooms()[:3])
for mask_f, col_name in zip(
    [putamen_lh_mask, putamen_rh_mask, subcollasal_mask],
    ["putamen_lh", "putamen_rh", "subcollasal"]
):
    mask = nib.load(mask_f)
    mask_data = mask.get_fdata().astype(bool)
    masked_data = data[mask_data, :]  
    volumes = masked_data.sum(axis=0) * voxel_vol
    df_main[col_name] = volumes

## prepare dfs
df_main.drop(columns=["group", "subject_ID", "distance", "weights", "subclass"], errors="ignore", inplace=True)
df_main.rename(columns={"tinception_id": "subjects"}, inplace=True)
region_cols = [
    "Accessory-Basal-nucleus-lh", "Basal-nucleus-lh", "Lateral-nucleus-lh",
    "GC-ML-DG-body-rh", "presubiculum-head-rh",
    "MDm-lh", "MDm-rh", "MDl-lh",
    "CeM-rh", "CeM-lh", "PuA-lh"
]
df_regions = df_regions[region_cols + ["subjects", "group"]]
df_main = df_main.merge(df_regions, on="subjects", how="inner")

covars = ["sex", "age", "site",	"PTA"]
mask_vols = ["putamen_lh", "putamen_rh", "subcollasal"]
reorder_cols = ["subjects", "group", "THI"] + covars + mask_vols + region_cols
df_main = df_main[reorder_cols]
df = df_main.copy()



regions = [
    'putamen_lh', 'putamen_rh', 'subcollasal', 'Accessory-Basal-nucleus-lh', 
    'Basal-nucleus-lh', 'Lateral-nucleus-lh', 'GC-ML-DG-body-rh', 
    'presubiculum-head-rh', 'MDm-lh', 'MDm-rh', 'MDl-lh', 'CeM-rh', 
    'CeM-lh', 'PuA-lh'
]

df_ti = df[df['group'].isin(['TI_0', 'TI_1'])].copy()

def prepare_for_pg(data, cols):
    temp_df = data.copy()
    for col in cols:
        if temp_df[col].dtype == 'object':
            temp_df[col] = pd.Categorical(temp_df[col]).codes
    return temp_df


covs_thi = ['sex', 'age', 'site', 'PTA']
df_pg_thi = prepare_for_pg(df_ti, covs_thi + ['group'])

# 1. Compare THI between TI_1 and TI_2 (ANCOVA)
thi_comp = pg.ancova(data=df_pg_thi, dv='THI', covar=covs_thi, between='group')


# 2. Partial Correlations (THI vs Regions)
thi_corr_list = []
for reg in regions:
    res = pg.partial_corr(data=df_pg_thi, x='THI', y=reg, covar=covs_thi)
    res['region'] = reg
    thi_corr_list.append(res)

df_thi_results = pd.concat(thi_corr_list).reset_index(drop=True)
df_thi_results['p-adj'] = multipletests(df_thi_results['p-val'], method='bonferroni')[1]
df_thi_results.sort_values(by="p-val", inplace=True)

covs_pta = ['sex', 'age', 'site']
df_pg_pta = prepare_for_pg(df_ti, covs_pta + ['group'])

# 3. Compare PTA between TI_1 and TI_2
pta_comp = pg.ancova(data=df_pg_pta, dv='PTA', covar=covs_pta, between='group')

# 4. Partial Correlations (THI vs Regions)
pta_corr_list = []
for reg in regions:
    res = pg.partial_corr(data=df_pg_pta, x='PTA', y=reg, covar=covs_pta)
    res['region'] = reg
    pta_corr_list.append(res)

df_pta_results = pd.concat(pta_corr_list).reset_index(drop=True)
df_pta_results['p-adj'] = multipletests(df_pta_results['p-val'], method='bonferroni')[1]

# --- Display Results ---
print("--- THI Group Comparison ---")
print(thi_comp)
print("\n--- THI vs Regions (Top 5) ---")
print(df_thi_results[['region', 'r', 'p-val', 'p-adj']].sort_values('p-val').head())

print("\n--- PTA Group Comparison ---")
print(pta_comp)
print("\n--- PTA vs Regions (Top 5) ---")
print(df_pta_results[['region', 'r', 'p-val', 'p-adj']].sort_values('p-val').head())



columns_to_keep = ['group', 'THI', 'putamen_lh', 'putamen_rh']
df_subset = df_ti[columns_to_keep].copy()

df_long = df_subset.melt(id_vars=['group', 'THI'], 
                        var_name='Region', 
                        value_name='Volume_mm3')
df_long["hemi"] = df_long["Region"].str.split("_").str[-1]
df_long["Region"] = df_long["Region"].str.split("_").str[0]

point_palette = sns.color_palette("viridis", n_colors=2)[:1]
line_palette = ["green"]

# Map hemi to colors
hemi_colors = dict(zip(["lh", "rh"], point_palette * 2))
hemi_lines = dict(zip(["lh", "rh"], line_palette * 2))

# Create FacetGrid
g = sns.FacetGrid(
    data=df_long,
    col="hemi",
    col_order=["lh", "rh"],
    height=3.3,
    aspect=1.1,
    sharex=True,
    sharey=True
)

for ax, hemi in zip(g.axes.flat, ["lh", "rh"]):
    subset = df_long[df_long['hemi'] == hemi]
    sns.regplot(
        data=subset,
        x='Volume_mm3',
        y='THI',
        scatter_kws={'s': 60, 'alpha': 0.8, 'color': hemi_colors[hemi], 'edgecolor':'k'},
        line_kws={'color': hemi_lines[hemi], 'linewidth': 2.5},
        #ci=None,
        ax=ax
    )
g.set_axis_labels("Volume (mm³)", "THI Score")
g.tight_layout()
plt.show()
g.savefig(tinception_dir / "plots" / "correlations" / "putamen_vs_THI.pdf", dpi=600, bbox_inches='tight')



df_thi = pd.read_csv(tinception_dir / "VBM_design" / "thi.csv")
df_main = df_main.merge(
            df_thi[["subjects", "UMAP1", "UMAP2"]],
            on="subjects",
            how="left"
            )
df = df_main.copy()

regions = ['UMAP1', 'UMAP2']

df_ti = df[df['group'].isin(['TI_0', 'TI_1'])].copy()

def prepare_for_pg(data, cols):
    temp_df = data.copy()
    for col in cols:
        if temp_df[col].dtype == 'object':
            temp_df[col] = pd.Categorical(temp_df[col]).codes
    return temp_df


covs_thi = ['sex', 'age', 'site', 'PTA']
df_pg_thi = prepare_for_pg(df_ti, covs_thi + ['group'])

# 2. Partial Correlations (THI vs Regions)
thi_corr_list = []
for reg in regions:
    res = pg.partial_corr(data=df_pg_thi, x='THI', y=reg, covar=covs_thi)
    res['region'] = reg
    thi_corr_list.append(res)

df_thi_results = pd.concat(thi_corr_list).reset_index(drop=True)
df_thi_results['p-adj'] = multipletests(df_thi_results['p-val'], method='bonferroni')[1]
df_thi_results.sort_values(by="p-val", inplace=True)

covs_pta = ['sex', 'age', 'site']
df_pg_pta = prepare_for_pg(df_ti, covs_pta + ['group'])

# 4. Partial Correlations (PTA vs Regions)
pta_corr_list = []
for reg in regions:
    res = pg.partial_corr(data=df_pg_pta, x='PTA', y=reg, covar=covs_pta)
    res['region'] = reg
    pta_corr_list.append(res)

df_pta_results = pd.concat(pta_corr_list).reset_index(drop=True)
df_pta_results['p-adj'] = multipletests(df_pta_results['p-val'], method='bonferroni')[1]

# --- Display Results ---
print("\n--- THI vs Regions (Top 5) ---")
print(df_thi_results[['region', 'r', 'p-val', 'p-adj']].sort_values('p-val').head())
print("\n--- PTA vs Regions (Top 5) ---")
print(df_pta_results[['region', 'r', 'p-val', 'p-adj']].sort_values('p-val').head())



columns_to_keep = ['group', 'PTA', 'UMAP1', 'UMAP2']
df_subset = df_ti[columns_to_keep].copy()

df_long = df_subset.melt(id_vars=['group', 'PTA'], 
                        var_name='umap', 
                        value_name='value')

df_long["axis"] = df_long["umap"].str.split("P").str[-1]
df_long["umap"] = df_long["umap"].str[:-1]

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="whitegrid", rc=custom_params)

g = sns.FacetGrid(
    data=df_long,
    col="axis",
    col_order=["1", "2"],
    height=4,
    aspect=1,
    sharex=False,
    sharey=True
)
g.map_dataframe(
            sns.regplot,
            scatter=False,
            x='value',
            y='PTA',
            line_kws={'color': "k", 'linewidth': 2.5}
            )
g.map_dataframe(
        sns.scatterplot,
        palette= ['#440154', '#21908C'],
        x='value',
        y='PTA',
        hue="group",
        edgecolor='k',
        s=60
        )
g.set_axis_labels("Value", "PTA Score")
g.tight_layout()

for ax in g.axes[0]:
    ax.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth=0.8, color='#d1d1d1', alpha=0.7, zorder=0)
    ax.grid(which='minor', linestyle=':', linewidth=0.5, color='#e5e5e5', alpha=0.5, zorder=0)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=15))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=12))

    for spine in ax.spines.values():
        spine.set_edgecolor('k')
        spine.set_linewidth(2)

plt.show()
g.savefig(tinception_dir / "plots" / "correlations" / "PTA_vs_UMAP.pdf", dpi=600, bbox_inches='tight')