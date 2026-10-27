from pathlib import Path
import pandas as pd
import abagen
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, spearmanr

from neuromaps import nulls
from neuromaps import nulls, datasets
import nibabel as nib
from neuromaps import nulls, stats as nmstats
from neuromaps.images import annot_to_gifti


## paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
ssa_results_dir = tinception_dir / "ssa_results"
df_grad = pd.read_csv(ssa_results_dir / "grads.csv")


################ initialization ################
def get_t_map(group_a, group_b):
    t_stats, _ = ttest_ind(group_a, group_b, axis=0, nan_policy='omit')
    return np.nan_to_num(t_stats)

def get_surface_centroids(coords, labels):
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]
    centroids = np.array([coords[labels == l].mean(axis=0) for l in unique_labels])
    centroids /= np.linalg.norm(centroids, axis=1, keepdims=True)
    return centroids

lh_annot = str(ssa_results_dir / "lh.Schaefer2018_1000Parcels_7Networks_order.annot")
rh_annot = str(ssa_results_dir /"rh.Schaefer2018_1000Parcels_7Networks_order.annot")
lh_parc = annot_to_gifti(lh_annot)
rh_parc = annot_to_gifti(rh_annot)

lh_labels, lh_ctab, lh_names = nib.freesurfer.read_annot(str(lh_annot))
rh_labels, rh_ctab, rh_names = nib.freesurfer.read_annot(str(rh_annot))
surfs = datasets.fetch_fsaverage(density='164k')

lh_sphere_coords, _ = nib.load(surfs['sphere'][0]).agg_data()
rh_sphere_coords, _ = nib.load(surfs['sphere'][1]).agg_data()

lh_centroids = get_surface_centroids(lh_sphere_coords, lh_labels)
rh_centroids = get_surface_centroids(rh_sphere_coords, rh_labels)

all_centroids = np.vstack([lh_centroids, rh_centroids])
hemi_labels = np.hstack([np.zeros(len(lh_centroids)), np.ones(len(rh_centroids))])

################ load gene data ################
# 1. Fix set_axis compatibility
_old_set_axis = pd.DataFrame.set_axis
def _new_set_axis(self, labels, axis=0, inplace=None, **kwargs):
    if inplace is not None:
        if inplace:
            if axis in [0, 'index']: self.index = labels
            else: self.columns = labels
            return None
    return _old_set_axis(self, labels, axis=axis, **kwargs)
pd.DataFrame.set_axis = _new_set_axis

# 2. Fix the missing .append() method
def _patched_append(self, other, ignore_index=False, verify_integrity=False, sort=False):
    # This mimics the old df.append behavior using pd.concat
    if isinstance(other, (pd.Series, dict)):
        if isinstance(other, dict):
            other = pd.Series(other)
        if other.name is None and not ignore_index:
            raise TypeError("Can only append a Series if its name is set")
        return pd.concat([self, other.to_frame().T], ignore_index=ignore_index)
    return pd.concat([self, other], ignore_index=ignore_index)

pd.DataFrame.append = _patched_append

# fetch abagen
atlas_path = ssa_results_dir / "Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.nii.gz"
expression = abagen.get_expression_data(atlas_path, 
                                        region_agg='donors', 
                                        verbose=1)


################ load deviation scores ################
df_devs = []
for n_grad in range(1, 4):
    fname = ssa_results_dir / "norm_models"/ f"gradient_{n_grad}" / "results" / "Z_test.csv"
    df_dev = pd.read_csv(fname)
    df_g = df_grad.query(f'component == {n_grad}')

    df_g.rename(columns={"subjects": "subject_ids"}, inplace=True)
    df_dev = df_dev.merge(
                    df_g[["subject_ids", "group"]],
                    on="subject_ids",
                    how="inner"
                    )
    df_dev["component"] = n_grad
    df_dev.drop(columns="observations", inplace=True)

    cols = df_dev.columns.tolist()
    last_two = cols[-2:] 
    remaining_cols = cols[:-2]
    new_order = [remaining_cols[0]] + last_two + remaining_cols[1:]
    df_devs.append(df_dev[new_order])

df_zs = pd.concat(df_devs)
feature_cols = list(df_zs.columns[-1000:])

stats_list = []
for n_grad in range(1, 4):
    df_z = df_zs.query(f'component == {n_grad}')
    df_ti = df_z.query('group == "TI"')
    df_co = df_z.query('group == "CO"')
    # df_bio1 =
    # df_bio2 = 

    ########################## compute stats ##########################
    t_map = get_t_map(df_ti[feature_cols].values, df_co[feature_cols].values)

    brain_map_series = pd.Series(t_map, index=range(1, 1001))
    clean_expression = expression.dropna(how='all')
    common_labels = clean_expression.index.intersection(brain_map_series.index)
    final_expression = clean_expression.loc[common_labels]
    final_brain_map = brain_map_series.loc[common_labels]

    print(f"Original parcels: 1000")
    print(f"Parcels with genetic data: {final_expression.shape[0]}")

    gene_corrs = []
    for gene_name in expression.columns:
        rho, p = spearmanr(final_brain_map, final_expression[gene_name])
        gene_corrs.append(
            {'gradient': n_grad,
            'gene': gene_name,
            'rho': rho,
            'p': p})

    df_results = pd.DataFrame(gene_corrs).sort_values('rho', ascending=False)
    top_15_genes = df_results.head(15)['gene'].tolist()

    for gene in top_15_genes:
        x_map = expression[gene].values
        spins = nulls.alexander_bloch(
                                    x_map, 
                                    atlas='fsaverage', 
                                    density='164k',
                                    parcellation=(lh_parc[0], rh_parc[0]),
                                    n_perm=1000, 
                                    seed=42
                                    )
    
        corr, p_spatial = nmstats.compare_images(
                                                x_map, 
                                                brain_map_series.values, 
                                                nulls=spins, 
                                                metric='spearmanr'
                                                )
        stats_list.append({
            'Gradient': n_grad,
            'Rho': rho,
            'p': p,
            'Gene': gene,
            'Spatial_Rho': corr,
            'P_Spatial': p_spatial,
            'Significant': p_spatial < 0.05
        })
        print(f"Processed {gene}: p_spatial = {p_spatial:.4f}")

df_top_stats = pd.DataFrame(stats_list)
df_top_stats = df_top_stats.sort_values('P_Spatial')
df_top_stats.to_csv(ssa_results_dir / "genes.csv", index=False)