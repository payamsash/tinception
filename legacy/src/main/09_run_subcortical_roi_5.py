from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

from math import pi
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import StandardScaler
import umap


tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
plots_dir = tinception_dir / "plots" / "rois"
rois_dir = tinception_dir / "subcortical_roi"

atlas_names_main = ["Amygdala", "Hippo", "Thalamic-nuclei"]
atlas_names_ukb = ["amygdalar_nuclei", "hippo_subfields", "thalamic_nuclei"]

## read main dataset
dfs = []
for atlas_name in atlas_names_main:
    fname = rois_dir / "stats" / f"{atlas_name}.csv"
    dfs.append(pd.read_csv(fname))

df_main = pd.concat(dfs)
df_main["Region"] = (
    df_main["Region"]
    .str.replace("_", "-", regex=False)
)
df_main["Region"] = df_main["Region"].str.replace(
    r"^(Left|Right)-(.*)$",
    lambda m: f"{m.group(2)}-{'lh' if m.group(1) == 'Left' else 'rh'}",
    regex=True
)
df_main["Region"] = df_main["Region"].str.replace(r"[()]", "", regex=True)
to_drop_main = ["Paralaminar-nucleus-rh"]
df_main = df_main[~df_main["Region"].isin(to_drop_main)]

## read ukb dataset
dfs = []
for atlas_name in atlas_names_ukb:
    fname = rois_dir / "ukb" / f"{atlas_name}_stat_bonferroni_results.csv"
    dfs.append(pd.read_csv(fname))

df_ukb = pd.concat(dfs)
df_ukb["brain_label"] = (
    df_ukb["brain_label"]
    .str.replace("_", "-", regex=False)                 # "_" → "-"
    .str.replace(r"^Volume-of-", "", regex=True)        # remove prefix at start
    .str.replace("left-hemisphere", "lh", regex=False)  # left → lh
    .str.replace("right-hemisphere", "rh", regex=False) # right → rh
)
to_drop_ukb = [
    "Whole-hippocampal-body-rh",
    "Whole-hippocampal-body-lh",
    "Whole-hippocampal-head-lh",
]
df_ukb = df_ukb[~df_ukb["brain_label"].isin(to_drop_ukb)]

cols_to_keep = ["label", "n1", "n2", "d", "pval", "pval_adj", "category"]
df_main.rename(columns={"Region": "label", "Cohen_d": "d", "N_CO": "n1", "N_TI": "n2", "p_val": "pval", "structure": "category", "p_fdr": "pval_adj"}, inplace=True)
df_ukb.rename(columns={"brain_label": "label", "cohen_d": "d", "n_control": "n1", "n_tinnitus": "n2"}, inplace=True)
df_main = df_main[cols_to_keep]
df_ukb = df_ukb[cols_to_keep[:-1]]
df_ukb = df_ukb.merge(
                    df_main[["label", "category"]],
                    on="label",
                    how="inner"
                    )

df_main["dataset"] = "main"
df_ukb["dataset"] = "ukb"
df_plot = pd.concat([df_main, df_ukb], axis=0)

## compute se
df_plot['se'] = np.sqrt(((df_plot['n1'] + df_plot['n2']) / (df_plot['n1'] * df_plot['n2'])) + 
                        (df_plot['d']**2 / (2 * (df_plot['n1'] + df_plot['n2']))))

df_plot['ci_upper'] = df_plot['d'] + (1.96 * df_plot['se'])
df_plot['ci_lower'] = df_plot['d'] - (1.96 * df_plot['se'])

## sorting and pval extraction
df_plot = df_plot.sort_values(['category', 'd'], ascending=[True, False])
df_plot = df_plot.reset_index(drop=True)
df_plot['log_p'] = -np.log10(df_plot['pval'].replace(0, 1e-10))



### plotting
g = sns.FacetGrid(
    data=df_plot, 
    col="category", 
    col_wrap=3, 
    height=6, 
    aspect=0.8, 
    sharey=False, 
    sharex=True
)

def draw_forest_plot(data, **kwargs):
    ax = plt.gca()
    colors = {'ukb': '#2c3e50', 'main': '#d35400'}
    
    # 1. Establish a rigid order based ONLY on 'main' effect sizes
    # Ascending=True puts the largest effects at the top of the plot
    ref_df = data[data['dataset'] == 'main'].sort_values('d', ascending=True)
    ref_order = ref_df['label'].tolist()
    
    # 2. Map every region label to a fixed integer Y-coordinate
    y_map = {label: i for i, label in enumerate(ref_order)}
    
    # 3. Separate data for clean plotting
    main_subset = data[data['dataset'] == 'main']
    ukb_subset = data[data['dataset'] == 'ukb']
    
    # 4. Plot Main
    if not main_subset.empty:
        # Get Y-coords for Main from the map
        y_main = [y_map[l] for l in main_subset['label']]
        ax.errorbar(
            x=main_subset['d'], y=np.array(y_main) - 0.15, 
            xerr=1.96 * main_subset['se'], 
            fmt='o', color=colors['main'], 
            label='main', markersize=5, capsize=3, elinewidth=2
        )
        
    # 5. Plot UKB (Forces squares to sit on the same line as the circles)
    if not ukb_subset.empty:
        # Ensure we only plot UKB labels that exist in our Main reference
        ukb_filtered = ukb_subset[ukb_subset['label'].isin(y_map.keys())]
        y_ukb = [y_map[l] for l in ukb_filtered['label']]
        ax.errorbar(
            x=ukb_filtered['d'], y=np.array(y_ukb) + 0.15, 
            xerr=1.96 * ukb_filtered['se'], 
            fmt='s', color=colors['ukb'], 
            label='ukb', markersize=5, capsize=3, elinewidth=2
        )

    # 6. Set Y-axis labels strictly following the reference order
    ax.set_yticks(range(len(ref_order)))
    ax.set_yticklabels(ref_order, fontdict={"size": 7, "style": "italic"})
    
    # Add the zero-line
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.3, zorder=0)


g.map_dataframe(draw_forest_plot)
g.set_axis_labels("Effect Size (Cohen's $d$)", "")
g.set_titles(col_template="{col_name}")
g.add_legend(title="Cohort")
plt.subplots_adjust(top=0.9, wspace=0.4)
g.fig.suptitle("Subcortical Tinnitus Phenotypes", fontsize=16)
g.savefig(
        plots_dir / "cohen.pdf",
        format="pdf",
        dpi=600,
        bbox_inches="tight"
        )



### heatmap plot
# 1. Metric Calculation & Tier Assignment (Same as before)
df_reliable = df_plot.pivot_table(
    index=['label', 'category'], 
    columns='dataset', 
    values=['d', 'pval', 'pval_adj', 'ci_lower', 'ci_upper']
)
df_reliable.columns = [f"{col[0]}_{col[1]}" for col in df_reliable.columns]
df_reliable = df_reliable.reset_index()

# Convergence Metrics
df_reliable['dir_agree'] = np.sign(df_reliable['d_main']) == np.sign(df_reliable['d_ukb'])
df_reliable['ovl_min'] = np.maximum(df_reliable['ci_lower_main'], df_reliable['ci_lower_ukb'])
df_reliable['ovl_max'] = np.minimum(df_reliable['ci_upper_main'], df_reliable['ci_upper_ukb'])
df_reliable['overlap_dist'] = (df_reliable['ovl_max'] - df_reliable['ovl_min']).clip(lower=0)
df_reliable['total_span'] = np.maximum(df_reliable['ci_upper_main'], df_reliable['ci_upper_ukb']) - \
                            np.minimum(df_reliable['ci_lower_main'], df_reliable['ci_lower_ukb'])
df_reliable['overlap_ratio'] = df_reliable['overlap_dist'] / df_reliable['total_span']

# 2. Assign Numeric Tiers to avoid KeyErrors
def assign_tier_num(row):
    if row['pval_adj_ukb'] < 0.05 and row['pval_main'] < 0.05 and row['dir_agree']:
        return 1
    elif row['pval_ukb'] < 0.05 and row['pval_main'] < 0.05 and row['dir_agree']:
        return 2
    elif (row['pval_ukb'] < 0.05 or row['pval_main'] < 0.05):
        return 3
    return 4

df_reliable['tier_num'] = df_reliable.apply(assign_tier_num, axis=1)
df_reliable['consensus_score'] = df_reliable['overlap_ratio'] * \
                                 (-np.log10(df_reliable['pval_ukb'] + 1e-10)) * \
                                 (-np.log10(df_reliable['pval_main'] + 1e-10))

# Sorting
df_reliable = df_reliable.sort_values(['tier_num', 'consensus_score'], ascending=[True, False])


fig, (ax_main, ax_tier) = plt.subplots(2, 1, figsize=(18, 2.5), 
                                    gridspec_kw={'height_ratios': [6, 0.5], 'hspace': 0.05}, layout="tight")

viz_data_h = df_reliable.set_index('label')[['d_ukb', 'd_main']].T
viz_data_h.index = ['UKB', 'Main']

tier_colors = {1: '#1E8449', 2: '#D4AC0D', 3: '#CA6F1E', 4: '#7F8C8D'} # Professional deep shades
current_tiers = sorted(df_reliable['tier_num'].unique())
tier_cmap = mcolors.ListedColormap([tier_colors[t] for t in current_tiers])


tier_data_h = df_reliable[['tier_num']].T


sns.heatmap(
    viz_data_h, 
    ax=ax_main, 
    annot=False, 
    fmt=".3f", 
    cmap="vlag",
    center=0, 
    cbar=False,
    cbar_kws={'label': "Cohen's $d$", 'pad': 0.01, 'shrink': 0.8},
    linewidths=0.7,
    annot_kws={"size": 9, "weight": "bold"}
)

ax_main.set_xticklabels([])
ax_main.set_xlabel("")
ax_main.set_xticks([])

sns.heatmap(
    tier_data_h, 
    ax=ax_tier, 
    cmap=tier_cmap, 
    cbar=False, 
    xticklabels=df_reliable['label'], 
    annot=False, 
)

plt.setp(ax_tier.get_xticklabels(), rotation=25, ha='right', rotation_mode='anchor', fontsize=5, style='italic')
ax_tier.set_yticklabels(ax_tier.get_yticklabels(), rotation=0)
plt.subplots_adjust(bottom=0.25)
tier_changes = np.where(df_reliable['tier_num'].values[:-1] != df_reliable['tier_num'].values[1:])[0] + 1
for change in tier_changes:
    ax_main.axvline(x=change, color='black', linewidth=2)
    ax_tier.axvline(x=change, color='black', linewidth=2)

plt.tight_layout()
plt.show()
fig.savefig(
    plots_dir / "heatmap.pdf",
    format="pdf",
    dpi=1200,
    bbox_inches="tight"
)


## umap plot
tier_nums = [2, 3]
regions = list(df_reliable.query('tier_num == @tier_nums')["label"].values)

def get_tinnitus_features(dataset="ukb", base_path="/Volumes/Extreme_SSD/payam_data/Tinception"):
    base_path = Path(base_path)
    
    # 1. Define dataset-specific parameters
    if dataset == "ukb":
        config = {
            "atlases": ["amygdalar_nuclei", "hippo_subfields", "thalamic_nuclei"],
            "model_path": base_path / "subcortical_roi" / "norm_models_ukb",
            "file": "Z_test.csv",
            "id_col": "subject_ids",
            "thalamus_atlas": "thalamic_nuclei"
        }
    else: # dataset == "main"
        config = {
            "atlases": ["Amygdala", "Hippo", "Thalamic-nuclei"],
            "model_path": base_path / "subcortical_roi" / "new_norm_models",
            "file": "Z_main.csv",
            "id_col": "tinception_id",
            "thalamus_atlas": "Thalamic-nuclei"
        }

    df_tis = []
    
    # 2. Process each atlas
    for atlas in config["atlases"]:
        df = pd.read_csv(config["model_path"] / atlas / "results" / config["file"])
        
        # Identification Logic
        if dataset == "ukb":
            df_tr = pd.read_csv(config["model_path"] / atlas / "results" / "Z_train.csv")
            df["group"] = df["subject_ids"].isin(df_tr["subject_ids"]).map({True: "CO", False: "TI"})
        else:
            df_master = pd.read_csv(base_path / "VBM_design" / "covars.csv")
            df = df.rename(columns={"subject_ids": "tinception_id"}).merge(
                df_master[["tinception_id", "group"]], on="tinception_id"
            )

        # Filter for Tinnitus group and drop non-feature metadata
        ti_subset = df.query('group == "TI"').copy()
        meta_cols = ["observations", "group", "subject_ids", "tinception_id"]
        ti_subset.drop(columns=[c for c in meta_cols if c in ti_subset.columns], inplace=True)

        # Atlas-specific cleanup
        if atlas == config["thalamus_atlas"]:
            bad_cols = ["Volume_of_Whole_thalamus", "Volume_of_Pf_right_hemisphere"]
            ti_subset.drop(columns=[c for c in bad_cols if c in ti_subset.columns], inplace=True)

        df_tis.append(ti_subset)

    # 3. Horizontal Concatenation and Column Prettifying
    full_df = pd.concat(df_tis, axis=1)
    full_df.columns = (
        full_df.columns
        .str.removeprefix("Volume_of_")
        .str.replace("_right_hemisphere$", "-rh", regex=True)
        .str.replace("_left_hemisphere$", "-lh", regex=True)
        .str.replace("_", "-", regex=False)
    )
    
    return full_df

df_dev_ti = get_tinnitus_features(dataset="ukb")
df_dev_ti.columns = (
    df_dev_ti.columns
    .str.removeprefix("Volume_of_")  # remove prefix
    .str.replace("_right_hemisphere$", "-rh", regex=True)
    .str.replace("_left_hemisphere$", "-lh", regex=True)
    .str.replace("_", "-", regex=False)  # "_" -> "-"
)

feature_cols = [col for col in df_dev_ti.columns if col in regions]

X = df_dev_ti[feature_cols].to_numpy()
X_scaled = StandardScaler().fit_transform(X)


reducer = umap.UMAP(n_neighbors=10, n_components=2, metric='cosine', min_dist=0.001, random_state=42)
embedding = reducer.fit_transform(X_scaled)

df_dev_ti['UMAP1'] = embedding[:, 0]
df_dev_ti['UMAP2'] = embedding[:, 1]


clusterer = HDBSCAN(min_cluster_size=20)
df_dev_ti['Biotype_UMAP'] = clusterer.fit_predict(embedding)
df_dev_ti['Global_Burden'] = np.sqrt((df_dev_ti[feature_cols]**2).sum(axis=1))


custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="whitegrid", rc=custom_params)

fig, ax = plt.subplots(figsize=(7, 4), layout="tight")
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth='0.8', color='#d1d1d1', alpha=0.7, zorder=0)
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='#e5e5e5', alpha=0.5, zorder=0)

# Adjust tick density - increase N to make it even more fine-grained
ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=15))
ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=12))

scatter = sns.scatterplot(
    data=df_dev_ti, 
    x='UMAP1', 
    y='UMAP2', 
    hue='Biotype_UMAP', 
    size='Global_Burden',
    sizes=(40, 200),       # Range of point sizes
    palette=['#440154', '#21908C'],     
    alpha=0.7, 
    edgecolor='black',
    linewidth=0.5,
    ax=ax
)

plt.xlabel("UMAP Dimension 1")
plt.ylabel("UMAP Dimension 2")

handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles=handles, 
    labels=labels, 
    title="Biotype & Structural Load", 
    bbox_to_anchor=(1.02, 1), 
    loc='upper left', 
    frameon=False,
    fontsize=10
)
for spine in ax.spines.values():
    spine.set_edgecolor('k')
    spine.set_linewidth(2)

plt.tight_layout()
plt.show()
fig.savefig(plots_dir / "rois" / f"umap.pdf", format="pdf",
    dpi=600, bbox_inches="tight")


## spider plot
df_profile = df_dev_ti[df_dev_ti['Biotype_UMAP'] != -1].groupby('Biotype_UMAP')[feature_cols].mean().T

# 2. Select top 8 regions with highest variance between clusters
variance = df_profile.var(axis=1)
top_regions = variance.nlargest(10).index.tolist()
df_radar = df_profile.loc[top_regions]

categories = df_radar.index
N = len(categories)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]
fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True), layout="tight")
ax.set_yticklabels([])
plt.xticks(angles[:-1], categories, color='grey', size=10)

# Plot Biotype 0
values = df_radar[0].values.tolist()
values += values[:1]
ax.plot(angles, values, linewidth=2, linestyle='solid', label='Biotype 0', color='#440154')
ax.fill(angles, values, '#440154', alpha=0.25)

# Plot Biotype 1
values = df_radar[1].values.tolist()
values += values[:1]
ax.plot(angles, values, linewidth=2, linestyle='solid', label='Biotype 1', color='#21908C')
ax.fill(angles, values, '#21908C', alpha=0.25)

plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.2), frameon=False)
plt.tight_layout()
fig.savefig(plots_dir / "rois" / "biotype_profiles.pdf", bbox_inches="tight", format="pdf",
    dpi=600)