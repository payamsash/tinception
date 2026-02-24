from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as stats
import pingouin as pg

## plot for number of extremes
atlas_name = "thalamic_nuclei"
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
plots_dir = tinception_dir / "plots"
vbm_design = tinception_dir / "VBM_design"

threshold = 1.96
models_dir = tinception_dir / "subcortical_roi" / "norm_models_ukb"
fname = models_dir / atlas_name / "results" / f"Z_test.csv"
fname_train = models_dir / atlas_name / "results" / f"Z_train.csv"
df_dev = pd.read_csv(fname)
df_tr = pd.read_csv(fname_train)

df_dev["group"] = df_dev["subject_ids"].isin(df_tr["subject_ids"]).map({True: "CO", False: "TI"})
df_dev_ti = df_dev.query('group == "TI"')
df_dev_cols = df_dev.columns[2:-1] # exclude group and observation and tinception_id
df_dev_ti = df_dev_ti[df_dev_cols]

results = []
for region in df_dev_cols:
    is_extreme = (df_dev[region].abs() > threshold) #.abs()
    contingency = pd.crosstab(df_dev['group'], is_extreme)
    chi2, p, dof, ex = stats.chi2_contingency(contingency)
    
    # Calculate percentages for the plot
    counts = contingency.get(True, pd.Series({'CO': 0, 'TI': 0}))
    total = df_dev['group'].value_counts()
    
    results.append({
        'region': region,
        'perc_CO': (counts['CO'] / total['CO']) * 100,
        'perc_TI': (counts['TI'] / total['TI']) * 100,
        'chi2_p': p
    })

df_stats = pd.DataFrame(results)
df_stats['p_fdr'] = pg.multicomp(df_stats['chi2_p'].values, method='fdr_bh')[1]
df_stats.sort_values(by="p_fdr", ascending=True, inplace=True)
#df_stats.reset_index(drop=True).to_csv(tinception_dir / "subcortical_roi" / "stats" / f"{atlas_name}_extreme.csv", index=False)
df_stats.sort_values(by="p_fdr", ascending=False, inplace=True)

## plot 1
extreme_counts = (
    ((df_dev_ti > threshold) | (df_dev_ti < -threshold))
    .sum()
    .reset_index(name="n_extreme")
    .rename(columns={"index": "region"})
)
extreme_counts.sort_values(by="n_extreme", inplace=True)
extreme_counts["region"] = (
    extreme_counts["region"]
    .str.replace("^Volume_of_", "", regex=True)
    .str.replace("_left_hemisphere$", " (lh)", regex=True)
    .str.replace("_right_hemisphere$", " (rh)", regex=True)
    .str.replace("_", " ")
)
order = extreme_counts["region"].values

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="white", rc=custom_params)
fig, ax = plt.subplots(1, 1, figsize=(10, 3), layout="tight")
ax.set_facecolor("white")
fig.patch.set_facecolor("white")    
sns.barplot(
            data=extreme_counts,
            x="region",
            y="n_extreme",
            palette="rocket",
            fill=True,
            dodge=False,
            gap=0.05,
            ax=ax,
            )

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.tick_params(axis="y", which="minor", length=8, width=5.5)
ax.set_xlabel("")
fig.savefig(
            plots_dir / "rois" / f"{atlas_name}_extreme_count_ukb.pdf",
            format="pdf",
            dpi=300,
            bbox_inches="tight"
            )

df_plot = df_stats.melt(
    id_vars=['region', 'p_fdr'], 
    value_vars=['perc_CO', 'perc_TI'],
    var_name='Group', 
    value_name='Percentage'
)
df_plot['Group'] = df_plot['Group'].str.replace('perc_', '')
df_plot['region'] = (
    df_plot['region']
    .str.replace("^Volume_of_", "", regex=True)
    .str.replace("_left_hemisphere$", " (lh)", regex=True)
    .str.replace("_right_hemisphere$", " (rh)", regex=True)
    .str.replace("_", " ")
)

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="white", rc=custom_params)
fig, ax = plt.subplots(1, 1, figsize=(10, 3), layout="tight")
ax.set_facecolor("white")
fig.patch.set_facecolor("white")    

palette = {"CO": '#1f77b4', "TI": '#d62728'}
sns.barplot(
    data=df_plot, 
    x="region", 
    y="Percentage", 
    hue="Group", 
    palette=palette,
    ax=ax,
    order=order,
    edgecolor='black',
    alpha=0.8
)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.tick_params(axis="y", which="minor", length=8, width=5.5)
ax.set_xlabel("")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.01, 1.2))
fig.savefig(
            plots_dir / "rois" / f"{atlas_name}_extreme_stats_ukb.pdf",
            format="pdf",
            dpi=300,
            bbox_inches="tight"
            )

## manhattan plot
df_stats['-log10p'] = -np.log10(df_stats['chi2_p'])
df_stats['is_significant'] = df_stats['p_fdr'] < 0.05

# Clean names for plotting
df_stats['clean_region'] = (
    df_stats['region']
    .str.replace("^Volume_of_", "", regex=True)
    .str.replace("_left_hemisphere$", " (lh)", regex=True)
    .str.replace("_right_hemisphere$", " (rh)", regex=True)
    .str.replace("_", " ")
)

fig, ax = plt.subplots(1, 1, figsize=(9, 3.5), layout="tight")

sns.scatterplot(data=df_stats[~df_stats['is_significant']], 
                x='clean_region', y='-log10p', color='grey', alpha=0.5, s=100, ax=ax)
sns.scatterplot(data=df_stats[df_stats['is_significant']], 
                x='clean_region', y='-log10p', color='#d62728', s=200, edgecolor='black', ax=ax)


fdr_thresh = df_stats[df_stats['p_fdr'] <= 0.05]['chi2_p'].max()
if not np.isnan(fdr_thresh):
    ax.axhline(-np.log10(fdr_thresh), color='black', linestyle='--', alpha=0.7, label='FDR < 0.05')


plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_ylabel("-log10(p-value)", fontweight='bold')
ax.set_title(f"Thalamic Extremes: Significance Map ({atlas_name})", fontsize=14, pad=20)
ax.set_ylim([-0.2, 3.5])
sns.despine()
fig.savefig(plots_dir / "rois" / f"{atlas_name}_manhattan_ukb.pdf", bbox_inches="tight")

## rainfall plot
df_stats['Ratio'] = df_stats['perc_TI'] / (df_stats['perc_CO'] + 0.1) 

fig, ax = plt.subplots(1, 1, figsize=(9, 3.5), layout="tight")
df_stats = df_stats.sort_values('Ratio', ascending=False)

sns.pointplot(data=df_stats, x='clean_region', y='Ratio', join=False, 
              hue='is_significant', palette={True: '#d62728', False: '#1f77b4'}, ax=ax)
ax.axhline(1, color='grey', linestyle='--') # 1 = Equality

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Odds Ratio (TI % / CO %)", fontweight='bold')
ax.set_title("Relative Prevalence of Structural Extremes", fontsize=14)
ax.legend(title="FDR Significant", frameon=False, bbox_to_anchor=(1, 0.9))
sns.despine()
fig.savefig(plots_dir / "rois" / f"{atlas_name}_ratio_ukb.pdf", bbox_inches="tight")