from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as stats
import pingouin as pg

## compute if deviations from norm is significant for tinnitus subjects?
def main(atlas_name):

    tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
    plots_dir = tinception_dir / "plots"
    vbm_design = tinception_dir / "VBM_design"

    threshold = 1.96
    models_dir = tinception_dir / "subcortical_roi" / "norm_models"
    fname = models_dir / atlas_name / "results" / f"Z_test.csv"
    df_dev = pd.read_csv(fname)
    df_master = pd.read_csv(vbm_design / "covars.csv")

    if atlas_name == "Hippo":
        df_dev.drop(columns=['Whole_hippocampus-lh', 'Whole_hippocampus-rh'], inplace=True)
    if atlas_name == "Amygdala":
        df_dev.drop(columns=['Whole_amygdala-lh', 'Whole_amygdala-rh'], inplace=True)
    if atlas_name == "Thalamic-nuclei":
        df_dev.drop(columns=['Right-Whole_thalamus', 'Left-Whole_thalamus'], inplace=True)   

    df_dev.rename(columns={"subject_ids" : "tinception_id"}, inplace=True)
    df_dev = df_dev.merge(
                            df_master[["tinception_id", "group"]],
                            on="tinception_id",
                            how="inner"
                        )
    df_dev_ti = df_dev.query('group == "TI"')
    df_dev_cols = df_dev.columns[2:-1] # exclude group and observation and tinception_id
    df_dev_ti = df_dev_ti[df_dev_cols]

    results = []
    for region in df_dev_cols:
        is_extreme = (df_dev[region].abs() > threshold)
        contingency = pd.crosstab(df_dev['group'], is_extreme)
        chi2, p, dof, ex = stats.chi2_contingency(contingency)
        
        # Calculate percentages for the plot
        counts = contingency.get(True, pd.Series({'CO': 0, 'TI': 0}))
        total = df_dev['group'].value_counts()
        
        results.append({
            'region': region,
            'perc_CO': (counts['CO'] / total['CO']) * 100,
            'perc_TI': (counts['TI'] / total['TI']) * 100,
            'chi2': chi2,
            'dof': dof,
            'chi2_p': p
        })

    df_stats = pd.DataFrame(results)
    df_stats['p_fdr'] = pg.multicomp(df_stats['chi2_p'].values, method='fdr_bh')[1]
    df_stats.sort_values(by="p_fdr", ascending=True, inplace=True)
    df_stats.reset_index(drop=True).to_csv(tinception_dir / "subcortical_roi" / "stats" / f"{atlas_name}_extreme.csv", index=False)
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
                ax=ax
                )

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.tick_params(axis="y", which="minor", length=8, width=5.5)
    ax.set_xlabel("")
    fig.savefig(
                plots_dir / "rois" / f"{atlas_name}_extreme_count.pdf",
                format="pdf",
                dpi=300,
                bbox_inches="tight"
                )

    ## plot 2
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
        edgecolor='black',
        order=order,
        alpha=0.8
    )

    for i, region in enumerate(df_stats['region']):
        p_val = df_stats.loc[df_stats['region'] == region, 'p_fdr'].values[0]
        if p_val < 0.05:
            max_h = df_stats.loc[df_stats['region'] == region, ['perc_CO', 'perc_TI']].max(axis=1).values[0]
            ax.text(i, max_h + 0.5, '*', ha='center', va='bottom', color='black', fontsize=20, fontweight='bold')

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.tick_params(axis="y", which="minor", length=8, width=5.5)
    ax.set_xlabel("")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.01, 1.2))
    fig.savefig(
                plots_dir / "rois" / f"{atlas_name}_extreme_stats.pdf",
                format="pdf",
                dpi=300,
                bbox_inches="tight"
                )

if __name__ == "__main__":
    atlas_names = ["Hippo", "Amygdala", "Thalamic-nuclei"]
    for atlas in atlas_names:
        main(atlas)