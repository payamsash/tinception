from pathlib import Path
import pandas as pd
import seaborn as sns
from nilearn import image, plotting
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img


tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
demographics_dir = Path("../../master_files")
df_master = pd.read_csv(demographics_dir / "matched" / "matched_optimal.csv")

## read site informations
sites = [
        "andes", "antinomics", "cogtail",
        "erlangen", "iccac", "inditms", 
        "london_1", "london_2", "neuropren",
        "new-castle", "tinmeg", "tinspect",
        "triple"
]
dfs = []
for site in sites:
    fname = demographics_dir / f"{site}_master.xlsx"
    dfs.append(pd.read_excel(fname, dtype={"subject_ID": str}))

df_sites = pd.concat(dfs, axis=0)

## create audio df
audio_cols = list(df_sites.columns[7:])
audio_cols_sorted = sorted(
                        audio_cols,
                        key=lambda x: (x.split('-')[0] != 'LH', int(x.split('-')[1]))
                    )

df_audio = df_master.merge(
                df_sites[["subject_ID"] + audio_cols_sorted],
                on="subject_ID",
                how="left"
)
cols_to_keep = ["subject_ID", "tinception_id", "group", "site", "PTA"] + audio_cols_sorted
df_audio = df_audio[cols_to_keep]

## add biotype information
dfs = []
for mode in ["CO", "TI"]:
    df_bio = pd.read_csv(tinception_dir / "biotypes" / f"main_{mode}.csv")
    dfs.append(df_bio[["subjects", "Biotype"]])
df_bio = pd.concat(dfs, axis=0)
df_bio.rename(columns={"subjects": "tinception_id"}, inplace=True)
df_audio = df_audio.merge(df_bio, on="tinception_id", how="left")

## make long format
df_long = df_audio.melt(
    id_vars=["subject_ID", "tinception_id", "group", "site", "PTA", "Biotype"],
    value_vars=audio_cols,
    var_name="hemi_freq",
    value_name="threshold"
)

df_long[["hemi", "freq"]] = df_long["hemi_freq"].str.split("-", expand=True)
df_long["freq"] = pd.to_numeric(df_long["freq"], errors="coerce")
df_long = df_long.drop(columns="hemi_freq")
df_long = df_long.query('site != "triple"')

## create biotypes df
biotypes = [0, 1]
df_b = df_long.query('Biotype == @biotypes')


def plot_audiometry(data, hue, row=None, col="hemi", palette=None, 
                    ylim=(70, 0), filename="plot.pdf", markersize=6, 
                    err_lw=2.7, show_xticklabels=True):
    
    # 1. Initialize Grid
    g = sns.FacetGrid(
        data=data, row=row, col=col, 
        col_order=["LH", "RH"] if col=="hemi" else None,
        height=3, aspect=1.4, legend_out=True, despine=False
    )
    
    # 2. Map Pointplot
    g.map_dataframe(
        sns.pointplot, x="freq", y="threshold", hue=hue,
        errorbar="se", palette=palette, markers=["o", "o"],
        markersize=markersize, linestyles="--", linewidth=2,
        capsize=0.35, err_kws={"linewidth": err_lw}
    )

    # 3. Frequency Label Logic
    tick_pos = sorted(data["freq"].unique())
    tick_labs = [f"{v/1000:g}k" if v >= 1000 else str(v) for v in tick_pos]

    # 4. Standardize Axes Styling
    for ax in g.axes.flat:
        ax.invert_yaxis()
        ax.set_ylim(ylim)
        ax.grid(True, axis="both", linewidth=0.7, alpha=0.3)
        
        # X-axis on top
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.tick_params(axis='x', rotation=45)
        
        if show_xticklabels:
            ax.set_xticklabels(tick_labs, ha="left", fontsize=10)
        else:
            ax.set_xticklabels([])

        # Spine styling
        ax.spines[["right", "bottom"]].set_visible(False)
        for s in ["top", "left"]: ax.spines[s].set_linewidth(1.5)

    # 5. Final Layout & Save
    g.set_xlabels("")
    g.add_legend()
    g.tight_layout()
    g.fig.subplots_adjust(top=0.85)
    
    save_path = tinception_dir / "plots" / "audiometry" / filename
    g.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")
    return g

# 1. Main Audio Plot
plot_audiometry(df_long, hue="group", palette={"CO": "#1f77b4", "TI": "#d62728"}, 
                filename="audio.pdf")

# 2. Sites Plot (Specific tweaks for row/col and markers)
plot_audiometry(df_long, hue="group", row="hemi", col="site", 
                palette={"CO": "#1f77b4", "TI": "#d62728"}, 
                ylim=(70, -10), markersize=3, err_lw=3.2, 
                show_xticklabels=False, filename="sites.pdf")

# 3. Biotypes Plot
plot_audiometry(df_b, hue="Biotype", row="group", 
                palette={0: "#440154", 1: "#21908C"}, 
                filename="biotypes.pdf")


## first add PTA_HF column
all_hf_freqs = ['6000', '8000', '9000', '10000', '11200', '12500', '14000', '16000']
hf_freqs = ['6000', '8000', '12500']
hf_cols = [col for col in audio_cols if col[3:] in hf_freqs]

cols_to_keep = ["subject_ID", "tinception_id", "group", "site", "PTA"] + hf_cols
df_hf = df_audio[cols_to_keep].copy()
df_hf[hf_cols] = df_hf[hf_cols].apply(pd.to_numeric, errors='coerce')
df_hf = df_hf.dropna(subset=hf_cols, how="any")

keep_cols = [f"{hemi}-{freq}" for freq in hf_freqs for hemi in ["LH", "RH"]]
df_hf[keep_cols] = df_hf[keep_cols].apply(pd.to_numeric, errors='coerce')
df_hf["PTA_HF"] = df_hf[keep_cols].mean(axis=1)
df_hf.drop(columns=keep_cols + ["subject_ID"], inplace=True)

## now add other columns
df_main = pd.read_csv(tinception_dir / "VBM_design" / "covars.csv")
df_hf = df_hf.merge(
            df_main[["tinception_id", "sex", "age", "TIV"]],
            on="tinception_id",
            how="left"
            )

## fix indexes and get GMs
df_model = df_main[df_main["tinception_id"].isin(df_hf["tinception_id"])].copy()
subset_idx = df_model.index.to_list()
subset_idx = df_model.index.to_list()

gms_fname = tinception_dir / "vbm_norm" / "GM_mod_merg_s3.nii.gz"
img_4d_subset = image.index_img(gms_fname, subset_idx)

## create covars and design matrix
design_df = df_hf[["group", "age", "sex", "site", "PTA", "PTA_HF", "TIV"]].copy()
for col in ["age", "PTA", "PTA_HF", "TIV"]:
    design_df[col] = pd.to_numeric(design_df[col], errors="coerce")
    design_df[col] = (design_df[col] - design_df[col].mean()) / design_df[col].std()


design_df["group"] = design_df["group"].astype("category")
design_df["sex"] = design_df["sex"].astype("category")
design_df["site"] = design_df["site"].astype("category")

design_matrix = pd.get_dummies(
    design_df,
    columns=["group", "sex", "site"],
    drop_first=True
)
design_matrix["intercept"] = 1.0
design_matrix = design_matrix.astype(float)

## run GLMs
group_cols = [c for c in design_matrix.columns if c.startswith("group_")]
group_col = group_cols[0]
print("Using group column:", group_col)
slm = SecondLevelModel(
    smoothing_fwhm=3.0,
    mask_img=None
)
slm = slm.fit(img_4d_subset, design_matrix=design_matrix)
z_map_ti_gt_co = slm.compute_contrast(group_col, output_type="z_score")
thresholded_ti_gt_co, threshold_ti_gt_co = threshold_stats_img(
                                                        z_map_ti_gt_co,
                                                        alpha=0.05,
                                                        height_control="fdr"
                                                    )

z_map_co_gt_ti = slm.compute_contrast(f"-{group_col}", output_type="z_score")
thresholded_co_gt_ti, threshold_co_gt_ti = threshold_stats_img(
                                                        z_map_co_gt_ti,
                                                        alpha=0.05,
                                                        height_control="fdr"
                                                    )

fsl_dir = Path("/Users/payamsadeghishabestari/fsl")
img_bg = fsl_dir / "data" / "standard" / "MNI152_T1_0.5mm.nii.gz"

for sel_image, thr, coord, title in zip(
                                [thresholded_co_gt_ti, thresholded_ti_gt_co],
                                [threshold_co_gt_ti, threshold_ti_gt_co],
                                [(-13.35, 24.13, -16.7), (-24.95, 3.01, -10.14)],
                                ["co_gt_ti", "ti_gt_co"]
                                ):
    kwargs = {
                "colorbar": False,
                "cbar_tick_format": "%.2g",
                "annotate": False,
                "draw_cross": False,
                "radiological": False,
                "cmap": 'Reds',
                "threshold": thr,
                "symmetric_cbar": False,
                "vmin": None,
                "vmax": None,
                "dim": -0.3,
                "black_bg": True,
                "cut_coords": coord
            }
    fig = plotting.plot_stat_map(
                stat_map_img=sel_image,
                bg_img=img_bg,
                display_mode="ortho",
                **kwargs
                )
    fig.savefig(tinception_dir / "plots" / "vbm_group" / f"{title}_pta_hf.pdf", dpi=600, bbox_inches='tight')