from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import mne
from mne.viz import Brain


tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
results_dir = tinception_dir / "SBM" / "results"
plots_dir = tinception_dir / "plots" / "sbm"
df_master = pd.read_csv(tinception_dir / "vbm_norm" / "covars.csv")

## fix these 3 subjects later
subjects_to_drop = ["AG_e1", "CK_C_e1", "DW_C_e1"]
df_master = df_master[~df_master['subject_ID'].isin(subjects_to_drop)]

## find significant clusters and create a df
smoothing = 10
measure = "thickness"
mode = "ti_gt_co"
hemis = ["lh", "rh"]
path_dict = dict(zip(hemis, [[], []]))
dfs = []
for hemi in hemis:
    for idx in range(1, 4):
        txt_dir = results_dir / f"{hemi}.{measure}.{smoothing}.glmdir" / mode / f"{hemi}.cluster{idx}.txt"
        summery_dir = results_dir / f"{hemi}.{measure}.{smoothing}.glmdir" / mode / "cache.th20.pos.sig.cluster.summary"
        if txt_dir.is_file():
            path_dict[hemi].append(txt_dir)
            df = pd.read_csv(txt_dir, header=None)
            df.rename(columns={0: measure}, inplace=True)
            df["hemi"] = len(df) * [hemi]
            df["cluster_idx"] = len(df) * [idx]
            df["group"] = df_master["group"]
            dfs.append(df)

df_thickness = pd.concat(dfs, axis=0)
df_thickness["group"] = df_thickness["group"].map({"CO": "Controls", "TI": "Tinnitus"})

## boxplots
pal = ['#1f77b4', '#d62728']
order = ["Controls", "Tinnitus"]

g = sns.FacetGrid(
                data=df_thickness,
                row="cluster_idx",
                col="hemi",
                #xlim=[0, 1],
                height=1.8,
                aspect=3
                )
g.map_dataframe(
            sns.stripplot,
            x="thickness",
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
            x="thickness",
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
g.axes[2, 0].remove()
g.tight_layout()
g.savefig(
        plots_dir / "boxplots.pdf",
        format="pdf",
        dpi=300,
        bbox_inches="tight"
        )

## plot brain
view = "lateral"
brain_kwargs = dict(background="white", surf="inflated", cortex=["#b8b4ac", "#b8b4ac"])
alpha = 0.95
parc = "aparc.a2009s"

brain_scrs = []
for hemi in ["rh", "lh"]:
    ## cluster corrected
    fname = results_dir / f"{hemi}.{measure}.{smoothing}.glmdir" / mode / f"cache.th20.pos.sig.ocn.annot"
    brain = Brain("fsaverage", subjects_dir=None, hemi=hemi, views=view, **brain_kwargs)
    brain.add_annotation(str(fname), hemi=hemi, borders=False, color="#C5A059", alpha=alpha)
    brain.add_annotation(parc, borders=True, color="white")
    brain_scr = brain.screenshot()
    brain_scrs.append(brain_scr)

fig, axes = plt.subplots(1, 2, figsize=(7, 5), layout="constrained")
fig.subplots_adjust(hspace=-0.1, wspace=-0.1)
for ax, brain in zip(axes, brain_scrs):
    ax.imshow(brain)
    ax.axis("off")
        
fig.savefig(
            plots_dir / f"brain.pdf",
            format="pdf",
            dpi=300,
            bbox_inches="tight")


## find the labels
def get_detailed_label(hemi, vertex_index):
    # Load the Destrieux annotation
    labels = mne.read_labels_from_annot(
        subject="fsaverage", 
        parc='aparc.a2009s', 
        hemi=hemi, 
        subjects_dir=None,
        verbose=False
    )
    for label in labels:
        if vertex_index in label.vertices:
            return label.name
    return "Label not found"

peaks = {
    'lh': [86151, 40295],
    'rh': [87022, 112526, 114571]
}
print(f"{'Hemi':<5} | {'Vertex':<8} | {'Destrieux Label'}")
print("-" * 40)
for hemi, vtx_list in peaks.items():
    for v in vtx_list:
        region = get_detailed_label(hemi, v)
        print(f"{hemi:<5} | {v:<8} | {region}")


## compute stats
for hemi in ["lh", "rh"]:
    for cluster_idx in [1, 2, 3]:
        df_sub = df_thickness.query(f'hemi == "{hemi}" & cluster_idx == {cluster_idx}')
        
        group_co = df_sub[df_sub["group"] == "Controls"]["thickness"]
        group_ti = df_sub[df_sub["group"] == "Tinnitus"]["thickness"]
        t_stat, p_val = stats.ttest_ind(group_ti, group_co)
            
        n1, n2 = len(group_ti), len(group_co)
        var1, var2 = np.var(group_ti, ddof=1), np.var(group_co, ddof=1)
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        d = (np.mean(group_ti) - np.mean(group_co)) / pooled_std
            
        print(f"Results for {hemi} and {cluster_idx}:")
        print(f"  Mean (TI): {np.mean(group_ti):.3f}mm | Mean (CO): {np.mean(group_co):.3f}mm")
        print(f"  t-value: {t_stat:.3f}")
        print(f"  Cohen's d: {d:.3f}\n")