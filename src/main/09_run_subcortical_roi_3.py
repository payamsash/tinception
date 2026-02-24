from pathlib import Path
import pandas as pd
import seaborn as sns

regions_dict = {
    "Amygdala": {
        "regions": [
            "Central-nucleus",
            "Cortical-nucleus",
            "Cortical-nucleus",
            "Medial-nucleus"
        ],
        "hemis": ["rh", "rh", "lh", "rh"]
    },
    "Thalamic-nuclei": {
        "regions": ["Right-LP", "Left-LD"],
        "hemis": ["rh", "lh"]
    }
}

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
plots_dir = tinception_dir / "plots" / "rois"
df = pd.read_csv(tinception_dir / "subcortical_roi" / "master.csv")

df_plots = []
for atlas in regions_dict.keys():
    sub_structures = regions_dict[atlas]["regions"]
    hemis = regions_dict[atlas]["hemis"]
    
    for sub_structure, hemi in zip(sub_structures, hemis):
        if atlas == "Thalamic-nuclei":
            df_plot = df.query(f'structure == "{atlas}" & region == "{sub_structure}"').copy()
        else:
            df_plot = df.query(f'structure == "{atlas}" & region == "{sub_structure}" & hemi == "{hemi}"').copy()
        
        # KEY FIX: Create a unique label for the facet (e.g., "rh: Cortical-nucleus")
        df_plot["plot_label"] = f"{hemi}: {sub_structure}"
        df_plots.append(df_plot)

df_plot = pd.concat(df_plots, axis=0)
df_plot["group"] = df_plot["group"].map({"CO": "Control", "TI": "Tinnitus"})

pal = ['#1f77b4', '#d62728']
order = ["Control", "Tinnitus"]

# Define order based on the unique labels created above
for atlas in list(regions_dict.keys())[:1]:
    df_atlas = df_plot.query(f'structure == "{atlas}"')
    
    # Generate the order list dynamically from your dictionary
    col_order = [f"{h}: {r}" for r, h in zip(regions_dict[atlas]["regions"], regions_dict[atlas]["hemis"])]

    g = sns.FacetGrid(
                data=df_atlas,
                col="plot_label",
                col_wrap=2,
                col_order=col_order,
                height=1.8,
                sharex=False,
                aspect=3
                )
    g.map_dataframe(
                sns.stripplot,
                x="volume",
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
                x="volume",
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
    
    g.savefig(
        plots_dir / f"boxplots_{atlas}.pdf",
        format="pdf",
        dpi=300,
        bbox_inches="tight"
        )
                