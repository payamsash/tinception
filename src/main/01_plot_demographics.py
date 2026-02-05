
from pathlib import Path
import pandas as pd
import seaborn as sns

## read the file
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
demographics_dir = Path("../../master_files")
saving_dir = tinception_dir / "plots" / "demographics"
saving_dir.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(demographics_dir / "master.csv", dtype={"subject_ID": str})

## fix values
hues = ["group", "sex"]
df_plot = df[["subject_ID", "site", "sex", "age", "group"]]

df_plot["sex"] = df_plot["sex"].astype(int).map({0: "Female", 1: "Male"})
df_plot["group"] = df_plot["group"].map({"CO": "Control", "TI": "Tinnitus"})

site_names = df_plot["site"].unique()
print(site_names)
pal = [
        sns.cubehelix_palette(3, rot=-.1, light=.2).as_hex()[1],
        sns.color_palette("ch:s=-.2,r=.8", as_cmap=False).as_hex()[2]
]

## start plotting
xlim = [-5, 100]
bw_adjust = 1
for hue in hues:
        g = sns.FacetGrid(
                    df_plot,
                    row="site",
                    hue=hue,
                    aspect=3.5,
                    height=1.6,
                    palette=pal,
                    row_order=site_names,
                    xlim=xlim,
                    sharex=True,
                    sharey=False
                    )
        g.map(
                sns.kdeplot,
                "age",
                bw_adjust=bw_adjust,
                clip_on=False,
                clip=xlim,
                fill=True,
                alpha=0.7,
                linewidth=1.5
                )
        g.map(
                sns.kdeplot,
                "age",
                clip_on=False,
                color="w",
                clip=xlim,
                lw=1.5,
                bw_adjust=bw_adjust
                )
        
        g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)
        g.figure.subplots_adjust(hspace=.15, top=0.72)
        g.set_titles("")
        g.add_legend(title="")
        g.set(yticks=[], ylabel="", xlabel=r"Age")
        g.despine(bottom=True, left=True)
        g.figure.savefig(saving_dir / f"{hue}_distribution.pdf", 
                        format="pdf",       
                        dpi=300,            
                        bbox_inches="tight"
                        )