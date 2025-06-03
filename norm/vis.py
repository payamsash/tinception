import os
from pathlib import Path
import pandas as pd
import seaborn as sns
import joypy

################# plot metrics of norm brain

directory = Path.cwd().parent / "data" / "norm" / "roi_models"
parcs = sorted(os.listdir(directory))[1:]

dfs_list = []
for parc in parcs:
    fname = directory / parc / "metrics_on_controls.csv"
    df = pd.read_csv(fname, index_col=0)
    df["parc"] = parc
    dfs_list.append(df)

df_all = pd.concat(dfs_list)
cmap1 = sns.cubehelix_palette(10, rot=-.25, light=.7, as_cmap=True)
cmap2 = sns.cubehelix_palette(10, dark=.15, light=.75, as_cmap=True)
kwargs = {
            "background": "k",
            "bins": 100,
            "overlap": 1,
            "linecolor": "k",
            }
    
fig, axs = joypy.joyplot(df_all[["SMSE", "parc"]], by="parc", x_range=[0, 1.7], colormap=cmap1, **kwargs)
fig, axs = joypy.joyplot(df_all[["MSLL", "parc"]], by="parc", x_range=[-1.5, 0.3], colormap=cmap2, **kwargs)


