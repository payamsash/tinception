import os
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import perform_roi_stats


def plot_roi_stats(roi, roi_sel):
    
    ## select what to read
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / f"{roi}.csv"
    df = pd.read_csv(fname, index_col=0)

    ## select subjects
    fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
    df_cov = pd.read_csv(fname_cov, index_col=0, delimiter="\t")
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["Unnamed: 1"])
    df = df.query('subjects == @subjects')
    df["group"] = df_cov["Unnamed: 2"].values

    ## some renaming and fixes
    if "MV(Re)" in df.columns.to_list():
        df = df.rename(columns={"MV(Re)": "MV"})
    df.columns = df.columns.str.replace('-', '_')

    palette_color = ['#1f77b4', '#d62728']
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    sns.boxplot(data=df, x="group", y=roi_sel,
                palette=palette_color, fill=False, width=0.6,
                gap=0, linewidth=1.8, ax=ax)
    sns.stripplot(data=df, x="group", y=roi_sel,
                palette=palette_color, linewidth=0, size=3.5,
                edgecolor=None, ax=ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(roi_sel)

    plt.rcParams['figure.facecolor'] = 'black'
    plt.rcParams['axes.facecolor'] = 'black'
    plt.rcParams['axes.edgecolor'] = 'white'
    plt.rcParams['xtick.color'] = 'white'
    plt.rcParams['ytick.color'] = 'white'
    plt.rcParams['axes.labelcolor'] = 'white'
    plt.rcParams['text.color'] = 'white'
    plt.rcParams['legend.facecolor'] = 'black'
    plt.rcParams['legend.edgecolor'] = 'white'

    return fig


roi_dir = Path.cwd().parent / "data" / "roi_stats"
rois = [item[:-4] for item in sorted(os.listdir(roi_dir)) if item.endswith(".csv")]
rois.remove("aseg")

for roi in rois:
    print(f"working on {roi} ...")
    df = perform_roi_stats(roi, method="fdr_bh", add_site=True)
    df.to_csv(Path.cwd().parent / "data" / "roi_stats" / "results" / "stats" / f"{roi}.csv")
    df_sub = df.query('reject == True')

    ## plotting
    for roi_sel in df_sub["Region"].values:
        fig = plot_roi_stats(roi, roi_sel)
        if roi.endswith("rh"):
            fig.savefig(Path.cwd().parent / "data" / "roi_stats" / "results" / "figures" / f"{roi_sel}_rh.pdf")
        elif roi.endswith("lh"):
            fig.savefig(Path.cwd().parent / "data" / "roi_stats" / "results" / "figures" / f"{roi_sel}_lh.pdf")
        else:
            fig.savefig(Path.cwd().parent / "data" / "roi_stats" / "results" / "figures" / f"{roi_sel}.pdf")


if __name__ == "__main__":
    plot_roi_stats()