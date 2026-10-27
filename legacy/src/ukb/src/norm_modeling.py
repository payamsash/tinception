from pathlib import Path
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
from pcntoolkit import (
    NormData,
    BLR,
    BsplineBasisFunction,
    NormativeModel,
    plot_centiles_advanced
)
import seaborn as sns



def run_norm_model(atlas_name):

    fname = f"./results/atlases/{atlas_name}.csv"
    subprocess.run(["dx", "download", fname], check=True)
    models_dir = Path.cwd() / "norm_models"
    df = pd.read_csv(f"{atlas_name}.csv")
    
    print(df.columns)

    match_fname = f"./results/ukb_vol_harmonized.csv"
    subprocess.run(["dx", "download", match_fname], check=True)
    df_matched = pd.read_csv("ukb_vol_harmonized.csv")

    df = df.merge(
                df_matched[["subject_id", "SITE"]],
                on="subject_id",
                how="inner"
                )
    df.drop(columns='Unnamed: 0', inplace=True)

    ## first run normative modeling
    demographic_cols = ["age", "sex", "srt"]
    response_cols = df.columns[7:-1].tolist()

    kwargs = {
                "covariates": demographic_cols,
                "batch_effects": ["SITE"],
                "response_vars": response_cols, 
                "subject_ids": "subject_id"
                }

    norm_train_all = NormData.from_dataframe(
                                            name="train",
                                            dataframe=df.query('tin_status == "Control"'),
                                            **kwargs
                                            )
    norm_test_tinnitus = NormData.from_dataframe(
                                            name="test",
                                            dataframe=df.query('tin_status == "Tinnitus"'),
                                            **kwargs
                                            )

    cfg = {
            "save_dir": models_dir / atlas_name,
            "train": norm_train_all,
            "test": norm_test_tinnitus,
            "savemodel": True,
            }

    template_blr = BLR(
                    name="payam_blr",
                    basis_function_mean=BsplineBasisFunction(degree=3, nknots=5),
                    fixed_effect=True,
                    heteroskedastic=True,
                    warp_name="warpsinharcsinh"
                    )

    model = NormativeModel(
                    template_regression_model=template_blr,
                    savemodel=cfg["savemodel"],
                    evaluate_model=True,
                    saveresults=True,
                    saveplots=False,
                    save_dir=str(cfg["save_dir"]),
                    inscaler="standardize",
                    outscaler="standardize",
                    )
    model.fit_predict(cfg["train"], cfg["test"])
    del model
    del template_blr


    
def plot_centiles(atlas_name, rois):
    
    

    models_dir = Path("./norm_models")
    df = pd.read_csv(f"{atlas_name}.csv")
    
    df_matched = pd.read_csv("ukb_vol_harmonized.csv")
    df = df.merge(
                df_matched[["subject_id", "SITE"]],
                on="subject_id",
                how="inner"
                )
    df.drop(columns='Unnamed: 0', inplace=True)

    ## first run normative modeling
    demographic_cols = ["age", "sex", "srt"]
    response_cols = df.columns[7:-1].tolist()

    kwargs = {
                "covariates": demographic_cols,
                "batch_effects": ["SITE"],
                "response_vars": response_cols, 
                "subject_ids": "subject_id"
                }

    norm_train_all = NormData.from_dataframe(
                                            name="train",
                                            dataframe=df.query('tin_status == "Control"'),
                                            **kwargs
                                            )
    norm_test_tinnitus = NormData.from_dataframe(
                                            name="test",
                                            dataframe=df.query('tin_status == "Tinnitus"'),
                                            **kwargs
                                            )

    model_dir = models_dir / atlas_name
    model = NormativeModel.load(str(model_dir))

    data_mode = "train"

    if data_mode == "train":
        scatter_data = norm_train_all
    if data_mode == "test":
        scatter_data = norm_test_tinnitus

    covars = ["age", "srt"]
    for covar in covars:
        figs, re_vars = plot_centiles_advanced(
                                model,
                                centiles=[0.05, 0.5, 0.95],
                                covariate=covar,
                                scatter_data=scatter_data,
                                show_other_data=False,
                                harmonize_data=True,
                                show_yhat=True,
                                brain_labels=rois
                                )
        for fig, roi in zip(figs, rois):
            figname = f"{roi}_{covar}_centile.pdf"
            fig.savefig(
                        figname,
                        format="pdf",
                        dpi=300,
                        bbox_inches="tight"
                        )
            subprocess.run(["dx", "upload", figname, "--dest", "/figures/atlases/"], check=True)
            
    
    
def plot_extremes(atlas_name):
    
    threshold = 1.96
    models_dir = Path.cwd() / "norm_models"
    fname = models_dir / atlas_name / "results" / f"Z_test.csv"
    df_dev = pd.read_csv(fname)
    df_dev = df_dev[df_dev.columns[2:]]

    extreme_counts = (
        ((df_dev > threshold) | (df_dev < -threshold))
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

    ax.set_yticks([10, 50, 90])
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.spines.left.set_bounds((10, 90))
    ax.spines["left"].set_position(("outward", 5))
    ax.tick_params(axis="y", which="minor", length=8, width=5.5)
    
    figname = f"{atlas_name}_extreme_count.pdf"
    fig.savefig(
                figname,
                format="pdf",
                dpi=300,
                bbox_inches="tight"
                )
    subprocess.run(["dx", "upload", figname, "--dest", "/figures/atlases/"], check=True)
    
    
    

if __name__ == "__main__":
    atlas_name = "thalamic_nuclei"
    run_model = False
    
    if run_model:
        run_norm_model(atlas_name)
    
    if atlas_name == "hippo_subfields":
        rois = [
            "Volume_of_presubiculum_body_right_hemisphere",
            "Volume_of_GC_ML_DG_body_right_hemisphere"
                ]
    if atlas_name == "amygdalar_nuclei":
        rois = [
                "Volume_of_Basal_nucleus_left_hemisphere"
               ]
    if atlas_name == "thalamic_nuclei":
        rois = [
               "Volume_of_CeM_left_hemisphere",
               "Volume_of_CeM_right_hemisphere",
               "Volume_of_MDm_left_hemisphere",
               "Volume_of_MDm_right_hemisphere",
               "Volume_of_MVRe_left_hemisphere",
               "Volume_of_MVRe_right_hemisphere"
               ]
        
    plot_centiles(atlas_name, rois=rois)
    plot_extremes(atlas_name)