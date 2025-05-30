import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd

from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis


def estimate_norm_on_rois(roi_fname):    

    ###################### prepare data

    ## load files done
    norm_dir = Path.cwd().parent / "data" / "norm"
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / roi_fname
    df = pd.read_csv(fname, index_col=0)
    site_names = ["karolinska", "cogtail", "talaska", "tinspect", "neuropren"]

    ## add covariates
    fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
    df_cov = pd.read_csv(fname_cov, delimiter="\t")
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["Unnamed: 1"])
    df = df.query('subjects == @subjects')

    ## add covariates
    df["sex"] = df_cov["Unnamed: 3"].values
    df["age"] = df_cov["Unnamed: 4"].values
    df["TIV"] = df_cov["Unnamed: 5"].values
    df["site"] = df_cov["Unnamed: 6"].values
    df["group"] = df_cov["Unnamed: 2"].values

    df.reset_index(inplace=True, drop=True)

    ## add variable to model site/scanner effects
    if "MV(Re)" in df.columns.to_list():
        df = df.rename(columns={"MV(Re)": "MV"})
    df.columns = df.columns.str.replace('-', '_')
    df.columns = df.columns.str.replace('(', '')
    df.columns = df.columns.str.replace(')', '')

    df['group'] = df['group'].astype('category')
    df['sex'] = df['sex'].astype('category')
    df['site'] = df['site'].astype('category')

    # fix this part
    df.drop(columns=[
                    "lh_WhiteSurfArea_area",
                    "rh_WhiteSurfArea_area",
                    "BrainSegVolNotVent",
                    "eTIV",
                    "lh_MeanThickness_thickness",
                    "rh_WhiteSurfArea_area",
                    "5th_Ventricle",
                    "Left_WM_hypointensities",
                    "Right_WM_hypointensities",
                    "non_WM_hypointensities",
                    "Left_non_WM_hypointensities",
                    "Right_non_WM_hypointensities",
                    "anterior_horn_of_lateral_ventricle"
                    ], errors="ignore", inplace=True)

    rois = df.columns[1:-5]
    df.dropna(inplace=True)
    df["age_demean"] = df["age"] - df["age"].mean()
    df_controls = df.query('group == "CO"')
    df_controls.reset_index(inplace=True, drop=True)

    ## get site indexes
    site_indices = {}
    for site_val in df["site"].unique():
        site_indices[site_val] = df.index[df["site"] == site_val].tolist()

    ## create feature dfs
    df_features = df[rois]
    df_con_features = df_controls[rois]

    df_cov = df[["age_demean", "sex", "TIV", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)

    df_con_cov = df_controls[["age_demean", "sex", "TIV", "site"]]
    df_con_cov = pd.get_dummies(df_con_cov, columns=['site'], dtype=float)

    rois = rois.to_list()
    roi_models_dir = norm_dir / "roi_models" / roi_fname[:-4]
    os.makedirs(roi_models_dir, exist_ok=True)

    ###################### bsplines

    xmin, xmax = (10, 95) 
    B = create_bspline_basis(xmin, xmax)

    X_tr = df_con_cov.to_numpy()
    X_te = df_cov.to_numpy()
    X_tr = np.concatenate((X_tr, np.ones((X_tr.shape[0], 1))), axis=1)
    X_te = np.concatenate((X_te, np.ones((X_te.shape[0], 1))), axis=1)

    Phi = np.array([B(i) for i in X_tr[:, 0]])
    Phis = np.array([B(i) for i in X_te[:, 0]])
    X_tr = np.concatenate((X_tr, Phi), axis=1)
    X_te = np.concatenate((X_te, Phis), axis=1)
    np.savetxt(roi_models_dir / 'cov_bspline_tr.txt', X_tr)
    np.savetxt(roi_models_dir / 'cov_bspline_te.txt', X_te)

    for roi in df_features.columns:
        df_con_features[roi].to_csv(roi_models_dir / f'resp_tr_{roi}.txt', header=False, index=False)
        df_features[roi].to_csv(roi_models_dir / f'resp_te_{roi}.txt', header=False, index=False)

    cov_file_tr = str(roi_models_dir / 'cov_bspline_tr.txt')
    cov_file_te = str(roi_models_dir / 'cov_bspline_te.txt')

    os.makedirs(roi_models_dir / "results", exist_ok=True)

    ###################### estimate norm brain
    df_metrics_list = []
    df_site_metrics = pd.DataFrame(columns = ['roi', 'site', 'MSLL', 'EXPV', 'SMSE', 'RMSE', 'Rho'])

    for roi in rois:
        resp_file_tr = str(roi_models_dir / f'resp_tr_{roi}.txt')
        resp_file_te = str(roi_models_dir / f'resp_te_{roi}.txt')
        
        # run a basic model
        yhat_te, s2_te, nm, z_scores, metric_te = estimate(
                                                            cov_file_tr,
                                                            resp_file_tr,
                                                            testresp=resp_file_te,
                                                            testcov=cov_file_te,
                                                            alg="blr",
                                                            optimizer="powell",
                                                            savemodel=False,
                                                            saveoutput=False,
                                                            standardize=False
                                                            )
        
        df_norm = df[["subjects", roi, "age", "sex", "site", "group"]]
        df_norm["yhat_te"] = yhat_te
        df_norm["s2_te"] = s2_te
        df_norm["z_scores"] = z_scores
        df_norm.to_csv(roi_models_dir / "results" / f"{roi}.csv")

        metric_te_flattened = {k: v.flatten()[0] for k, v in metric_te.items()}
        df_metric = pd.DataFrame([metric_te_flattened])
        df_metric["ROI"] = roi
        df_metrics_list.append(df_metric)


        for key, val in site_indices.items():
            y_mean_te_site = np.array([df_features[roi].iloc[val].values.mean()])
            y_var_te_site = np.array([df_features[roi].iloc[val].values.var()])


            metrics_te_site = evaluate(np.array([df_features[roi].iloc[val].values]).T,
                                        yhat_te[val], s2_te[val], y_mean_te_site, y_var_te_site)
            
            df_site_metrics.loc[len(df_site_metrics)] = [
                                                            roi,
                                                            site_names[int(key)],
                                                            metrics_te_site['MSLL'][0],
                                                            metrics_te_site['EXPV'][0],
                                                            metrics_te_site['SMSE'][0],
                                                            metrics_te_site['RMSE'][0],
                                                            metrics_te_site['Rho'][0]
                                                            ]

    df_metrics = pd.concat(df_metrics_list)
    df_metrics.to_csv(roi_models_dir / "results" / "metrics.csv")
    df_site_metrics.to_csv(roi_models_dir / "results" / "site_metrics.csv")

    ## cleaning
    pattern = os.path.join(roi_models_dir, 'resp*.txt')
    files_to_delete = glob.glob(pattern)
    for file_path in files_to_delete:
        os.remove(file_path)


if __name__ == "__main__":
    
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    roi_fnames = [fname for fname in sorted(os.listdir(roi_dir)) if fname.endswith(".csv")]

    for roi_fname in roi_fnames:
        estimate_norm_on_rois(roi_fname)