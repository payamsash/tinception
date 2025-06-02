import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis


def estimate_norm_on_rois(roi_fname, method, test_size=0.2, random_state=42):    

    ###################### prepare data

    ## load files done
    norm_dir = Path.cwd().parent / "data" / "norm"
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / roi_fname
    df = pd.read_csv(fname, index_col=0)
    site_names = ["karolinska", "cogtail", "tinspect", "neuropren"] # no talaska

    ## add covariates
    fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
    df_cov = pd.read_csv(fname_cov, delimiter="\t")
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["Unnamed: 1"])
    df = df.query('subjects == @subjects')

    ## add covariates
    df["sex"] = df_cov["Unnamed: 3"].values
    df["age"] = df_cov["Unnamed: 4"].values / 100
    df["TIV"] = df_cov["Unnamed: 5"].values
    df["site"] = df_cov["Unnamed: 6"].values
    df["group"] = df_cov["Unnamed: 2"].values

    df = df[df["site"] != 2] # drop talaske
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
    df = df.query('group == "CO"')
    df.reset_index(inplace=True, drop=True)

    df_features = df[rois]
    df_cov = df[["age_demean", "sex", "TIV", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)

    df_cov = df[["age_demean", "sex", "TIV", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)

    rois = rois.to_list()
    roi_models_dir = norm_dir / "roi_models" / f"{roi_fname[:-4]}"
    os.makedirs(roi_models_dir, exist_ok=True)


    X_tr, X_te, y_train, y_test, z_train, z_test = train_test_split(
                                                        df_cov,
                                                        df_features,
                                                        df,
                                                        stratify=df["site"],
                                                        test_size=test_size,
                                                        random_state=random_state
                                                        )

    z_test.reset_index(drop="index", inplace=True)
    sites_dict = dict(zip([0, 1, 3, 4], ["karolinska", "cogtail", "tinspect", "neuropren"]))

    np.savetxt(roi_models_dir / 'cov_tr.txt', X_tr)
    np.savetxt(roi_models_dir / 'cov_te.txt', X_te)

    ###################### bsplines

    X_tr = np.concatenate((X_tr, np.ones((X_tr.shape[0], 1))), axis=1)
    X_te = np.concatenate((X_te, np.ones((X_te.shape[0], 1))), axis=1)

    B = create_bspline_basis(xmin=10, xmax=95, p=3, nknots=5)
    Phi = np.array([B(i) for i in X_tr[:, 0]])
    Phis = np.array([B(i) for i in X_te[:, 0]])
    X_tr = np.concatenate((X_tr, Phi), axis=1)
    X_te = np.concatenate((X_te, Phis), axis=1)
    np.savetxt(roi_models_dir / 'cov_bspline_tr.txt', X_tr)
    np.savetxt(roi_models_dir / 'cov_bspline_te.txt', X_te)

    np.savetxt(roi_models_dir / 'batch_tr.txt', np.array(z_train["site"].values))
    np.savetxt(roi_models_dir / 'batch_te.txt', np.array(z_test["site"].values))

    for roi in df_features.columns:
        y_train[roi].to_csv(roi_models_dir / f'resp_tr_{roi}.txt', header=False, index=False)
        y_test[roi].to_csv(roi_models_dir / f'resp_te_{roi}.txt', header=False, index=False)

    os.makedirs(roi_models_dir / "results", exist_ok=True)

    ###################### estimate norm brain
    df_metrics_list = []
    df_site_metrics = pd.DataFrame(columns = ['roi', 'site', 'MSLL', 'EXPV', 'SMSE', 'RMSE', 'Rho'])

    for roi in rois:
        resp_file_tr = str(roi_models_dir / f'resp_tr_{roi}.txt')
        resp_file_te = str(roi_models_dir / f'resp_te_{roi}.txt')
        
        # run a basic model
        if method == "blr_spline":
            cov_file_tr = str(roi_models_dir / 'cov_bspline_tr.txt')
            cov_file_te = str(roi_models_dir / 'cov_bspline_te.txt')
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
        if method == "blr":
            cov_file_tr = str(roi_models_dir / 'cov_tr.txt')
            cov_file_te = str(roi_models_dir / 'cov_te.txt')
            yhat_te, s2_te, nm, z_scores, metric_te = estimate(
                                                                cov_file_tr,
                                                                resp_file_tr,
                                                                testresp=resp_file_te,
                                                                testcov=cov_file_te,
                                                                alg="hbr",
                                                                binary=True,
                                                                # optimizer="powell",
                                                                savemodel=False,
                                                                saveoutput=False,
                                                                standardize=False
                                                                )
        

        if method == "hbr":
            trbefile = str(roi_models_dir / 'batch_tr.txt')
            tsbefile = str(roi_models_dir / 'batch_te.txt')
            yhat_te, s2_te, nm, z_scores, metric_te = estimate(
                                                                covfile=cov_file_tr,
                                                                respfile=resp_file_tr,
                                                                tsbefile=tsbefile,
                                                                trbefile=trbefile,
                                                                inscaler="standardize",
                                                                outscaler="standardize",
                                                                linear_mu="True",
                                                                random_intercept_mu="True",
                                                                centered_intercept_mu="True",
                                                                alg="hbr",
                                                                binary=True,
                                                                testcov=cov_file_te,
                                                                testresp=resp_file_te,
                                                                savemodel=False,
                                                                nuts_sampler="nutpie",
                                                                )

        df_norm = z_test[["subjects", roi, "age", "sex", "site", "group"]]
        df_norm["yhat_te"] = yhat_te
        df_norm["s2_te"] = s2_te
        df_norm["z_scores"] = z_scores
        df_norm.to_csv(roi_models_dir / "results" / f"{roi}.csv")

        metric_te_flattened = {k: v.flatten()[0] for k, v in metric_te.items()}
        df_metric = pd.DataFrame([metric_te_flattened])
        df_metric["ROI"] = roi
        df_metrics_list.append(df_metric)


        for key, val in sites_dict.items():
            z_sub = z_test[z_test["site"] == key]
            idxs = np.array(z_test[z_test["site"] == key].index)
            y_mean_te_site = np.array([z_sub[roi].values.mean()])
            y_var_te_site = np.array([z_sub[roi].values.var()])


            metrics_te_site = evaluate(np.array([z_sub[roi].values]).T,
                                        yhat_te[idxs], s2_te[idxs], y_mean_te_site, y_var_te_site)
            
            df_site_metrics.loc[len(df_site_metrics)] = [
                                                            roi,
                                                            sites_dict[int(key)],
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
    methods = ["hbr"]
    for roi_fname in ["aparc_thickness_lh.csv"]:
        for method in methods:
            estimate_norm_on_rois(roi_fname, method, test_size=0.2, random_state=42)
