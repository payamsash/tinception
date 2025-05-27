import os
import shutil
import glob
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis


def estimate_norm_1(roi_fname):

    ## load files
    norm_dir = Path.cwd().parent / "data" / "norm"
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / roi_fname
    df = pd.read_csv(fname, index_col=0)

    ## add covariates
    fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
    df_cov = pd.read_csv(fname_cov, index_col=0, delimiter="\t")
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["Unnamed: 1"])
    df = df.query('subjects == @subjects')

    ## add covariates
    df["sex"] = df_cov["Unnamed: 3"].values
    df["age"] = df_cov["Unnamed: 4"].values
    df["TIV"] = df_cov["Unnamed: 5"].values
    df["site"] = df_cov["Unnamed: 6"].values
    df["group"] = df_cov["Unnamed: 2"].values

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

    regions = df.columns[1:-5]
    df.dropna(inplace=True)

    ## train–test split
    df_features = df[regions]
    df_cov = df[["age", "sex", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)
    regions = regions.to_list()

    X_train, X_test, y_train, y_test = train_test_split(
                                                        df_cov,
                                                        df_features,
                                                        stratify=df["site"],
                                                        test_size=0.2,
                                                        random_state=42
                                                        )
    X_train.reset_index(drop=True, inplace=True) 
    X_test.reset_index(drop=True, inplace=True) 
    y_train.reset_index(drop=True, inplace=True)
    y_test.reset_index(drop=True, inplace=True)

    k_idx = X_test.index[X_test['site_0.0'] == 1].to_list() 
    c_idx = X_test.index[X_test['site_1.0'] == 1].to_list()
    ta_idx = X_test.index[X_test['site_2.0'] == 1].to_list()
    ti_idx = X_test.index[X_test['site_3.0'] == 1].to_list()
    n_idx = X_test.index[X_test['site_4.0'] == 1].to_list()

    sites = [k_idx, c_idx, ta_idx, ti_idx, n_idx]
    site_names = ["karolinska", "cogtail", "talaska", "tinspect", "neuropren"]

    ## set up output directories
    roi_models_dir = norm_dir / "roi_models" / roi_fname[:-4]
    os.makedirs(roi_models_dir, exist_ok=True)

    for c in y_train.columns:
        y_train[c].to_csv(roi_models_dir / f'resp_tr_{c}.txt', header=False, index=False)
        X_train.to_csv(roi_models_dir / 'cov_tr.txt', sep = '\t', header=False, index=False)
        y_train.to_csv(roi_models_dir / 'resp_tr.txt', sep = '\t', header=False, index=False)
    for c in y_test.columns:
        y_test[c].to_csv(roi_models_dir / f'resp_te_{c}.txt', header=False, index=False)
        X_test.to_csv(roi_models_dir / 'cov_te.txt', sep = '\t', header=False, index=False)
        y_test.to_csv(roi_models_dir / 'resp_te.txt', sep = '\t', header=False, index=False)

    for roi in regions:
        resp_tr_file = roi_models_dir / f"resp_tr_{roi}.txt"
        resp_te_file = roi_models_dir / f"resp_te_{roi}.txt"
        if not os.path.exists(resp_tr_file):
            raise FileExistsError(f"{resp_tr_file} does not exist.")
        
        roi_dir = os.path.join(roi_models_dir, roi)
        os.makedirs(roi_dir, exist_ok=True)

        shutil.copy(resp_tr_file, os.path.join(roi_dir, "resp_tr.txt"))
        shutil.copy(resp_te_file, os.path.join(roi_dir, "resp_te.txt"))
        shutil.copy(roi_models_dir / "cov_tr.txt", os.path.join(roi_dir, "cov_tr.txt"))
        shutil.copy(roi_models_dir / "cov_te.txt", os.path.join(roi_dir, "cov_te.txt"))

    [os.remove(file) for file in glob.glob(str(roi_models_dir/"resp_*.txt"))]
    [os.remove(file) for file in glob.glob(str(roi_models_dir/"cov_*.txt"))]

    ## basis expansion using B-splines
    xmin, xmax = (10, 95) 
    B = create_bspline_basis(xmin, xmax)
    for roi in regions:
        roi_dir = roi_models_dir / roi
        os.makedirs(roi_dir / 'blr', exist_ok=True)
        
        X_tr = np.loadtxt(roi_dir / 'cov_tr.txt')
        X_te = np.loadtxt(roi_dir / 'cov_te.txt')

        X_tr = np.concatenate((X_tr, np.ones((X_tr.shape[0], 1))), axis=1)
        X_te = np.concatenate((X_te, np.ones((X_te.shape[0], 1))), axis=1)
        np.savetxt(roi_dir / 'cov_int_tr.txt', X_tr)
        np.savetxt(roi_dir / 'cov_int_te.txt', X_te)

        Phi = np.array([B(i) for i in X_tr[:, 0]])
        Phis = np.array([B(i) for i in X_te[:, 0]])
        X_tr = np.concatenate((X_tr, Phi), axis=1)
        X_te = np.concatenate((X_te, Phis), axis=1)
        np.savetxt(roi_dir / 'cov_bspline_tr.txt', X_tr)
        np.savetxt(roi_dir / 'cov_bspline_te.txt', X_te)


    ## estimate normative model
    blr_metrics = pd.DataFrame(columns=['ROI', 'MSLL', 'EV', 'SMSE', 'RMSE', 'Rho'])
    blr_site_metrics = pd.DataFrame(columns=['ROI', 'site', 'y_mean', 'y_var', 'yhat_mean',
                                                'yhat_var', 'MSLL', 'EV', 'SMSE', 'RMSE', 'Rho'])
    df_z = pd.DataFrame(columns=['ROI', 'z_value', 'site'])
    
    for roi in regions:
        # configure the covariates to use
        roi_dir = roi_models_dir / roi
        cov_file_tr = roi_dir / 'cov_bspline_tr.txt'
        cov_file_te = roi_dir / 'cov_bspline_te.txt'

        # load train & test response files
        resp_file_tr = roi_dir / 'resp_tr.txt'
        resp_file_te = roi_dir / 'resp_te.txt'

        # run a basic model
        yhat_te, s2_te, nm, z_scores, metrics_te = estimate(
                                                            str(cov_file_tr),
                                                            str(resp_file_tr),
                                                            testresp=str(resp_file_te),
                                                            testcov=str(cov_file_te),
                                                            alg="blr",
                                                            optimizer="powell",
                                                            savemodel=False,
                                                            saveoutput=False,
                                                            standardize=False
                                                            )
        # create a dataframe for z_score
        for num, site in enumerate(sites):
            zs = np.array(z_scores[site])[:, 0]
            site_name = site_names[num]
            for z in zs:
                df_z.loc[len(df_z)] = [roi, z, site_name]

        # save metrics
        blr_metrics.loc[len(blr_metrics)] = [
                                            roi,
                                            metrics_te['MSLL'][0],
                                            metrics_te['EXPV'][0],
                                            metrics_te['SMSE'][0],
                                            metrics_te['RMSE'][0],
                                            metrics_te['Rho'][0]
                                            ]
        # Compute metrics per site in test set
        X_te = np.loadtxt(cov_file_te)
        y_te = np.loadtxt(resp_file_te)

        for num, site in enumerate(sites):
            y_mean_te_site = np.array([[np.mean(y_te[site])]])
            y_var_te_site = np.array([[np.var(y_te[site])]])
            yhat_mean_te_site = np.array([[np.mean(yhat_te[site])]])
            yhat_var_te_site = np.array([[np.var(yhat_te[site])]])

            if y_te.ndim != 2: y_te = np.expand_dims(y_te, axis=-1)
            if yhat_te.ndim != 2: yhat_te = np.expand_dims(yhat_te, axis=-1)

            metrics_te_site = evaluate(y_te[site], yhat_te[site], s2_te[site], y_mean_te_site, y_var_te_site)
            site_name = site_names[num]
            blr_site_metrics.loc[len(blr_site_metrics)] = [
                                                            roi,
                                                            site_names[num],
                                                            y_mean_te_site[0][0],
                                                            y_var_te_site[0][0],
                                                            yhat_mean_te_site[0][0],
                                                            yhat_var_te_site[0][0],
                                                            metrics_te_site ['MSLL'][0],
                                                            metrics_te_site ['EXPV'][0],
                                                            metrics_te_site ['SMSE'][0],
                                                            metrics_te_site ['RMSE'][0],
                                                            metrics_te_site ['Rho'][0]
                                                            ]
    
    ## saving
    saving_dir = norm_dir / "results" / roi_fname[:-4]
    os.makedirs(saving_dir, exist_ok=True)
    blr_metrics.to_csv(saving_dir / f"blr_metrics.csv")
    blr_site_metrics.to_csv(saving_dir / f"blr_site_metrics.csv")
    df_z.to_csv(saving_dir / f"z_scores.csv")



if __name__ == "__main__":
    
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    roi_fnames = [fname for fname in sorted(os.listdir(roi_dir)) if fname.endswith(".csv")]

    for roi_fname in roi_fnames:
        estimate_norm_1(roi_fname)