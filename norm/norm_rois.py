import os
import glob
from pathlib import Path
import shutil
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
import pcntoolkit as ptk
import numpy as np
import pickle


def estimate_hbr_norm_rois(roi_fname):
    
    ######## part 1
    os.chdir("/Users/payamsadeghishabestari/tinception/norm")
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
    df["age"] = df_cov["Unnamed: 4"].values / 100 # normalize
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

    # drop to avoid problem
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
    df = df.query('group == "CO"')
    df.reset_index(inplace=True, drop=True)

    df_features = df[rois]
    df_cov = df[["age", "sex", "TIV", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)

    rois = rois.to_list()
    roi_models_dir = norm_dir / "roi_models" / f"{roi_fname[:-4]}"
    os.makedirs(roi_models_dir, exist_ok=True)

    ######## part 2
    X_tr, X_te, y_train, y_test, z_train, z_test = \
                                    train_test_split(
                                                    df_cov,
                                                    df_features,
                                                    df,
                                                    stratify=df["site"],
                                                    test_size=0.2,
                                                    random_state=42
                                                    )

    z_test.reset_index(drop="index", inplace=True)
    sites_dict = dict(zip([0, 1, 3, 4], site_names))

    X_tr = X_tr[["age", "sex", "TIV"]]
    X_te = df_cov[["age", "sex", "TIV"]]
    y_test = df_features
    be_tr = z_train["site"].astype(int)
    be_te = df["site"].astype(int)

    datas = [X_tr, y_train, be_tr, X_te, y_test, be_te]
    file_names = ["X_train", "Y_train", "trbefile", "X_test", "Y_test", "tsbefile"]
    
    os.chdir(roi_models_dir)
    for data, file_name in zip(datas, file_names):
        with open(f"{file_name}.pkl", "wb") as file:
            pickle.dump(pd.DataFrame(data), file)

    respfile = os.path.join(roi_models_dir, "Y_train.pkl")
    covfile = os.path.join(roi_models_dir, "X_train.pkl")  
    testrespfile_path = os.path.join(roi_models_dir, "Y_test.pkl") 
    testcovfile_path = os.path.join(roi_models_dir, "X_test.pkl") 
    trbefile = os.path.join(roi_models_dir, "trbefile.pkl") 
    tsbefile = os.path.join(roi_models_dir, "tsbefile.pkl")

    output_path = os.path.join(roi_models_dir, "Models/")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)


    ############### on controls only to evaluate the model
    ptk.normative.estimate(
                            covfile=covfile,
                            respfile=respfile,
                            tsbefile=tsbefile,
                            trbefile=trbefile,
                            inscaler="standardize",
                            outscaler="standardize",
                            linear_mu="True",
                            random_intercept_mu="True",
                            centered_intercept_mu="True",
                            alg="hbr",
                            binary=True,
                            output_path=output_path,
                            testcov=testcovfile_path,
                            testresp=testrespfile_path,
                            outputsuffix="_estimate",
                            savemodel=True,
                            nuts_sampler="nutpie"
                        )
    
    ## create a dataframe
    metrics = ["EXPV", "MSLL", "NLL", "pRho", "Rho", "RMSE", "SMSE", "yhat", "ys2", "Z"]
    lbs = list(df_features.columns)
    dfs_list_1, dfs_list_2 = [], []
    
    for metric in metrics:
        df_metric = pd.read_pickle(f"{metric}_estimate.pkl")
        if metric in metrics[:-3]:
            df_metric.rename(columns={df_metric.columns[0]: metric}, inplace=True)
            dfs_list_1.append(df_metric)
        else:
            lbs_up = [f"{lb}__{metric}" for lb in lbs]
            df_metric.rename(columns=dict(zip(df_metric.columns, lbs_up)), inplace=True)
            dfs_list_2.append(df_metric)

    df_1 = pd.concat(dfs_list_1, axis=1)
    df_1["lbs"] = df_features.columns.to_list()
    df_2 = pd.concat(dfs_list_2, axis=1)
    cols_to_add = ['subjects', 'sex', 'age', 'TIV', 'site', 'group']
    df_2[cols_to_add] = df[cols_to_add]
    df_2['site'] = df_2['site'].map(sites_dict)

    df_1.to_csv("metrics_on_controls.csv")
    df_2.to_csv("scores_on_controls.csv")

    ## remove unnecessary stuff
    pkl_files = glob.glob(os.path.join(roi_models_dir, '*.pkl'))
    for file in pkl_files:
        os.remove(file)




def predict_hbr_norm_rois(roi_fname):

    ######## part 1
    os.chdir("/Users/payamsadeghishabestari/tinception/norm")
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
    df["age"] = df_cov["Unnamed: 4"].values / 100 # normalize
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
    df = df.query('group == "TI"')
    df.reset_index(inplace=True, drop=True)

    df_features = df[rois]
    df_cov = df[["age", "sex", "TIV", "site"]]
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)

    rois = rois.to_list()
    roi_models_dir = norm_dir / "roi_models" / f"{roi_fname[:-4]}"
    os.makedirs(roi_models_dir, exist_ok=True)

    ######## part 2
    X_tr, X_te, y_train, y_test, z_train, z_test = \
                                    train_test_split(
                                                    df_cov,
                                                    df_features,
                                                    df,
                                                    stratify=df["site"],
                                                    test_size=0.2,
                                                    random_state=42
                                                    )

    z_test.reset_index(drop="index", inplace=True)
    sites_dict = dict(zip([0, 1, 3, 4], site_names))

    X_te = df_cov[["age", "sex", "TIV"]]
    be_te = df["site"].astype(int)

    datas = [X_te, df_features, be_te]
    file_names = ["X_test_fit", "Y_test_fit", "tsbefile_fit"]

    os.chdir(roi_models_dir)
    for data, file_name in zip(datas, file_names):
        with open(f"{file_name}.pkl", "wb") as file:
            pickle.dump(pd.DataFrame(data), file)

    testrespfile_path = os.path.join(roi_models_dir, "Y_test_fit.pkl") 
    testcovfile_path = os.path.join(roi_models_dir, "X_test_fit.pkl") 
    tsbefile = os.path.join(roi_models_dir, "tsbefile_fit.pkl")

    ############### on all subjects from sites (except talaska)
    model_path = str(roi_models_dir / "Models")
    ptk.normative.predict(
                        covfile=testcovfile_path,
                        respfile=testrespfile_path,
                        model_path=model_path,
                        outputsuffix="fit",
                        inputsuffix="estimate",
                        return_y=True,
                        alg='hbr',
                        tsbefile=tsbefile,
                        inscaler="standardize",
                        outscaler="standardize",
                        linear_mu="True",
                        random_intercept_mu="True",
                        centered_intercept_mu="True",
                        )
    
    ## create a dataframe
    metrics = ["EXPV", "MSLL", "pRho", "Rho", "RMSE", "SMSE", "yhat", "ys2", "Z"]
    lbs = list(df_features.columns)
    dfs_list_1, dfs_list_2 = [], []
    
    for metric in metrics:
        df_metric = pd.read_pickle(f"{metric}_fit.pkl")
        if metric in metrics[:-3]:
            df_metric.rename(columns={df_metric.columns[0]: metric}, inplace=True)
            dfs_list_1.append(df_metric)
        else:
            lbs_up = [f"{lb}__{metric}" for lb in lbs]
            df_metric.rename(columns=dict(zip(df_metric.columns, lbs_up)), inplace=True)
            dfs_list_2.append(df_metric)

    df_1 = pd.concat(dfs_list_1, axis=1)
    df_1["lbs"] = df_features.columns.to_list()
    df_2 = pd.concat(dfs_list_2, axis=1)
    cols_to_add = ['subjects', 'sex', 'age', 'TIV', 'site', 'group']
    df_2[cols_to_add] = df[cols_to_add]
    df_2['site'] = df_2['site'].map(sites_dict)

    df_1.to_csv("metrics_on_tinnitus.csv")
    df_2.to_csv("scores_on_tinnitus.csv")

    ## remove unnecessary stuff
    pkl_files = glob.glob(os.path.join(roi_models_dir, '*.pkl'))
    for file in pkl_files:
        os.remove(file)



if __name__ == "__main__":
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    roi_fnames = [fname for fname in sorted(os.listdir(roi_dir)) if fname.endswith(".csv")]
    for roi_fname in roi_fnames[:1]:
        estimate_hbr_norm_rois(roi_fname)
        predict_hbr_norm_rois(roi_fname)


