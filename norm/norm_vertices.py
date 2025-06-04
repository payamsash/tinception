import os
import glob
from pathlib import Path
import re
import shutil
import pickle
import numpy as np
import pandas as pd

import nibabel as nib
import pcntoolkit as ptk
from sklearn.model_selection import train_test_split


def estimate_hbr_norm_rois(hemi):
    
    os.chdir("/home/ubuntu/data/src_codes/norm")
    surfs_dir = Path("/home/ubuntu/volume/NORM/surfs")
    vertex_norm_dir = Path("/home/ubuntu/volume/NORM/vertex_hbr")
    fname_cov = "/home/ubuntu/volume/SBM/FSGD/tinception_tiv_fsgd.txt"

    # os.chdir("/Users/payamsadeghishabestari/tinception/norm")
    # vertex_norm_dir = Path("/Users/payamsadeghishabestari/tinception/data/norm/vertex_models")
    # vertex_norm_dir = Path("/home/ubuntu/volume/NORM/vertex_hbr")
    # fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"

    df_cov = pd.read_csv(fname_cov, delimiter="\t")
    df_cov.dropna(inplace=True)

    mapper = {
                "Unnamed: 1": "subjects",
                "Unnamed: 2": "group",
                "Unnamed: 3": "sex",
                "Unnamed: 4": "age",
                "Unnamed: 5": "TIV",
                "Unnamed: 6": "site"
            }
    
    df_cov.rename(columns=mapper, inplace=True)
    df_cov.drop(columns="GroupDescriptorFile 1", inplace=True)
    df_cov.reset_index(drop="index", inplace=True)
    df_cov["age"] = df_cov["age"].values / 100 # normalize

    df_cov['group'] = df_cov['group'].astype('category')
    df_cov['sex'] = df_cov['sex'].astype('category')
    df_cov['site'] = df_cov['site'].astype('category')

    df_cov_controls = df_cov[df_cov["group"] == "CO"]
    subjects_control = df_cov_controls["subjects"]

    ########### all controls resp
    data_te = []
    for subject in subjects_control:
        fname_thickness = surfs_dir / f"{subject}_{hemi}.thickness.fsaverage5.nii"
        data = nib.load(fname_thickness).get_fdata().squeeze()
        data_te.append(data)

    data_te = np.vstack(data_te)

    ############ now split them
    X_tr, X_te, y_train, y_test = train_test_split(
                                                        df_cov_controls,
                                                        data_te,
                                                        stratify=df_cov_controls["site"],
                                                        test_size=0.2,
                                                        random_state=42
                                                        )

    be_tr = X_tr["site"].astype(int)
    be_te = df_cov_controls["site"].astype(int)
    X_tr.reset_index(drop="index", inplace=True)
    X_tr = X_tr[["age", "sex", "TIV"]]
    X_te = df_cov_controls[["age", "sex", "TIV"]]
    y_test = data_te

    datas = [X_tr, y_train, be_tr, X_te, y_test, be_te]
    file_names = ["X_train", "Y_train", "trbefile", "X_test", "Y_test", "tsbefile"]

    os.chdir(vertex_norm_dir)
    for data, file_name in zip(datas, file_names):
        with open(f"{file_name}.pkl", "wb") as file:
            pickle.dump(pd.DataFrame(data), file)

    respfile = os.path.join(vertex_norm_dir, "Y_train.pkl")
    covfile = os.path.join(vertex_norm_dir, "X_train.pkl")  
    testrespfile_path = os.path.join(vertex_norm_dir, "Y_test.pkl") 
    testcovfile_path = os.path.join(vertex_norm_dir, "X_test.pkl") 
    trbefile = os.path.join(vertex_norm_dir, "trbefile.pkl") 
    tsbefile = os.path.join(vertex_norm_dir, "tsbefile.pkl")

    output_path = os.path.join(vertex_norm_dir, "Models/")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    site_names = ["karolinska", "cogtail", "tinspect", "neuropren"] # no talaska
    sites_dict = dict(zip([0, 1, 3, 4], site_names))

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
    lbs = range(data_te.shape[1])
    dfs_list_1 = []

    for metric in metrics[:-3]:
        df_metric = pd.read_pickle(f"{metric}_estimate.pkl")
        df_metric.rename(columns={df_metric.columns[0]: metric}, inplace=True)
        dfs_list_1.append(df_metric)

    df_1 = pd.concat(dfs_list_1, axis=1)
    df_1.to_csv(f"metrics_on_controls_{hemi}.csv")

    cols_to_add = ['subjects', 'sex', 'age', 'TIV', 'site', 'group']
    for metric in metrics[-3:]:
        df_metric = pd.read_pickle(f"{metric}_estimate.pkl")
        df_metric.rename(columns=dict(zip(df_metric.columns, lbs)), inplace=True)
        df_metric[cols_to_add] = df_cov[cols_to_add]
        df_metric['site'] = df_metric['site'].map(sites_dict)
        df_metric = df_metric.T
        df_metric.to_csv(f"scores_on_controls_{metric}_{hemi}.csv")

    ## remove unnecessary stuff
    pkl_files = glob.glob(os.path.join(vertex_norm_dir, '*.pkl'))
    for file in pkl_files:
        os.remove(file)


def predict_hbr_norm_rois(hemi):

    
    os.chdir("/home/ubuntu/data/src_codes/norm")
    surfs_dir = Path("/home/ubuntu/volume/NORM/surfs")
    vertex_norm_dir = Path("/home/ubuntu/volume/NORM/vertex_hbr")
    fname_cov = "/home/ubuntu/volume/SBM/FSGD/tinception_tiv_fsgd.txt"

    # os.chdir("/Users/payamsadeghishabestari/tinception/norm")
    # vertex_norm_dir = Path("/Users/payamsadeghishabestari/tinception/data/norm/vertex_models")
    #fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"

    df_cov = pd.read_csv(fname_cov, delimiter="\t")
    df_cov.dropna(inplace=True)

    mapper = {
                "Unnamed: 1": "subjects",
                "Unnamed: 2": "group",
                "Unnamed: 3": "sex",
                "Unnamed: 4": "age",
                "Unnamed: 5": "TIV",
                "Unnamed: 6": "site"
            }
    df_cov.rename(columns=mapper, inplace=True)
    df_cov.drop(columns="GroupDescriptorFile 1", inplace=True)
    df_cov.reset_index(drop="index", inplace=True)
    df_cov["age"] = df_cov["age"].values / 100 # normalize
    df_cov = df_cov[df_cov["site"] != 2] # drop talaske

    df_cov['group'] = df_cov['group'].astype('category')
    df_cov['sex'] = df_cov['sex'].astype('category')
    df_cov['site'] = df_cov['site'].astype('category')
    subjects = df_cov["subjects"]

    df_cov_tinnitus = df_cov[df_cov["group"] == "TI"]
    subjects_tinnitus = df_cov_tinnitus["subjects"]

    site_names = ["karolinska", "cogtail", "tinspect", "neuropren"] # no talaska
    sites_dict = dict(zip([0, 1, 3, 4], site_names))

    ########### all controls resp
    data_te = []
    for subject in subjects_tinnitus:
        fname_thickness = surfs_dir / f"{subject}_{hemi}.thickness.fsaverage5.nii"
        # fname_thickness = f"/Users/payamsadeghishabestari/antinomics_clean_codes/dvob_processed/sMRI/surfs/{subject}_{hemi}.thickness.fsaverage.mgh"
        data = nib.load(fname_thickness).get_fdata().squeeze()
        data_te.append(data)

    data_te = np.vstack(data_te)

    ############ now split them
    X_te = df_cov_tinnitus[["sex", "age", "TIV"]]
    be_te = df_cov_tinnitus["site"].astype(int)

    ## fix missing vertex models
    model_path = str(vertex_norm_dir / "Models")
    model_id_mapping, data_te = fix_model_idxs(model_path, data_te, move=True)


    datas = [X_te, data_te, be_te]
    file_names = ["X_test_fit", "Y_test_fit", "tsbefile_fit"]

    os.chdir(vertex_norm_dir)
    for data, file_name in zip(datas, file_names):
        with open(f"{file_name}.pkl", "wb") as file:
            pickle.dump(pd.DataFrame(data), file)

    testrespfile_path = os.path.join(vertex_norm_dir, "Y_test_fit.pkl") 
    testcovfile_path = os.path.join(vertex_norm_dir, "X_test_fit.pkl") 
    tsbefile = os.path.join(vertex_norm_dir, "tsbefile_fit.pkl")

    ############### on all subjects from sites (except talaska)

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
    lbs = list(model_id_mapping.keys())
    dfs_list_1 = []

    for metric in metrics[:-3]:
        df_metric = pd.read_pickle(f"{metric}_fit.pkl")
        df_metric.rename(columns={df_metric.columns[0]: metric}, inplace=True)
        dfs_list_1.append(df_metric)

    df_1 = pd.concat(dfs_list_1, axis=1)
    df_1["lbs"] = list(model_id_mapping.keys()) 
    df_1.to_csv(f"metrics_on_tinnitus_{hemi}.csv")

    cols_to_add = ['subjects', 'sex', 'age', 'TIV', 'site', 'group']
    for metric in metrics[-3:]:
        df_metric = pd.read_pickle(f"{metric}_fit.pkl")
        df_metric.rename(columns=dict(zip(df_metric.columns, lbs)), inplace=True)
        df_metric[cols_to_add] = df_cov[cols_to_add]
        df_metric['site'] = df_metric['site'].map(sites_dict)
        df_metric = df_metric.T
        df_metric.to_csv(f"scores_on_tinnitus_{metric}_{hemi}.csv")

    ## remove unnecessary stuff
    pkl_files = glob.glob(os.path.join(vertex_norm_dir, '*.pkl'))
    for file in pkl_files:
        os.remove(file)




def fix_model_idxs(model_path, data_te, move=True):
    models_list = sorted(glob.glob(os.path.join(model_path, '*.pkl')))

    ## find missing ones
    found_idxs = set(
        int(re.search(r'NM_0_(\d+)_estimate\.pkl', f).group(1))
        for f in models_list
        if re.search(r'NM_0_(\d+)_estimate\.pkl', f)
    )
    old_idxs = list(found_idxs)
    new_idxs = range(len(old_idxs))
    model_id_mapping = dict(zip(old_idxs, new_idxs))
    data_te = data_te[:, old_idxs]

    if move:
        for old_idx, new_idx in zip(old_idxs, new_idxs):
            shutil.move(
                        Path(model_path) / f"NM_0_{old_idx}_estimate.pkl",
                        Path(model_path) / f"NM_0_{new_idx}_estimate.pkl"
                        )

    return model_id_mapping, data_te

if __name__ == "__main__":

    hemi = "lh"
    estimate_hbr_norm_rois(hemi)
    predict_hbr_norm_rois(hemi)


