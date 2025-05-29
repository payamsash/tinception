import os
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
import nibabel as nib

from pcntoolkit.normative import estimate
from pcntoolkit.util.utils import create_bspline_basis

def estimate_norm_on_vertices(hemi):

    ###################### prepare data
    subjects_dir = Path("/Users/payamsadeghishabestari/antinomics_clean_codes/dvob_processed/sMRI")
    vertex_model_dir = Path("/Users/payamsadeghishabestari/antinomics_clean_codes/vertex_norm")
    os.makedirs(vertex_model_dir, exist_ok=True)

    fname = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
    df = pd.read_csv(fname, index_col=0, delimiter="\t")
    df.dropna(inplace=True)
    subjects = list(df["Unnamed: 1"])

    data_te = []
    for subject in subjects:
        fname_thickness_lh = subjects_dir / subject / "surf" / f"{hemi}.thickness.fsaverage.mgh"
        data = nib.load(fname_thickness_lh).get_fdata().squeeze()
        data_te.append(data)

    data_te = np.vstack(data_te)

    ####################### create cov files for train and test
    df_cov = pd.DataFrame()
    groups = df["Unnamed: 2"].values
    df_cov["group"] = groups
    df_cov["sex"] = df["Unnamed: 3"].values
    df_cov["age"] = df["Unnamed: 4"].values
    df_cov["TIV"] = df["Unnamed: 5"].values
    df_cov["site"] = df["Unnamed: 6"].values

    df_cov['group'] = df_cov['group'].astype('category')
    df_cov['sex'] = df_cov['sex'].astype('category')
    df_cov['site'] = df_cov['site'].astype('category')

    df_cov["age_demean"] = df_cov["age"] - df_cov["age"].mean()
    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)
    df_con_cov = df_cov.query('group == "CO"')
    control_idx = df_cov[df_cov['group'] == 'CO'].index
    [df_sub.drop(columns=["group", "age"], inplace=True) for df_sub in [df_cov, df_con_cov]]
    data_tr = data_te[control_idx]
    
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
    np.savetxt(vertex_model_dir / 'cov_bspline_tr.txt', X_tr)
    np.savetxt(vertex_model_dir / 'cov_bspline_te.txt', X_te)

    cov_file_tr = str(vertex_model_dir / 'cov_bspline_tr.txt')
    cov_file_te = str(vertex_model_dir / 'cov_bspline_te.txt')

    ###################### estimate norm brain
    os.makedirs(vertex_model_dir / "results", exist_ok=True)

    ver_idxs = range(data_tr.shape[1])

    df_norm = pd.DataFrame(columns=subjects)

    for ver_idx in tqdm(ver_idxs):
        np.savetxt(vertex_model_dir / f'resp_tr_{ver_idx}.txt', data_tr[:, ver_idx])
        np.savetxt(vertex_model_dir / f'resp_te_{ver_idx}.txt', data_te[:, ver_idx])

        resp_file_tr = str(vertex_model_dir / f'resp_tr_{ver_idx}.txt')
        resp_file_te = str(vertex_model_dir / f'resp_te_{ver_idx}.txt')

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
        
        df_norm.loc[len(df_norm)] = z_scores

        os.remove(vertex_model_dir / f'resp_tr_{ver_idx}.txt')
        os.remove(vertex_model_dir / f'resp_te_{ver_idx}.txt')

    df_norm.to_csv(vertex_model_dir / "results" / f"df_{hemi}_zscores.csv")

if __name__ == "__main__":
    for hemi in ["lh", "rh"]:
        estimate_norm_on_vertices(hemi)