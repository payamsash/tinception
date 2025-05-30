import os
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
import nibabel as nib

from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis

def estimate_norm_on_vertices(hemi):

    ###################### prepare data
    subjects_dir = Path("/home/ubuntu/volume/subjects_fs_dir")
    vertex_model_dir = Path("/Users/payamsadeghishabestari/tinception/data/norm/vertex_models")
    site_names = ["karolinska", "cogtail", "talaska", "tinspect", "neuropren"]
    os.makedirs(vertex_model_dir, exist_ok=True)

    # fname_fsgd = "/home/ubuntu/volume/SBM/FSGD/tinception_tiv_fsgd.txt"
    fname_fsgd = "/Users/payamsadeghishabestari/tinception/data/tinception_tiv_fsgd.txt"
    df = pd.read_csv(fname_fsgd, index_col=0, delimiter="\t")
    df.dropna(inplace=True)
    subjects = list(df["Unnamed: 1"])

    data_te = []
    for subject in subjects:
        fname_thickness = subjects_dir / subject / "surf" / f"{hemi}.thickness.fsaverage5.mgh"
        # fname_thickness = f"/Users/payamsadeghishabestari/antinomics_clean_codes/dvob_processed/sMRI/surfs/{subject}_{hemi}.thickness.fsaverage.mgh"
        data = nib.load(fname_thickness).get_fdata().squeeze()
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

    site_indices = {}
    for site_val in df_cov["site"].unique():
        site_indices[site_val] = df_cov.index[df_cov["site"] == site_val].tolist()


    df_cov = pd.get_dummies(df_cov, columns=['site'], dtype=float)
    df_con_cov = df_cov.query('group == "CO"')
    control_idx = df_cov[df_cov['group'] == 'CO'].index
    [df_sub.drop(columns=["group", "age"], inplace=True) for df_sub in [df_cov, df_con_cov]]
    df_con_cov.reset_index(inplace=True, drop=True)
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

    df_norm = pd.DataFrame(columns=subjects + ["vertex_idx"])
    df_metrics_list = []
    df_site_metrics = pd.DataFrame(columns = ['vertex_idx', 'site', 'MSLL', 'EXPV', 'SMSE', 'RMSE', 'Rho'])

    for ver_idx in tqdm(ver_idxs[:1000]):
        np.savetxt(vertex_model_dir / f'resp_tr_{ver_idx}.txt', data_tr[:, ver_idx])
        np.savetxt(vertex_model_dir / f'resp_te_{ver_idx}.txt', data_te[:, ver_idx])

        resp_file_tr = str(vertex_model_dir / f'resp_tr_{ver_idx}.txt')
        resp_file_te = str(vertex_model_dir / f'resp_te_{ver_idx}.txt')

        # run a basic model
        try:
            yhat_te, s2_te, _, z_scores, metric_te = estimate(
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
        
            z_scores = np.round(z_scores, 4)
            df_norm.loc[len(df_norm)] = list(z_scores[:,0]) + [ver_idx]

            os.remove(vertex_model_dir / f'resp_tr_{ver_idx}.txt')
            os.remove(vertex_model_dir / f'resp_te_{ver_idx}.txt')


            metric_te_flattened = {k: v.flatten()[0] for k, v in metric_te.items()}
            df_metric = pd.DataFrame([metric_te_flattened])
            df_metric["vertex_idx"] = ver_idx
            df_metrics_list.append(df_metric)


            for key, val in site_indices.items():
                y_mean_te_site = np.array([data_te[val, ver_idx].mean()])
                y_var_te_site = np.array([data_te[val, ver_idx].var()])


                metrics_te_site = evaluate(data_te[val, ver_idx][:, np.newaxis],
                                            yhat_te[val], s2_te[val], y_mean_te_site, y_var_te_site)
                
                df_site_metrics.loc[len(df_site_metrics)] = [
                                                                ver_idx,
                                                                site_names[int(key)],
                                                                metrics_te_site['MSLL'][0],
                                                                metrics_te_site['EXPV'][0],
                                                                metrics_te_site['SMSE'][0],
                                                                metrics_te_site['RMSE'][0],
                                                                metrics_te_site['Rho'][0]
                                                                ]
        
        except:
            os.remove(vertex_model_dir / f'resp_tr_{ver_idx}.txt')
            os.remove(vertex_model_dir / f'resp_te_{ver_idx}.txt')

            continue

    df_norm.to_csv(vertex_model_dir / "results" / f"df_{hemi}_zscores.csv")
    df_metrics = pd.concat(df_metrics_list)
    df_metrics.to_csv(vertex_model_dir / "results" / f"metrics_{hemi}.csv")
    df_site_metrics.to_csv(vertex_model_dir / "results" / f"site_metrics_{hemi}.csv")

if __name__ == "__main__":
    for hemi in ["lh"]:
        estimate_norm_on_vertices(hemi)