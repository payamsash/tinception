from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps


def compute_gradients(n_parcels, approach, n_components):

    ## get subjects list
    fname_subjects = "/Users/payamsadeghishabestari/tinception/material/SBM/FSGD/tinception_fsgd_PTA.txt"
    ssa_directory = Path("")
    kernel = None


    df = pd.read_csv(fname_subjects, delimiter="\t")
    df.dropna(inplace=True)
    df.rename(columns={"Unnamed: 1": "subjects", "Unnamed: 2": "group"}, inplace=True)
    subjects_to_drop = [
                        "3xo04l_TI-HA", "S16_hwEa_MRT",
                        "2kyuoc_TI-HL", "3s9x2q_NT-HL",
                        "3ucbt4_NT-HL", "4a71bn_TI-HA",
                        "6ddsfu_NT-HL", "71rpqm_NT-HL",
                        "7jindx_TI-HA", "8i29bt_TI-HA"
                        ]
    
    df = df.query('subjects != @subjects_to_drop')
    df = df[["subjects", "group"]]
    subjects = df["subjects"].values.tolist()

    df.reset_index(drop="index", inplace=True)
    co_idxs = df[df["group"] == "CO"].index.to_list()
    ti_idxs = df[df["group"] == "TI"].index.to_list()

    ## create list of connectome
    connectomes = []
    for subject in tqdm(subjects):
        fname = ssa_directory / f"{subject}_Schaefer2018_{n_parcels}Parcels_7Networks_order.csv"
        df = pd.read_csv(fname)
        df.drop(columns="Unnamed: 0", inplace=True)
        connectome = df.to_numpy()
        np.fill_diagonal(connectome, 1)
        connectomes.append(connectome)


    ## compute and align gradients
    gm = GradientMaps(
                        n_components=n_components,
                        approach=approach,
                        kernel=kernel,
                        alignment='procrustes'
                        )
    gm.fit(connectomes)

    ## compute mean connectome per group
    avg_connectome_co = np.array(connectomes)[co_idxs].mean(axis=0)
    avg_connectome_ti = np.array(connectomes)[ti_idxs].mean(axis=0)

    ## save the lambdas and gradients
    np.save(ssa_directory / f"lambda_{n_parcels}_{approach}.npy", gm.lambdas_)
    np.save(ssa_directory / f"gradients_{n_parcels}_{approach}.npy", gm.gradients_)
    np.save(ssa_directory / f"avg_connectome_co.npy", avg_connectome_co)
    np.save(ssa_directory / f"avg_connectome_ti.npy", avg_connectome_ti)

## maybe add normalizations
if __name__ == "__main__":
    for approach in ["dm", "pca"]:
        compute_gradients(
                            1000,
                            approach,
                            n_components=10
                            )
