import os
from pathlib import Path
import subprocess
import pandas as pd

def get_roi_volumes():

    ## set paths and rois
    fs_dir = "/usr/local/freesurfer/8.0.0"
    subjects_dir = "/home/ubuntu/volume/subjects_fs_dir"
    tables_dir = f"{subjects_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)
    
    subjects = []
    for subject in os.listdir(subjects_dir):
        if os.path.isdir(f'{subjects_dir}/{subject}/hist'): 
            subjects.append(f"{subject}")

    roi_names = [
                    "amygdalar-nuclei.lh.T1.v22.stats",
                    "amygdalar-nuclei.rh.T1.v22.stats",
                    "arousalNetworkVolumes.v10.stats",
                    "aseg.stats",
                    "brainstem.v13.stats",
                    "hipposubfields.lh.T1.v22.stats",
                    "hipposubfields.rh.T1.v22.stats",
                    "hypothalamic_subunits_volumes.v1.stats",
                    "lh.aparc.a2009s.stats",
                    "lh.aparc.DKTatlas.stats",
                    "rh.aparc.a2009s.stats",
                    "rh.aparc.DKTatlas.stats",
                    "thalamic-nuclei.lh.v13.T1.stats",
                    "thalamic-nuclei.rh.v13.T1.stats"
                ]

    ## extract values
    for roi_name in roi_names:
        parts = roi_name.split(".")
        roi = parts[0]
        hemi = parts[1]
        if hemi in ["lh", "rh"]:
            table_fname = f"{roi}_{hemi}.txt"
        elif roi in ["lh", "rh"]:
            table_fname = f"{parts[2]}.txt"
        else:
            table_fname = f"{roi}.txt" 

        command = [
            f"{fs_dir}/bin/asegstats2table",
            "--subjects", *subjects,
            "--sd", subjects_dir,
            "--meas", "volume", 
            "--statsfile", "hipposubfields.lh.T1.v22.stats", 
            "--tablefile", f"{tables_dir}/{table_fname}"
        ]
        subprocess.run(command, check=True)

        ## modify and convert to .csv
        df = pd.read_csv(f"{tables_dir}/{table_fname}", sep="\t")
        df = df.rename(columns={'Measure:volume': 'subjects'})
        df.to_csv(f"{tables_dir}/{table_fname[:-4]}.csv")
        os.remove(f"{tables_dir}/{table_fname}")

    ## Next Brain atlas
    dfs_lh, dfs_rh = [], []
    for subject in subjects:
        fname_lh = Path(subjects_dir) / subject / "hist" / f"vols.left.csv"
        fname_rh = Path(subjects_dir) / subject / "hist" / f"vols.right.csv"
        dfs_lh.append(pd.read_csv(fname_lh))
        dfs_rh.append(pd.read_csv(fname_rh))
    
    df_lh = pd.concat(dfs_lh)
    df_rh = pd.concat(dfs_rh)

    df_lh["subjects"] = subjects
    df_rh["subjects"] = subjects

    col_list = list(df_lh.columns)
    col_list = [col_list[-1]] + col_list[:-1]
    df_lh = df_lh[col_list]
    col_list = list(df_rh.columns)
    col_list = [col_list[-1]] + col_list[:-1]
    df_rh = df_rh[col_list]

    df_lh.to_csv(f"{tables_dir}/hist_lh.csv")
    df_rh.to_csv(f"{tables_dir}/hist_rh.csv")

if __name__ == "__main__":
    get_roi_volumes()