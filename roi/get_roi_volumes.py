import os
from pathlib import Path
import subprocess
import pandas as pd

def get_roi_volumes():

    ## set paths and rois
    subjects_dir = "/home/ubuntu/volume/subjects_fs_dir"
    tables_dir = f"{subjects_dir}/tables"
    os.makedirs(tables_dir, exist_ok=True)
    
    subjects = []
    for subject in sorted(os.listdir(subjects_dir)):
        if os.path.isdir(f'{subjects_dir}/{subject}/hist'): 
            subjects.append(f"{subject}")

    roi_names = [
                    "amygdalar-nuclei.lh.T1.v22.stats",
                    "amygdalar-nuclei.rh.T1.v22.stats",
                    "hipposubfields.lh.T1.v22.stats",
                    "hipposubfields.rh.T1.v22.stats",
                    "thalamic-nuclei.lh.v13.T1.stats",
                    "thalamic-nuclei.rh.v13.T1.stats",
                    "aseg.stats",
                    "brainstem.v13.stats"
                ]
    
    ## aseg
    for roi_name in roi_names:
        parts = roi_name.split(".")
        roi = parts[0]
        hemi = parts[1]
        if hemi in ["lh", "rh"]: table_fname = f"{roi}_{hemi}.txt"
        else: table_fname = f"{roi}.txt" 
        
        command = [
                    "asegstats2table",
                    "--subjects", *subjects,
                    "--sd", subjects_dir,
                    "--meas", "volume", 
                    "--statsfile", roi_name, 
                    "--tablefile", f"{tables_dir}/{table_fname}"
                ]
        subprocess.run(command, check=True)

        ## modify and convert to .csv
        df = pd.read_csv(f"{tables_dir}/{table_fname}", sep="\t")
        df = df.rename(columns={'Measure:volume': 'subjects'})
        df.to_csv(f"{tables_dir}/{table_fname[:-4]}.csv")
        os.remove(f"{tables_dir}/{table_fname}")
        print(f"{roi_name} is ready ...!")


    ## arousal network
    df_aan = pd.DataFrame(columns=['subjects', 'DR', 'PBC_R',
                            'LC_L', 'LC_R', 'PAG',
                            'MnR', 'LDTg_L', 'LDTg_R',
                            'PBC_L', 'PnO_L', 'PnO_R',
                            'mRt_L', 'mRt_R', 'PTg_R',
                            'PTg_L', 'VTA', 'Whole_arousal_network',
                            ])
    for idx, subject in enumerate(subjects):
        print(f"working on subjet {subject} aan data ...!")
        aan_fname = Path(subjects_dir) / subject / "stats" / "arousalNetworkVolumes.v10.stats"
        if aan_fname.exists():
            with open(aan_fname) as fyle:
                lines = fyle.readlines()

            vols = [subject]
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    name, vol = parts[0], float(parts[1])
                    vols.append(vol)
        else:
            raise FileNotFoundError(f"aan stats not found for subject {subject}.")
        df_aan.loc[idx+1] = vols
        
    print(f"AAN is ready ...!")
    df_aan.to_csv(f"{tables_dir}/AAN.csv")

    ## apacs
    for parc in ["aparc", "aparc.a2009s"]:
        for measure in ["area", "volume", "thickness"]:
            for hemi in ["lh", "rh"]:
                table_fname = f"{parc}_{measure}_{hemi}.txt"
                command = [
                            "aparcstats2table",
                            "--subjects", *subjects,
                            "--meas", measure, 
                            "--hemi", hemi,
                            "--parc", parc, 
                            "--tablefile", f"{tables_dir}/{table_fname}"
                        ]
                subprocess.run(command, check=True)

                ## modify and convert to .csv
                df = pd.read_csv(f"{tables_dir}/{table_fname}", sep="\t")
                df = df.rename(columns={f'{hemi}.{parc}.{measure}': 'subjects'})
                df.to_csv(f"{tables_dir}/{table_fname[:-4]}.csv")
                os.remove(f"{tables_dir}/{table_fname}")
        
        print(f"{parc} is ready ...!")


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