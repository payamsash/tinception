import os
from pathlib import Path
from MIND import compute_MIND

def compute_ssa():
    subjects_dir = "/home/ubuntu/volume/Tinception/subjects_fs_dir"
    ssa_dir = Path(f"/home/ubuntu/volume/Tinception/SSA")
    os.makedirs(ssa_dir, exist_ok=True)

    features = ["CT", "MC", "Vol", "SD", "SA"]
    parc = "Schaefer2018_1000Parcels_7Networks_order"
    template_subjects = ["MNI152", "fsaverage", "fsaverage5"]

    for subject in sorted(os.listdir(subjects_dir)):
        ssa_fname = ssa_dir / f"{subject}_{parc}.csv"
            
        if ssa_fname.exists() or subject in template_subjects:
            print(f"SSA already computed on {subject}!")
            continue
            
        else:
            try: 
                df = compute_MIND(
                            surf_dir=f"{subjects_dir}/{subject}",
                            features=features,
                            parcellation=parc,
                            filter_vertices=False,
                            resample=False,
                            n_samples=4000
                            ) 
                fname_save = Path(ssa_dir) / f"{subject}_{parc}.csv"
                df.to_csv(fname_save)
                print(f"SSA computed on {subject} with {parc} atlas!")
            except:
                print(f"---------- {subject} failed -----------")
            

if __name__ == "__main__":
    compute_ssa()

