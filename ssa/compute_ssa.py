import os
from pathlib import Path
from MIND import compute_MIND

def compute_ssa():
    subjects_dir = "/home/ubuntu/volume/subjects_fs_dir"
    ssa_dir = f"{subjects_dir}/SSA"
    os.makedirs(ssa_dir, exist_ok=True)

    features = ["CT", "MC", "Vol", "SD", "SA"] 
    parcs = [
            "aparc", "aparc.a2009s",
            "Schaefer2018_400Parcels_17Networks_order", 
            "Schaefer2018_600Parcels_17Networks_order",
            "Schaefer2018_800Parcels_17Networks_order"
            ]

    subjects = []
    for subject in sorted(os.listdir(subjects_dir)):
        if os.path.isdir(f'{subjects_dir}/{subject}/hist'):
            subjects.append(f"{subject}")

    for subject in subjects:
        for parc in parcs:
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

if __name__ == "__main__":
    compute_ssa()