import os
from pathlib import Path
import pandas as pd
import numpy as np

def extract_graphs():
    ssa_dir = "/home/ubuntu/volume/subjects_fs_dir/SSA"
    parc = "aparc"

    graphs, subjects = [], []
    for fname in sorted(os.listdir(ssa_dir)):
        if fname.endswith(f"{parc}.csv"):
            
            ## extract graphs
            print(f"working on {fname[:-4]} ...")
            df = pd.read_csv(Path(ssa_dir) / fname, index_col=0)
            graph = df.to_numpy()
            graphs.append(graph)

            l_idx = fname.find(f"_{parc}.csv")
            subject = fname[:l_idx]
            subjects.append(subject)

    graphs = np.array(graphs)
    subjects = np.array(subjects)
    np.save(f"graphs_{parc}.npy", graphs)
    np.save(f"subjects_{parc}.npy", subjects)

if __name__ == "__main__":
    extract_graphs()