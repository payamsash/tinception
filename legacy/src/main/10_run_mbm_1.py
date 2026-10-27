from pathlib import Path
import pandas as pd
import numpy as np


# Paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
subjects_dir = tinception_dir / "subjects_fs_dir"
fsgd_fname = tinception_dir / "SBM" / "FSGD" / "tinception.fsgd"
out_dir = tinception_dir / "MBM"
out_dir.mkdir(parents=True, exist_ok=True)

# Settings
mode = "thickness"
hemis = ["lh", "rh"]
group_map = {
    "CO": 1,
    "TI": 2,
}


# Read FSGD
df = pd.read_csv(
    fsgd_fname,
    skiprows=5,
    sep=r"\s+",
    names=["input", "subjects", "group", "sex", "age", "PTA", "TIV"]
).drop(columns="input")

# Clean text columns
df["subjects"] = df["subjects"].astype(str).str.strip()
df["group"] = df["group"].astype(str).str.strip().str.upper()

# Convert numeric columns
df["sex"] = pd.to_numeric(df["sex"], errors="raise")   # already 0/1
df["age"] = pd.to_numeric(df["age"], errors="raise")
df["PTA"] = pd.to_numeric(df["PTA"], errors="raise")
df["TIV"] = pd.to_numeric(df["TIV"], errors="raise")

# Map group labels to numeric codes
df["group_code"] = df["group"].map(group_map)

if df["group_code"].isna().any():
    bad_groups = df.loc[df["group_code"].isna(), "group"].unique().tolist()
    raise ValueError(
        f"These group labels are not in group_map: {bad_groups}\n"
        f"Found groups in file: {df['group'].unique().tolist()}"
    )


# Create shared design matrices
# Two-sample design: one-hot encoding
groups_sorted = sorted(group_map.items(), key=lambda x: x[1])
g1_name, g1_code = groups_sorted[0]
g2_name, g2_code = groups_sorted[1]

G_two_sample = pd.DataFrame({
    g1_name: (df["group_code"] == g1_code).astype(int),
    g2_name: (df["group_code"] == g2_code).astype(int),
})

two_sample_file = out_dir / f"{mode}_design_two_sample.txt"
np.savetxt(two_sample_file, G_two_sample.values, fmt="%d")

# ANCOVA design:
# group + age + sex + PTA + TIV
G_ancova = pd.DataFrame({
    "group_code": df["group_code"].astype(int),
    "age": df["age"] - df["age"].mean(),
    "sex": df["sex"].astype(int),
    "PTA": df["PTA"] - df["PTA"].mean(),
    "TIV": df["TIV"] - df["TIV"].mean(),
})

ancova_file = out_dir / f"{mode}_design_ancova.txt"
np.savetxt(ancova_file, G_ancova.values, fmt=["%d", "%.6f", "%d", "%.6f", "%.6f"])


# Create hemisphere-specific map lists
for hemi in hemis:
    map_paths = []
    missing = []

    for subject in df["subjects"]:
        fname = subjects_dir / subject / "surf" / f"{hemi}.{mode}.fsaverage.mgh"
        if fname.exists():
            map_paths.append(str(fname))
        else:
            missing.append(str(fname))

    if missing:
        print(f"\nMissing files for {hemi}:")
        for m in missing[:20]:
            print(m)
        if len(missing) > 20:
            print(f"... and {len(missing) - 20} more")
        raise FileNotFoundError(f"Missing {len(missing)} files for {hemi}")

    map_list_file = out_dir / f"{hemi}_{mode}_map_list.txt"
    with open(map_list_file, "w") as f:
        for p in map_paths:
            f.write(p + "\n")

    print(f"Saved: {map_list_file}")

# Summary
print("\nDone.")
print(f"Saved: {two_sample_file}")
print(f"Saved: {ancova_file}")
print("Group mapping used:", group_map)
print("Two-sample design columns:", list(G_two_sample.columns))
print("ANCOVA columns:", list(G_ancova.columns))
print("N subjects:", len(df))
print("Groups found:", df["group"].unique().tolist())