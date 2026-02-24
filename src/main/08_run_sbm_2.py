import os
from pathlib import Path
import pandas as pd

## paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
subjects_dir = tinception_dir / "subjects_fs_dir"
df = pd.read_csv(tinception_dir / "vbm_norm" / "covars.csv")
fsgd_path = tinception_dir / "SBM" / "FSGD" / "tinception.fsgd"

## fix these 3 subjects later
subjects_to_drop = ["AG_e1", "CK_C_e1", "DW_C_e1"]
df = df[~df['subject_ID'].isin(subjects_to_drop)]

## check missing subjects
folder_subjects = {f for f in os.listdir(subjects_dir) if os.path.isdir(os.path.join(subjects_dir, f))}
df_subjects = set(df['subject_ID'].astype(str))
missing_folders = df_subjects - folder_subjects
not_in_df = folder_subjects - df_subjects
print(f"Total in DF: {len(df_subjects)}")
print(f"Total in Folder: {len(folder_subjects)}")
print("-" * 30)

if missing_folders:
    print(f"CRITICAL: {len(missing_folders)} subjects missing folders (mris_preproc WILL fail):")
    print(list(missing_folders))
else:
    print("All subjects in DataFrame have matching folders.")

if not_in_df:
    print(f"Note: {len(not_in_df)} folders exist that are not in your DF (Safe to ignore).")

## check missing qcache files
fwhm = 10
target = "fsaverage"

lh_pattern = f"lh.thickness.fwhm{fwhm}.{target}.mgh"
rh_pattern = f"rh.thickness.fwhm{fwhm}.{target}.mgh"

missing_qcache = []
print(f"Checking {len(df)} subjects for qcache files...")
for sub in df['subject_ID']:
    lh_path = os.path.join(subjects_dir, sub, 'surf', lh_pattern)
    rh_path = os.path.join(subjects_dir, sub, 'surf', rh_pattern)
    if not os.path.exists(lh_path) or not os.path.exists(rh_path):
        missing_qcache.append(sub)

if not missing_qcache:
    print(f"Success! All subjects have {fwhm}mm qcache files.")
else:
    print(f"Found {len(missing_qcache)} subjects with missing qcache files:")
    print(missing_qcache)


## make it clean and demean
df.drop(columns=["THI", "distance",	"weights", "subclass", "site"], inplace=True)
for col in ['age', 'PTA', 'TIV']:
    df[f'{col}'] = df[col] - df[col].mean()

## write fsgd file
with open(fsgd_path, "w") as f:
    f.write("GroupDescriptorFile 1\n")
    f.write("Title Tinnitus_vs_Controls_SBM\n")
    
    # Define Classes
    f.write("Class CO\n")
    f.write("Class TI\n")
    
    # Define Variables (must match the order in the data lines)
    f.write("Variables sex age PTA TIV\n")
    
    for _, row in df.iterrows():
        line = f"Input {row['subject_ID']} {row['group']} {row['sex']} {row['age']} {row['PTA']} {row['TIV']}\n"
        f.write(line)

print(f"FSGD file saved to {fsgd_path}")

## write contrast files
# Define the contrast vectors
# Order: CO, TI, Sex, Age, PTA, TIV
contrasts = {
    "ti_gt_co": [ -1,  1,  0,  0,  0,  0,  0,  0,  0,  0],  # TI > CO (Hypertrophy)
    "co_gt_ti": [  1, -1,  0,  0,  0,  0,  0,  0,  0,  0],  # CO > TI (Atrophy)
}

for name, vector in contrasts.items():
    with open(tinception_dir / "SBM" / "FSGD" / f"{name}.mtx", "w") as f:
        # Join the numbers with spaces and write as a single line
        line = " ".join(map(str, vector))
        f.write(line + "\n")
        print(f"Created {name}.mtx with: {line}")