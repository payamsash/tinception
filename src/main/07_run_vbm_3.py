from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

# --- Configuration ---
match_method = "optimal"
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_sites_dir = tinception_dir / "vbm_per_site"
vbm_dir = tinception_dir / "VBM"
design_dir = tinception_dir / "VBM_design"
vbm_struc_dir = vbm_dir / "struc"


vbm_struc_dir.mkdir(parents=True, exist_ok=True)
dry_run = True

df = pd.read_csv(f"../../master_files/matched/matched_{match_method}.csv")
sites = df["site"].unique()
expected_suffixes = ["_struc_brain.nii.gz", "_struc_cut.nii.gz", "_struc.nii.gz"]
id_map = dict(zip(df['subject_ID'].astype(str), df['tinception_id'].astype(str)))

print("\n********** Computing TIV per subject ***********\n")

## compute TIV
def get_tissue_volume(file_path):
    if not file_path.exists():
        print(f"***** Sth is wrong with {file_path.name}*****")
        return 0.0
    
    cmd = ["fslstats", str(file_path), "-m", "-v"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout.split()
    if len(output) < 3:
        print(f"***** Sth is wrong with {file_path.name}*****")
        return 0.0
        
    mean_val = float(output[0])
    total_vol = float(output[2])
    
    return mean_val * total_vol

tiv_list = []
print(f"{'Subject ID':<15} | {'TIV (mm3)':<15}")
print("-" * 35)

for tid in df['tinception_id'].astype(str):
    # Mapping the files based on your directory listing
    csf_file = vbm_struc_dir / f"{tid}_struc_brain_pve_0.nii.gz"
    gm_file  = vbm_struc_dir / f"{tid}_struc_GM.nii.gz"
    wm_file  = vbm_struc_dir / f"{tid}_struc_brain_pve_2.nii.gz"
    
    # Calculate each tissue volume
    vol_csf = get_tissue_volume(csf_file)
    vol_gm  = get_tissue_volume(gm_file)
    vol_wm  = get_tissue_volume(wm_file)
    
    tiv = vol_csf + vol_gm + vol_wm
    tiv_list.append(tiv)
    print(f"{tid:<15} | {tiv:>15.2f}")

# Update CSV
df['TIV'] = tiv_list
df.to_csv(design_dir / "covars.csv", index=False)
print(f"\nDone! TIV added to {design_dir}")


print("\n********** Creating design and contrast files ***********\n")

## create the design.mat and design.con file
df = pd.read_csv(design_dir / "covars.csv")
existing_ids = [f.name.split('_')[0] for f in vbm_struc_dir.glob("*_struc.nii.gz")]
df = df[df['tinception_id'].astype(str).isin(existing_ids)]
df = df.sort_values("tinception_id").reset_index(drop=True)
print(f"Generating design for {len(df)} subjects...")

## creating covars
age = df['age'] - df['age'].mean()
sex = df['sex']
site = df['site']
site_dummies = pd.get_dummies(df['site'], prefix='site', drop_first=True).astype(int)
pta = df['PTA'] - df['PTA'].mean()
tiv = df['TIV'] - df['TIV'].mean()
group1 = (df['group'] == "CO").astype(int)
group2 = (df['group'] == "TI").astype(int)

# create design matrix: Group1, Group2, Age, Sex, Site, PTA, TIV
design_matrix = np.column_stack([
    group1, group2, age, sex, site_dummies, pta, tiv
])

def write_fsl_mat(matrix, output_path):
    with open(output_path, 'w') as f:
        f.write(f"/NumWaves {matrix.shape[1]}\n")
        f.write(f"/NumPoints {matrix.shape[0]}\n")
        f.write("/Matrix\n")
        np.savetxt(f, matrix, fmt='%.6f', delimiter='\t')

# create contrasts
contrasts = np.array([
    [1, -1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # CO > TI
    [-1, 1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # TI > CO
])

def write_fsl_con(con_matrix, output_path):
    with open(output_path, 'w') as f:
        f.write(f"/NumWaves {con_matrix.shape[1]}\n")
        f.write(f"/NumContrasts {con_matrix.shape[0]}\n")
        f.write("/Matrix\n")
        np.savetxt(f, con_matrix, fmt='%.6f', delimiter='\t')

# run it
write_fsl_mat(design_matrix, design_dir / "design.mat")
write_fsl_con(contrasts, design_dir / "design.con")
print(f"Success! Files saved to: {design_dir.resolve()}")