from pathlib import Path
import shutil
import pandas as pd

# --- Configuration ---
match_method = "optimal"
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_sites_dir = tinception_dir / "vbm_per_site"
vbm_dir = tinception_dir / "VBM"
vbm_struc_dir = vbm_dir / "struc"


vbm_struc_dir.mkdir(parents=True, exist_ok=True)
dry_run = False

df = pd.read_csv(f"../../master_files/matched/matched_{match_method}.csv")
sites = df["site"].unique()
expected_suffixes = ["_struc_brain.nii.gz", "_struc_cut.nii.gz", "_struc.nii.gz"]
id_map = dict(zip(df['subject_ID'].astype(str), df['tinception_id'].astype(str)))

print(f"Starting {'DUMMY RUN' if dry_run else 'REAL RUN'} for {len(sites)} sites...")

## real part (task 1)
for site_dir in vbm_sites_dir.iterdir():
    if not site_dir.is_dir() or site_dir.name not in sites:
        continue
    
    print(f"\nProcessing site: {site_dir.name}")
    
    for nii_file in site_dir.glob("*.nii*"):
        subject_id = nii_file.name.split('.')[0]
        match = df[df["subject_ID"] == subject_id]

        if match.empty:
            print(f"  [Skip] No ID match found for: {subject_id}")
            continue
            
        t_id = match["tinception_id"].values[0]
        dest_file = vbm_dir / f"{t_id}.nii.gz"
        
        if dry_run:
            print(f"  [Ready] {nii_file.name} -> {dest_file.name}")
        else:
            if not dest_file.exists():
                shutil.copy2(nii_file, dest_file)
                print(f"  [Copied] {subject_id} as {t_id}")
            else:
                print(f"  [Exists] Skipping {t_id}")

print("\n********** Task 1 complete ***********\n")
print(f"--- Starting {'DUMMY' if dry_run else 'REAL'} Run ---")

for site_dir in vbm_sites_dir.iterdir():
    if not site_dir.is_dir() or site_dir.name not in sites:
        continue

    site_struc_dir = site_dir / "struc"
    if not site_struc_dir.exists():
        continue

    print(f"\nProcessing Site: {site_dir.name}")

    for nii_file in site_struc_dir.glob("*.nii.gz"):
        filename = nii_file.name
        matched_suffix = None
        
        # 1. Identify which suffix this file has
        for suffix in expected_suffixes:
            if filename.endswith(suffix):
                matched_suffix = suffix
                break
        
        if not matched_suffix:
            # Skip files that don't match our 3 specific types
            continue

        # 2. Isolate the subject_id by removing the suffix
        subject_id = filename.replace(matched_suffix, "")
        
        # 3. Check if this subject_id is in our mapping
        if subject_id in id_map:
            t_id = id_map[subject_id]
            new_filename = f"{t_id}{matched_suffix}"
            dest_path = vbm_struc_dir / new_filename

            if dry_run:
                print(f"  [Ready] {filename} -> {new_filename}")
            else:
                if not dest_path.exists():
                    shutil.copy2(nii_file, dest_path)
                    print(f"  [Copied] {new_filename}")
                else:
                    print(f"  [Exists] Skipping {new_filename}")
        else:
            print(f"  [No ID Match] Subject '{subject_id}' not found in dataframe: {filename})")

print("\n********** Task 2 complete ***********\n")