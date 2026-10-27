from pathlib import Path
import shutil
import pandas as pd
import json

# 1. Setup paths
data_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
qc_dir = data_dir / "QC"
vbm_dir = data_dir / "vbm"
fname_master = "../../master_files/master.csv"

# 2. Site list from your prompt
sites = [
        "andes",
        "antinomics",
        "cogtail",
        "erlangen",
        "iccac",
        "inditms",
        "london_1",
        "london_2", 
        "neuropren",
        "new-castle",
        "tinmeg",
        "tinspect",
        "triple"
        ]

# 3. Load Mapping (assuming 'subject_ID' matches the filename without .nii.gz)
df = pd.read_csv(fname_master)
df['subject_ID'] = df['subject_ID'].astype(str)
id_map = dict(zip(df['subject_ID'], df['tinception_id']))

# 4. Initialize BIDS structure
qc_dir.mkdir(parents=True, exist_ok=True)
desc = {"Name": "Tinception_QC", "BIDSVersion": "1.4.0", "DatasetType": "raw"}
with open(qc_dir / "dataset_description.json", "w") as f:
    json.dump(desc, f, indent=4)

# 5. Loop over each site and copy
success_count = 0
missing_in_csv = 0

for site in sites:
    site_path = vbm_dir / site
    
    if not site_path.exists():
        print(f"Skipping site {site}: Folder not found at {site_path}")
        continue
        
    print(f"Processing site: {site}...")
    
    # Grab only .nii.gz files in the top level of the site folder
    for fpath in site_path.glob("*.nii.gz"):
        raw_id = fpath.name.replace(".nii.gz", "")
        
        if raw_id in id_map:
            bids_id = id_map[raw_id]
            
            # Ensure the ID starts with 'sub-'
            sub_label = bids_id if bids_id.startswith("sub-") else f"sub-{bids_id}"
            
            # Create BIDS destination: QC/sub-XXXX/anat/sub-XXXX_T1w.nii.gz
            dest_dir = qc_dir / sub_label / "anat"
            dest_file = dest_dir / f"{sub_label}_T1w.nii.gz"
            
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(fpath, dest_file)
            success_count += 1
        else:
            missing_in_csv += 1

print(f"\n--- Migration Complete ---")
print(f"Successfully copied: {success_count} files")
print(f"Files found but missing in CSV: {missing_in_csv}")