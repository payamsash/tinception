from pathlib import Path
import pandas as pd

demographics_dir = Path("../../master_files")
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
vbm_dir = tinception_dir / "VBM"

## add missing sites (erlangen, andes, antinomics)
site_base = {
            "andes":      1000,
            "antinomics": 2000,
            "cogtail":    3000,
            "iccac":      4000,
            "inditms":    5000,
            "london_1":   6000,
            "london_2":   7000,
            "neuropren":  8000,
            "new-castle": 9000,
            "tinmeg":     10000,
            "tinspect":   11000,
            "triple":     12000,
            }

## read and sort files
dfs = []
for fname in sorted(demographics_dir.iterdir()):
    if fname.stem in ["andes_master", "antinomics_master"]:
        continue
    if fname.suffix == ".xlsx":
        dfs.append(pd.read_excel(fname, dtype={"subject_ID": str}))

df = pd.concat(dfs, axis=0)
df["subject_ID"] = df["subject_ID"].astype(str)
df.sort_values(by=["site", "subject_ID"], inplace=True)
df["tinception_id"] = (
    "sub-"
    + (
        df["site"].map(site_base)
        + df.groupby("site").cumcount()
        + 1
    )
    .astype(int)
    .astype(str)
    .str.zfill(5)
)

## drop missing vars and number checks
df = df.dropna(subset=["age", "sex", "PTA"])
def extract_subject_id(fname: str) -> str:
    name = fname.replace(".nii.gz", "")
    name = name.replace("sub-", "")
    return name

for site in sorted(df["site"].unique()):
    print(f"\n{'=' * 60}")
    print(f"Site: {site}")
    print(f"{'=' * 60}")

    df_site = df[df["site"] == site]
    df_ids = set(df_site["subject_ID"])
    mri_dir = vbm_dir / site

    if not mri_dir.exists():
        print(f"MRI directory missing: {mri_dir}")
        continue

    nii_files = list(mri_dir.glob("*.nii.gz"))
    file_ids = set(extract_subject_id(p.name) for p in nii_files)
    in_files_not_df = sorted(file_ids - df_ids)
    in_df_not_files = sorted(df_ids - file_ids)

    print(f"Subjects in MRI files: {len(file_ids)}")
    print(f"Subjects in demographics: {len(df_ids)}")

    if in_files_not_df:
        print("\n In MRI files but NOT in demographics:")
        for sid in in_files_not_df:
            print(f"  - {sid}")
    else:
        print("\n No subjects missing from demographics")

    if in_df_not_files:
        print("\n In demographics but NOT in MRI files:")
        for sid in in_df_not_files:
            print(f"  - {sid}")
    else:
        print("\n No subjects missing from MRI files")

## 3. add a final check to see order of happening in df and vbm, sbm dir are same
## 4. drop some from sbm_dir as well

counts = (
        pd.crosstab(df["site"], df["group"])
        .rename(columns={0: "control", 1: "patient"})
        )

print("\n************** SITE × GROUP COUNTS **************")
print(f"{counts}")
print("*************************************************\n")