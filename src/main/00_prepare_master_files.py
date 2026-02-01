from pathlib import Path
import pandas as pd

## 1. check sex if 0=female and 1=male
## 2. add andes and antinomics, and erlangen

demographics_dir = Path("../../master_files")

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

    if fname.name in ["andes_master", "antinomics_master"]:
        continue

    if fname.suffix == ".xlsx":
        dfs.append(pd.read_excel(fname))

df = pd.concat(dfs, axis=0)
df.sort_values(by=["site", "subject_ID"], inplace=True)
df.drop(columns="subject ID", inplace=True)

## create tinception_id
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

## check missing files and drop the ones with age, sex, PTA missing

## cp the files in VBM dir
## mv the files in SBM dir