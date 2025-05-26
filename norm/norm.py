## load files
norm_dir = Path.cwd().parent / "data" / "norm"
roi_dir = Path.cwd().parent / "data" / "roi_stats"
roi_fname = "aparc_area_lh.csv"
fname = roi_dir / roi_fname
df = pd.read_csv(fname, index_col=0)


## add covariates
fname_cov = Path.cwd().parent / "data" / "tinception_tiv_fsgd.txt"
df_cov = pd.read_csv(fname_cov, index_col=0, delimiter="\t")
df_cov.dropna(inplace=True)
subjects = list(df_cov["Unnamed: 1"])
df = df.query('subjects == @subjects')

## add covariates
df["sex"] = df_cov["Unnamed: 3"].values
df["age"] = df_cov["Unnamed: 4"].values
df["TIV"] = df_cov["Unnamed: 5"].values
df["site"] = df_cov["Unnamed: 6"].values
df["group"] = df_cov["Unnamed: 2"].values

## some renaming and fixes
if "MV(Re)" in df.columns.to_list():
    df = df.rename(columns={"MV(Re)": "MV"})
df.columns = df.columns.str.replace('-', '_')
df.columns = df.columns.str.replace('(', '')
df.columns = df.columns.str.replace(')', '')
df = df.drop(columns=["anterior_horn_of_lateral_ventricle"], errors='ignore')

df['group'] = df['group'].astype('category')
df['sex'] = df['sex'].astype('category')
df['site'] = df['site'].astype('category')

# fix this part
df.drop(columns=["lh_WhiteSurfArea_area", "BrainSegVolNotVent",	"eTIV"], errors="ignore", inplace=True)
regions = df.columns[1:-5]

df = pd.get_dummies(df, columns=['site'])
df.dropna(inplace=True)

## split the data
df_fera = df[[subset=roi_ids]]
df_cov = df[["age", "sex", "site"]]