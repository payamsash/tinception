from pathlib import Path
import pandas as pd
import gseapy as gp


tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
pheno_dir = tinception_dir / "phenotypes"
dfs = []
for fname in pheno_dir.iterdir():
    if fname.name.startswith("._"):
        continue
    df_pheno = pd.read_csv(fname, sep="\t")
    df_pheno["structure"] = fname.stem.split("-")[1]
    dfs.append(df_pheno)

df = pd.concat(dfs, axis=0)
df_sig = df.query('pval < 1e-6')
print(f"significant gene is {df_sig['nearest_genes'].values}")
df_candidates = df.query('pval < 1e-5')
genes = df_candidates["nearest_genes"].unique()
print(f"candidate genes are {genes}")

genes = ["CAMTA1", "WNT4"] # removing non coding one
enr = gp.enrichr(
    gene_list=genes,
    organism="Human",
    gene_sets=[
        "GO_Biological_Process_2021",
        "KEGG_2021_Human",
        "Reactome_2022"
    ],
    cutoff=0.05
)
print(enr.results.head())