import os
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
import pingouin as pg

def analyze_ssa():

    ## set param
    ssa_dir = "/home/ubuntu/volume/subjects_fs_dir/SSA"
    parc = "aparc"
    method = "fdr_bh"
    add_site = True

    dfs_list = []
    for fname in sorted(os.listdir(ssa_dir)):
        if fname.endswith(f"{parc}.csv"):
            
            ## extract graphs
            print(f"working on {fname[:-4]} ...")
            df = pd.read_csv(Path(ssa_dir) / fname, index_col=0)

            ## create a flat df
            lbs = df.columns.to_list()
            triu_indices = np.triu_indices_from(df, k=1)
            values = df.values[triu_indices]
            label_pairs = [f"{lbs[i]}_vs_{lbs[j]}" for i, j in zip(*triu_indices)]
            flat_df = pd.DataFrame([values], columns=label_pairs)

            # removing very small connections
            threshold = 0.01 * flat_df.values.max()
            flat_df[flat_df < threshold] = 0

            ## add subject
            l_idx = fname.find(f"_{parc}.csv")
            subject = fname[:l_idx]
            flat_df["subjects"] = subject
            dfs_list.append(flat_df)

    df_glob = pd.concat(dfs_list)
    df_glob.to_csv("df_glob.csv")

    ## now I need to add some columns here
    fname_cov = "/home/ubuntu/volume/SBM/FSGD/tinception_tiv_fsgd.txt"
    df_cov = pd.read_csv(fname_cov, index_col=0, delimiter="\t")
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["Unnamed: 1"])
    df_glob = df_glob.query('subjects == @subjects')
    df_glob["sex"] = df_cov["Unnamed: 3"].values
    df_glob["age"] = df_cov["Unnamed: 4"].values
    df_glob["TIV"] = df_cov["Unnamed: 5"].values
    df_glob["site"] = df_cov["Unnamed: 6"].values
    df_glob["group"] = df_cov["Unnamed: 2"].values

    ## fix some naming in columns
    df_glob.columns = df_glob.columns.str.replace('-', '_')
    df_glob.columns = df_glob.columns.str.replace('(', '')
    df_glob.columns = df_glob.columns.str.replace(')', '')

    df_glob['group'] = df_glob['group'].astype('category')
    df_glob['sex'] = df_glob['sex'].astype('category')
    df_glob['site'] = df_glob['site'].astype('category')
    edges = df_glob.columns[:-6]

    ## now the real part
    f_vals, df_groups, df_resids, p_vals, effect_sizes = ([] for _ in range(5))
    for edge in tqdm(edges):
        if add_site:
            formula = f"{edge} ~ group + sex + age + TIV + site"
        else:
            formula = f"{edge} ~ group + sex + age + TIV"

        model = ols(formula, data=df_glob).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        df_glob['sex_num'] = df_glob['sex'].astype('category').cat.codes
        df_glob['site_num'] = df_glob['site'].astype('category').cat.codes

        if add_site:
            ancova = pg.ancova(data=df_glob, dv=edge,
                                covar=["sex_num", "age", "TIV", "site_num"],
                                between='group')
        else:
            ancova = pg.ancova(data=df_glob, dv=edge,
                                covar=["sex_num", "age", "TIV"],
                                between='group')

        f_vals.append(anova_table['F']['group'])
        df_groups.append(int(anova_table['df']['group']))
        df_resids.append(int(anova_table['df']['Residual']))
        p_vals.append(anova_table['PR(>F)']['group'])
        effect_sizes.append(ancova["np2"][0])

    reject, pvals_corrected, _, _ = multipletests(p_vals, method=method) 
    df_results = pd.DataFrame(columns=["Region", "F_val", "df_group",
                                        "df_resid", "p_val", 
                                        "corrected_p_val",
                                        "np2", "reject"])

    counter = 0
    for elements in zip(edges, f_vals, df_groups, df_resids, \
                        p_vals, pvals_corrected, effect_sizes, reject):
        df_results.loc[counter] = [*elements]
        counter += 1

    df_results.to_csv(f"ssa_df_{parc}_site_{add_site}.csv")

if __name__ == "__main__":
    analyze_ssa()