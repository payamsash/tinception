from pathlib import Path
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import pingouin as pg

def perform_roi_stats(roi, method="bonferroni", add_site=True):

    ## select what to read
    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / f"{roi}.csv"
    df = pd.read_csv(fname, index_col=0)

    ## select subjects
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

    if "area" in roi or "thickness" in roi:
        regions = df.columns[1:-8]
    if "volume" in roi:
        regions = df.columns[1:-7]
    else:
        regions = df.columns[1:-5]
    
    ## now the real part
    f_vals, df_groups, df_resids, p_vals, effect_sizes = ([] for _ in range(5))
    for region in regions:
        if add_site:
            formula = f"{region} ~ group + sex + age + TIV + site"
        else:
            formula = f"{region} ~ group + sex + age + TIV"

        model = ols(formula, data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        df['sex_num'] = df['sex'].astype('category').cat.codes
        df['site_num'] = df['site'].astype('category').cat.codes

        if add_site:
            ancova = pg.ancova(data=df, dv=region,
                                covar=["sex_num", "age", "TIV", "site_num"],
                                between='group')
        else:
            ancova = pg.ancova(data=df, dv=region,
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
    for elements in zip(regions, f_vals, df_groups, df_resids, \
                        p_vals, pvals_corrected, effect_sizes, reject):
        df_results.loc[counter] = [*elements]
        counter += 1
    
    return df_results


if __name__ == "__main__":
    perform_roi_stats()