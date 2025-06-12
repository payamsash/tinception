from pathlib import Path
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import pingouin as pg
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant

def perform_roi_stats(roi, method="bonferroni", normalize=True):

    roi_dir = Path.cwd().parent / "data" / "roi_stats"
    fname = roi_dir / f"{roi}.csv"
    df = pd.read_csv(fname, index_col=0)

    ## select subjects
    fname_cov = Path.cwd().parent / "material" / "no_talaska" / "behavioural" / "tinception_matched_optimal.xlsx"
    df_cov = pd.read_excel(fname_cov)
    df_cov.dropna(inplace=True)
    subjects = list(df_cov["subject ID"])
    df = df.query('subjects == @subjects')

    fname_design = Path.cwd().parent / "material" / "no_talaska" / "VBM" / "design.txt"
    df_design = pd.read_csv(fname_design, sep=" ", names=["A1", "A2", "sex", "age", "TIV", "PTA", "site_1", "site_2", "site_3"])

    ## add covariates
    df["sex"] = df_cov["sex"].values
    df["age"] = (df_cov["age"] - df_cov["age"].mean()).values
    df["site"] = df_cov["site"].values
    df["TIV"] = (df_design["TIV"] - df_design["TIV"].mean()).values
    df["PTA"] = (df_design["PTA"] - df_design["PTA"].mean()).values
    df["group"] = df_cov["group"].values

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

    if roi.startswith("hippo"):
        regions = df.columns[1:-9] 
    else:
        regions = df.columns[1:-7]

    roi_vol = df.columns[-7]
    print(roi_vol)

    ## check VIF
    X = add_constant(df[['age', 'PTA']])
    df_vif_age_pta = pd.DataFrame({
        'variable': X.columns,
        'VIF': [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    })
    X = add_constant(df[[roi_vol, 'TIV']])
    df_vif_tiv = pd.DataFrame({
        'variable': X.columns,
        'VIF': [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    })


    if normalize:
        df[regions] = df[regions].div(df[roi_vol], axis=0)

    ## now the real part
    f_vals, df_groups, df_resids, p_vals, effect_sizes = ([] for _ in range(5))
    for region in regions:
        formula = f"{region} ~ group + sex + age + site + TIV + PTA + {roi_vol}"
        
        model = ols(formula, data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        df['sex_num'] = df['sex'].astype('category').cat.codes
        df['site_num'] = df['site'].astype('category').cat.codes

        ancova = pg.ancova(
                            data=df,
                            dv=region,
                            covar=["sex_num", "age", "site_num", "TIV", "PTA", roi_vol],
                            between='group'
                            )

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

    df_results = df_results.sort_values(by='corrected_p_val', ascending=True)

    return  df_results, df_vif_age_pta, df_vif_tiv

if __name__ == "__main__":
    perform_roi_stats()