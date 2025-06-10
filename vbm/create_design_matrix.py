from pathlib import Path
import pandas as pd

def main():
    fname_master = Path.cwd().parent / "material" / "behavioural" / "tinception_matched_optimal.xlsx"
    df = pd.read_excel(fname_master)

    fname_tiv = Path.cwd().parent / "material" / "VBM" / "tiv_results.txt"
    df_tiv = pd.read_csv(fname_tiv, delimiter='\t')

    df.drop(columns=['...1', 'distance', 'weights', 'subclass'], inplace=True)
    group_dummies = pd.get_dummies(df['group'], prefix='group')
    df = pd.concat([df, group_dummies], axis=1)
    df['age'] = df['age'] - df['age'].mean()
    df['PTA'] = df['PTA'] - df['PTA'].mean()
    site_dummies = pd.get_dummies(df['site'], prefix='site')
    df = pd.concat([df, site_dummies], axis=1)
    for col in ["group_CO", "group_TI", "site_cogtail", "site_karolinska", "site_neuropren", "site_talaska"]:
        df[col] = df[col].astype(int)
    df.drop(columns=['site', 'site_tinspect', 'group', 'subject ID'], inplace=True)
    df["TIV"] = df_tiv["TIV_mm3"]
    df['TIV'] = df['TIV'] - df['TIV'].mean()

    df = df[['group_CO', 'group_TI', 'sex', 'age', 'TIV', 'PTA', 'site_cogtail', 'site_karolinska', 'site_neuropren', 'site_talaska']]
    df.to_csv(Path.cwd().parent / "material" / "VBM" / "design_with_PTA.txt", sep=' ', header=False, index=False)

    df = df[['group_CO', 'group_TI', 'sex', 'age', 'TIV', 'site_cogtail', 'site_karolinska', 'site_neuropren', 'site_talaska']]
    df.to_csv(Path.cwd().parent / "material" / "VBM" / "design_without_PTA.txt", sep=' ', header=False, index=False)

if __name__ == "__main__":
    main()