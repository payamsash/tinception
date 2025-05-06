from pathlib import Path
import pandas as pd

def main():
    fname_master = Path.cwd().parent / "data" / "tinception_matched_optimal.xlsx"
    df = pd.read_excel(fname_master)

    df.drop(columns=['...1', 'distance', 'weights', 'subclass'], inplace=True)
    group_dummies = pd.get_dummies(df['group'], prefix='group')
    df = pd.concat([df, group_dummies], axis=1)
    df['age'] = df['age'] - df['age'].mean()
    site_dummies = pd.get_dummies(df['site'], prefix='site')
    df = pd.concat([df, site_dummies], axis=1)
    for col in ["group_CO", "group_TI", "site_cogtail", "site_karolinska", "site_neuropren", "site_talaska"]:
        df[col] = df[col].astype(int)
    df.drop(columns=['site', 'site_tinspect', 'group', 'subject ID'], inplace=True)
    df = df[['group_CO', 'group_TI', 'age', 'sex', 'site_cogtail', 'site_karolinska', 'site_neuropren', 'site_talaska']]
    df.to_csv(Path.cwd().parent / "data" / "design.txt", sep=' ', header=False, index=False)


if __name__ == "__main__":
    main()