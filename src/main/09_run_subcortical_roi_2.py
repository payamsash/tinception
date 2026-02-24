from pathlib import Path
import pandas as pd
import numpy as np
import pingouin as pg
import seaborn as sns
from statsmodels.stats.outliers_influence import variance_inflation_factor


def compute_stats(df, structure, covariates, fdr_method="fdr_bh", separate_hemi=None):
    
    df = df.query(f'structure == "{structure}"').copy()

    if structure == "Hippo":
        df["region"] = df["region"] + "-" + df["hemi"]
        hippo_totals = df[df['region'].isin(['Whole_hippocampus-lh', 'Whole_hippocampus-rh'])]
        total_vol_map = hippo_totals.groupby('tinception_id')['volume'].sum()
        exclude = [
                    'Whole_hippocampal_body-lh', 'Whole_hippocampal_head-lh', 'Whole_hippocampus-lh',
                    'Whole_hippocampal_body-rh', 'Whole_hippocampal_head-rh', 'Whole_hippocampus-rh'
                    ]
        regions = [r for r in df['region'].unique() if r not in exclude]
        
    elif structure == "Amygdala":
        df["region"] = df["region"] + "-" + df["hemi"]
        amygdala_totals = df[df['region'].isin(['Whole_amygdala-lh', 'Whole_amygdala-rh'])]
        total_vol_map = amygdala_totals.groupby('tinception_id')['volume'].sum()
        exclude = ['Whole_amygdala-lh', 'Whole_amygdala-rh']
        regions = [r for r in df['region'].unique() if r not in exclude]

    elif structure == "Thalamic-nuclei":
        thalamus_totals = df[df['region'].isin(['Left-Whole_thalamus', 'Right-Whole_thalamus'])]
        total_vol_map = thalamus_totals.groupby('tinception_id')['volume'].sum()
        exclude = ['Left-Whole_thalamus', 'Right-Whole_thalamus']
        regions = [r for r in df['region'].unique() if r not in exclude]
        df.drop(columns="hemi", inplace=True)
    
    else:
        regions = list(df['region'].unique())
        total_vol_map = None

    # Map the total volume for use as a potential covariate
    if total_vol_map is not None:
        df['total_struct_vol'] = df['tinception_id'].map(total_vol_map)

    for col in ['site', 'sex']:
        if col in df.columns:
            df[col] = df[col].astype('category').cat.codes

    # VIF Check
    vif_cols = ['age', 'PTA', 'TIV']
    if 'total_struct_vol' in df.columns:
        vif_cols.append('total_struct_vol')
        
    vif_temp = df[vif_cols + ['volume']].dropna()
    if not vif_temp.empty:
        vif_df = pd.DataFrame()
        vif_df["feature"] = vif_cols
        vif_df["VIF"] = [variance_inflation_factor(vif_temp[vif_cols].values, i) for i in range(len(vif_cols))]

    ## compute ratio
    df["ratio"] = df["total_struct_vol"] / df["TIV"]

    results = []
    for reg in regions:
        sub_df = df[df['region'] == reg].dropna(subset=['volume', 'group'] + covariates)

        if separate_hemi is None:
            pass
        else:
            sub_df = sub_df.query(f'hemi == "{separate_hemi}"')

        if sub_df.empty: continue
        # ANCOVA
        ancova = pg.ancova(data=sub_df, dv='volume', covar=covariates, between='group')
        res_row = ancova[ancova['Source'] == 'group'].iloc[0]
        
        co_vols = sub_df[sub_df['group'] == 'CO']['volume']
        ti_vols = sub_df[sub_df['group'] == 'TI']['volume']
        
        n1, n2 = len(co_vols), len(ti_vols)
        s1, s2 = co_vols.std(), ti_vols.std()
        pooled_std = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
        d = (ti_vols.mean() - co_vols.mean()) / pooled_std
        
        results.append({
            'Region': reg,
            'N_CO': n1,
            'N_TI': n2,
            'Mean_CO': co_vols.mean(),
            'Mean_TI': ti_vols.mean(),
            'Diff': ti_vols.mean() - co_vols.mean(),
            'Cohen_d': d,
            'F_stat': res_row['F'],
            'p_val': res_row['p-unc'],
            'partial_eta_sq': res_row['np2']
        })

    # finalizing
    final_df = pd.DataFrame(results)
    if not final_df.empty:
        _, final_df['p_fdr'] = pg.multicomp(final_df['p_val'].values, method=fdr_method)
        final_df['Significant'] = final_df['p_fdr'] < 0.05
        final_df["structure"] = [structure] * len(final_df)

        ## rename and order cols

        return final_df.sort_values('p_val'), vif_df
    
    return "No results found."

if __name__ == "__main__":

    ## check missing files
    tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
    subjects_dir = tinception_dir / "subjects_fs_dir"
    vbm_norm_dir = tinception_dir / "vbm_norm"
    df_master = pd.read_csv(vbm_norm_dir / "covars.csv")
    cols_to_keep = ["subject_ID", "group", "tinception_id", "age", "sex", "site", "PTA", "THI", "TIV"]
    df_master = df_master[cols_to_keep]
    subject_ids = df_master["subject_ID"].values

    missing_subjects = []
    completed_ids = []
    df_hs, df_as, df_ts = [], [], []
    hemis = ["lh", "rh"]
    for folder in sorted(subjects_dir.iterdir()):
        if folder.is_dir():
            subject_id = folder.name
            
            if subject_id not in list(subject_ids):
                continue
            else:
                tin_id = df_master.loc[df_master['subject_ID'] == subject_id, 'tinception_id'].values[0]
                site = df_master.loc[df_master['subject_ID'] == subject_id, 'site'].values[0]
            
            if site == "antinomics":
                mode_hippo = "T1-T2"
            else:
                mode_hippo = "T1"

            mri_dir = subjects_dir / subject_id / "mri"
            hippo_file = mri_dir / f"lh.hippoSfVolumes-{mode_hippo}.v22.txt"
            amygdala_file = mri_dir / f"lh.amygNucVolumes-{mode_hippo}.v22.txt"
            thalamic_file = mri_dir / "ThalamicNuclei.v13.T1.volumes.txt"

            if not hippo_file.is_file():
                print(f"missing hippo file for subject {subject_id} -> {tin_id}")
            else:
                for hemi in hemis:
                    hippo_file = mri_dir / f"{hemi}.hippoSfVolumes-{mode_hippo}.v22.txt"
                    df_h = pd.read_csv(hippo_file, sep=" ", names=["region", "volume"], header=None)
                    df_h["subject_ID"] = [subject_id] * len(df_h)
                    df_h["hemi"] = [hemi] * len(df_h)
                    df_hs.append(df_h)

            if not amygdala_file.is_file():
                print(f"missing amygdala file for subject {subject_id} -> {tin_id}")
            else:
                for hemi in hemis:
                    amygdala_file = mri_dir / f"{hemi}.amygNucVolumes-{mode_hippo}.v22.txt"
                    df_a = pd.read_csv(amygdala_file, sep=" ", names=["region", "volume"], header=None)
                    df_a["subject_ID"] = [subject_id] * len(df_a)
                    df_a["hemi"] = [hemi] * len(df_a)
                    df_as.append(df_a)

            if not thalamic_file.is_file():
                print(f"missing thalamic file for subject {subject_id} -> {tin_id}")
                missing_subjects.append(subject_id)
            else:
                df_t = pd.read_csv(thalamic_file, sep=" ", names=["region", "volume"], header=None)
                df_t["subject_ID"] = [subject_id] * len(df_t)
                df_ts.append(df_t)


    ## compute stats
    df_h = pd.concat(df_hs, axis=0)
    df_a = pd.concat(df_as, axis=0)
    df_t = pd.concat(df_ts, axis=0)

    dfs = [df_h, df_a, df_t]    
    df_h["structure"] = ["Hippo"] * len(df_h)
    df_a["structure"] = ["Amygdala"] * len(df_a)
    df_t["structure"] = ["Thalamic-nuclei"] * len(df_t)
    df = pd.concat(dfs, axis=0)
    df = pd.merge(df, df_master, on='subject_ID', how='left')
    subcortical_dir = tinception_dir / "subcortical_roi"
    df.to_csv(subcortical_dir / "master.csv", index=False)


    ## possible covariates are: age, sex, PTA, site, TIV, total_struct_vol, ratio
    covariates = ['age', 'sex', 'site', 'PTA', 'TIV']
    
    df["TIV"] = df["TIV"] - df["TIV"].mean()
    df["PTA"] = df["PTA"] - df["PTA"].mean()
    for structure in ["Hippo", "Amygdala", "Thalamic-nuclei"]:
        df_stat, df_vif = compute_stats(
                                        df,
                                        structure,
                                        covariates=covariates,
                                        separate_hemi=None
                                        )
        df_stat.reset_index(drop=True).to_csv(subcortical_dir / "stats" / f"{structure}.csv")
        df_vif.reset_index(drop=True).to_csv(subcortical_dir / "stats" / f"{structure}_vif.csv")
    
    print("******* Completed sucessfully! ******")