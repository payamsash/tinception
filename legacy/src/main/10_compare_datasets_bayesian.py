from pathlib import Path
import pandas as pd
import pymc as pm
import arviz as az
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

# ── 0. Sample sizes ────────────────────────────────────────────────────────────
N_MAIN_TI, N_MAIN_CO = 285, 284
N_UKB_TI,  N_UKB_CO  = 2180, 2176
D_MIN = 0.05   # minimally relevant effect size

# ── 1. Helper functions ────────────────────────────────────────────────────────

def cohend_se(d, n1, n2):
    """Standard error of Cohen's d (two independent groups)."""
    return np.sqrt((n1 + n2) / (n1 * n2) + d**2 / (2 * (n1 + n2)))

def fstat_to_cohend(F, n1, n2):
    """
    Convert F-statistic (df1=1) to Cohen's d.
    Only valid for one df-numerator F (i.e. two-group comparison).
    Sign is unknown from F alone — we recover it from mean difference.
    """
    t = np.sqrt(F)
    return t * np.sqrt((n1 + n2) / (n1 * n2))   # unsigned; sign added below

def tval_to_cohend(t, n1, n2):
    """Convert t-value to Cohen's d (signed)."""
    return t * np.sqrt((n1 + n2) / (n1 * n2))

# ── 2. Compute Cohen's d and SE for each dataset ───────────────────────────────

def prepare_df(df, n_ti, n_co, d_col=None, t_col=None, f_col=None):
    """
    Add cohen_d (signed) and se_d columns.
    Priority: use existing cohen_d → t-value → F-statistic.
    Sign is always taken from (mean_tinnitus - mean_control).
    """
    df = df.copy()
    sign = np.sign(df["mean_tinnitus"] - df["mean_control"])

    if d_col and d_col in df.columns:
        df["cohen_d_calc"] = df[d_col].abs() * sign   # enforce correct sign

    elif t_col and t_col in df.columns:
        df["cohen_d_calc"] = tval_to_cohend(df[t_col], n_ti, n_co)

    elif f_col and f_col in df.columns:
        df["cohen_d_calc"] = fstat_to_cohend(df[f_col], n_ti, n_co) * sign

    else:
        raise ValueError("No usable effect-size column found.")

    df["se_d"] = cohend_se(df["cohen_d_calc"], n_ti, n_co)
    return df

def meta_fixed(d1, se1, d2, se2):
    w1, w2   = 1/se1**2, 1/se2**2
    d_pool   = (w1*d1 + w2*d2) / (w1 + w2)
    se_pool  = np.sqrt(1 / (w1 + w2))
    z        = d_pool / se_pool
    pval     = 2 * (1 - stats.norm.cdf(abs(z)))
    ci_lo    = d_pool - 1.96*se_pool
    ci_hi    = d_pool + 1.96*se_pool
    return d_pool, se_pool, z, pval, ci_lo, ci_hi

def assign_tier_meta(row, alpha=0.05):
    sig_ukb  = row["padj_ukb"]      < alpha
    sig_main = row["padj_main"]     < alpha
    sig_meta = row["pval_meta_fdr"] < alpha
    big      = abs(row["d_pooled"]) > D_MIN
    same     = row["same_sign"]

    if sig_meta and big and same and sig_ukb and sig_main:
        return "Tier 1: Robust Replication"
    elif sig_meta and big and same:
        return "Tier 2: Consistent (Main Underpowered)"
    elif sig_meta and same:
        return "Tier 3: Consistent but Small Effect"
    elif sig_ukb or sig_main:
        return "Tier 4: Dataset-Specific"
    else:
        return "Tier 5: Non-Significant"
    
def run_bayesian(df, d_min=D_MIN, n_samples=2000, tune=1000):
    records = []
    n_regions = len(df)

    for i, (_, row) in enumerate(df.iterrows()):
        region  = row["brain_label"]
        d_obs   = np.array([row["cohen_d_calc_main"], row["cohen_d_calc_ukb"]])
        se_obs  = np.array([row["se_d_main"],          row["se_d_ukb"]])

        with pm.Model():
            # Shared effect
            theta = pm.Normal("theta", mu=0, sigma=0.5)
            # Between-dataset heterogeneity
            tau   = pm.HalfNormal("tau", sigma=0.2)
            # Dataset-specific deviations
            delta = pm.Normal("delta", mu=0, sigma=tau, shape=2)
            # Likelihood
            pm.Normal("obs", mu=theta + delta, sigma=se_obs, observed=d_obs)

            trace = pm.sample(
                n_samples, tune=tune, chains=4,
                progressbar=False,
                return_inferencedata=True,
                nuts_sampler="numpyro"   # faster; falls back to default if missing
            )

        theta_s         = trace.posterior["theta"].values.flatten()
        p_exceeds_dmin  = float(np.mean(np.abs(theta_s) > d_min))
        p_positive      = float(np.mean(theta_s > 0))
        p_same_sign     = float(max(p_positive, 1 - p_positive))
        hdi             = az.hdi(theta_s, hdi_prob=0.95)

        records.append(dict(
            brain_label     = region,
            d_main          = row["cohen_d_calc_main"],
            d_ukb           = row["cohen_d_calc_ukb"],
            theta_mean      = float(theta_s.mean()),
            theta_hdi_lo    = float(hdi[0]),
            theta_hdi_hi    = float(hdi[1]),
            p_exceeds_dmin  = p_exceeds_dmin,
            p_same_sign     = p_same_sign,
        ))

        print(f"  [{i+1}/{n_regions}] {region}: "
              f"θ={theta_s.mean():.3f}, "
              f"P(|θ|>{d_min})={p_exceeds_dmin:.2f}, "
              f"P(same sign)={p_same_sign:.2f}")

    return pd.DataFrame(records)

def assign_tier_bayes(row, p_eff=0.95, p_eff_mod=0.80, p_sign=0.90):
    high  = row["p_exceeds_dmin"] >= p_eff
    mod   = row["p_exceeds_dmin"] >= p_eff_mod
    cons  = row["p_same_sign"]    >= p_sign

    if high and cons:
        return "Tier 1: High-Confidence Shared Effect"
    elif mod and cons:
        return "Tier 2: Probable Shared Effect"
    elif cons:
        return "Tier 3: Directionally Consistent, Small"
    else:
        return "Tier 4: Inconsistent / No Evidence"
    

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
main_dir = tinception_dir / "subcortical_roi" / "stats"
ukb_dir = tinception_dir / "subcortical_roi" / "ukb"

atlases_main = ["Amygdala", "Hippo", "Thalamic-nuclei"]
atlases_ukb = ["amygdalar_nuclei_stat_bonferroni_results", "hippo_subfields_stat_bonferroni_results", "thalamic_nuclei_stat_bonferroni_results"]

df_finals = []
for atlas_m, atlas_u in zip(atlases_main, atlases_ukb): 

        df1 = pd.read_csv(main_dir / f"{atlas_m}.csv")
        df2 = pd.read_csv(ukb_dir / f"{atlas_u}.csv")

        df1.drop(columns=["Unnamed: 0", "Diff", "structure"], inplace=True)
        df2.drop(columns="Unnamed: 0", inplace=True)

        mapping = {
                "Region": "brain_label",
                "N_CO": "n_control",
                "N_TI": "n_tinnitus",
                "Mean_CO": "mean_control",
                "Mean_TI": "mean_tinnitus",
                "Cohen_d": "cohen_d",
                "p_val": "pval",
                "p_fdr": "pval_adj",
                "Significant": "significant"
        }
        df1.rename(columns=mapping, inplace=True)
        df1.sort_values(by="brain_label", inplace=True)
        df2.sort_values(by="brain_label", inplace=True)
        df2["brain_label"] = (
                df2["brain_label"]
                .str.replace(r"^Volume_of_", "", regex=True)  # remove prefix
                .str.replace(r"_left_hemisphere$", "-lh", regex=True)  # left → -lh
                .str.replace(r"_right_hemisphere$", "-rh", regex=True)  # right → -rh
                .str.replace("_", "-")  # remaining underscores → hyphens
        )

        if atlas_m == "Thalamic-nuclei":
                is_left = df1["brain_label"].str.startswith("Left")
                is_right = df1["brain_label"].str.startswith("Right")

                df1["brain_label"] = (
                df1["brain_label"]
                .str.replace(r"^Left-\s*", "", regex=True)
                .str.replace(r"^Right-\s*", "", regex=True)
                .str.replace(r"\s*\(.*?\)", "", regex=True)  # remove parentheses + content
                .str.strip()
                )

                df1.loc[is_left, "brain_label"] += "-lh"
                df1.loc[is_right, "brain_label"] += "-rh"

        df1_prep = prepare_df(df1, N_MAIN_TI, N_MAIN_CO,
                        d_col="cohen_d", f_col="F_stat")
        df2_prep = prepare_df(df2, N_UKB_TI, N_UKB_CO,
                        d_col="cohen_d", t_col="t")
        df = df1_prep.merge(
        df2_prep,
        on="brain_label",
        suffixes=("_main", "_ukb"),
        how="inner"
        )
        df["same_sign"] = (
        np.sign(df["cohen_d_calc_main"]) == np.sign(df["cohen_d_calc_ukb"])
        )

        print(f"Regions in main : {len(df1_prep)}")
        print(f"Regions in UKB  : {len(df2_prep)}")
        print(f"Shared regions  : {len(df)}")

        meta_rows = []
        for _, row in df.iterrows():
                d_pool, se_pool, z, pval, ci_lo, ci_hi = meta_fixed(
                        row["cohen_d_calc_main"], row["se_d_main"],
                        row["cohen_d_calc_ukb"],  row["se_d_ukb"]
                )
                meta_rows.append(dict(
                        brain_label  = row["brain_label"],
                        d_main       = row["cohen_d_calc_main"],
                        d_ukb        = row["cohen_d_calc_ukb"],
                        pval_main    = row["pval_main"],
                        pval_ukb     = row["pval_ukb"],
                        padj_main    = row["pval_adj_main"],
                        padj_ukb     = row["pval_adj_ukb"],
                        same_sign    = row["same_sign"],
                        d_pooled     = d_pool,
                        se_pooled    = se_pool,
                        z_meta       = z,
                        pval_meta    = pval,
                        ci_lo        = ci_lo,
                        ci_hi        = ci_hi,
                ))

        df_meta = pd.DataFrame(meta_rows)

        # FDR on pooled p-values
        _, df_meta["pval_meta_fdr"], _, _ = multipletests(
        df_meta["pval_meta"], method="fdr_bh"
        )

        df_meta["tier_meta"] = df_meta.apply(assign_tier_meta, axis=1)
        print("\n── Meta-analysis tier summary ──")
        print(df_meta["tier_meta"].value_counts())

        print("\nRunning Bayesian hierarchical model …")
        df_bayes = run_bayesian(df)

        df_bayes["tier_bayes"] = df_bayes.apply(assign_tier_bayes, axis=1)
        print("\n── Bayesian tier summary ──")
        print(df_bayes["tier_bayes"].value_counts())


        df_final = df_meta.merge(
        df_bayes[["brain_label", "theta_mean", "theta_hdi_lo",
                "theta_hdi_hi", "p_exceeds_dmin", "p_same_sign", "tier_bayes"]],
        on="brain_label"
        )

        cols_ordered = [
        "brain_label",
        "d_main", "d_ukb",
        "pval_main", "padj_main",
        "pval_ukb",  "padj_ukb",
        "same_sign",
        "d_pooled", "pval_meta", "pval_meta_fdr",
        "theta_mean", "theta_hdi_lo", "theta_hdi_hi",
        "p_exceeds_dmin", "p_same_sign",
        "tier_meta", "tier_bayes"
        ]
        df_final = df_final[cols_ordered]
        df_finals.append(df_final)

df_final = pd.concat(df_finals, axis=0)

cols_ordered = [
    "brain_label",
    "d_main", "d_ukb",
    "pval_main", "padj_main",
    "pval_ukb",  "padj_ukb",
    "same_sign",
    "d_pooled", "pval_meta", "pval_meta_fdr",
    "theta_mean", "theta_hdi_lo", "theta_hdi_hi",
    "p_exceeds_dmin", "p_same_sign",
    "tier_meta", "tier_bayes"
]
df_final = df_final[cols_ordered]
df_final.to_csv(tinception_dir / "subcortical_roi" / "bayesian" / "bayesian.csv", index=False)