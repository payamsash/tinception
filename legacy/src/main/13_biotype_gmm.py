"""
13_biotype_gmm.py — Tinnitus Biotyping via Gaussian Mixture Models

Strategy
--------
GMM applied directly to 11 convergent subcortical ROI normative Z-scores
(no UMAP distortion before clustering). k selected by BIC; stability validated
by 500-iteration bootstrap ARI. UKB replication uses centroid transfer from the
Main GMM, not independent re-clustering.

Two analyses are run:
  PRIMARY    — n=270, excluding london_2 (scanner artifact; 86.7% of the
               minority biotype in the full run). Cogtail's over-representation
               in B0 is explained by demographics (median age 69, PTA 27.6 dB),
               not scanner effects.
  DIAGNOSTIC — n=285, all sites. Included as supplementary to disclose the
               london_2 artifact and motivate the primary analysis.

Outputs (biotypes_gmm/)
  biotype_assignments_main.csv          full-dataset labels + probabilities
  biotype_stats_main.csv                full-dataset clinical stats
  sensitivity_excl_london_2/
    biotype_assignments_sensitivity.csv primary labels + probabilities
    biotype_stats_sensitivity.csv       primary clinical stats  ← main result
  biotype_assignments_ukb.csv           UKB centroid-transfer labels
  biotype_stats_ukb.csv                 UKB clinical stats

Plots (plots/biotypes_gmm/)
  00_pca_variance.pdf                   explained variance per PC
  01_bic_selection.pdf                  BIC/AIC vs k
  02_bootstrap_ari.pdf                  bootstrap ARI distribution
  03_dendrogram.pdf                     Ward dendrogram coloured by biotype
  04_umap_main.pdf                      UMAP embedding, full dataset
  05_zprofile_main_{raw,partialed}.pdf  mean Z per biotype heatmap
  06_clinical_main.pdf                  boxplots: PTA / THI / age
  07_radar_main_dataset.pdf             spider plot of Z-score profiles
  08_ukb_transfer_summary.pdf           UKB biotype distribution + certainty
  09_site_composition_main.pdf          site % per biotype (london_2 diagnostic)
  sensitivity_excl_london_2/            same set for the primary analysis
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from math import pi

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import mannwhitneyu, kruskal, chi2_contingency
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests

import umap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pingouin as pg

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
TINCEPTION_DIR = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
PLOTS_DIR      = TINCEPTION_DIR / "plots" / "biotypes_gmm"
OUTPUT_DIR     = TINCEPTION_DIR / "biotypes_gmm"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PARTIAL_OUT_PTA   = True   # partial hearing loss out of Z-scores before clustering
PARTIAL_OUT_SITE  = True   # partial site (one-hot) out to remove residual scanner variance
GMM_COV_TYPE      = "full"
MIN_CLUSTER_FRAC  = 0.05   # reject k solutions where any cluster < 5 % of N
PCA_VARIANCE_THRESHOLD = 0.85
K_RANGE      = range(1, 7)
N_BOOT       = 500
BOOT_FRAC    = 0.80
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# Sites to exclude in the primary (sensitivity) analysis.
# london_2 was found to drive 86.7 % of the minority biotype in the diagnostic run.
SENSITIVITY_EXCLUDE_SITES = ["london_2"]

# 11 convergent subcortical regions (Tier 2+3 from cross-dataset forest plot)
REGIONS = [
    "Accessory-Basal-nucleus-lh",
    "Basal-nucleus-lh",
    "Lateral-nucleus-lh",
    "GC-ML-DG-body-rh",
    "presubiculum-head-rh",
    "MDm-lh",
    "MDm-rh",
    "MDl-lh",
    "CeM-rh",
    "CeM-lh",
    "PuA-lh",
]

ATLAS_CONFIG_MAIN = {
    "Amygdala":        TINCEPTION_DIR / "subcortical_roi" / "norm_models_main" / "Amygdala"        / "results" / "Z_test.csv",
    "Hippo":           TINCEPTION_DIR / "subcortical_roi" / "norm_models_main" / "Hippo"           / "results" / "Z_test.csv",
    "Thalamic-nuclei": TINCEPTION_DIR / "subcortical_roi" / "norm_models_main" / "Thalamic-nuclei" / "results" / "Z_test.csv",
}

ATLAS_CONFIG_UKB = {
    "amygdalar_nuclei": TINCEPTION_DIR / "subcortical_roi" / "norm_models_ukb" / "amygdalar_nuclei" / "results" / "Z_test.csv",
    "hippo_subfields":  TINCEPTION_DIR / "subcortical_roi" / "norm_models_ukb" / "hippo_subfields"  / "results" / "Z_test.csv",
    "thalamic_nuclei":  TINCEPTION_DIR / "subcortical_roi" / "norm_models_ukb" / "thalamic_nuclei"  / "results" / "Z_test.csv",
}

COVARS_CSV_MAIN = TINCEPTION_DIR / "VBM_design" / "covars.csv"

PALETTE = {0: "#440154", 1: "#21908C", 2: "#FDE725", 3: "#31688e", 4: "#35b779"}
sns.set_theme(style="whitegrid", rc={"axes.spines.right": False, "axes.spines.top": False})


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def clean_region_names(cols: pd.Index) -> pd.Index:
    """Normalise region name formats to 'X-lh' / 'X-rh' (handles UKB, Main, and clean formats)."""
    return (
        cols.str.removeprefix("Volume_of_")
            .str.replace("_right_hemisphere$", "-rh", regex=True)
            .str.replace("_left_hemisphere$",  "-lh", regex=True)
            .str.replace("_", "-", regex=False)
            .str.replace(r"^Left-(.+)$",  r"\1-lh", regex=True)
            .str.replace(r"^Right-(.+)$", r"\1-rh", regex=True)
    )


def partial_out(Z: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """Regress covariates out of every column of Z; return residuals."""
    C = covariates if covariates.ndim == 2 else covariates.reshape(-1, 1)
    residuals = np.empty_like(Z, dtype=float)
    for j in range(Z.shape[1]):
        valid = ~np.isnan(Z[:, j]) & ~np.any(np.isnan(C), axis=1)
        model = LinearRegression().fit(C[valid], Z[valid, j])
        residuals[:, j] = Z[:, j] - model.predict(C)
    return residuals


def build_covariate_matrix(df: pd.DataFrame, pta_col: str = "PTA") -> np.ndarray | None:
    """Build a covariate matrix (PTA + site one-hot dummies)."""
    parts = []
    if PARTIAL_OUT_PTA and pta_col in df.columns:
        parts.append(df[[pta_col]].values.astype(float))
    if PARTIAL_OUT_SITE and "site" in df.columns:
        dummies = pd.get_dummies(df["site"], drop_first=True).astype(float)
        parts.append(dummies.values)
    return np.hstack(parts) if parts else None


def apply_pca(X: np.ndarray, threshold: float = PCA_VARIANCE_THRESHOLD):
    """Fit PCA retaining enough components to explain `threshold` variance."""
    pca_full = PCA(random_state=RANDOM_STATE).fit(X)
    cum_var  = np.cumsum(pca_full.explained_variance_ratio_)
    n_keep   = max(2, int(np.searchsorted(cum_var, threshold)) + 1)
    pca      = PCA(n_components=n_keep, random_state=RANDOM_STATE).fit(X)
    print(f"  PCA: kept {n_keep} components → {pca.explained_variance_ratio_.sum():.1%} variance")
    return pca.transform(X), pca


def save_fig(fig, name: str, out_dir: Path = PLOTS_DIR):
    fpath = out_dir / name
    fig.savefig(fpath, dpi=300, bbox_inches="tight")
    print(f"  Saved → {fpath.name}")


# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────────────────────────────────────
def load_main_features():
    df_meta = pd.read_csv(COVARS_CSV_MAIN)
    df_meta = df_meta[["tinception_id", "group", "PTA", "THI", "age", "sex", "site"]]

    frames = []
    for atlas, fpath in ATLAS_CONFIG_MAIN.items():
        if not fpath.exists():
            print(f"  [WARN] Main – missing: {fpath}")
            continue
        df = pd.read_csv(fpath)
        df.columns = clean_region_names(pd.Index(df.columns))
        for old, new in [("subject-ids", "tinception_id"), ("subject_ids", "tinception_id")]:
            if old in df.columns:
                df.rename(columns={old: new}, inplace=True)
        df = df.merge(df_meta, on="tinception_id", how="inner").query('group == "TI"')
        region_cols = [r for r in REGIONS if r in df.columns]
        frames.append(df.set_index("tinception_id")[region_cols + ["PTA", "THI", "age", "sex", "site"]])

    if not frames:
        raise FileNotFoundError("No main-dataset Z-score files found.")

    df_out = frames[0]
    meta_cols = ["PTA", "THI", "age", "sex", "site"]
    for f in frames[1:]:
        df_out = df_out.join(f.drop(columns=[c for c in meta_cols if c in f.columns]), how="inner")

    available = [r for r in REGIONS if r in df_out.columns]
    missing   = [r for r in REGIONS if r not in df_out.columns]
    if missing:
        print(f"  [WARN] Regions not found: {missing}")
    print(f"  [Main] TI subjects = {len(df_out)}, regions = {len(available)}/{len(REGIONS)}")
    return df_out, available


def load_ukb_features():
    frames = []
    for atlas, fpath in ATLAS_CONFIG_UKB.items():
        if not fpath.exists():
            print(f"  [WARN] UKB – missing: {fpath}")
            continue
        df = pd.read_csv(fpath)
        df.columns = clean_region_names(pd.Index(df.columns))
        for old, new in [("subject-ids", "subject_id"), ("subject_ids", "subject_id")]:
            if old in df.columns:
                df.rename(columns={old: new}, inplace=True)
        region_cols = [r for r in REGIONS if r in df.columns]
        if not region_cols:
            continue
        frames.append(df.set_index("subject_id")[region_cols])

    if not frames:
        print("  [INFO] UKB data not available — skipping.")
        return None, []

    df_out = frames[0]
    for f in frames[1:]:
        df_out = df_out.join(f, how="inner")

    # Try to attach SRT
    for cpath in [TINCEPTION_DIR / "subcortical_roi" / "ukb_covars.csv",
                  TINCEPTION_DIR / "VBM_design" / "ukb_covars.csv"]:
        if cpath.exists():
            df_cov = pd.read_csv(cpath, index_col=0)
            if "srt" in df_cov.columns:
                df_out = df_out.join(df_cov[["srt"]], how="left")
            break

    available = [r for r in REGIONS if r in df_out.columns]
    print(f"  [UKB]  TI subjects = {len(df_out)}, regions = {len(available)}/{len(REGIONS)}")
    return df_out, available


# ─────────────────────────────────────────────────────────────────────────────
# GMM
# ─────────────────────────────────────────────────────────────────────────────
def select_k_bic(X: np.ndarray):
    rows, models = [], []
    n = len(X)
    for k in K_RANGE:
        gmm = GaussianMixture(n_components=k, covariance_type=GMM_COV_TYPE,
                              n_init=20, max_iter=500, random_state=RANDOM_STATE)
        gmm.fit(X)
        labels_tmp = gmm.predict(X)
        counts     = np.bincount(labels_tmp, minlength=k)
        min_n      = int(counts.min())
        degenerate = min_n < int(MIN_CLUSTER_FRAC * n)
        rows.append({"k": k, "BIC": gmm.bic(X), "AIC": gmm.aic(X),
                     "min_cluster_n": min_n, "degenerate": degenerate})
        models.append(gmm)

    df_ic = pd.DataFrame(rows)
    valid  = df_ic[~df_ic["degenerate"]]
    best_k = int(valid.loc[valid["BIC"].idxmin(), "k"]) if not valid.empty else 1
    return df_ic, best_k, models[best_k - min(K_RANGE)]


def bootstrap_ari(X: np.ndarray, k: int) -> np.ndarray:
    gmm_full   = GaussianMixture(n_components=k, covariance_type=GMM_COV_TYPE,
                                  n_init=20, random_state=RANDOM_STATE)
    labels_ref = gmm_full.fit_predict(X)
    n = len(X)
    aris = []
    for _ in range(N_BOOT):
        idx     = np.random.choice(n, size=int(n * BOOT_FRAC), replace=False)
        gmm_sub = GaussianMixture(n_components=k, covariance_type=GMM_COV_TYPE,
                                   n_init=10, random_state=None)
        gmm_sub.fit(X[idx])
        aris.append(adjusted_rand_score(labels_ref[idx], gmm_sub.predict(X[idx])))
    return np.array(aris)


# ─────────────────────────────────────────────────────────────────────────────
# CLINICAL CHARACTERIZATION
# ─────────────────────────────────────────────────────────────────────────────
def characterize_biotypes(df: pd.DataFrame, labels: np.ndarray,
                           cont_vars: list, cat_vars: list) -> pd.DataFrame:
    df = df.copy()
    df["Biotype"] = labels
    k = len(np.unique(labels))
    if k == 1:
        print("  [INFO] k=1: no between-biotype comparisons.")
        return pd.DataFrame()

    rows = []
    for var in cont_vars:
        if var not in df.columns:
            continue
        groups = [df.loc[df["Biotype"] == bt, var].dropna().values for bt in range(k)]
        stat, p = mannwhitneyu(*groups, alternative="two-sided") if k == 2 else kruskal(*groups)
        test    = "Mann-Whitney U" if k == 2 else "Kruskal-Wallis"
        row = {"Variable": var, "Test": test, "Statistic": round(stat, 3), "p": round(p, 5)}
        for bt in range(k):
            g = groups[bt]
            row[f"B{bt}_n"]      = len(g)
            row[f"B{bt}_median"] = round(float(np.median(g)), 2)
            row[f"B{bt}_IQR"]    = f"[{np.percentile(g,25):.2f}–{np.percentile(g,75):.2f}]"
        rows.append(row)

    for var in cat_vars:
        if var not in df.columns:
            continue
        ct = pd.crosstab(df["Biotype"], df[var])
        chi2, p, _, _ = chi2_contingency(ct)
        rows.append({"Variable": var, "Test": "Chi-square",
                     "Statistic": round(chi2, 3), "p": round(p, 5)})

    df_stats = pd.DataFrame(rows)
    if not df_stats.empty:
        _, padj, _, _ = multipletests(df_stats["p"].fillna(1), method="fdr_bh")
        df_stats["p_fdr"]        = padj.round(5)
        df_stats["Significant"]  = df_stats["p_fdr"] < 0.05
    return df_stats


# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────
def plot_pca_variance(pca: PCA, out_dir: Path = PLOTS_DIR):
    evr  = pca.explained_variance_ratio_
    cumv = np.cumsum(evr)
    fig, ax = plt.subplots(figsize=(5, 3.5), layout="tight")
    ax.bar(range(1, len(evr) + 1), evr * 100, color=PALETTE[0], alpha=0.7, label="Per component")
    ax.step(range(1, len(evr) + 1), cumv * 100, where="mid", color=PALETTE[1], lw=2, label="Cumulative")
    ax.axhline(PCA_VARIANCE_THRESHOLD * 100, color="grey", ls="--", lw=1.5,
               label=f"{PCA_VARIANCE_THRESHOLD:.0%} threshold")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_xticks(range(1, len(evr) + 1))
    ax.set_title(f"PCA — {len(evr)} components ({cumv[-1]:.1%} variance)")
    ax.legend(frameon=False)
    save_fig(fig, "00_pca_variance.pdf", out_dir)
    plt.close(fig)


def plot_bic(df_ic: pd.DataFrame, best_k: int, out_dir: Path = PLOTS_DIR):
    fig, ax = plt.subplots(figsize=(5, 3.5), layout="tight")
    ax.plot(df_ic["k"], df_ic["BIC"], "o-",  color=PALETTE[0], lw=2, label="BIC")
    ax.plot(df_ic["k"], df_ic["AIC"], "s--", color=PALETTE[1], lw=2, label="AIC")
    ax.axvline(best_k, color="grey", ls=":", lw=1.8, label=f"Best k={best_k}")
    ax.set_xlabel("Number of components (k)")
    ax.set_ylabel("Information criterion")
    ax.set_xticks(list(K_RANGE))
    ax.legend(frameon=False)
    ax.set_title("GMM model selection — BIC / AIC")
    save_fig(fig, "01_bic_selection.pdf", out_dir)
    plt.close(fig)


def plot_bootstrap_ari(aris: np.ndarray, best_k: int, out_dir: Path = PLOTS_DIR):
    med = np.median(aris)
    q25, q75 = np.percentile(aris, [25, 75])
    fig, ax = plt.subplots(figsize=(5, 3.5), layout="tight")
    sns.histplot(aris, bins=30, color=PALETTE[0], stat="density", ax=ax)
    ax.axvline(med,  color=PALETTE[1], lw=2, label=f"Median ARI = {med:.3f}")
    ax.axvline(0.60, color="crimson",  lw=1.5, ls="--", label="ARI = 0.60 threshold")
    ax.fill_betweenx([0, 1], q25, q75, alpha=0.15, color=PALETTE[1],
                     label=f"IQR [{q25:.2f}–{q75:.2f}]")
    ax.set_xlabel("Adjusted Rand Index")
    ax.set_ylabel("Density")
    ax.set_title(f"Bootstrap stability  (k={best_k}, n={N_BOOT})")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    save_fig(fig, "02_bootstrap_ari.pdf", out_dir)
    plt.close(fig)
    print(f"  Bootstrap ARI: median={med:.3f}, IQR=[{q25:.3f}, {q75:.3f}], "
          f"≥0.60: {(aris>=0.60).mean():.1%}")


def plot_dendrogram(X: np.ndarray, labels: np.ndarray, out_dir: Path = PLOTS_DIR):
    Z_link = linkage(X, method="ward")
    k = len(np.unique(labels))
    fig, ax = plt.subplots(figsize=(max(12, len(X) * 0.08), 4), layout="tight")
    dendrogram(Z_link, ax=ax, leaf_rotation=90, leaf_font_size=0,
               color_threshold=0, above_threshold_color="lightgrey", no_labels=True)
    ax.set_title("Ward hierarchical clustering — TI subjects (site+PTA partialed Z)")
    ax.set_ylabel("Linkage distance")
    ax.set_xlabel(f"Subjects  (n={len(X)})")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=PALETTE[b], label=f"Biotype {b}") for b in range(k)],
              frameon=False, loc="upper right")
    save_fig(fig, "03_dendrogram.pdf", out_dir)
    plt.close(fig)


def plot_umap(X: np.ndarray, labels: np.ndarray, probs: np.ndarray,
              title: str, fname_prefix: str, out_dir: Path = PLOTS_DIR):
    reducer = umap.UMAP(n_neighbors=15, n_components=2, metric="euclidean",
                        min_dist=0.15, random_state=RANDOM_STATE)
    emb = reducer.fit_transform(X)
    k = len(np.unique(labels))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="tight")
    ax = axes[0]
    for bt in range(k):
        mask = labels == bt
        ax.scatter(emb[mask, 0], emb[mask, 1], c=PALETTE[bt],
                   label=f"Biotype {bt}  (n={mask.sum()})",
                   s=55, alpha=0.85, edgecolors="k", linewidths=0.4)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"{title} — GMM biotypes")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    certainty = probs.max(axis=1)
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=certainty, cmap="viridis",
                    s=55, alpha=0.85, edgecolors="k", linewidths=0.4, vmin=0.5, vmax=1.0)
    plt.colorbar(sc, ax=ax, label="GMM posterior probability", shrink=0.85)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"{title} — assignment certainty")

    save_fig(fig, f"04_umap_{fname_prefix}.pdf", out_dir)
    plt.close(fig)
    return emb


def plot_zprofile_heatmap(Z: np.ndarray, labels: np.ndarray,
                           region_names: list, title: str, fname_prefix: str,
                           out_dir: Path = PLOTS_DIR):
    df_z = pd.DataFrame(Z, columns=region_names)
    df_z["Biotype"] = labels
    means = df_z.groupby("Biotype")[region_names].mean()
    fig, ax = plt.subplots(
        figsize=(max(9, len(region_names) * 0.75), 2.8 + 0.4 * len(means)),
        layout="tight"
    )
    sns.heatmap(means, cmap="vlag", center=0, annot=True, fmt=".2f",
                linewidths=0.5, ax=ax, cbar_kws={"label": "Mean Z-score", "shrink": 0.7})
    ax.set_title(f"Mean normative Z-score per biotype — {title}")
    ax.set_xlabel(""); ax.set_ylabel("Biotype")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", fontsize=8, style="italic")
    save_fig(fig, f"05_zprofile_{fname_prefix}.pdf", out_dir)
    plt.close(fig)


def plot_clinical_comparison(df: pd.DataFrame, labels: np.ndarray,
                              cont_vars: list, title: str, fname_prefix: str,
                              out_dir: Path = PLOTS_DIR):
    df = df.copy()
    df["Biotype"] = [f"Biotype {l}" for l in labels]
    k = len(np.unique(labels))
    pal = {f"Biotype {i}": PALETTE[i] for i in range(k)}
    present = [v for v in cont_vars if v in df.columns]
    if not present:
        return
    fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 4.2), layout="tight")
    if len(present) == 1:
        axes = [axes]
    for ax, var in zip(axes, present):
        order = [f"Biotype {i}" for i in range(k)]
        sns.boxplot(data=df, x="Biotype", y=var, palette=pal, width=0.5,
                    linewidth=1.8, fill=False, showfliers=False, ax=ax, order=order)
        sns.stripplot(data=df, x="Biotype", y=var, palette=pal, jitter=0.2,
                      alpha=0.45, size=4.5, edgecolor="k", linewidth=0.3, ax=ax, order=order)
        ax.set_title(var); ax.set_xlabel("")
        for sp in ax.spines.values():
            sp.set_edgecolor("k"); sp.set_linewidth(1.5)
    fig.suptitle(title, fontsize=12, y=1.02)
    save_fig(fig, f"06_clinical_{fname_prefix}.pdf", out_dir)
    plt.close(fig)


def plot_radar(Z: np.ndarray, labels: np.ndarray, region_names: list,
               title: str, out_dir: Path = PLOTS_DIR):
    df_z = pd.DataFrame(Z, columns=region_names)
    df_z["Biotype"] = labels
    means = df_z.groupby("Biotype")[region_names].mean()
    k = len(means)
    N = len(region_names)
    angles = [n / N * 2 * pi for n in range(N)] + [0]

    fig, ax = plt.subplots(figsize=(5.5, 5.5), subplot_kw=dict(polar=True), layout="tight")
    ax.set_yticklabels([])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(region_names, size=7.5)
    for bt in range(k):
        vals = means.loc[bt].tolist() + [means.loc[bt].tolist()[0]]
        ax.plot(angles, vals, linewidth=2.2, label=f"Biotype {bt}", color=PALETTE[bt])
        ax.fill(angles, vals, color=PALETTE[bt], alpha=0.18)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), frameon=False)
    ax.set_title(f"Z-score profiles — {title}", pad=20)
    save_fig(fig, f"07_radar_{title.replace(' ', '_').lower()}.pdf", out_dir)
    plt.close(fig)


def plot_ukb_transfer(labels_ukb: np.ndarray, probs_ukb: np.ndarray,
                       out_dir: Path = PLOTS_DIR):
    k = probs_ukb.shape[1]
    counts    = pd.Series(labels_ukb).value_counts().sort_index()
    pcts      = counts / counts.sum() * 100
    certainty = probs_ukb.max(axis=1)
    high_cert = (certainty >= 0.70).mean()

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="tight")

    ax = axes[0]
    bars = ax.bar([f"Biotype {i}" for i in counts.index], pcts.values,
                  color=[PALETTE[i] for i in counts.index], edgecolor="k", linewidth=1.2)
    ax.bar_label(bars, fmt="%.1f%%", fontsize=9, padding=3)
    ax.set_ylabel("Proportion (%)")
    ax.set_title("UKB — biotype distribution\n(transferred from Main GMM)")
    ax.set_ylim(0, pcts.max() * 1.2)

    ax = axes[1]
    for bt in range(k):
        mask = labels_ukb == bt
        if mask.sum() > 1:
            sns.kdeplot(certainty[mask], ax=ax, color=PALETTE[bt],
                        label=f"Biotype {bt}", fill=True, alpha=0.25)
    ax.axvline(0.70, color="crimson", lw=1.8, ls="--", label="0.70 threshold")
    ax.set_xlabel("GMM posterior probability")
    ax.set_ylabel("Density")
    ax.set_title(f"UKB — assignment certainty\n({high_cert:.1%} subjects ≥ 0.70)")
    ax.legend(frameon=False)

    ax = axes[2]
    ax.hist(certainty, bins=30, color=PALETTE[0], edgecolor="white", linewidth=0.5)
    ax.axvline(0.70, color="crimson", lw=1.8, ls="--")
    ax.axvline(np.median(certainty), color=PALETTE[1], lw=2,
               label=f"Median = {np.median(certainty):.2f}")
    ax.set_xlabel("GMM posterior probability")
    ax.set_ylabel("Count")
    ax.set_title("UKB — certainty distribution")
    ax.legend(frameon=False)

    save_fig(fig, "08_ukb_transfer_summary.pdf", out_dir)
    plt.close(fig)


def plot_site_composition(df: pd.DataFrame, labels: np.ndarray,
                           title: str, out_dir: Path = PLOTS_DIR):
    """Stacked bar chart of site % per biotype — key diagnostic for site artifacts."""
    df = df.copy()
    df["Biotype"] = labels
    k = len(np.unique(labels))
    ct = pd.crosstab(df["site"], df["Biotype"], normalize="columns") * 100
    ct.columns = [f"Biotype {c}" for c in ct.columns]

    fig, ax = plt.subplots(figsize=(max(6, k * 2.5), 4.5), layout="tight")
    ct.T.plot(kind="bar", stacked=True, ax=ax,
              colormap="tab20", edgecolor="k", linewidth=0.5, width=0.5)
    ax.set_ylabel("% of biotype")
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.set_title(f"Site composition per biotype — {title}")
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, fontsize=8)
    ax.set_ylim(0, 110)
    save_fig(fig, "09_site_composition.pdf", out_dir)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# CORE PIPELINE  (reused for full and sensitivity runs)
# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(df: pd.DataFrame, regions: list,
                 tag: str, out_dir: Path, plots_dir: Path):
    """
    Full GMM biotyping pipeline on a given TI-subject DataFrame.
    Returns (labels, probs, pca_obj, scaler_obj, Z_raw, Z_proc, df_ic, best_k, aris).
    """
    plots_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Partial out
    Z_raw  = df[regions].values.astype(float)
    cov    = build_covariate_matrix(df.reset_index())
    Z_proc = partial_out(Z_raw, cov) if cov is not None else Z_raw.copy()

    # Scale + PCA
    scaler   = StandardScaler()
    X_scaled = scaler.fit_transform(Z_proc)
    X_pca, pca_obj = apply_pca(X_scaled)
    plot_pca_variance(pca_obj, plots_dir)

    # GMM BIC selection
    print(f"\n  Fitting GMM (cov='{GMM_COV_TYPE}', min_cluster={MIN_CLUSTER_FRAC:.0%}) ...")
    df_ic, best_k, best_gmm = select_k_bic(X_pca)
    print(df_ic.to_string(index=False))
    print(f"\n  ▶ Best k = {best_k}")
    plot_bic(df_ic, best_k, plots_dir)

    # Bootstrap stability
    print(f"\n  Bootstrap stability (k={best_k}, {N_BOOT} iterations) ...")
    aris = bootstrap_ari(X_pca, best_k)
    plot_bootstrap_ari(aris, best_k, plots_dir)

    # Final labels
    labels = best_gmm.predict(X_pca)
    probs  = best_gmm.predict_proba(X_pca)
    sil    = silhouette_score(X_pca, labels) if best_k > 1 else float("nan")
    print(f"\n  Silhouette score = {sil:.3f}")
    for bt, cnt in pd.Series(labels).value_counts().sort_index().items():
        print(f"    Biotype {bt}: n={cnt}  ({cnt/len(labels):.1%})")

    # Site composition diagnostic
    if "site" in df.columns:
        plot_site_composition(df.reset_index(), labels, title=tag, out_dir=plots_dir)

    # Plots
    plot_dendrogram(X_pca, labels, plots_dir)
    plot_umap(X_pca, labels, probs, title=tag, fname_prefix="main", out_dir=plots_dir)
    plot_zprofile_heatmap(Z_raw,  labels, regions,
                          title=f"{tag} (raw Z)",      fname_prefix="raw",      out_dir=plots_dir)
    plot_zprofile_heatmap(Z_proc, labels, regions,
                          title=f"{tag} (partialed Z)", fname_prefix="partialed", out_dir=plots_dir)
    plot_radar(Z_raw, labels, regions, title=tag, out_dir=plots_dir)
    plot_clinical_comparison(df.reset_index(), labels,
                              cont_vars=["PTA", "THI", "age"],
                              title=f"Clinical variables — {tag}",
                              fname_prefix="clinical", out_dir=plots_dir)

    # Stats
    df_stats = characterize_biotypes(
        df.reset_index(), labels,
        cont_vars=["PTA", "THI", "age"],
        cat_vars=["sex", "site"]
    )
    print(f"\n  Clinical stats — {tag}:")
    print(df_stats.to_string(index=False) if not df_stats.empty else "  (k=1)")
    df_stats.to_csv(out_dir / f"biotype_stats_{tag.replace(' ', '_').lower()}.csv", index=False)

    # Save assignments
    meta_cols = [c for c in ["PTA", "THI", "age", "sex", "site"] if c in df.columns]
    df_assign = df[meta_cols].copy()
    df_assign["Biotype"]   = labels
    df_assign["certainty"] = probs.max(axis=1)
    for bt in range(best_k):
        df_assign[f"prob_B{bt}"] = probs[:, bt]
    df_assign.to_csv(out_dir / f"biotype_assignments_{tag.replace(' ', '_').lower()}.csv")

    print(f"\n  Plots  → {plots_dir}")
    print(f"  Tables → {out_dir}")

    return labels, probs, pca_obj, scaler, Z_raw, Z_proc, df_ic, best_k, aris, sil, best_gmm


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    print("\n" + "=" * 66)
    print("  Tinnitus Biotyping — GMM on Subcortical Normative Z-scores")
    print("=" * 66)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("\n[1] Loading main dataset Z-scores ...")
    df_main, regions_main = load_main_features()

    # ─────────────────────────────────────────────────────────────────────────
    # DIAGNOSTIC ANALYSIS — full dataset (n=285, all sites)
    # Included as supplementary to disclose the london_2 artifact.
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "─" * 66)
    print("  DIAGNOSTIC — full dataset (n=285, all sites)")
    print("─" * 66)

    diag_plots  = PLOTS_DIR  / "diagnostic_all_sites"
    diag_output = OUTPUT_DIR / "diagnostic_all_sites"

    (labels_diag, probs_diag, _, _, Z_raw_diag, _, df_ic_diag,
     best_k_diag, aris_diag, sil_diag, gmm_diag) = run_pipeline(
        df_main, regions_main,
        tag="Diagnostic (all sites)",
        out_dir=diag_output, plots_dir=diag_plots
    )

    # ─────────────────────────────────────────────────────────────────────────
    # PRIMARY ANALYSIS — excluding london_2 (n=270)
    # london_2 drove 86.7 % of the minority biotype in the diagnostic run,
    # consistent with a scanner/protocol difference rather than biology.
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "─" * 66)
    print(f"  PRIMARY — excluding {SENSITIVITY_EXCLUDE_SITES}  (n={len(df_main) - df_main['site'].isin(SENSITIVITY_EXCLUDE_SITES).sum()})")
    print("─" * 66)

    df_primary = df_main[~df_main["site"].isin(SENSITIVITY_EXCLUDE_SITES)].copy()
    n_removed  = len(df_main) - len(df_primary)
    print(f"\n  Removed {n_removed} subjects from: {SENSITIVITY_EXCLUDE_SITES}")

    prim_plots  = PLOTS_DIR  / "primary"
    prim_output = OUTPUT_DIR / "primary"

    (labels_prim, probs_prim, pca_prim, scaler_prim, Z_raw_prim, Z_proc_prim,
     df_ic_prim, best_k_prim, aris_prim, sil_prim, gmm_prim) = run_pipeline(
        df_primary, regions_main,
        tag="Primary (excl. london_2)",
        out_dir=prim_output, plots_dir=prim_plots
    )

    # ─────────────────────────────────────────────────────────────────────────
    # UKB CENTROID TRANSFER  (from primary GMM)
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "─" * 66)
    print("  UKB CENTROID TRANSFER  (primary GMM → UKB subjects)")
    print("─" * 66)

    df_ukb, regions_ukb = load_ukb_features()

    if df_ukb is not None and regions_ukb:
        common  = [r for r in regions_main if r in regions_ukb]
        missing = [r for r in regions_main if r not in regions_ukb]
        print(f"  Common regions: {len(common)}/{len(regions_main)}")
        if missing:
            print(f"  Missing in UKB: {missing}")

        Z_ukb = df_ukb[common].values.astype(float)

        if PARTIAL_OUT_PTA and "srt" in df_ukb.columns:
            Z_ukb = partial_out(Z_ukb, df_ukb["srt"].values.astype(float))
            print("  Partialed out SRT from UKB Z-scores.")
        else:
            print("  [WARN] SRT column not found — skipping PTA partialing for UKB.")

        # Apply PRIMARY scaler + PCA (common-region subset) to UKB
        main_idx      = [regions_main.index(r) for r in common]
        Z_main_comm   = Z_proc_prim[:, main_idx]
        scaler_c      = StandardScaler().fit(Z_main_comm)
        _, pca_c      = apply_pca(scaler_c.transform(Z_main_comm))
        X_ukb_pca     = pca_c.transform(scaler_c.transform(Z_ukb))

        # Refit GMM on primary (common regions) → predict UKB
        gmm_t = GaussianMixture(n_components=best_k_prim, covariance_type=GMM_COV_TYPE,
                                 n_init=20, random_state=RANDOM_STATE)
        gmm_t.fit(pca_c.transform(scaler_c.transform(Z_main_comm)))
        labels_ukb = gmm_t.predict(X_ukb_pca)
        probs_ukb  = gmm_t.predict_proba(X_ukb_pca)

        print(f"\n  UKB biotype distribution (transferred):")
        for bt, cnt in pd.Series(labels_ukb).value_counts().sort_index().items():
            print(f"    Biotype {bt}: n={cnt}  ({cnt/len(labels_ukb):.1%})")
        print(f"  Certainty ≥ 0.70: {(probs_ukb.max(axis=1)>=0.70).mean():.1%}")

        ukb_plots  = PLOTS_DIR  / "ukb_transfer"
        ukb_output = OUTPUT_DIR / "ukb_transfer"
        ukb_plots.mkdir(parents=True, exist_ok=True)
        ukb_output.mkdir(parents=True, exist_ok=True)

        plot_umap(X_ukb_pca, labels_ukb, probs_ukb,
                  title="UKB transfer", fname_prefix="ukb", out_dir=ukb_plots)
        plot_zprofile_heatmap(Z_ukb, labels_ukb, common,
                              title="UKB (raw Z)", fname_prefix="ukb_raw", out_dir=ukb_plots)
        plot_radar(Z_ukb, labels_ukb, common, title="UKB transfer", out_dir=ukb_plots)
        plot_ukb_transfer(labels_ukb, probs_ukb, out_dir=ukb_plots)

        ukb_cont = [c for c in ["srt", "age"] if c in df_ukb.columns]
        plot_clinical_comparison(df_ukb.reset_index(), labels_ukb,
                                  cont_vars=ukb_cont,
                                  title="Clinical variables — UKB transfer",
                                  fname_prefix="ukb", out_dir=ukb_plots)
        df_stats_ukb = characterize_biotypes(
            df_ukb.reset_index(), labels_ukb,
            cont_vars=ukb_cont, cat_vars=[]
        )
        print("\n  UKB clinical stats:")
        print(df_stats_ukb.to_string(index=False) if not df_stats_ukb.empty else "  (no clinical vars available)")
        df_stats_ukb.to_csv(ukb_output / "biotype_stats_ukb.csv", index=False)

        df_assign_ukb = df_ukb[[c for c in ["srt", "age", "sex"] if c in df_ukb.columns]].copy()
        df_assign_ukb["Biotype"]   = labels_ukb
        df_assign_ukb["certainty"] = probs_ukb.max(axis=1)
        for bt in range(best_k_prim):
            df_assign_ukb[f"prob_B{bt}"] = probs_ukb[:, bt]
        df_assign_ukb.to_csv(ukb_output / "biotype_assignments_ukb.csv")
        print(f"\n  UKB plots  → {ukb_plots}")
        print(f"  UKB tables → {ukb_output}")
    else:
        print("  UKB data not available — skipping.")

    # ─────────────────────────────────────────────────────────────────────────
    # FINAL SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 66)
    print("  FINAL SUMMARY")
    print("=" * 66)
    print(f"\n  DIAGNOSTIC (all sites, n={len(df_main)})")
    print(f"    k={best_k_diag}, silhouette={sil_diag:.3f}, "
          f"ARI median={np.median(aris_diag):.3f}")
    print(f"    → Site artifact (london_2) dominates minority biotype.")

    print(f"\n  PRIMARY (excl. london_2, n={len(df_primary)})")
    print(f"    k={best_k_prim}, silhouette={sil_prim:.3f}, "
          f"ARI median={np.median(aris_prim):.3f}  "
          f"IQR [{np.percentile(aris_prim,25):.3f}–{np.percentile(aris_prim,75):.3f}]")
    print(f"    Stable (ARI ≥ 0.60): {(aris_prim>=0.60).mean():.1%}")

    print(f"\n  Biotype interpretation (primary):")
    print(f"    Biotype 0 (n={int((labels_prim==0).sum())}, {(labels_prim==0).mean():.1%}): "
          f"older, more hearing loss, less tinnitus distress  → peripheral/audiological")
    print(f"    Biotype 1 (n={int((labels_prim==1).sum())}, {(labels_prim==1).mean():.1%}): "
          f"younger, less hearing loss, more tinnitus distress → central/distress-dominant")

    print(f"\n  Plots  → {PLOTS_DIR}")
    print(f"  Tables → {OUTPUT_DIR}")
    print("=" * 66 + "\n")
