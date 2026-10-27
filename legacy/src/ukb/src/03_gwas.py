from pathlib import Path
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
gwas_dir = tinception_dir / "GWAS"
RAW_DIR = gwas_dir / "raw"
PROCESSED_DIR = gwas_dir / "processed"
PLOTS_DIR = gwas_dir / "plots"
HITS_DIR = gwas_dir / "significant_hits"
CLUMP_DIR = gwas_dir / "clumping"

os.makedirs(PROCESSED_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(HITS_DIR, exist_ok=True)
os.makedirs(CLUMP_DIR, exist_ok=True)

GENOME_WIDE_P = 5e-8
SUGGESTIVE_P = 1e-5

REQUIRED_COLS = [
    "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
    "A1FREQ", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P"
]

def compute_lambda_gc(pvals):
    pvals = np.asarray(pvals)
    pvals = pvals[np.isfinite(pvals)]
    pvals = pvals[(pvals > 0) & (pvals <= 1)]
    if len(pvals) == 0:
        return np.nan
    chisq = chi2.isf(pvals, df=1)
    return np.median(chisq) / 0.4549364

def qq_plot(pvals, title, outpath):
    pvals = np.asarray(pvals)
    pvals = pvals[np.isfinite(pvals)]
    pvals = pvals[(pvals > 0) & (pvals <= 1)]
    pvals = np.sort(pvals)

    n = len(pvals)
    if n == 0:
        return

    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(pvals)

    plt.figure(figsize=(6, 6))
    plt.scatter(exp, obs, s=8, alpha=0.6)
    maxv = max(exp.max(), obs.max())
    plt.plot([0, maxv], [0, maxv], linestyle="--", linewidth=1)
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def manhattan_plot(df, title, outpath):
    d = df.copy()
    d = d[np.isfinite(d["P"])]
    d = d[(d["P"] > 0) & (d["P"] <= 1)]
    d["minuslog10p"] = -np.log10(d["P"])

    d["CHR"] = pd.to_numeric(d["CHR"], errors="coerce")
    d["BP"] = pd.to_numeric(d["BP"], errors="coerce")
    d = d.dropna(subset=["CHR", "BP"])

    d = d.sort_values(["CHR", "BP"]).reset_index(drop=True)

    # cumulative position
    chroms = sorted(d["CHR"].unique())
    x_labels = []
    x_labels_pos = []
    offset = 0
    cumulative_positions = []

    for chrom in chroms:
        idx = d["CHR"] == chrom
        positions = d.loc[idx, "BP"].values
        cumulative = positions + offset
        cumulative_positions.extend(cumulative)

        x_labels.append(str(int(chrom)))
        x_labels_pos.append(cumulative.mean())

        offset = cumulative.max()

    d["cumulative_pos"] = cumulative_positions

    plt.figure(figsize=(14, 6))
    for i, chrom in enumerate(chroms):
        idx = d["CHR"] == chrom
        plt.scatter(
            d.loc[idx, "cumulative_pos"],
            d.loc[idx, "minuslog10p"],
            s=8,
            alpha=0.7,
            label=str(int(chrom)) if i < 2 else None
        )

    plt.axhline(-np.log10(5e-8), linestyle="--", linewidth=1)
    plt.axhline(-np.log10(1e-5), linestyle=":", linewidth=1)

    plt.xticks(x_labels_pos, x_labels, rotation=0)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def standardize_file(filepath):
    trait = os.path.basename(filepath).replace("_assoc_merged.txt", "")

    df = pd.read_csv(filepath, sep="\t", low_memory=False)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{filepath} missing columns: {missing}")

    df = df.copy()

    # Keep only additive model if multiple TEST categories exist
    if "TEST" in df.columns:
        additive_mask = df["TEST"].astype(str).str.upper().isin(["ADD", "ADDITIVE", "DOM", "GENO"])
        if additive_mask.any():
            # Prefer ADD if present
            add_mask = df["TEST"].astype(str).str.upper() == "ADD"
            if add_mask.any():
                df = df.loc[add_mask].copy()
            else:
                df = df.loc[additive_mask].copy()

    df["CHR"] = pd.to_numeric(df["CHROM"], errors="coerce")
    df["BP"] = pd.to_numeric(df["GENPOS"], errors="coerce")
    df["SNP"] = df["ID"].astype(str)
    df["A1"] = df["ALLELE1"].astype(str)
    df["A2"] = df["ALLELE0"].astype(str)
    df["FRQ"] = pd.to_numeric(df["A1FREQ"], errors="coerce")
    df["BETA"] = pd.to_numeric(df["BETA"], errors="coerce")
    df["SE"] = pd.to_numeric(df["SE"], errors="coerce")
    df["CHISQ"] = pd.to_numeric(df["CHISQ"], errors="coerce")
    df["LOG10P"] = pd.to_numeric(df["LOG10P"], errors="coerce")
    df["N"] = pd.to_numeric(df["N"], errors="coerce")

    df["P"] = np.power(10.0, -df["LOG10P"])
    df.loc[~np.isfinite(df["P"]), "P"] = np.nan
    df.loc[df["P"] == 0, "P"] = 1e-300

    df["Z"] = df["BETA"] / df["SE"]

    keep_cols = [
        "CHR", "BP", "SNP", "A1", "A2", "FRQ", "N",
        "BETA", "SE", "Z", "CHISQ", "LOG10P", "P", "TEST"
    ]
    out = df[keep_cols].copy()
    out = out.dropna(subset=["CHR", "BP", "SNP", "P"])
    out = out[np.isfinite(out["P"])]
    out = out[(out["P"] > 0) & (out["P"] <= 1)]

    out["CHR"] = out["CHR"].astype(int)
    out["BP"] = out["BP"].astype(int)

    # save processed summary stats
    processed_path = os.path.join(PROCESSED_DIR, f"{trait}.sumstats.tsv")
    out.to_csv(processed_path, sep="\t", index=False)

    # save clumping input
    clump_input = out[["SNP", "P"]].drop_duplicates()
    clump_input.to_csv(
        os.path.join(CLUMP_DIR, f"{trait}.clump_input.tsv"),
        sep="\t", index=False
    )

    # significant hits
    gws = out[out["P"] < GENOME_WIDE_P].sort_values("P")
    sugg = out[out["P"] < SUGGESTIVE_P].sort_values("P")

    gws.to_csv(os.path.join(HITS_DIR, f"{trait}.genomewide_hits.tsv"), sep="\t", index=False)
    sugg.to_csv(os.path.join(HITS_DIR, f"{trait}.suggestive_hits.tsv"), sep="\t", index=False)

    # plots
    manhattan_plot(out, f"{trait} Manhattan", os.path.join(PLOTS_DIR, f"{trait}.manhattan.png"))
    qq_plot(out["P"].values, f"{trait} QQ", os.path.join(PLOTS_DIR, f"{trait}.qq.png"))

    lambda_gc = compute_lambda_gc(out["P"].values)

    summary = {
        "trait": trait,
        "n_snps": len(out),
        "min_p": out["P"].min(),
        "n_gws": len(gws),
        "n_suggestive": len(sugg),
        "lambda_gc": lambda_gc
    }
    return summary

def main_1():
    files = sorted(glob.glob(os.path.join(RAW_DIR, "*.txt")))
    if not files:
        raise FileNotFoundError(f"No .txt files found in {RAW_DIR}")

    summaries = []
    for f in files:
        print(f"Processing {f}")
        summaries.append(standardize_file(f))

    summary_df = pd.DataFrame(summaries).sort_values("min_p")
    summary_df.to_csv(os.path.join(PROCESSED_DIR, "gwas_summary_overview.tsv"), sep="\t", index=False)
    print(summary_df)

if __name__ == "__main__":
    main_1()



###########
HITS_DIR = gwas_dir / "significant_hits"
OUT_DIR = gwas_dir / "overlap"
os.makedirs(OUT_DIR, exist_ok=True)

def load_hits(pattern):
    files = sorted(glob.glob(os.path.join(HITS_DIR, pattern)))
    data = {}
    for f in files:
        trait = os.path.basename(f).replace(pattern.replace("*", ""), "").replace(".tsv", "")
        trait = os.path.basename(f).replace(".genomewide_hits.tsv", "").replace(".suggestive_hits.tsv", "")
        df = pd.read_csv(f, sep="\t")
        if "SNP" in df.columns:
            data[trait] = set(df["SNP"].astype(str))
        else:
            data[trait] = set()
    return data

def pairwise_overlap(data, label):
    traits = sorted(data.keys())
    rows = []
    for i, t1 in enumerate(traits):
        for t2 in traits[i+1:]:
            s1 = data[t1]
            s2 = data[t2]
            inter = s1 & s2
            rows.append({
                "trait1": t1,
                "trait2": t2,
                "n_trait1": len(s1),
                "n_trait2": len(s2),
                "n_overlap": len(inter),
                "overlap_snps": ";".join(list(sorted(inter))[:50])
            })
    out = pd.DataFrame(rows).sort_values("n_overlap", ascending=False)
    out.to_csv(os.path.join(OUT_DIR, f"pairwise_overlap_{label}.tsv"), sep="\t", index=False)
    return out

def main_2():
    gws = load_hits("*.genomewide_hits.tsv")
    sugg = load_hits("*.suggestive_hits.tsv")

    gws_out = pairwise_overlap(gws, "genomewide")
    sugg_out = pairwise_overlap(sugg, "suggestive")

    print("Top genome-wide overlaps:")
    print(gws_out.head(20))
    print("\nTop suggestive overlaps:")
    print(sugg_out.head(20))

if __name__ == "__main__":
    main_2()

########
IN_DIR = gwas_dir / "processed"
OUT_DIR = gwas_dir / "processed/fuma_ready"
os.makedirs(OUT_DIR, exist_ok=True)

for f in glob.glob(os.path.join(IN_DIR, "*.sumstats.tsv")):
    trait = os.path.basename(f).replace(".sumstats.tsv", "")
    df = pd.read_csv(f, sep="\t")

    out = df[["SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE", "N"]].copy()
    out = out.dropna(subset=["SNP", "CHR", "BP", "P"])
    out.to_csv(os.path.join(OUT_DIR, f"{trait}.fuma.tsv"), sep="\t", index=False)

######## final stage
from pyliftover import LiftOver
import pandas as pd

lo = LiftOver('hg38', 'hg19')

fnames = ["PC2_brain", "tinnitus_subtype", "MDm-lh", "Lateral-nucleus-lh", "CeM-lh"]

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
gwas_dir = tinception_dir / "GWAS"

def normalize_chr(ch):
    ch = str(ch).strip()
    if ch.lower().startswith("chr"):
        ch = ch[3:]
    if ch == "23":
        ch = "X"
    return f"chr{ch}"

def convert(row):
    chrom = normalize_chr(row["CHR"])
    pos = int(row["BP"])
    res = lo.convert_coordinate(chrom, pos)
    if not res:
        return None
    return res[0][1]


for fname in fnames:
    df = pd.read_csv(gwas_dir / "processed" / f"{fname}.sumstats.tsv", sep="\t")
    for _, r in df.head(5).iterrows():
        print(r["CHR"], r["BP"], normalize_chr(r["CHR"]), convert(r))

    df["BP_hg19"] = df.apply(convert, axis=1)

    print(f"{fname} Before drop:", len(df))
    print("Mapped:", df["BP_hg19"].notna().sum())
    df2 = df.dropna(subset=["BP_hg19"]).copy()
    df2["BP"] = df2["BP_hg19"].astype(int)
    df2 = df2.drop(columns=["BP_hg19"])
    print(f"{fname} After drop:", len(df2))
    df2.to_csv(gwas_dir / "processed" / "fuma_ready" / f"{fname}.hg19.tsv", sep="\t", index=False, compression="zip")