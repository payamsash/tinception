from pathlib import Path
import pandas as pd
import math

tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
gwas_dir = tinception_dir / "GWAS"

GENOME_WIDE_LOG10P = -math.log10(5e-8)
SUGGESTIVE_LOG10P = -math.log10(1e-5)

results = {}
top_variants = []
suggestive_variants = []

for fname in sorted(gwas_dir.iterdir()):
    gwas_fname = fname.name

    if gwas_fname.startswith("."):
        continue
    if not gwas_fname.endswith(".txt"):
        continue
    if not gwas_fname.endswith("_assoc_merged.txt"):
        continue

    df = pd.read_csv(fname, sep="\t")

    # Basic cleanup
    if "LOG10P" not in df.columns:
        print(f"Skipping {gwas_fname}: no LOG10P column")
        continue

    df["LOG10P"] = pd.to_numeric(df["LOG10P"], errors="coerce")
    df = df.dropna(subset=["LOG10P"]).copy()

    if df.empty:
        print(f"Skipping {gwas_fname}: empty after cleaning")
        continue

    # Create a variant key
    required_cols = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Skipping {gwas_fname}: missing columns {missing}")
        continue

    df["variant_key"] = (
        df["CHROM"].astype(str) + ":" +
        df["GENPOS"].astype(str) + ":" +
        df["ALLELE0"].astype(str) + ":" +
        df["ALLELE1"].astype(str)
    )

    # Sort strongest first
    df = df.sort_values("LOG10P", ascending=False).reset_index(drop=True)

    top_row = df.iloc[0]
    n_gws = (df["LOG10P"] >= GENOME_WIDE_LOG10P).sum()
    n_sugg = (df["LOG10P"] >= SUGGESTIVE_LOG10P).sum()

    results[gwas_fname] = {
        "n_variants": len(df),
        "max_log10p": top_row["LOG10P"],
        "top_variant_key": top_row["variant_key"],
        "top_id": top_row["ID"],
        "top_chr": top_row["CHROM"],
        "top_pos": top_row["GENPOS"],
        "n_genome_wide": int(n_gws),
        "n_suggestive": int(n_sugg),
    }

    top_variants.append({
        "file": gwas_fname,
        "variant_key": top_row["variant_key"],
        "ID": top_row["ID"],
        "CHROM": top_row["CHROM"],
        "GENPOS": top_row["GENPOS"],
        "LOG10P": top_row["LOG10P"],
    })

    sugg_df = df.loc[df["LOG10P"] >= SUGGESTIVE_LOG10P,
                     ["variant_key", "ID", "CHROM", "GENPOS", "LOG10P"]].copy()
    sugg_df["file"] = gwas_fname
    suggestive_variants.append(sugg_df)

# -------------------------
# Summary per file
# -------------------------
summary_df = pd.DataFrame(results).T.sort_values("max_log10p", ascending=False)

print("\n=== Summary per file ===")
print(summary_df.to_string())

# -------------------------
# Top hit in each file
# -------------------------
top_df = pd.DataFrame(top_variants).sort_values("LOG10P", ascending=False)

print("\n=== Top hit per file ===")
print(top_df.to_string(index=False))

# -------------------------
# Shared top variants across files
# -------------------------
if not top_df.empty:
    shared_top_df = (
        top_df.groupby("variant_key")
        .agg(
            n_files=("file", "nunique"),
            files=("file", lambda x: ", ".join(sorted(x))),
            max_LOG10P=("LOG10P", "max"),
            ID=("ID", "first"),
            CHROM=("CHROM", "first"),
            GENPOS=("GENPOS", "first"),
        )
        .reset_index()
        .sort_values(["n_files", "max_LOG10P"], ascending=[False, False])
    )

    shared_top_df = shared_top_df[shared_top_df["n_files"] > 1]

    print("\n=== Variants that are top hits in more than one file ===")
    if shared_top_df.empty:
        print("None")
    else:
        print(shared_top_df.to_string(index=False))

# -------------------------
# Shared suggestive variants across files
# -------------------------
if suggestive_variants:
    sugg_all = pd.concat(suggestive_variants, ignore_index=True)

    shared_sugg_df = (
        sugg_all.groupby("variant_key")
        .agg(
            n_files=("file", "nunique"),
            files=("file", lambda x: ", ".join(sorted(set(x)))),
            max_LOG10P=("LOG10P", "max"),
            ID=("ID", "first"),
            CHROM=("CHROM", "first"),
            GENPOS=("GENPOS", "first"),
        )
        .reset_index()
        .sort_values(["n_files", "max_LOG10P"], ascending=[False, False])
    )

    shared_sugg_df = shared_sugg_df[shared_sugg_df["n_files"] > 1]

    print("\n=== Suggestive variants shared across files (-log10P >= 5) ===")
    if shared_sugg_df.empty:
        print("None")
    else:
        print(shared_sugg_df.to_string(index=False))

# -------------------------
# Save outputs
# -------------------------
# summary_df.to_csv(gwas_dir / "gwas_summary_per_file.tsv", sep="\t")
# top_df.to_csv(gwas_dir / "gwas_top_hit_per_file.tsv", sep="\t", index=False)

if 'shared_top_df' in locals():
    shared_top_df.to_csv(gwas_dir / "gwas_shared_top_hits.tsv", sep="\t", index=False)

if 'shared_sugg_df' in locals():
    shared_sugg_df.to_csv(gwas_dir / "gwas_shared_suggestive_hits.tsv", sep="\t", index=False)

print("\nSaved:")
print(gwas_dir / "gwas_summary_per_file.tsv")
print(gwas_dir / "gwas_top_hit_per_file.tsv")
print(gwas_dir / "gwas_shared_top_hits.tsv")
print(gwas_dir / "gwas_shared_suggestive_hits.tsv")