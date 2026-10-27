# Tinception — Project Overview

## What This Project Is

**Tinception** is a neuroimaging study investigating structural brain differences in chronic tinnitus using a **multiscale, multi-dataset approach**. The core question: are tinnitus-related brain changes spatially patterned, biologically meaningful, and replicable?

Two datasets are used:
- **Main dataset** — a multi-site clinical cohort assembled from 13 sites (ANDES, ANTINOMICS, COGTAIL, ERLANGEN, ICCAC, INDITMS, LONDON_1, LONDON_2, NEUROPREN, NEW-CASTLE, TINMEG, TINSPECT, TRIPLE), with Tinnitus (TI) and Control (CO) subjects matched on sex, age, site, and PTA.
- **UK Biobank (UKB)** — a large population cohort (~40,000 subjects), used as a replication and normative reference dataset.

---

## Pipeline Overview (step by step)

### Step 1 — Data Preparation (`01_prepare_master_files.py`)
- Reads per-site Excel demographic files from `master_files/`
- Assigns unique `tinception_id` per subject (e.g., `sub-01001`)
- Drops subjects missing age, sex, or PTA
- Saves merged `master.csv`

### Step 2 — FreeSurfer Reconstruction (`02_reconstruct.sh`, `02_subsegment.sh`)
- Runs `recon-all` for cortical surface reconstruction
- Runs sub-segmentations: hippocampal subfields, amygdalar nuclei, thalamic nuclei

### Step 3 — Quality Control (`03_run_qc_1/2/3.py/.sh`)
- Visual and automated QC of FreeSurfer outputs

### Step 4 — Propensity Score Matching (`04_match.Rmd`)
- Uses R `MatchIt` (nearest-neighbor, 1:1, caliper=0.2)
- Matching variables: sex, age, site, PTA
- Outputs `matched_nearest.csv` and HTML summary tables in `master_files/matched/`

### Step 5 — VBM Pipeline (`07_run_vbm_1` through `07_run_vbm_10`)
- **1**: Copies NIfTI files from per-site directories to central structure, renamed by `tinception_id`
- **2**: FSL VBM preprocessing (`fslvbm_1_bet`, `fslvbm_2_template`, `fslvbm_3_proc`)
- **3**: Computes TIV from tissue probability maps; builds FSL design matrix + contrast files (TI>CO, CO>TI); covariates: age, sex, site dummies, PTA, TIV
- **4**: Runs FSL randomise with TFCE for whole-brain group comparison
- **5**: Normative modeling at **voxel level** (subcortical mask) using PCNToolkit BLR with B-spline basis (degree=3, nknots=5), heteroskedastic + SinhArcSinh warp; trained on controls only, applied to tinnitus subjects
- **6**: Reads VBM normative outputs; compares proportions of extreme deviations (|Z|>1.96) between CO and TI per voxel using chi-square + FDR correction; saves significant voxel mask as NIfTI
- **7**: Visualizes significant voxels over MNI152 template
- **8**: Evaluates normative model quality (Rho, MSLL, EXPV) with glass-brain plots
- **9**: Extracts cluster masks (putamen L/R, subcallosal) from VBM TFCE results for downstream correlations
- **10**: Normative modeling of VBM for UKB dataset (separate script)

### Step 6 — SBM Pipeline (`08_run_sbm_1` through `08_run_sbm_7`)
- **2**: Prepares FSGD file + contrast matrices for FreeSurfer's mris_preproc/mri_glmfit
- **3/4**: Runs mris_preproc and mri_glmfit for cortical thickness group comparison
- **5**: Normative modeling on surface (same BLR framework, vertex-wise)
- **7**: Evaluates SBM normative model quality on inflated surface (fsaverage5)

### Step 7 — Subcortical ROI Analysis (`09_run_subcortical_roi_2` through `_7`)
- **2**: Reads hippocampal subfields, amygdalar nuclei, thalamic nuclei, and cortical parcellation (aparc, aparc.a2009s, aseg) volumes from FreeSurfer stats files for all subjects
- Runs ANCOVA (pingouin) per region controlling for age, sex, site, PTA, TIV; computes Cohen's d and FDR-corrected p-values
- **5**: Normative modeling per atlas (Hippo, Amygdala, Thalamus) for both the main dataset and UKB; identifies subjects with extreme deviations per region
- **6**: Forest plot comparing Cohen's d effect sizes between main dataset and UKB, with a tier system (Tier 1 = significant after correction in both; Tier 2 = nominally significant in both; Tier 3 = one only)

### Step 8 — Biotyping (`10_compare_datasets.py`, `10_compare_biotypes_1/2.py`)
- Takes normative Z-scores of subcortical regions (Tier 2+3 regions from both datasets)
- Applies UMAP (cosine metric) + HDBSCAN clustering to define **tinnitus biotypes** based on the pattern of structural deviations
- Saves biotype assignments for each subject
- Compares biotypes using whole-brain VBM contrast (nilearn SecondLevelModel, FDR-corrected)
- Plots UMAP state-space colored by biotype and global deviation burden; spider/radar plots of biotype profiles

### Step 9 — Correlations with Clinical Variables (`12_run_correlations.py`)
- Partial correlations (controlling for sex, age, site, PTA) between THI (tinnitus handicap inventory) and brain region volumes in tinnitus subjects
- Same for PTA
- Also correlates THI/PTA with UMAP axes
- Plots: scatter + regression for putamen vs THI; FacetGrid for PTA vs UMAP dimensions

### Step 10 — Structural Similarity / Gradient Analysis (`11_compute_ssa_0` through `_4`)
- **1**: Loads individual structural covariance connectomes (Schaefer-1000 parcellation), computes group-averaged reference, runs diffusion map embedding (BrainSpace `GradientMaps`) with Procrustes alignment; saves per-subject gradient scores and eigenvalues
- **4**: Computes T-maps (TI vs CO) per gradient component; correlates T-maps with Allen Human Brain Atlas gene expression (via `abagen`) across 1000 parcels; runs spatial null models (Alexander-Bloch spin test via `neuromaps`) for top genes; saves significant gene associations

### Step 11 — UKB-Specific Pipeline (`src/ukb/`)
- **00**: Downloads matched UKB demographics + FreeSurfer derivative features via DNAnexus (`dx` CLI); neuroHarmonize batch correction (Combat) across UKB imaging sites trained on controls
- **01**: Mass ANCOVA (OLS) per atlas (aparc, aparc.a2009s, amygdalar nuclei, hippo subfields, thalamic nuclei, aseg) controlling for age, sex, SRT; Bonferroni correction; partial correlation with tinnitus severity and duration
- Normative modeling (`norm_modeling.py`, `sbm_norm_modeling.py`, `vbm_norm.py`): same BLR framework, trained on UKB controls; covariates: age, sex, SRT (speech reception threshold), SITE as batch effect
- **GWAS** (`02_gwas.sh`, `03_gwas.py`): GWAS results for structural measures, enrichment analysis with `gseapy` (GO, KEGG, Reactome); candidate genes include CAMTA1, WNT4

---

## Key Results / Findings (from scripts and output files)

### VBM (Whole-Brain)
- Main classical analysis: CO > TI cluster in **putamen (bilateral)** and **subcallosal cortex** (TFCE-corrected)
- Normative analysis (voxel-wise chi-square of extreme deviations): significant voxels in subcortical regions after FDR correction

### Subcortical ROI — Main Dataset vs UKB (tier system)
Regions with convergent evidence across both cohorts (from `10_compare_datasets.py`):

| Region | Direction | Tier |
|---|---|---|
| GC-ML-DG body (rh) | reduced in TI | significant in both |
| Presubiculum body/head (rh/lh) | reduced | convergent |
| Molecular layer HP body (rh) | reduced | convergent |
| CA3 head (lh) | reduced | convergent |
| Fimbria (lh) | reduced | convergent |
| Basal nucleus (lh), Accessory Basal, Lateral nucleus (lh) | reduced | amygdala convergence |
| Medial nucleus amygdala (lh/rh) | altered | amygdala convergence |
| CeM (lh/rh) | altered | thalamic convergence |
| MDm (lh/rh) | altered | thalamic convergence |
| MVRe (lh/rh) | altered | thalamic convergence |
| Accumbens area (lh/rh), Inf Lat Vent | altered | aseg convergence |

Thalamic nuclei with Bonferroni-significant results in UKB: CeM, MDm, MVRe (bilateral).

### Normative Models (UKB)
- Centile plots exist for hippocampal subfields (presubiculum, GC-ML-DG), amygdalar nuclei (Basal), and thalamic nuclei (CeM, MDm, MVRe) — saved as PDFs in `src/ukb/data/plots/`

### Biotypes
- Two stable biotypes identified from UMAP+HDBSCAN on normative Z-score profiles across converging regions
- Biotypes differ in whole-brain VBM maps (separate for TI and CO groups)
- Biotypes differ in THI score and correlate with PTA

### Gene Expression (SSA/Gradients)
- Top genes correlated with gradient T-maps (TI>CO) identified via spin-test-validated Spearman correlations against Allen Human Brain Atlas
- Results saved in `ssa_results/genes.csv`

---

## File/Directory Layout

```
Tinception/
├── master_files/             # Per-site demographics (xlsx) + merged master.csv + matched CSVs
├── src/
│   ├── main/                 # Main pipeline scripts (01_–12_)
│   │   ├── 01_prepare_master_files.py
│   │   ├── 04_match.Rmd
│   │   ├── 07_run_vbm_{1-10}.py/.sh
│   │   ├── 08_run_sbm_{1-7}.py/.sh
│   │   ├── 09_run_subcortical_roi_{2-7}.py
│   │   ├── 10_compare_datasets.py       # Forest plot, tier assignment
│   │   ├── 10_compare_biotypes_1/2.py   # Biotype VBM contrast
│   │   ├── 11_compute_ssa_{0-4}.py      # Gradient + gene analysis
│   │   └── 12_run_correlations.py       # THI/PTA vs brain regions
│   └── ukb/
│       ├── src/
│       │   ├── 00_get_and_match_data.ipynb
│       │   ├── 01_harmonize_and_stats.ipynb
│       │   ├── norm_modeling.py          # UKB normative modeling (subcortical)
│       │   ├── sbm_norm_modeling.py      # UKB SBM normative modeling
│       │   ├── vbm_norm.py               # UKB VBM normative modeling
│       │   ├── 02_gwas.sh / 03_gwas.py   # GWAS + enrichment
│       │   └── run_sbm.py / run_vbm.py
│       └── data/
│           ├── norm_models/hippo_subfields/  # Fitted normative models
│           ├── plots/                        # Centile + extreme count PDFs
│           └── results/                      # Bonferroni stat tables per atlas
```

---

## Key Technical Choices

| Choice | Detail |
|---|---|
| Normative modeling | PCNToolkit BLR, B-spline basis (deg=3, 5 knots), heteroskedastic, SinhArcSinh warp |
| Batch correction (UKB) | neuroHarmonize (Combat), trained on controls |
| Matching | MatchIt nearest-neighbor 1:1, caliper=0.2, on sex/age/site/PTA |
| Multiple comparison | FDR-BH for main dataset; Bonferroni for UKB; chi-square + FDR for voxel-wise normative |
| Clustering | UMAP (cosine, n_neighbors=10) + HDBSCAN (min_cluster_size=20) |
| Gradient analysis | Diffusion map, Schaefer-1000, Procrustes alignment to control reference |
| Gene spatial null | Alexander-Bloch spin test (neuromaps, 1000 permutations) |
| VBM stats | FSL randomise + TFCE (non-parametric permutation) |
| SBM stats | FreeSurfer mri_glmfit |

---

## Open Questions / Things to Discuss

1. **Which results to highlight?** The tier system in `10_compare_datasets.py` gives a principled ranking — are Tier 1 regions (significant after correction in both datasets) the headline finding?
2. **Biotyping stability**: are the UMAP/HDBSCAN biotypes stable across random seeds and parameter choices?
3. **Gene results**: spin-test validated gene correlations from SSA — are any of these publishable as a "molecular virtual biopsy"?
4. **GWAS integration**: the GWAS script is minimal (`pval < 1e-6` threshold, two candidate genes: CAMTA1, WNT4) — how far has this analysis progressed?
5. **VBM normative vs classical**: both approaches were run — how do they compare? Do the normative voxels overlap with the TFCE clusters?
