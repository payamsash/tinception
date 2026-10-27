# Tinception: Neuroanatomical and Genetic Characterization of Chronic Tinnitus

*A comprehensive project summary — methodology, results, discussion points, and open issues.*
*Compiled from: repository source code (`src/main/`, `src/ukb/`), result files on `/Volumes/Extreme_SSD/payam_data/Tinception`, and `master_files/wtinception_v2.pptx` (67-slide master deck, as of the latest revision).*

---

## 1. Executive Summary

**Tinception** asks whether chronic subjective tinnitus is associated with reproducible structural brain differences, using a **discovery-replication design**: a **multi-site clinical cohort** (13 recruiting sites, ~285 tinnitus vs. ~285 matched controls after QC/matching) as the discovery sample, and the **UK Biobank (UKB)** (2,180 tinnitus vs. 2,176 controls) as a large, independently-ascertained population replication sample. The pipeline spans whole-brain voxel- and surface-based morphometry, fine-grained subcortical nuclei segmentation (hippocampal subfields, amygdalar nuclei, thalamic nuclei), subject-level normative modeling, unsupervised biotyping of tinnitus patients, structural-covariance gradient/transcriptomic analysis, mode-based (eigenmode/spectral) morphometry, and a full genome-wide association pipeline with polygenic risk scoring.

**Headline pattern:** the most consistent finding across both datasets and multiple analysis strategies is a **small, bilateral increase in putamen grey matter volume in tinnitus** (not a decrease), alongside a family of small-effect increases in specific **amygdalar, hippocampal, and thalamic nuclei**, and a possible **subcallosal cortex** decrease that is inconsistently significant. Effect sizes throughout are small (Cohen's d mostly 0.1–0.3), which is itself an important finding: this is a subtle, distributed structural signature rather than a large focal lesion. A two-cluster **biotype structure** (~20–25% "peripheral/audiological" vs. ~75–80% "central/distress-dominant") replicates closely between the clinical cohort and UKB. Genetics (GWAS) is well-powered and clean (no inflation) but has **not yet produced a genome-wide-significant hit**; suggestive loci have been annotated via FUMA with biologically plausible candidate genes (HPA-axis, calcium signaling, synaptic vesicle machinery).

**Where the project currently stands:** most core analyses are complete and have been run on both datasets; the newest and most rigorous layers (GMM-based biotyping, Bayesian hierarchical meta-analysis of cross-dataset effects) have only partially been folded back into the main narrative and slide deck. There are several internal inconsistencies (four parallel "which regions replicate" scoring systems; a confirmed bug in one results file; an untracked script for the PRS analysis) that should be resolved before this becomes a manuscript.

---

## 2. Background & Rationale

Chronic subjective tinnitus (phantom auditory perception without an external sound source) affects a substantial fraction of the population and its neurobiology remains contested. Prior single-site, small-N structural MRI studies have produced conflicting results, plausibly because (a) sample sizes were too small to detect the small effect sizes typical of psychiatric/audiological structural phenotypes, (b) tinnitus is neurobiologically heterogeneous (auditory-driven vs. limbic/distress-driven presentations), and (c) hearing loss (a near-universal comorbidity) is a strong confound that is not always adequately controlled.

Tinception addresses all three: **pooling 13 independent sites** for power, **explicitly modeling and matching on hearing threshold (PTA)** as a covariate throughout, and **explicitly testing for neuroanatomical subtypes** (biotypes) rather than assuming a single homogeneous "tinnitus brain." The UK Biobank arm adds a population-representative replication sample roughly 8x larger than the clinical cohort, at the cost of coarser tinnitus phenotyping (self-report severity/duration rather than clinical instruments) and a healthier, older, non-clinically-selected population.

---

## 3. Study Design & Datasets

### 3.1 Main multi-site clinical cohort

13 contributing sites/consortium partners, each with local ethics approval and MRI acquisition, harmonized centrally:

| Site | PI / Partner | Location | Sequence | Field | Loudness measure | Controls (raw) | Tinnitus (raw) |
|---|---|---|---|---|---|---|---|
| Andes | Paul Delano | Chile | MPRAGE | 3T | — | 48 | 58 |
| Antinomics | Patrick Neff | USZ | TFE | 3T | THI | 38 | 41 |
| Cogtail | Martin Meyer | USZ | TFE | 3T | THI | 79 | 64 |
| Erlangen | Nadia Müller | Erlangen | MPRAGE | 3T | — | 0 | 22 |
| ICCAC | Tobias Kleinjung | PUK | TFE | 3T | THI | 0 | 24 |
| IndiTMS | Stefan Schoisswohl | Regensburg | MPRAGE | 1.5T | THI | 0 | 13 |
| London_1 | William Sedley | London | MPM | 3T | THI | 13 | 16 |
| London_2 | William Sedley | London | MPRAGE | 3T | — | 14 | 15 |
| Neuropren | Patrick Neff | Regensburg | MPRAGE | 1.5T | THI | 12 | 28 |
| New Castle | William Sedley | New Castle | TFE | 3T | THI | 14 | 14 |
| Tinmeg | Christopher Cederroth | Karolinska Huddinge | MPRAGE | 3T | THI | 21 | 18 |
| Tinspect | Patrick Neff | USZ | TFE | 3T | THI | 22 | 20 |
| Triple | Stefan Schoisswohl | Regensburg | TFE | 3T | THI | 25 | 43 |

Covariates recorded uniformly: age, sex, PTA (pure-tone average, hearing threshold), TIV (total intracranial volume). Raw pooled N before QC/matching ≈ 376 controls / 286 tinnitus (per the deck's metadata slide).

**After QC and 1:1 propensity-score matching** (nearest-neighbor, caliper=0.2, on sex/age/site/PTA), the *analysis-ready* matched cohort used for the subcortical ROI, VBM design matrix, and MBM pipelines is **N=569** (284 controls / 285 tinnitus) across **10 sites** — Erlangen, ICCAC, and IndiTMS (all tinnitus-only recruitment, no local controls) are present in the raw VBM-per-site data but **drop out of the final matched analysis set**, since 1:1 matching requires control availability at each site. *(This is a real discrepancy from the "13 sites" framing used elsewhere in project documentation — worth stating explicitly as "13 recruiting sites, 10 contribute to the matched analysis" going forward.)*

### 3.2 UK Biobank (UKB) validation cohort

Drawn from >500,000 UKB participants with imaging. Matched cohort: **2,180 self-reported tinnitus vs. 2,176/2,180 controls** (age 56.7 ± 7.5 vs. 56.8 ± 7.1, ~40% female in both groups; no significant age/sex difference, p>0.5–0.9).

Matching/adjustment variables: sex, age, MRI site, frequency of loud-music exposure, noisy-workplace history, hearing-aid use, hearing difficulty in background noise, alcohol intake frequency, handedness, average income, and **speech-reception threshold (SRT)** left and right (UKB's proxy for PTA, used throughout as the hearing-loss covariate in place of clinical PTA).

UKB structural derivatives were harmonized across imaging sites/scanner batches with **neuroHarmonize (ComBat)**, trained on controls only, before any group comparison or normative modeling.

---

## 4. Methodological Pipeline

### 4.1 Data preparation & matching
Per-site Excel demographics → merged `master.csv` with a unique `tinception_id` per subject → subjects missing age/sex/PTA dropped → **R `MatchIt`**, nearest-neighbor 1:1, caliper 0.2, matched on sex/age/site/PTA.

### 4.2 Structural MRI processing
FreeSurfer `recon-all` (cortical surface reconstruction) plus dedicated sub-segmentation streams for **hippocampal subfields**, **amygdalar nuclei**, and **thalamic nuclei** (FreeSurfer's Bayesian nucleus/subfield segmentation tools). Visual + automated QC gates the pipeline before any group-level analysis (`QC/` directory: 664 subject-level BIDS anatomical folders on the SSD; `qc_dir/` holds ~3,620 MRIQC-style run logs — no single consolidated pass/fail table currently exists, which is worth generating for a Methods section).

### 4.3 Voxel-Based Morphometry (VBM)
FSL-VBM stream: brain extraction → study-specific template building → GM modulation/registration → **`randomise` with TFCE** (threshold-free cluster enhancement, non-parametric permutation), design matrix covariates = age, sex, site (dummy-coded), PTA, TIV; contrasts TI>CO and CO>TI. Run at two significance conventions throughout the deck: p<0.005 (primary) and p<0.05 (exploratory/uncorrected-adjacent). A secondary "97-subject subsample" analysis adds **high-frequency PTA (PTA_HF)** as an extra covariate to test whether findings survive stricter hearing-loss control. On UKB, the equivalent whole-brain contrast is implemented independently in **nilearn** (`SecondLevelModel` + `non_parametric_inference`, TFCE, 5,000 permutations, one-sided directional tests), covariates age/sex/SRT/eTIV + site dummies — this UKB classical group-level VBM stream exists in the codebase (`vbm_group.py`) but is not documented in the original project overview.

### 4.4 Surface-Based Morphometry (SBM)
FreeSurfer `mris_preproc` + `mri_glmfit` on fsaverage-resampled cortical thickness (and area/volume) maps, two design variants: a simple two-sample model and a **DOSS** (different-offset, same-slope) model with site as a fixed effect. Cluster-forming thresholds tested at both th=2.0/2.3/3.0 (≈p<0.01/0.005/0.001), cluster-wise significance via `mri_glmfit-sim` Monte Carlo simulation (`--cache`, `--cwp 0.05`, `--2spaces` correcting across both hemispheres/contrasts). Subject-level normative modeling is run on the same vertex-wise thickness maps (see §4.6). The UKB analog (`run_sbm_2.ipynb`) mirrors this with the same `mri_glmfit`/`mri_glmfit-sim` approach, fed by a dedicated FreeSurfer `-qcache` re-run (`run_qcache.py`) on UKB's bulk FreeSurfer derivatives.

### 4.5 Subcortical nuclei segmentation & ROI analysis
Region volumes extracted from FreeSurfer stats files across six atlases: hippocampal subfields, amygdalar nuclei, thalamic nuclei, `aparc` (Desikan-Killiany, 68 labels), `aparc.a2009s` (Destrieux, 148 labels), `aseg` (39 subcortical/ventricular labels). Statistical test: **ANCOVA (pingouin)** per region, covariates age/sex/site/PTA/TIV, effect size Cohen's d, multiple-comparison correction **FDR-BH** for the main dataset and **Bonferroni** for UKB (UKB's larger sample and single-cohort nature make the stricter correction tractable).

### 4.6 Normative modeling framework
Used at voxel level (VBM), vertex level (SBM), and ROI level (subcortical nuclei), consistently implemented with **PCNToolkit Bayesian Linear Regression (BLR)**: B-spline basis (degree=3, 5 knots), heteroskedastic noise model, **SinhArcSinh warp** (to handle non-Gaussian residuals). Models are trained on **controls only** and applied to predict expected volume/thickness/intensity for tinnitus subjects, yielding a per-subject, per-feature deviation Z-score. Extreme deviations (|Z|>1.96) are then compared between groups via **chi-square + FDR**, testing whether tinnitus subjects are *disproportionately* represented in the tails of the normative distribution — a complementary lens to the mean-shift ANCOVA/TFCE tests, sensitive to heterogeneous/bidirectional individual effects that a group-mean test would miss.

Two additional refinements beyond the original design:
- **Model-extension / transfer learning for the main dataset** (`run_nm_subcortical_ukb_main.py`): rather than fitting an independent BLR model from scratch on the ~284 main-dataset controls, the *UKB-trained* normative model (built on ~2,176 controls) is extended with the main-dataset controls via PCNToolkit's `extend_predict`, borrowing statistical strength from the much larger UKB reference population. **This is a meaningful modeling choice that should be stated explicitly in any Methods section**, since it changes what "deviation from norm" means for the clinical cohort (deviation from a UKB-anchored population norm, not a small site-specific one).
- **Quality control of the normative models themselves** (`10_evaluate_norm_models.py` + UKB counterparts): Rho (predicted-observed correlation), MSLL (mean standardized log-loss), and EXPV (explained variance) are computed per voxel/vertex and visualized as glass-brain / inflated-surface maps. Per the deck's Conclusion 2, **main-dataset-trained models substantially outperform UKB-trained models on these metrics**, attributed to the clinical cohort's protocol homogeneity vs. UKB's greater scanner/site heterogeneity — meaning individual-level deviation estimates from the main dataset are considered more trustworthy.

### 4.7 Biotyping / subtype discovery
Three distinct clustering strategies were tried over the project's course (an important detail for a Methods section, since only one should be reported as primary):

1. **UMAP (cosine, n_neighbors=10) + HDBSCAN (min_cluster_size=20)** on normative Z-scores of convergent subcortical regions — the original approach (`10_compare_biotypes_1/2.py`, `10_compare_datasets.py`). Used to generate the whole-brain VBM biotype-contrast maps still referenced in the deck.
2. **PLS-residual "Resilient vs. Vulnerable" split** (`test_discovery_replication.ipynb`, exploratory only): PLS regression of brain Z-scores against clinical/demographic variables, with the *residual* of brain-score-on-clinical-score used to define a "True Sensitivity Index," clustered or tertile-split. Not used downstream; useful mainly as the design precursor to (3) below (it first identified that PTA explains part of the cluster separation).
3. **Gaussian Mixture Model (GMM), current/most rigorous approach** (`13_biotype_gmm.py`): fit directly on normative Z-scores (11 convergent regions) after optionally partialing out PTA and site, standardizing, and PCA-reducing to 85% variance. Model order selected by **BIC** across k=1–6 with a minimum-cluster-size guard (5%); stability assessed via **500-iteration bootstrap** with Adjusted Rand Index (median ARI=0.747, 66.8% of resamples ≥0.60) and a Ward hierarchical dendrogram as an independent check. Critically, this script **discovered and corrected a site artifact**: the naive full-sample clustering (N=285) was 86.7% driven by a single site (London_2); the reported "primary" result excludes that site (N=270). UKB replication is done by **posterior-probability transfer** (predicting UKB subjects' cluster membership from the main-dataset-fitted GMM), not independent re-clustering — a stronger replication test than re-clustering UKB separately, since it asks "do UKB subjects fall into the *same* structure," not merely "does UKB also have two clusters."

**GMM biotyping is the analysis that should be treated as canonical going forward.** It is directly reflected in the current slide deck (the 21.5%/78.5% main and 24.9%/75.1% UKB proportions on Slides 41–42 match this script's output exactly), superseding the UMAP/HDBSCAN-derived whole-brain contrast maps still shown in Slides 45–48, which were not checked for the same site artifact.

### 4.8 Structural similarity network & gradient / transcriptomic analysis (SSA)
Individual **structural covariance connectomes** built from Schaefer-1000-parcel (7-network) cortical thickness/morphometric similarity, subject-level. Group-averaged reference connectome → **diffusion map embedding** (BrainSpace `GradientMaps`, `approach="dm"`, 10 components) with **Procrustes alignment** to the reference, yielding per-subject gradient scores. Group (TI vs. CO) T-maps computed per gradient component, then spatially correlated (Spearman) against **Allen Human Brain Atlas** regional gene expression (via `abagen`, same Schaefer-1000 parcellation), with statistical significance assessed via the **Alexander-Bloch spin test** (spatial autocorrelation-preserving null model, `neuromaps`, 1,000 permutations) — this null-model correction is essential because naive spatial correlations between two smooth cortical maps are almost always "significant" by chance.

Individual-subject **normative models were also fit on the gradient scores themselves** (parcel-wise BLR per gradient component; `test_ssa_grads.ipynb`), testing whether tinnitus subjects show excess extreme deviations in structural-covariance gradient position, independent of the group T-map/gene-expression analysis.

Network-level (Yeo-7 system) block-averaged connectome comparisons (TI vs. CO, and Biotype1 vs. Biotype0) via max-statistic permutation testing exist as an exploratory convergent-validity check (`test_ssa_conn.ipynb`) but are not yet in a production script.

### 4.9 Mode-based morphometry (MBM)
A previously undocumented analysis axis using the published **Mode-Based Morphometry** method (Cao, Pang, Segal, Chen, Aquino, Breakspear, Fornito — *Human Brain Mapping*, 2024; toolbox at `~/MBM`, MATLAB). Cortical thickness maps are decomposed onto the **Laplace-Beltrami eigenmodes** of cortical surface geometry — a spatial-frequency ("shape harmonic") basis analogous to a Fourier decomposition adapted to brain anatomy — allowing group differences to be characterized by *which spatial scales* (coarse global patterns vs. fine local patterns) carry the effect, in addition to *where*.

Pipeline: FreeSurfer `fsaverage` white-matter surfaces converted to VTK; ANCOVA design matrix built (age-centered, sex, PTA-centered, TIV-centered, N=569, matching the main matched cohort); `mbm_main()` run per hemisphere with **ANCOVA_F test, 5,000 permutations**, tail-refined via Generalized Pareto Distribution approximation below p=0.1 (PALM-style), **200 eigenmodes**, FDR-corrected beta spectrum, threshold 0.05. **Results exist on disk (`.mat` files + rendered PNGs per hemisphere) confirming the pipeline completed, but there is no downstream script interpreting the beta-spectrum output into the results narrative yet** — this is a completed-but-unreported analysis.

### 4.10 Genetics: GWAS, FUMA annotation, PRS
**GWAS (UKB, DNAnexus Research Analysis Platform).** Genotype-array QC (plink2: MAF≥0.01, MAC≥20, genotyping/sample missingness≤0.1, HWE p>1e-15) followed by **REGENIE two-step whole-genome regression**: Step 1 (null-model LOCO predictions) on genotyped SNPs with age, genotype batch, SRT, eTIV, PC1–10 as covariates; Step 2 association testing on **whole-exome-sequencing (WES) OQFE** variant calls, run separately for:
- Two **brain principal components** (PC1_brain, PC2_brain — the major axes of covariance across the 11 convergent subcortical regions), quantitative trait, linear model.
- **Tinnitus biotype/subtype** as a binary phenotype (Firth logistic regression with SE approximation, for the imbalanced ~25%/75% class split).
- Each of the **11 individual convergent subcortical regions** separately, as quantitative traits.

Post-processing (`03_gwas.py`): standardizes REGENIE output, computes genomic control λ, produces Manhattan/QQ plots, genome-wide (p<5e-8) and suggestive (p<1e-5) hit tables, **cross-phenotype SNP-overlap tables**, FUMA-ready export files, and hg38→hg19 liftover for a subset of traits.

**FUMA annotation** ("Functional Mapping and Annotation of GWAS," SNP→gene mapping using GTEx eQTL and chromatin-interaction/PsychENCODE data) has been run on the suggestive loci for several traits (see §5.6 for the actual gene list and biological narrative), producing gene-level and pathway-level (MAGMA) hypotheses even though no locus reached genome-wide significance.

**PRS (polygenic risk scores).** Five conditions tested against three subcortical atlases (amygdalar nuclei, hippocampal subfields, thalamic nuclei): **Alzheimer's disease (AD), bipolar disorder (BD), schizophrenia (SCZ), hypertension (HT), cardiovascular disease (CVD)**. **No script generating these results currently exists anywhere in the repository or on the results SSD** — only the 15 output CSVs are present (columns: `feature, r, pval, n, pval_adj, significant`, N=1,232). The exact PRS source (UKB's own released PRS fields vs. a custom LDpred2/PRSice2 computation) cannot be confirmed from the outputs alone and should be tracked down (a different machine, an unsaved notebook, or a colleague's script) before this is written up, since a reviewer will ask.

### 4.11 Clinical correlations
Partial correlations (pingouin `partial_corr`, covariates sex/age/site/PTA, Bonferroni-corrected) between **THI** (Tinnitus Handicap Inventory) and brain region volumes/putamen/subcallosal cluster masks in tinnitus subjects only; same for **PTA**. Also tested: THI/PTA vs. UMAP/GMM biotype axis position, and THI/PTA differences *between* biotypes (Mann-Whitney/Kruskal-Wallis). UKB equivalent uses self-reported tinnitus **severity and duration** in place of THI/PTA (UKB has no clinical THI instrument).

### 4.12 Statistical framework — multiple comparison and "which regions replicate" scoring systems
Because this point recurs throughout Results (§5.3) and is a genuine point of internal inconsistency, it's worth stating plainly here: **four different schemes exist in the codebase for deciding which subcortical regions "count" as cross-dataset replicated**, developed in sequence as the project matured:

| Scheme | Script | Logic | Regions flagged as top tier |
|---|---|---|---|
| 1. Nominal-p tier system | `10_compare_datasets.py` | Tier 1 = FDR-sig in UKB **and** nominal-p<0.05 in Main, same direction; Tier 2 = nominal-p<0.05 in both | Tier 1/2 combined: ~44 regions (see §5.3) |
| 2. OLS covariate-adjusted tiers | `test_discovery_replication.ipynb` (exploratory) | Same idea, but region effect from full OLS model instead of ANCOVA/Cohen's-d table | Not finalized/scratch quality |
| 3. Frequentist meta-analysis tiers | `10_compare_datasets_bayesian.py` | Inverse-variance-weighted pooled Cohen's d across datasets, Wald z-test, FDR-BH on pooled p | "Tier 2: Consistent (Main Underpowered)" — 44/97 regions |
| 4. Bayesian hierarchical tiers | `10_compare_datasets_bayesian.py` (PyMC) | Per-region random-effects model (`theta` = shared true effect, `tau` = between-dataset heterogeneity), tiered by posterior P(&#124;theta&#124;>0.05)≥0.95/0.80 and P(same sign)≥0.90 | **Only 3/97 regions** reach "probable shared effect" (see §5.3) |

**Scheme 4 (Bayesian) is the most statistically principled** — it properly propagates uncertainty from both datasets' standard errors rather than treating a nominal p<0.05 in an underpowered 284-vs-285 sample as equivalent evidence to p<0.05 in a 2,180-vs-2,176 sample — and should probably be adopted as the primary/headline result, with the nominal-tier table retained as a supplementary "broader consistency" view.

---

## 5. Results

### 5.1 VBM (whole-brain voxel-based morphometry)

**Main dataset (TFCE-corrected, p<0.005):**
- **TI > CO:** bilateral **putamen**, right uvula (cerebellar vermis) — the putamen finding is the project's most consistent structural result (see §6).
- **CO > TI:** nothing survives correction at p<0.005; **subcallosal cortex** shows a trend (Cohen's d≈0.3, T=2.5) that becomes p<0.05-significant only in the 97-subject PTA_HF-adjusted subsample.
- At the looser p<0.05 threshold only (unreplicated): supracalcarine cortex.
- In the PTA_HF-adjusted 97-subject subsample: putamen persists, the uvula finding **disappears** (suggesting it may be driven by hearing-loss severity rather than tinnitus per se), subcallosal becomes p<0.05-significant.

**UKB dataset (TFCE-corrected):**
- **TI > CO:** cerebellar lobules VIIa/VIIb and VI (p<0.005); putamen re-appears but only at p<0.05.
- **CO > TI:** nothing significant (max effect size 0.016).

Four cluster masks were extracted from the TFCE results for downstream use: **putamen-lh, putamen-rh, subcallosal, and uvula-rh** (this fourth mask is not mentioned in earlier project documentation).

*Note on evidence quality:* the classical TFCE cluster maps are stored only as NIfTI probability maps on the SSD — no FSL `cluster` text report (peak MNI coordinates, exact cluster extents, corrected p-values as a table) currently exists. Generating one would make the VBM section reviewer-ready.

**Individual-level / normative VBM:** no voxels survive FDR correction for extreme-deviation excess in tinnitus in either direction (positive-only or negative-only) in either dataset — consistent with the group-level putamen effect being a small, population-wide, unidirectional shift too subtle to produce excess *individual* extreme deviations, rather than a bimodal "some patients affected, some not" pattern. Main-dataset-trained normative models substantially outperform UKB-trained models on calibration metrics (Rho/MSLL/EXPV), so main-dataset individual deviation maps are considered more reliable.

### 5.2 SBM (cortical thickness)

**Main dataset**, two-sample design (cluster-forming th=2.0, p<0.005-equivalent):
- Cortical thickness, **TI > CO only**: lh — lateraloccipital (1020.8 mm², CWP=0.0002), fusiform (798.7 mm², CWP=0.003); rh — superiorparietal (1525.9 mm², CWP=0.0002), inferiorparietal (752.5 mm², CWP=0.0062), middletemporal (557.7 mm², CWP=0.036).
- No clusters survive at the stricter th=3.0 threshold; no CO>TI clusters anywhere; no significant clusters for cortical area or volume.
- **DOSS design (site as fixed effect, th=2.3, the deck's headline figure)**: only **one cluster survives** — right hemisphere thickness, TI>CO, in a region spanning inferior parietal cortex / superior temporal sulcus (396.13 mm³, T=2.403, CWP=0.047, Cohen's d=0.202).

This single surviving cluster is discussed at length in the deck's Conclusion 3: **d=0.202 is a small effect and CWP=0.047 sits just inside the conventional threshold**; power analysis suggests ~400/group would be needed for 80% power at this effect size, meaning the main cohort (285/285) is marginally underpowered. It **does not replicate in UKB**, plausibly due to phenotype ascertainment differences (self-report vs. clinical diagnosis) or acquisition-protocol effects on sulcal-region thickness estimation. The deck explicitly labels this a **hypothesis-generating, not established, finding**.

**Individual-level / normative SBM:** no vertices survive correction in either dataset — reinforcing the "small, spatially-restricted, non-bimodal effect" interpretation, or alternatively suggesting the single cortical-thickness cluster is a sampling artifact of the primary dataset.

### 5.3 Subcortical nuclei ROI analysis (ANCOVA, main dataset)

**Important finding on effect direction:** across both the main dataset and UKB, virtually every significant/near-significant subcortical effect is in the **TI > CO (larger volume in tinnitus)** direction — e.g., Central-nucleus-rh (amygdala): CO=45.33 mm³ vs. TI=47.74 mm³ (d=+0.21); UKB Basal-nucleus-lh (amygdala): d=+0.114. *(This corrects an earlier framing in project notes that described these effects as "reduced in tinnitus" — the primary source data show the opposite direction throughout.)*

**Main dataset, raw ANCOVA:** **zero regions reach FDR-significance** across all six atlases combined (367 regions tested). Best nominal (uncorrected) hits: Right-Putamen (aseg, d=+0.216, p=0.018); Right-VA thalamic nucleus reaches **FDR-significance only via the normative extreme-deviation chi-square test** (p_fdr=0.028), not the mean-shift ANCOVA — consistent with the "bidirectional individual deviations cancel out at the group mean" interpretation discussed in the deck's Conclusion 6.

**UKB, Bonferroni-corrected (much larger sample, so smaller effects reach significance):** **22 of 135 tested regions** significant, **all positive-direction except one**:

| Atlas | Sig./Tested | Notable regions (d, p_adj) |
|---|---|---|
| Thalamic nuclei | 8/50 | MVRe-rh (+0.120, 3.6e-3), MVRe-lh (+0.117, 5.6e-3), MDm-lh (+0.115, 7.8e-3), MDm-rh (+0.107, 2.2e-2), CeM-rh (+0.106, 2.4e-2), CeM-lh (+0.105, 2.8e-2), MDl-lh (+0.103, 3.3e-2), PuA-lh (+0.101, 4.3e-2) |
| Hippocampal subfields | 5/41 | GC-ML-DG-body-rh (+0.114, 7.3e-3), presubiculum-body-rh (+0.109, 1.3e-2), whole-hippocampal-body-rh/lh (+0.105/+0.102), presubiculum-head-rh (+0.101, 3.4e-2) |
| aseg | 6/27 | Accumbens-lh (+0.117, 3.0e-3), Accumbens-rh (+0.108, 1.0e-2), **Inf-Lat-Vent-lh (−0.107, 1.1e-2 — the only negative-direction hit)**, Thalamus-Proper-rh/lh (+0.103/+0.102), Amygdala-lh (+0.099, 3.0e-2) |
| Amygdalar nuclei | 3/17 | Basal-nucleus-lh (+0.114, 3.1e-3), Lateral-nucleus-lh (+0.102, 1.4e-2), Accessory-Basal-nucleus-lh (+0.097, 2.3e-2) |

**Cross-dataset synthesis (Bayesian hierarchical model, 97 regions tested, `bayesian_dmin_0.05.csv`):**
- Frequentist meta-analytic tiering: 44/97 regions "Tier 2: Consistent (Main Underpowered)" (i.e., pooled-effect direction agrees and is nominally supported, but the main-dataset arm alone is underpowered); 52/97 "Tier 5: Non-Significant"; 1/97 "Tier 4: Dataset-Specific."
- **Bayesian posterior tiering (stricter, accounts for between-dataset heterogeneity): only 3 regions reach "Tier 2: Probable Shared Effect"** — **Accessory-Basal-nucleus-lh** (d_main=0.153, d_ukb=0.097), **Medial-nucleus-rh** (d_main=0.199, d_ukb=0.082), **MDm-rh** (d_main=0.123, d_ukb=0.107). The remaining 94/97 regions fall into "Tier 4: Inconsistent/No Evidence" under this stricter criterion.

This is the single most important nuance for how the team frames the "converging regions" story: **the breadth of the finding (dozens of nominally consistent regions) is much larger than its statistical depth (3 regions with rigorously-quantified high-confidence shared effect)**. Both framings are legitimate but serve different purposes — the broad list is useful for hypothesis generation (e.g., feeding GWAS phenotypes, biotyping features), the narrow list is what should be defended as a "finding" in a paper's abstract.

**Individual-level thalamic nuclei (normative deviation, not group ANCOVA):** two dataset-specific, non-overlapping findings — **Right-VA (ventral anterior)** significant in the main dataset only (p=0.0283), **Left-PuL (pulvinar)** significant in UKB only (p=0.0273). Interpreted (Conclusion 6) as reflecting genuine dissociations: VA is the thalamic relay of the basal-ganglia/putaminal output circuit (mechanistically linking to the putamen VBM finding) and may specifically track clinically-recruited distress/chronicity not well captured in UKB's self-report; PuL is implicated in **thalamocortical dysrhythmia**, a leading tinnitus pathophysiology framework, and its detection may require UKB's much larger sample because only a minority subgroup (~15–20%, by the deck's own power-reasoning) shows the deviation.

### 5.4 Biotypes (GMM, current/canonical analysis)

**Main dataset (primary run, N=270, London_2 excluded as a site artifact):**
- **Biotype 0**: n=58 (21.5%) — older (mean 68 vs. 59.5 years, p=0.0002), **greater hearing loss** (PTA 27.7 vs. 18.1 dB, p=0.0001), **lower distress** (THI 14 vs. 28, p=0.023), sex-imbalanced (p=0.021). Interpreted as a **"peripheral/audiological" subtype**.
- **Biotype 1**: n=212 (78.5%) — younger, less hearing loss, higher THI distress. Interpreted as a **"central/distress-dominant" subtype**.
- Cluster quality: silhouette=0.121 (modest but expected for real biological subtypes with overlapping boundaries, not artificial well-separated clusters), bootstrap ARI median=0.747 (66.8% of 500 resamples ≥0.60 — reasonably stable).
- The diagnostic run including all sites (N=285) selected the same k=2 by BIC, but 86.7% of the minority cluster was attributable to a single site (London_2) — a confound the script explicitly detected and excluded from the primary reported result.

**UKB replication (posterior-probability transfer from the main-dataset-fitted model, N=4,356 tinnitus subjects):** Biotype 0 n=1,084 (24.9%), Biotype 1 n=3,272 (75.1%) — **closely replicates the main dataset's 21.5%/78.5% split**, with 84.3% of UKB subjects classified with high certainty (posterior probability ≥0.70).

**Whole-brain contrast between biotypes (from the earlier UMAP/HDBSCAN-based pipeline, not yet re-run under the GMM labels):** Biotype 0 (labeled in the deck as showing "limbic-salience network dysregulation") differs from controls in whole-brain VBM contrasts (localized around anterior cingulate gyrus); Biotype 1 ("low-deviation subtype") does not differ from controls. *Because the GMM analysis and the whole-brain-contrast analysis currently use different clustering methods, this specific whole-brain-contrast result should be re-verified under the GMM (site-artifact-corrected) labels before it's presented as final — the clinical characterization (peripheral/audiological vs. central/distress) and the neuroanatomical characterization (limbic-salience dysregulation vs. low-deviation) are complementary framings of what is very likely the same underlying two-cluster structure, given the near-identical group proportions, but this hasn't been formally confirmed.*

### 5.5 Clinical correlations (THI, PTA)

**Main dataset:** **PTA correlates with brain region volumes** (putamen, subcallosal, and/or the 11 convergent regions — partial correlation controlling for sex/age/site); **THI does not** survive correction (same result under FDR-BH and Bonferroni). This is a biologically sensible pattern (hearing-loss severity tracks structure) but means the *distress* dimension of tinnitus (THI) is currently structurally silent in this cohort at the whole-region level — despite THI differing significantly *between* biotypes (§5.4), which is a distinct, cluster-level (not linear-correlation) relationship.

**UKB dataset:** none of the three questions tested (subtype vs. severity/duration; region volume vs. severity/duration; UMAP axis vs. severity/duration) reached significance — likely reflecting the coarser self-report severity/duration phenotype relative to the clinical THI instrument.

### 5.6 Genetics: GWAS, FUMA, PRS

**GWAS summary (15 phenotypes, ~180K SNPs each, λ_GC 1.00–1.08 — well-controlled, no inflation):** **zero genome-wide-significant (p<5e-8) hits for any phenotype.** Suggestive-threshold (p<1e-5) hits only:

| Phenotype | Suggestive hits | Lead SNP (p) |
|---|---|---|
| CeM-lh | 5 | chr16:28843229 (β=+0.316, p=2.4e-6) |
| MDl-lh | 4 | chr13:21146817 (p=4.1e-6) |
| tinnitus_subtype (biotype GWAS) | 2 | chr16:3118828 (β=−1.086, p=1.9e-6; n=1,497) |
| tinnitus_vs_controls | 1 | chr22:27799344 (β=−1.071, p=3.1e-6; n=4,030) |
| PuA-lh, MDm-rh, Accessory-Basal-nucleus-lh | 0 | — (weakest signals) |

**Zero SNP overlap** between any pair of phenotypes at either genome-wide or suggestive thresholds — the individual-region, subtype, and brain-PC GWAS signals are genetically independent of one another so far.

**FUMA gene mapping (suggestive loci only, biologically interpretive, not confirmatory):**
- **Tinnitus-subtype loci:** chr9 locus (140,267,710–140,277,794) → **EXD3** exclusively; associated SNPs are hippocampal eQTLs (GTEx v8, FDR≈0.02) with subtype-risk alleles → **reduced hippocampal EXD3 expression** — notable given the parahippocampal deviations seen in normative modeling. chr16 locus → **ZNF205-AS1** + adjacent non-coding RNA (regulatory, not protein-coding). Nominal (non-Bonferroni-surviving) MAGMA pathway hits: calcium-ion-transport regulation, Parkinson's-disease gene sets (echoing the basal-ganglia/putamen structural theme).
- **PC2 (brain PCA) loci:** chr12 locus (31 SNPs) → **HRK** (pro-apoptotic BH3-only protein; implicated in cochlear hair-cell and neuronal survival under noise/neurotrophic stress) — MAGMA gene-based test reaches genome-wide significance here, the **single strongest genetic signal in the project**; nominal pathway enrichment: voltage-gated calcium channel regulation, sensory pain perception (echoing tinnitus/chronic-pain central-sensitization overlap hypotheses). Also eQTL hits for **TESC** (frontal cortex) and **RNFT2** (cerebellum). chr17 locus → **CAMTA2/MINK1** (both near-completely loss-of-function-intolerant, pLI≈1.0).
- **CeM-lh (left central medial thalamic nucleus) loci — richest architecture, 4 suggestive loci:** chr17q21.31 (a well-known inversion-polymorphism, pleiotropic region) → **CRHR1** (corticotropin-releasing-hormone receptor 1, primary HPA-axis stress mediator), hippocampal eQTL, risk alleles → **reduced hippocampal CRHR1 expression** — directly supports a stress-axis mechanistic link; also **KANSL1**, **NSF** (synaptic vesicle recycling) at the same locus, previously linked to neuroticism/worry/cortical morphology GWAS. Plus chr8, chr9 (**ZNF462**), chr16 (**ATP2A1**, extending the calcium-signaling theme) loci.
- **Lateral-nucleus-lh:** one chr4 locus → **DCTD** (pyrimidine metabolism/DNA repair), no brain eQTLs, no prior neuro/psych GWAS catalog hits — flagged in the deck itself as a lower-priority candidate, reported for completeness.
- **Gene2Func tissue-expression profiling:** EXD3 low overall brain expression but strong hippocampal eQTL regulation; **ZSCAN10** near-absent in adult brain (possible neurodevelopmental-only role); ZNF205/DPP7 broader regulatory/proteostatic expression. No pathway enrichment survives correction with only 4 input genes (expected, given the small gene-set size).

**Overall genetic narrative:** no confirmed hits, but a coherent, biologically-plausible *hypothesis space* has emerged — recurring themes across independent loci are **HPA-axis/stress signaling (CRHR1)**, **calcium signaling/channel regulation** (multiple independent loci), and **synaptic vesicle machinery (NSF)** — worth stating explicitly as the project's genetic hypothesis-generating narrative, with the caveat that none of it is genome-wide-significant and all gene assignments come from suggestive-threshold loci.

**PRS:** only **1 of 15 tested associations significant after correction**: **Bipolar Disorder PRS ↔ CA1-body (right) hippocampal subfield volume** (r=+0.093, p=0.0011, p_adj=0.047, n=1,232). Alzheimer's, schizophrenia, hypertension, and cardiovascular-disease PRS showed no corrected-significant associations with any tested subcortical structure (closest near-miss: AD-PRS vs. thalamic Pc-lh, p_adj=0.065). *(Source script for this analysis could not be located — see §7.)*

### 5.7 Gene expression / structural-covariance gradients (SSA)

Whole-map (group T-map level) spatial correlations between the three structural-covariance gradient T-maps (TI vs. CO) and Allen Human Brain Atlas gene expression were modest and not map-level-significant (Gradient 1: Rho=−0.022, p=0.57; Gradient 2: Rho=0.008, p=0.84; Gradient 3: Rho=−0.071, p=0.070, nominally closest to significance). However, at the **individual-gene level**, spin-test-corrected (spatial-null, 1,000 permutations) analysis flagged **38 of 45 examined genes as significant** (P_Spatial<0.05) despite the weak whole-map correlation — this is statistically coherent (individual genes can show strong, spatially-specific, null-surviving correlation patterns even when the aggregate map-level correlation is weak) but means the "top genes" list should be interpreted as **spatially-localized associations**, not evidence of a genome-wide-coordinated expression shift.

Strongest signals by gradient: **Gradient 3** (spin-test spatial rho up to 0.16–0.17, p_spatial as low as 0.001) — MRRF, LRP8, CDADC1, SLC25A37, INPP5A, PCNT, VAV3, and others; **Gradient 1** — ZNF451, EYA3, PI4K2B, ABCA6, PECR, and others (spatial rho 0.11–0.15); **Gradient 2** — TMC8, MYH3, ENGASE, F10, C5, and others (spatial rho ~0.09–0.12).

**⚠ Known bug, must be fixed before quoting raw statistics from this analysis:** in `11_compute_ssa_4.py`, the raw (non-spin-test) Spearman Rho/p columns in the saved `ssa_results/genes.csv` are computed inside a loop but overwritten by a stale variable reuse bug — every gene within a given gradient shares the identical raw Rho/p value (that of whichever gene happens to be last-alphabetical in the expression matrix for that gradient), rather than its own value. **The spin-test-corrected columns (`Spatial_Rho`, `P_Spatial`) are computed correctly and are unaffected** — so the significance calls above are trustworthy, but the raw correlation magnitudes reported anywhere from this file need to be regenerated from a corrected script before publication.

**Individual-level normative modeling on gradient scores** (a complementary, subject-level test independent of the group T-map/gene analysis above) found **37, 19, and 125 nominally-different parcels for gradients 1, 2, 3 respectively, but zero parcels survived FDR correction for any gradient.** This null result should be reported alongside the gene-expression findings for calibration — it indicates that, at the level of *individual-subject deviation from norm* in gradient position, there is no robust group-level signal, even though the *group-average* T-map shows spatially-organized, gene-correlated structure. These are different questions (mean shift vs. excess individual deviation) and are not contradictory, but both need to be in the write-up.

### 5.8 Mode-Based Morphometry (MBM)

Eigenmode/spectral decomposition completed on the full matched cohort (N=569) for both hemispheres (ANCOVA_F, 5,000 permutations, 200 eigenmodes, FDR-corrected). Output files exist (`.mat` results, rendered PNG figures of the top-10 most influential modes per hemisphere) confirming successful completion, but **no script currently interprets or reports the resulting beta-spectrum** (which spatial frequencies carry the TI-vs-CO effect) into the project narrative. This is a "ready to analyze" result, not yet a "reported" one.

---

## 6. Synthesis: The Converging Narrative

Pulling the above together, the project's strongest, most defensible storyline is:

1. **A small, consistent, bilateral putamen grey-matter increase in tinnitus** — present in the main-dataset TFCE analysis (p<0.005), preserved under hearing-loss-severity adjustment (97-subject PTA_HF subsample), and re-emerging in UKB despite its very different phenotyping (p<0.05). This is the project's headline finding.
2. **A basal-ganglia–thalamic circuit theme**, not just an isolated putamen finding: the putamen result is mechanistically linked to the **VA thalamic nucleus** individual-deviation finding (VA is the primary thalamic relay of the putaminal/pallidal output circuit) and to the Bayesian-confirmed **MDm** (mediodorsal thalamus) and **Medial/Accessory-Basal amygdala nucleus** effects — together sketching a distributed basal-ganglia/limbic-thalamic circuit rather than a single focal region.
3. **A possible, less certain subcallosal/anterior-cingulate-adjacent decrease**, only significant under stricter hearing-loss adjustment and absent in UKB — plausibly reflecting a clinically-selected-sample-specific (distress, chronicity, psychiatric comorbidity) effect not captured by UKB's population-based, self-report ascertainment.
4. **A genuine two-cluster biotype structure** (~20–25% peripheral/audiological, ~75–80% central/distress-dominant) that replicates with unusually good quantitative agreement between the clinical cohort and UKB (21.5/78.5 vs. 24.9/75.1) — arguably the project's second-strongest finding, since replication of *cluster proportions* under an independent transfer-prediction procedure is a fairly stringent test.
5. **A coherent but unconfirmed genetic hypothesis space** centered on HPA-axis/stress signaling (CRHR1), calcium-channel regulation (recurring across 3+ independent suggestive loci), and synaptic vesicle machinery (NSF) — currently only suggestive-threshold evidence, but thematically consistent enough to be worth stating as a hypothesis for future, larger GWAS.
6. **Small effect sizes throughout (mostly d=0.1–0.3) are themselves a finding**: they argue against a single dramatic tinnitus "lesion" and for either (a) a distributed, low-magnitude structural signature consistent with a functional/network disorder more than a focal structural one, or (b) meaningful biological heterogeneity that only the biotyping analysis (item 4) is currently positioned to resolve.

---

## 7. Discussion Points, Inconsistencies, and Open Questions

These are things worth actively discussing as a team before finalizing a paper or a public-facing presentation — not polish issues, but substantive decisions that change what gets claimed.

1. **Which "replication" scoring system is authoritative?** Four exist (§4.12); they disagree by an order of magnitude on how many regions "replicate" (44 vs. 3). Recommend adopting the Bayesian hierarchical tiering as the primary claim, with the broader nominal-tier list explicitly labeled as a secondary/exploratory/hypothesis-generating list (e.g., for GWAS phenotype selection, which is how it's already being used).
2. **Effect-direction correction.** Confirm and lock in the "larger volume in tinnitus" direction across all downstream materials (slides, any prior manuscript drafts, internal notes) — the raw data throughout (both PPTX table and SSD result files) show increases, not decreases, and any older material describing decreases should be corrected.
3. **Reconcile the two biotype characterizations.** "Peripheral/audiological vs. central/distress-dominant" (clinical characterization, from the GMM script) and "limbic-salience dysregulation vs. low-deviation" (neuroanatomical characterization, from the older UMAP/HDBSCAN whole-brain contrast) are very likely the same two clusters given near-identical proportions, but this hasn't been formally re-verified under the corrected (site-artifact-excluded) GMM labels. Re-running the whole-brain VBM biotype contrast under GMM labels would let both framings be presented together with confidence.
4. **The GMM script's site-artifact finding raises a broader QC question**: were other group-level or individual-level analyses (not just biotyping) checked for site-driven confounds beyond the covariate adjustment already in the ANCOVA/TFCE design matrices? The GMM script is the only one in the codebase that performed this kind of explicit sensitivity check (leave-one-site-out).
5. **Normative-model provenance for the main dataset.** Multiple, inconsistently named output directories exist across scripts (`norm_models`, `norm_models_main`, `new_norm_models`, `norm_models_ukb_main`) for what should be a single canonical set of main-dataset normative Z-scores, and at least one script (`run_nm_subcortical_ukb_main.py`) uses UKB-model-extension rather than from-scratch fitting. Before writing a Methods section, it would be worth auditing which directory each downstream figure/analysis (GMM biotyping, tier comparison, THI/PTA correlation) actually reads, to confirm consistency.
6. **The `ssa_results/genes.csv` bug** (§5.7) needs a rerun of `11_compute_ssa_4.py` (a one-line fix — compute rho/p from the current `gene` in the top-15 loop rather than reusing the last loop-iteration value) before any raw correlation magnitude from that file is quoted in a talk or paper. The significance calls (spin-test) are unaffected.
7. **PRS analysis provenance.** No generating script exists in the repo or on the results SSD for `PRS/prs_*.csv` — worth tracking down before presenting these results, both to confirm methodology (which PRS release/algorithm) and because an untraceable-but-real analysis is a reproducibility risk.
8. **MBM (mode-based morphometry) and biotype-GWAS/brain-PC-GWAS are "dark data"**: both have completed, unreported results. MBM's beta-spectrum (which spatial scales carry the TI-vs-CO cortical thickness signal) could add a genuinely new dimension to the SBM story; the biotype/brain-PC GWAS results have no FUMA annotation at all yet (only the 11-region GWAS does), despite the biotype GWAS being arguably the more novel and mechanistically interesting phenotype (a genome-wide test of the neuroanatomical subtype itself).
9. **Sample-size framing.** "13 sites" is used loosely in prior documentation; be precise about **13 recruiting sites, 10 contributing to the final 1:1-matched analysis cohort (N=569)** — Erlangen, ICCAC, and IndiTMS (tinnitus-only recruitment) don't have local controls to match against.
10. **THI vs. biotype relationship.** THI correlates with biotype membership (categorically) but not with brain volumes (continuously) in the main dataset — worth discussing explicitly as evidence that the *structural biotype* captures something the linear THI-volume correlation misses, strengthening the case for biotyping as the more sensitive clinical-structural link.
11. **Slide 66 of the master deck contains what appear to be real, plaintext lab-equipment credentials** (EEG/experiment computer logins). This is unrelated to the neuroimaging findings but should be flagged: consider removing that slide before the deck is shared more widely, and rotating those credentials if they're still active.

---

## 8. Limitations

- **Small effect sizes** (mostly d=0.1–0.3) throughout mean most individual findings are only detectable with UKB-scale power; the main clinical cohort (285 vs. 285) is underpowered for many of the region-level effects it nominally flags, which is precisely why the Bayesian meta-analysis (§4.12, §5.3) is important for calibrating confidence.
- **Phenotyping asymmetry between cohorts**: THI (validated clinical instrument) in the main dataset vs. self-reported severity/duration in UKB is a substantive measurement-quality gap, likely explaining several main-dataset-only clinical correlations that fail to replicate in UKB.
- **No genome-wide-significant GWAS hit** for any phenotype — all genetic results are suggestive-threshold and should be presented as hypothesis-generating, not confirmatory, regardless of how mechanistically appealing the FUMA gene narrative is.
- **Multiple-testing burden across the whole project** is substantial (dozens of atlases × regions × datasets × normative/classical approaches); while each individual analysis applies appropriate within-analysis correction (FDR/Bonferroni/TFCE), there is no project-wide correction across the *number of analyses run*, which is a standard and accepted limitation of multi-modal discovery pipelines but should be stated explicitly.
- **Site heterogeneity** in the main dataset (13 sites, multiple field strengths 1.5T/3T, multiple sequences MPRAGE/TFE/MPM) is handled via covariate adjustment and ComBat-style harmonization (UKB side only — confirm whether the main dataset also uses harmonization or relies solely on covariate adjustment, since this wasn't confirmed in the current investigation) but remains a residual concern, as illustrated concretely by the GMM biotyping site artifact.
- **Cross-sectional design**: no ability to determine whether structural differences precede tinnitus onset (risk factor / vulnerability marker) or result from it (consequence of chronic altered auditory/limbic processing) — this is worth stating explicitly as a design limitation in any discussion section, especially given the CRHR1/stress-axis and putamen/basal-ganglia findings could plausibly run in either causal direction.

---

## 9. What's Missing / Recommended Next Steps

1. Fix the `genes.csv` bug and regenerate raw correlation values (§5.7, §7.6).
2. Locate or re-derive the PRS analysis source code (§4.10, §7.7).
3. Generate a text/tabular FSL `cluster` report for the VBM TFCE results (peak MNI coordinates, exact extents) — currently only NIfTI probability maps exist, which is a barrier to writing a precise Results section or figure legend.
4. Re-run the whole-brain VBM biotype contrast under the corrected GMM (site-artifact-excluded) cluster labels, to unify the clinical and neuroanatomical biotype characterizations (§7.3).
5. Interpret and report the MBM beta-spectrum results (§5.8, §7.8) — this is a completed analysis sitting unused.
6. Run FUMA annotation on the biotype-GWAS and brain-PC-GWAS suggestive loci (currently only the 11-region GWAS has gene annotation) (§7.8).
7. Decide on and document a single canonical main-dataset normative-model directory/pipeline (§7.5) before further downstream analyses are built on it.
8. Produce a consolidated QC pass/fail summary table (currently QC logs exist but no single summary was found) (§4.2).
9. Consider the two open questions the deck's own Appendix (Slide 65) poses but doesn't yet answer: **does hearing asymmetry explain any of the structural findings, and does it interact with tinnitus status?** — not yet investigated anywhere in the codebase.
10. Decide which of the four tier/replication systems (§4.12) is the headline framing for external presentation, and update the slide deck's tier heatmap/forest plot accordingly if the Bayesian version is adopted.
11. Remove or redact Slide 66 (lab credentials) before any further sharing of the master deck, and consider whether those credentials need to be rotated.

---

## 10. Appendix: Data & Code Map (for reference while building slides)

| What | Where |
|---|---|
| Pipeline source code | `/Users/payamsadeghishabestari/Tinception/src/main/` (main-dataset pipeline, scripts 01–13) and `/Users/payamsadeghishabestari/Tinception/src/ukb/` (UKB-specific pipeline + GWAS) |
| Master demographics | `/Users/payamsadeghishabestari/Tinception/master_files/master.csv`, matched cohort under `master_files/matched/` |
| Master slide deck | `/Users/payamsadeghishabestari/Tinception/master_files/wtinception_v2.pptx` (67 slides, current) |
| All quantitative results (stats CSVs, GWAS, PRS, MBM, biotypes) | `/Volumes/Extreme_SSD/payam_data/Tinception/` — subfolders: `subcortical_roi/`, `GWAS/`, `PRS/`, `MBM/`, `biotypes_gmm/`, `VBM/`, `vbm_per_site/`, `SBM/`, `ssa_results/`, `QC/`, `plots/` |
| Biotyping full write-up (already exists, very complete) | `/Volumes/Extreme_SSD/payam_data/Tinception/biotypes_gmm/BIOTYPING_RESULTS.md` |
| GWAS summary table | `/Volumes/Extreme_SSD/payam_data/Tinception/GWAS/processed/gwas_summary_overview.tsv` |
| Gene expression / gradient results (has the known bug, §5.7) | `/Volumes/Extreme_SSD/payam_data/Tinception/ssa_results/genes.csv` |
| Bayesian cross-dataset tier table | `/Volumes/Extreme_SSD/payam_data/Tinception/subcortical_roi/bayesian/bayesian_dmin_0.05.csv` |

---

*Document compiled 2026-08-13 by synthesizing the repository source code, the results archive on the external SSD, and the master presentation deck. Cross-checked for internal consistency across all three sources; discrepancies are flagged explicitly rather than silently resolved, since several (effect direction, tier-system choice, biotype naming) are decisions the team should make deliberately rather than have made for them.*
