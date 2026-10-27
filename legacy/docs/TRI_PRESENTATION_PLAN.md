# TRI Conference Talk — Slide-by-Slide Build Plan

**Project:** Tinception — neuroanatomy of chronic tinnitus across a multi-site clinical cohort and UK Biobank
**Format:** ~22 min talk + ~5 min Q&A · **26 slides**
**Framing:** rigor-forward — the honest evidence base, with heterogeneity as the real story
**Companion document:** `TINCEPTION_PROJECT_SUMMARY.md` (the analytical reference; this file is the presentational one)

---

## ⚠️ READ THIS FIRST — verify before you present

Eleven things were re-derived directly from the raw result files while building this plan, and they change what can honestly be claimed on stage. **Items 1–3 and 10–11 affect slides you would otherwise present as positive findings.**

| # | Finding | Consequence for the talk |
|---|---|---|
| **1** | **UKB biotype "replication" is computed on the wrong sample.** In `src/main/13_biotype_gmm.py`, `load_main_features()` filters `.query('group == "TI"')` (line 194); `load_ukb_features()` (line 214) has **no group filter** and returns all 4,356 UKB subjects (2,180 TI + 2,176 CO) while printing "TI subjects". Also, `biotype_stats_ukb.csv` is a **1-byte empty file** in both `biotypes_gmm/` and `biotypes_gmm/ukb_transfer/`. | The 24.9%/75.1% split is of the **whole UKB sample**, not UKB tinnitus cases — so it is **not** a replication of the main cohort's 21.5%/78.5%. Split properly, **UKB controls come out 23.8%/76.2% and UKB cases 26.0%/74.0% — the model does not discriminate.** Compounding it, `ukb_covars.csv` exists at neither path the script looks for, so UKB Z-scores were **never residualised** while the clinical ones were. Both need fixing. **See slide 25.** |
| **2** | **Zero of 361 regions survive FDR in the main-cohort ROI ANCOVA** (min p_fdr = 0.156, in `subcortical_roi/stats/*.csv`; every `Significant` field is `False`). | Never source a "significant subcortical difference" to the main cohort. It must come from UKB, the Bayesian meta, or the normative analysis. |
| **3** | **Every biotype whole-brain VBM contrast is null** — `vbm_results/Bio_TFCE/`: 0 voxels at p<0.05 for b0>b1, b1>b0 and the sub80 variants (max −log10 p = 1.109). **And those maps are dated 2026-04-08, a month before the GMM ran (2026-05-06); no GMM biotype design file exists and no script in the repo writes them** — so they are the superseded UMAP/HDBSCAN labels. | The "limbic–salience dysregulation subtype" (v2 slides 45–48) is not supported: **do not present it.** But state the null accurately — the whole-brain contrast under the *current* mixture-model labels **has never been run**. |
| **4** | **MBM is effectively null.** Its "FDR" threshold map is an uncorrected p≤0.05 mask; 5.56% (lh) / 6.04% (rh) of vertices survive — the chance rate. | Move to "what's next" (slide 27). Do not show as a result. |
| **5** | **Bayesian tier-2 regions are weaker than previously written.** P(\|θ\|>0.05) is 0.85–0.87 — the ≥0.80 tier, not ≥0.95 — and **all three 95% HDIs cross zero**. | State this explicitly on slide 15. It is a strength of the talk, not a weakness. |
| **6** | **`cluster_uvula-rh_mask.nii.gz` is empty (0 voxels)**; `cluster_subcollasal_mask` has mean d = +0.002 (null). | The uvula and subcallosal narratives from the old deck have no mask-level support. Drop both. |
| **7** | **10 sites, not 13** in every analysis file. N drifts: **569** (ROI/SBM/MBM) vs **572** (VBM) vs **560** (SSA/gradients). | Say "13 recruiting sites, 10 with local controls contributing to the matched analysis". Pick one N per slide and footnote the others. |
| **8** | **GWAS used ~180K SNPs** — a pruned panel, not full imputed density. | One honest clause on slide 26. |
| **9** | **PRS: 1 of 555 tests survives** (BD × CA1-body-rh, p_adj = 0.047). | Indistinguishable from chance. Mention only if asked. |
| **10** | **The biotype–THI difference is probably a site effect.** THI mean by site ranges 17.2 (cogtail) → 58.9 (triple) and was never collected at andes or london_2. **cogtail supplies 25 of Biotype 0's 58 patients (43%) and has the lowest THI of any site.** The test is an unadjusted Mann-Whitney on 43 vs 124 patients, with Biotype 0 also 8.5 years older. | **Do not lead with "the subtype with less hearing loss has more distress".** Redo it with site in the model first. Also fix the merge bug that drops antinomics THI (n = 38) from `subcortical_roi/master.csv` — it takes THI from 167 to 205. **See slide 22.** |
| **11** | **The biotype GWAS and the THI/putamen correlations use the *superseded* UMAP/HDBSCAN labels**, not the mixture-model biotypes. The GWAS phenotype comes from `biotypes/ukb.csv` (n = 1,512, the "n = 1,497" in the summary); `12_run_correlations.py` reads `VBM_design/gm.csv` with the 132/153 split. Nothing in the repo reads `biotypes_gmm/*`, and `13_biotype_gmm.py` is untracked in git. | Say "a data-driven subtype GWAS", not "a GWAS of the biotypes I just showed you". Same for `plots/correlations/putamen_vs_THI.pdf`. **See slide 26.** |

**Two integrity items before the deck leaves your laptop:**
- **`wtinception_v2.pptx` slide 66 contains plaintext lab-equipment credentials** (EEG/experiment machine logins). Delete slides 66–67 from any conference copy, and rotate those credentials if they are still live.
- `ssa_results/genes.csv` has a stale-variable bug in the raw Rho/p columns (every gene within a gradient shares one value). The spin-test columns are fine. Don't quote raw rho from that file.

---

## Title

> # The Tinnitus Brain Is Not One Brain
> ### Structural heterogeneity across a multi-site clinical cohort and UK Biobank

The memorable clause carries the title on its own; the cohort detail belongs in the subtitle, where it reads as credentials rather than as more title. Six words on the screen, and it states your actual thesis rather than a region name.

**Shorter variants, if you want the cohorts in the main line:**

| Variant | Note |
|---|---|
| **The Tinnitus Brain Is Not One Brain** *(subtitle carries the rest)* | ★ recommended — cleanest, most quotable, fits any slide template |
| The Tinnitus Brain Is Not One Brain: Structural Heterogeneity Across Two Cohorts | Keeps the claim in the title; still short enough to read at a glance |
| Not One Tinnitus Brain: Structural Heterogeneity in 4,900 Scans | Tightest version; the number does the credibility work |

**Subtitle line for the title slide:** *Multi-site clinical cohort (N = 569) · UK Biobank (N = 4,356)*

**A note on delivering it:** the title makes a strong claim, and slide 23 is where you qualify it — say early (slide 4) that heterogeneity is the *hypothesis* the talk tests, not a result you're asserting up front. That framing makes the title an argument rather than an overclaim. And note the title still works under the axis reading: not one brain does not require exactly two.

---

## The thesis, in one sentence

> Chronic tinnitus is associated with small, reproducible, consistently *positive* volumetric differences concentrated in basal-ganglia–thalamic–limbic circuitry — effects far too subtle to explain the disorder on their own, and best understood as the average over a patient population that varies continuously in how far its brains sit from normative expectation, along an axis organised by hearing rather than by symptom severity.

## Three things the audience should leave with

1. **Effects are real but small** (d ≈ 0.10–0.35) — and they point *up*, not down. The atrophy framing that dominates the literature is not what two independent cohorts show.
2. **Statistical rigor changes the answer.** Forty-four regions look "consistent across datasets" by nominal criteria; three survive a Bayesian hierarchical model, and even those have HDIs crossing zero.
3. **Heterogeneity is the story — as an axis, not a taxonomy.** Patients differ continuously in how far their subcortical structure sits from normative expectation, and that variation is organised by **hearing, not by distress**. The same axis falls out of a completely separate resting-state EEG cohort. That is where the new datasets should be aimed.

## Timing budget

| Section | Slides | Minutes |
|---|---|---|
| Opening & problem | 1–4 | 3 |
| Design & methods | 5–8 | 4 |
| Whole-brain results (VBM, SBM) | 9–13 | 5 |
| Regional & cross-dataset rigor | 14–16 | 4 |
| Normative modelling | 17–18 | 2 |
| **Heterogeneity** | **19–25** | **6** |
| Genetics & what's next | 26–27 | 1.5 |
| Conclusions & thanks | 28–29 | 1.5 |
| | | **~27 min** |

*29 slides is the full set. The main path is 28 — slide 23 is marked optional and its conclusion is carried by 22 and 24.*

*Trim path if you run long, in order: cut slide 23 (−1 min), slide 12 (SBM boxplots) and slide 18 (normative model quality) (−1.5 min), then compress 9–11 into two (−1.5 min). Cutting the first three lands you at ~24.5 minutes without losing an argument — the heterogeneity section is now the second-largest block in the talk and is where the title gets paid off, so protect 19–22 and 24.*

---

# SLIDES

---

## Slide 1 — Title

**On the slide**
- **The Tinnitus Brain Is Not One Brain**
- *Structural heterogeneity across a multi-site clinical cohort and UK Biobank*
- Multi-site clinical cohort (N = 569) · UK Biobank (N = 4,356)
- Payam Sadeghi-Shabestari — [institution]
- TRI 2026 · [date]

**Figure** — consortium/partner logo strip. `NEEDS MAKING` (or reuse the layout of v2 slide 1).

**Say** — "I'm going to show you the largest and most tightly controlled structural study of tinnitus I'm aware of, and I'm going to be honest with you about how small the effects are — because I think that smallness is itself the finding."

**Watch out** — Set the humble tone in the first 15 seconds. It buys you credibility for everything after.

---

## Slide 2 — The problem: a literature that doesn't converge

**On the slide**
- Decades of structural MRI in tinnitus — findings rarely replicate
- Reported: auditory cortex, ACC, insula, hippocampus, thalamus… in both directions
- Typical N per group: 20–40
- Three plausible reasons: **power**, **hearing loss**, **heterogeneity**

**Figure** — two panels, both built and ready to drop in. Use **B alone** if you only have room for one.

- **A** `plots/talk_tri/slide02a_literature_landscape.pdf` — dot matrix of the real published literature: **57 findings from 25 studies across 12 regions**, one dot per finding, coloured by direction, sized by sample size. Regions reported in both directions sort to the top and are flagged *conflicting*. Data in `src/figures/slide02_literature.csv`, every row carrying a PMID.
- **B** `plots/talk_tri/slide02b_power_curves.pdf` — statistical power vs. sample size. Legend is deliberately anatomy-free, since this slide runs before any results.

Regenerate either with:
```
venv/bin/python src/figures/make_slide02_figures.py
```

**On the slide**
- 25 structural-MRI studies, ~20 years, **median 28 participants per group**
- **7 of 12 regions have been reported in both directions** — auditory cortex, thalamus, hippocampus, frontal cortex, precuneus, cerebellum, basal ganglia
- Mühlau's flagship subcallosal finding: **3 designed replications, 3 nulls** (Landgrebe 2009, Melcher 2013, Allan 2016)
- At d = 0.35 a 30-per-group study finds the effect **27%** of the time. At d = 0.10, **7%**.
- Three explanations: **power**, **hearing loss**, **heterogeneity**

**Say** — "If you read this literature you get a different tinnitus brain in every paper. Before assuming everyone before us was wrong, look at what a typical study could actually detect. At the effect sizes that turn out to be realistic here, a thirty-per-group study finds a moderate effect about one time in four, and a small one about one time in fifteen. That's not a field making mistakes. That's a field sampling noise."

**Two quotes you can use verbatim, both from abstracts:**
- Adjamian et al. 2014 (review of 17 morphometric studies): *"the overall evidence for structural abnormalities specifically related to tinnitus is poor."*
- Vanneste et al. 2015 (n = 154): *"not one single region is replicated in studies of different centers."*

**Watch out** — Don't name-and-shame individual papers; frame it as a field-wide power problem. And expect the obvious counter: *"then why should we believe your small effects?"* The answer is slides 11 and 15 — independent replication in a cohort that actually has the power, and a Bayesian model that explicitly discounts the underpowered arm. Have that ready; it's the strongest moment in the talk if you land it.

**A subtlety worth one sentence if the room is technical:** direction often tracks the *comparator*, not the pathology. Boyen 2013 and Husain 2011 both flip sign depending on whether controls are normal-hearing or hearing-matched, and Makani's 2022 meta-analysis reverses sign between its normal-hearing and hearing-loss strata — concluding *"hearing loss is the driving force of changes in cortical gray matter."* That is a cleaner account of the contradictions than "inconsistent findings", and it sets up slide 3.

**Key numbers from panel B** (printed by the script when you run it):

| Effect | Power at n = 30/group | n/group for 80% power |
|---|---|---|
| d = 0.35 | **27%** | 130 |
| d = 0.20 | **12%** | 394 |
| d = 0.10 | **7%** | 1,571 |

The middle row independently confirms the "~400 per group for 80% power" claim on slide 13 — worth repeating there, since it is exactly why your own STS cluster at d = 0.202 is fragile in a 285-per-group sample.

---

## Slide 3 — Why hearing loss is the hard confound

**On the slide**
- Tinnitus and hearing loss co-occur in the large majority of patients
- Hearing loss alone changes brain structure
- Most prior studies matched loosely, or not at all
- **Our approach: PTA matched at recruitment, modelled as a covariate, and tested again with high-frequency PTA**

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/audiometry/audio.pdf` (audiograms) and, if room, `plots/audiometry/sites.pdf` to show cross-site audiometric comparability.

**Say** — "Any structural finding in tinnitus has to survive the question 'is this just hearing loss?' We built that question into the design rather than answering it in a limitations paragraph."

**Watch out** — This slide pre-empts the single most likely hostile question. Land it clearly.

---

## Slide 4 — Heterogeneity: the third explanation

**On the slide**
- Tinnitus is a symptom, not a disease
- Patients may not be one population — and heterogeneity comes in two flavours: **discrete subtypes**, or a **continuous axis** along which patients vary
- Either one degrades the group mean; they are not the same claim and they need different evidence
- **If patients are heterogeneous, the "average tinnitus brain" describes nobody**
- This is the **hypothesis the talk tests**, not a result being asserted — and by slide 23 I'll tell you which flavour the data support

**Figure** — `NEEDS MAKING`: a three-panel schematic — left, the group mean you actually observe; middle, the categorical account (two separate distributions); right, the dimensional account (one distribution with a shifted tail). Both middle and right produce the same weak group mean. This is the conceptual spine of the whole talk and is worth 20 minutes of drawing time. **Draw the right-hand panel carefully — it is the one the results end up supporting**, and slide 28 reuses this figure.

**Say** — "This is the idea I want you to hold for the rest of the talk. If a population is heterogeneous, the group mean describes none of it. But notice there are two versions of that claim — patients falling into kinds, or patients spread along a dimension — and they look identical from the outside. Distinguishing them is most of what the second half of this talk is about."

---

## Slide 5 — Design: discovery and replication

**On the slide**
- **Discovery:** multi-site clinical cohort — 13 recruiting sites, **10 with local controls**
- **N = 569** after QC and 1:1 matching (284 CO / 285 TI)
- **Replication:** UK Biobank — **2,180 TI vs 2,176 CO** (N = 4,356)
- Deliberately different: clinical phenotyping vs. population self-report

**Figure** — `NEEDS MAKING`: a two-box design diagram (clinical cohort → UKB), with N and phenotyping method under each. Or reuse the layout of **v2 slide 2** (site table) alongside **v2 slide 3** (UKB).

**Say** — "These two cohorts are not redundant. The clinical cohort has real audiology and validated distress instruments; UKB has eight times the N and a completely independent ascertainment. Anything that shows up in both is unlikely to be a scanner effect or a recruitment effect."

**Watch out** — Say **"13 recruiting sites, 10 contribute to the matched analysis"** — Erlangen, ICCAC and IndiTMS recruited tinnitus patients only, so 1:1 matching drops them.

---

## Slide 6 — The discovery cohort

**On the slide**
- 10 sites, 3T and 1.5T, MPRAGE / TFE / MPM
- CO: age 57.0 ± 17.9 · TI: age 56.7 ± 16.7
- PTA: 19.7 ± 14.3 vs 19.9 ± 14.2 dB — **matched on hearing**
- THI available for 167 tinnitus patients (mean 30.9 ± 23.5)
- 1:1 nearest-neighbour matching, caliper 0.2, on sex/age/site/PTA

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/demographics/group_distribution.pdf` + `plots/demographics/sex_distribution.pdf`. Per-site N table: reuse **v2 slide 2**.

**Say** — "Matching on PTA rather than covarying for it after the fact is the design choice I'd defend hardest here."

**Watch out** — Sex is coded 0/1 in `subcortical_roi/master.csv` with **no legend anywhere in the tree**. Confirm which code is male before putting a sex breakdown on a slide.

---

## Slide 7 — The UK Biobank replication cohort

**On the slide**
- 2,180 tinnitus vs 2,176 controls
- Matched on sex, age, site, noise exposure, hearing-aid use, income, handedness
- Hearing proxy: **speech-reception threshold (SRT)** in place of clinical PTA
- Scanner harmonisation: neuroHarmonize / ComBat, trained on controls only
- Trade-off: 8× the power, coarser phenotype

**Figure** — `plots/demographics/ukb_distributions.pdf` (and `ukb_distributions_2.pdf` if you want the second panel).

**Say** — "UKB buys us power and independence. It costs us phenotype quality — self-reported tinnitus rather than a clinical diagnosis, and a speech-in-noise measure rather than an audiogram. Remember that trade when we get to the findings that don't replicate."

---

## Slide 8 — Analysis pipeline

**On the slide**
- Whole brain: **VBM** (FSL randomise, TFCE, 5,000 perms) · **SBM** (FreeSurfer, cluster-wise MC)
- Regional: 6 atlases, 361 regions — hippocampal subfields, amygdalar & thalamic nuclei, aseg, aparc, aparc.a2009s
- Individual level: **normative modelling** (PCNToolkit BLR, trained on controls)
- Cross-dataset: frequentist meta-analysis + **Bayesian hierarchical model**
- Subtypes: **GMM** on normative Z-scores

**Figure** — `NEEDS MAKING`: one clean pipeline diagram, left-to-right, five boxes. This replaces the "Outlook" navigation slides in the current deck (v2 slides 5, 14, 19, 22, 26, 29, 37, 40, 50, 60) — you do not have time to show that navigator ten times.

**Say** — "Five levels of analysis, from whole-brain voxels down to individual deviation from a normative population model. The point of that redundancy is that a real effect should show up at more than one level."

**Watch out** — Keep this to 45 seconds. It's orientation, not content.

---

## Slide 9 — VBM: the primary whole-brain result

**On the slide**
- TFCE-corrected, covariates age / sex / site / PTA / TIV
- **TI > CO: right putamen** (peak MNI 26, 10, −8; 4,496 voxels; p = 0.0014)
- **TI > CO: right cerebellum VIIIa** (peak 36, −46, −56; 8,070 voxels)
- **CO > TI: nothing** (max 1−p = 0.62)
- Putamen effect: **|d| ≈ 0.35**

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/vbm_group/ti_gt_co_x.pdf` (plus `_y.pdf` / `_z.pdf` if you want three planes). Equivalent to **v2 slide 6**.

**Say** — "One direction only. Tinnitus patients have *more* grey matter in the putamen and cerebellum, and nothing anywhere survives in the opposite direction."

**Watch out** — The archived `d_map.nii.gz` is signed CO−TI, so it reads negative in these clusters. **Quote |d| and state the direction in words.** Also: the cluster table here was regenerated with FSL for this talk — `vbm_results/` archives no stats text. Save a proper cluster table before the manuscript.

---

## Slide 10 — VBM: does it survive hearing loss?

**On the slide**
- Re-run with **high-frequency PTA** as an additional covariate
- **Left putamen / pallidum persists** (1,459 voxels, peak −22, 2, −12)
- Cerebellar uvula finding **disappears** → was likely hearing-loss driven
- CO > TI emerges weakly: left frontal pole, subgenual/vmPFC

**Figure** — `plots/vbm_group/ti_gt_co_pta_hf.pdf` and `plots/vbm_group/co_gt_ti_pta_hf.pdf`. Equivalent to **v2 slide 11**.

**Say** — "This is the slide that tells you which findings are about tinnitus and which are about hearing. The putamen survives. The cerebellar vermis does not — and I think that's genuinely informative rather than disappointing."

**Watch out** — Do **not** carry over the deck's subcallosal claim: `cluster_subcollasal_mask` has mean d = +0.002, i.e. null. The vmPFC cluster here is a separate, weaker observation — present it as such or leave it off.

---

## Slide 11 — VBM: independent replication in UK Biobank

**On the slide**
- Independent nilearn/PALM implementation, TFCE, 5,000 permutations
- **TI > CO** survives: left cerebellum VIIb (4,790 vox), right lateral occipital, **left putamen (71 vox, peak −16, 8, −10)**
- **CO > TI: nothing** (max −log10 p = 0.48)
- Same direction, same structures, entirely different cohort and pipeline

**Figure** — `plots/vbm_group/ukb/ti_gt_co_x_thr_0.005.pdf` (or the `_0.05` variant, which is where the putamen appears). Equivalent to **v2 slide 8**.

**Say** — "Different country, different scanners, different phenotyping, different software — and the putamen and cerebellum come back, in the same direction. That convergence is the strongest evidence in the talk."

**Watch out** — Be precise: the UKB putamen cluster is **small and at the more liberal threshold**. Don't oversell. Also `UKB/TI_gt_CO_effect_size.nii` ranges ±0.017 — that is **not** Cohen's d; never quote it as one.

---

## Slide 12 — SBM: cortical thickness

**On the slide**
- FreeSurfer `mri_glmfit`, N = 569, DOF = 559, cluster-wise MC correction
- **All significant clusters are TI > CO — thicker cortex in tinnitus**
- rh: superiorparietal 1,526 mm² (CWP = 0.0002) · inferiorparietal 753 mm² (0.006) · middletemporal 558 mm² (0.036)
- lh: lateraloccipital 1,021 mm² (0.0002) · fusiform 799 mm² (0.003)
- Nothing survives at the stricter p<0.001 cluster-forming threshold

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/sbm/brain.pdf`. Equivalent to **v2 slide 20**.

**Say** — "The cortical picture matches the subcortical one in direction — more tissue, not less — but it's threshold-fragile."

**Watch out** — These clusters exist only at the liberal (p<0.01) cluster-forming threshold. Say so out loud; someone will ask.

---

## Slide 13 — SBM: the honest version

**On the slide**
- Sensitivity analysis with **site as a fixed effect (DOSS)**: one cluster survives
- Right inferior parietal / STS — 396 mm², **CWP = 0.047**, **d = 0.202**
- ~400 per group needed for 80% power at this effect size
- **Does not replicate in UKB**
- Presented as **hypothesis-generating, not established**

**Figure** — `plots/sbm/brain_doss.pdf`. Equivalent to **v2 slide 20** (DOSS panel).

**Say** — "A cluster-wise p of 0.047, an effect size of 0.2, in a sample marginally below the power needed to find it, that doesn't replicate in a cohort eight times larger. I'm showing you this because I don't think it should be in a paper as a finding."

**Watch out** — This slide will win you more credibility than any positive result on the deck. Deliver it without apology.

---

## Slide 14 — Regional analysis: what 361 regions show

**On the slide**
- 6 atlases, 361 regions, ANCOVA controlling age / sex / site / PTA / TIV
- **Discovery cohort: 0 regions survive FDR** (min p_fdr = 0.156)
- Strongest uncorrected trends, all TI > CO: right putamen (d = +0.216), right central amygdala nucleus (+0.211), right fusiform (+0.237)
- **UKB: 22 of 135 regions survive Bonferroni** — thalamic (8), hippocampal (5), aseg (6), amygdalar (3)
- All UKB effects TI > CO, **d = 0.097–0.120**, except smaller ventricles in tinnitus

**Figure** — `plots/rois/cohen.pdf` (Cohen's d across nuclei) or `plots/rois/heatmap.pdf`. UKB-specific: `plots/rois/thalamic_nuclei_manhattan_ukb.pdf`. Equivalent to **v2 slides 30–32**.

**Say** — "In the clinical cohort, nothing survives correction across 361 regions. In UKB, twenty-two regions do — and every one of them points the same way. That contrast is a power story, not a contradiction."

**Watch out** — **Never** describe a main-cohort regional result as significant. This is the single easiest place to overclaim.

---

## Slide 15 — Do the two datasets tell the same story?

**On the slide**
- Bayesian hierarchical model per region: shared effect θ, between-dataset heterogeneity τ
- 79 of 97 regions agree in direction
- Nominal / meta-analytic criteria: **44 regions "consistent"**
- **P(\|θ\| > 0.05) ≥ 0.80 and P(same sign) ≥ 0.90: only 3 regions**
- **Accessory-basal amygdala (L)** · **medial amygdala (R)** · **MDm thalamus (R)**
- **All three 95% HDIs still cross zero** — and so does every one of the 97 regions tested

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/talk_tri/slide15_bayesian_forest.pdf` (built; `.png` alongside it). Two panels on a shared Cohen's d scale: left, the 3 surviving regions shown three ways — clinical cohort alone, UK Biobank alone, and the pooled posterior with its 95% HDI; right, all 44 nominally-consistent regions as pooled θ with HDIs, none of which excludes zero. Generated by `src/figures/make_slide15_forest.py` from `subcortical_roi/bayesian/bayesian_dmin_0.05.csv`. Keep `plots/tier_heatmap_v2.pdf` as a backup if you want the full 97-region view.

*The left panel is the one to talk over: UK Biobank is tight and excludes zero, the clinical cohort is wide, and the pooled posterior is **wider than UKB alone** because the hierarchical model pays for between-dataset heterogeneity. That is the mechanism behind the 44-vs-3 gap, and it is visible rather than asserted.*

**Say** — "This is the slide I care most about. Depending on which criterion you pick, this project has found forty-four regions or three. I'm showing you both numbers because the gap between them is what most imaging papers hide."

**Watch out** — Be exact on three points. (i) The tier threshold met is **P(|θ|>0.05) ≥ 0.80** (values 0.85–0.87), *not* the ≥0.95 tier. (ii) The HDIs are **95%**, not 94% — `10_compare_datasets_bayesian.py:121` passes `hdi_prob=0.95` explicitly, overriding the ArviZ default. (iii) The **binding criterion is direction consistency, not effect size**: of the 44, fourteen reach P(|θ|>0.05) ≥ 0.80 but only three reach P(same sign) ≥ 0.90. Say the HDIs cross zero before anyone finds it — and you can go further, because **not one of the 97 regions has a 95% HDI excluding zero.**

---

## Slide 16 — The direction nobody expects

**On the slide**
- Across both cohorts, across VBM, SBM and 361 regions: **effects run TI > CO**
- UKB: 21 of 22 Bonferroni-significant regions are increases
- Only consistent decrease: **ventricular volume** (i.e. more tissue, not less)
- Not the atrophy pattern the literature predicts
- Candidate explanations: activity-dependent plasticity · neuroinflammation · selection effects · reverse causation

**Figure** — `plots/rois/boxplots_Amygdala.pdf` and `plots/rois/boxplots_Thalamic-nuclei.pdf` side by side.

**Say** — "I want to flag this rather than explain it away. Almost everything we find is bigger in tinnitus. With cross-sectional data I cannot tell you whether that's a consequence of chronic altered input, a pre-existing risk trait, or something about who ends up in a tinnitus clinic."

**Watch out** — This is your most interesting discussion point. Leave room for it; expect a question.

---

## Slide 17 — Beyond group means: normative modelling

**On the slide**
- PCNToolkit Bayesian linear regression, trained on **controls only**
- Predicts expected volume given age, sex, site, hearing, TIV → per-subject deviation Z
- Test: are tinnitus patients over-represented in the tails (\|Z\| > 1.96)?
- **Right VA thalamic nucleus: 9.1% of patients vs 2.1% of controls** — χ²(1) = 11.88, p = 5.7 × 10⁻⁴, **p_FDR = 0.028**
- The **only** FDR-significant normative finding in the entire project

**Figure** — `plots/rois_norm/main/Thalamic-nuclei_extreme_count.pdf` + `..._extreme_stats.pdf`. Equivalent to **v2 slide 38**.

**Say** — "VA didn't show up in any group-mean comparison, and that's the point. A minority of patients deviate hard; the majority don't; the mean doesn't move. And VA is the primary thalamic output target of the putaminal circuit we saw in VBM — so these two findings are the same circuit seen two ways."

**Watch out** — Near-misses worth having ready: left cortical amygdala nucleus (3.2% vs 9.5%, p_FDR = 0.064) and right presubiculum body (4.9% vs 11.9%). Mention only if pressed.

---

## Slide 18 — Are the normative models trustworthy? *(cut first if short on time)*

**On the slide**
- Per-feature Rho, MSLL, explained variance
- Clinical-cohort models outperform UKB-trained models
- Attributed to protocol homogeneity vs. UKB scanner heterogeneity
- → individual deviation estimates from the clinical cohort are the more reliable ones

**Figure** — `plots/VBM_Normative_Metrics_Final.pdf` or `plots/SBM_Normative_Metrics.pdf`. Equivalent to **v2 slides 17 / 24**.

**Say** — "Before I trust a deviation score I want to know the model predicts controls well. It does — better in the clinical cohort than in UKB, which is worth knowing when we weigh the two."

**Watch out** — `vbm_norm` / `sbm_norm` `statistics_test.csv` contain **numerically exploded R²/SMSE/EXPV** (values up to ~1e78). **Rho is fine (mean 0.41 / 0.38) — show Rho only.**

---

## Slide 19 — Looking for structure inside the patient group

**On the slide**
- The group mean is the wrong summary if patients are not one population
- Normative modelling already hinted at this: **9.1% of patients deviate hard in VA thalamus vs 2.1% of controls** — a minority moves, the mean barely does
- So: cluster patients on their **own deviation profiles**, not on the case/control contrast
- Features: **11 convergent subcortical Z-scores** (tier 2+3 amygdalar, hippocampal, thalamic nuclei)
- Hearing loss and site **residualised out before clustering** — a second time, on top of the normative model
- → StandardScaler → **PCA to 85% variance (6 components)** → Gaussian mixture, k scanned 1–6, BIC with a 5%-of-N minimum-cluster guard
- **k = 2.** Stability: **median ARI 0.747, IQR 0.25–0.84, 67% of resamples ≥ 0.60**

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/biotypes_gmm/primary/01_bic_selection.pdf` and `.../02_bootstrap_ari.pdf`. Optionally `00_pca_variance.pdf` if you want to show the 6-component choice.

**Say** — "We partial hearing out *before* clustering, precisely so nobody can say we've rediscovered the audiogram. And I'm giving you the interquartile range on the stability, not just the median, because the median on its own flatters this result."

**Watch out** — Two things the earlier version of this slide got wrong, now fixed. The GMM runs on **6 principal components**, not on the 11 raw Z-scores — say "PCA then mixture model" if anyone asks about dimensionality. And PTA and site are already covariates and batch effects **inside** the PCNToolkit normative model (`covariates = ['age','sex','PTA','TIV']`, `batch_effects = ['site']`), so this is the *second* time they are removed. That is a strength, but state it accurately.

---

## Slide 20 — The site artefact we caught

**On the slide**
- First run (all 285 patients, 10 sites): clean k = 2, **silhouette 0.52, bootstrap ARI 0.94**
- But the minority cluster was a scanner: **87% of london_2's patients fell into it**, and london_2 supplied **48% of the whole cluster** (13 of 27)
- Site χ² = 118.0, p < 0.001 — and **nothing else** was significant: not PTA, not THI, not age, not sex
- Re-run excluding that site: **N = 270, 9 sites**, k = 2, silhouette 0.121, ARI 0.747
- **The clean-looking solution was the wrong one**

**Figure** — `plots/biotypes_gmm/primary/09_site_composition.pdf`. Optionally show `diagnostic_all_sites/02_bootstrap_ari.pdf` beside `primary/02_bootstrap_ari.pdf` to make the "better statistics, worse science" point visually.

**Say** — "The version with a silhouette of 0.52 and an ARI of 0.94 was a scanner. The version I'm about to show you has worse numbers and is the real result. I'm including this slide because leave-one-site-out is not standard practice in this literature and it should be."

**One more line, if you want the strongest version of this point:** the same thing happened independently in our resting-state EEG cohort — the best-silhouette clustering solution there (0.415) turned out to be dominated by a single site, and those clusters differed on no clinical variable at all. Two datasets, two modalities, the same trap. That generalises the slide from "we caught a bug" to a field-level methods warning.

**Watch out** — This is arguably the best slide on the deck. Don't rush it. Note the corrected attribution: **87% is of london_2's patients, not of the cluster.** The earlier phrasing ("86.7% of the minority cluster came from a single site") was wrong and is checkable.

---

## Slide 21 — What actually separates the two groups in the brain

**On the slide**
- Biotype 0 (n = 58) sits **above** the normative mean on **10 of 11** regions; Biotype 1 (n = 212) sits just below
- Strongest: **MDl-lh** (d = 0.43) · **Accessory-Basal-lh** (0.42) · **CeM-lh** (0.38) · **Basal-lh** (0.36) · **Lateral-lh** (0.35)
- **8 of 11 regions are left-hemisphere**, and every strong one is
- The difference is a **uniform shift in magnitude, not a different spatial pattern**
- No region is "on" in one group and "off" in the other — one group deviates more, in the same direction

**Figure** — `plots/biotypes_gmm/primary/05_zprofile_partialed.pdf` (the Z-profile) and `.../07_radar_primary_(excl._london_2).pdf` (the radar signature). Both exist. Note the filename contains parentheses and a period.

*Optional and better:* a forest plot of the 11 Cohen's d values with CIs. The table below has never been written to disk — only the PDFs exist — so making it is a small job worth doing.

| region | B0 mean Z | B1 mean Z | Cohen d | p_FDR |
|---|---|---|---|---|
| MDl-lh | +0.247 | −0.068 | 0.429 | 0.067 |
| Accessory-Basal-nucleus-lh | +0.299 | −0.082 | 0.418 | **0.027** |
| CeM-lh | +0.202 | −0.055 | 0.378 | **0.027** |
| Basal-nucleus-lh | +0.216 | −0.059 | 0.357 | **0.027** |
| Lateral-nucleus-lh | +0.177 | −0.048 | 0.352 | **0.029** |
| PuA-lh | +0.230 | −0.063 | 0.345 | 0.150 |
| MDm-lh | +0.239 | −0.066 | 0.308 | 0.298 |
| presubiculum-head-rh | +0.135 | −0.037 | 0.244 | 0.101 |
| CeM-rh | +0.123 | −0.034 | 0.235 | 0.116 |
| MDm-rh | +0.088 | −0.024 | 0.143 | 0.574 |
| GC-ML-DG-body-rh | −0.061 | +0.017 | −0.099 | 0.503 |

**Say** — "This is the slide I'd want to see if I were in the audience, and it's the one that changed how I describe this result. The two groups don't have different brains. One group's limbic and thalamic nuclei sit further from normative expectation than the other's, in the same direction, on essentially every feature. That's a difference of degree."

**Watch out** — **Say the circularity out loud before anyone else does.** These are the same 11 features the mixture model clustered on, so the p-values are not inferential — they describe the profile, they don't test it. And in **raw, unpartialled Z-space nothing survives FDR at all**. If you show this table, show it as a description of the solution, not as evidence that the solution is real. Slide 23 is where the evidence question gets answered.

---

## Slide 22 — The axis is audiological

**On the slide**
- **PTA: 27.7 vs 18.1 dB** — p_FDR = 0.0001
- **Age: 68 vs 59.5 years** — p_FDR = 0.0002
- Sex differs too (p_FDR = 0.021)
- **The group whose brains deviate more is the group with more hearing loss** — even though hearing was removed twice
- THI: 14 vs 28 (p_FDR = 0.023) — **but see the caveat; do not lead with this**

**Figure** — `plots/biotypes_gmm/primary/06_clinical_clinical.pdf`.

**Say** — "Hearing was a covariate in the normative model and was regressed out again before clustering, and the clusters still separate on it at p = 0.0001. I read that as the audiogram not being the right hearing variable rather than as the brain secretly encoding PTA — and I'll come back to why on the EEG slide."

**Watch out — three things, and the third is the important one.**

1. **Only five variables were testable**: PTA, THI, age, sex, site. That is the entire clinical characterisation available in the merged cohort. "Clinically distinct" is more than the data support; say "separates on the audiological axis".
2. **Site is one of the five significant results** (χ² = 22.09, p_FDR = **0.0079**) — *after* excluding london_2 and *after* regressing site out of the features. Biotype-0 share runs from 5% at tinspect to 39% at cogtail. Disclose it; someone will ask why site is still there.
3. **The THI result is probably a site effect and should not be presented as a finding.** THI mean by site ranges from 17.2 (cogtail) to 58.9 (triple), and it was never collected at andes (58 patients) or london_2. **cogtail contributes 25 of Biotype 0's 58 patients — 43% — and has the lowest THI of any site.** The test is a plain Mann-Whitney with no site adjustment, on 43 vs 124 patients, with Biotype 0 also 8.5 years older. The previous version of this deck led with "the subtype with less hearing loss has more distress." That claim cannot be defended until the test is redone with site in the model. If asked, say exactly that.

---

## Slide 23 — Category or continuum? *(optional — cut this first if short on time; 22 and 24 both carry the conclusion)*

**On the slide**
- Silhouette **0.121**
- BIC prefers k = 2 over k = 1 by **ΔBIC = 33.9** — real, but modest
- Bootstrap ARI **IQR 0.25–0.84**; a third of resamples give a substantially different partition
- Membership varies **5% → 39% across nine sites**
- In UK Biobank the same model splits **controls (23.8%) almost exactly like cases (26.0%)**
- **No null model was ever run** — no dip test, no label-shuffling null, no simulated-unimodal comparison

**→ We report an axis of deviation, not two diseases.**

**Figure** — reuse `01_bic_selection.pdf` and `02_bootstrap_ari.pdf` from slide 19, or build a single "evidence for discreteness" summary panel. `NEEDS MAKING` if you want the summary panel.

**Say** — "I want to be straight about what this is. k = 2 is a convenience — it's where BIC lands once you forbid tiny clusters. What the data actually support is that patients differ continuously in how far their subcortical structure sits from normative expectation, and that this variation is organised by hearing. Calling that two biotypes would be dressing a dimension up as a taxonomy."

**What would change my mind** — a box worth putting on the slide:
- A Hartigan dip test rejecting unimodality on the mixture posterior
- Bootstrap ARI stable when the PCA is refitted inside each resample
- Membership predicting something the model never saw — treatment response, EEG signature, trajectory
- The split surviving in UK Biobank **tinnitus cases but not controls**

**Watch out** — The UKB number is preliminary. The transfer has two known defects (no group filter, and UKB Z-scores were never residualised because `ukb_covars.csv` does not exist), so present it as "the first look is not encouraging", not as a settled negative result. See slide 25.

---

## Slide 24 — The same axis appears in resting-state EEG

**On the slide**
- Separate cohort, separate modality, separate paradigm: **TIDE consortium, 7 sites, 268 tinnitus / 276 controls**, resting-state EEG connectivity
- A classifier (AUC 0.82) **misses 85 of 268 patients (32%)** — and these are not borderline cases: median predicted probability **0.17**
- The missed patients have **better high-frequency hearing** (HF-PTA 38.6 vs 48.8 dB, d = −0.34) and **identical distress** (THI d = 0.10, TFI d = 0.01)
- Match them to identical HF-PTA and the neural difference **grows** (d = −1.14 → −1.27) — so it is not a hearing confound
- Same conclusion about structure: **no discrete cluster** (silhouette 0.072), a continuum
- **And the dose-response runs the same way in both.** Proportion of patients whose brains deviate from norms, by hearing-loss quartile:
  - **MRI: 8% → 16% → 28% → 34%** (computed for this talk; a 4-fold gradient, monotonic)
  - **EEG: 63% → 57% → 69% → 85%** (i.e. miss rate 37% → 43% → 31% → 15%)

| | Tinception (MRI) | TIDE (EEG) |
|---|---|---|
| Deviates from norms | Biotype 0, n = 58 (21.5%) | detected, n = 183 (68.5%) |
| Looks normative | Biotype 1, n = 212 (78.5%) | missed, n = 85 (31.5%) |
| Hearing, deviant group | PTA 27.7 dB | HF-PTA 48.8 dB |
| Hearing, normative-looking group | PTA 18.1 dB | HF-PTA 38.6 dB |
| Discrete? | silhouette 0.121 | silhouette 0.072 |

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/plots/talk_tri/slide24_mri_eeg_convergence.pdf` (built; `.png` alongside it). Two columns, three rows: the deviation distributions, median hearing per group, and the quartile dose-response. Generated by `src/figures/make_slide24_convergence.py` — the MRI side is computed from the result files at run time, the EEG side is transcribed from `~/Tinnorm` with per-value provenance in the constants block. The two x-axes are deliberately **not** on a common scale, and the figure says so in its footer.

*Note the MRI row-1 curves are the observed distribution of all 270 patients; the EEG row-1 curves are normal approximations from the reported mean ± SD, which is stated under the panel. If anyone asks, say that plainly.*

**Say** — "This is from a different project, and I'm showing it because it answers the obvious objection to everything on the last three slides. Look at the bottom row: in both cohorts, the worse the hearing, the more likely a patient's brain is to look abnormal — eight percent to thirty-four in ours, and the same direction in the EEG data. Two multi-site cohorts, largely different people, one structural and one functional, one unsupervised and one supervised. Both find that the dominant axis of variation inside the tinnitus group is hearing, not severity. And both find it's a gradient, not a category. In the EEG cohort you can equate the two groups exactly on high-frequency hearing and the brain difference gets *bigger* — so this isn't the audiogram in disguise."

**And the payoff for slide 22:** the EEG effect lives in **high-frequency** hearing. Our structural analysis controls for PTA — a 0.5–4 kHz average that cannot see above 4 kHz at all. That is a concrete, testable reason why hearing keeps reappearing after we thought we had removed it, and it does not require anything mysterious.

**Watch out**
- **Clear this with the TIDE co-authors before you present it.** It is their unpublished data, from a different consortium.
- **Do not quote the "controls with worse hearing don't get higher scores" sensitivity check (SA1).** It conditions on the classifier's own outcome, and the misclassification analysis in the same project shows false-positive controls *do* have much higher HF-PTA than true-negative controls (37.5 vs 20.4 dB, p = 0.003). The matched-pairs result and the ANCOVA are the ones that carry the argument — use those.
- The two cohorts are **largely different people** (different consortia, different site sets), which makes the agreement worth more, not less — but **verify there is no subject overlap** before you say "independent".
- Two further parallels if you have time: both projects find tinnitus means *more*, not less (EEG hyperconnectivity 18:1 positive-to-negative deviation; MRI larger volumes), and **neither modality points to auditory cortex** (their top feature is right pars orbitalis; auditory cortex ranks 8th).

---

## Slide 25 — Does the axis exist in UK Biobank?

**On the slide**
- Main-cohort mixture model applied to UKB by posterior-probability transfer
- 84.3% of subjects classified with certainty ≥ 0.70 — the structure transfers technically
- **But it does not discriminate**: UKB tinnitus splits **26.0% / 74.0%**, UKB controls **23.8% / 76.2%**
- The proportion "match" to the clinical cohort's 21.5% / 78.5% is a **whole-sample** number, cases and controls pooled
- **Open question: is this a tinnitus axis, or normal population variation in subcortical morphometry?**

**Figure** — `plots/biotypes_gmm/ukb_transfer/08_ukb_transfer_summary.pdf` and `.../04_umap_ukb.pdf`. Equivalent to **v2 slide 42**. **Both must be regenerated once the transfer is fixed.**

**Say** — "I want to be careful here. In preparing this talk I found the UKB arm wasn't restricted to tinnitus cases, and when you split it properly the controls come out looking like the cases. That is not what you want to see. It's a fixable analysis and it's at the top of my list — but until it's fixed, I'm not going to tell you the axis is tinnitus-specific."

**Watch out**
- The transfer has **two** defects, not one. `load_ukb_features()` (`src/main/13_biotype_gmm.py:214`) has no group filter; and `ukb_covars.csv` does not exist at either path the script looks for, so **UKB Z-scores were never residualised for hearing or site** while the clinical ones were. The model is being applied across a preprocessing mismatch. Fix both before drawing a conclusion.
- Traps in the results folder: the top-level `biotypes_gmm/biotype_assignments_ukb.csv` says 4,019/337 and **contradicts** `ukb_transfer/` (1,084/3,272) — use only `ukb_transfer/`. `biotype_stats_ukb.csv` is a **1-byte empty file** in both locations.
- `VBM_design/thi.csv` and the `biotypes/` folder hold the **rejected UMAP/HDBSCAN** labels (PTA 9.4 vs 29.2 dB — the exact artefact you discarded). Don't mix them in.

**Do not present** the "limbic–salience dysregulation subtype" whole-brain contrast from v2 slides 45–48. Those maps (`vbm_results/Bio_TFCE/`) are null — 0 voxels survive FWE at p < 0.05, max −log10 p = 1.109 — **and they were computed from the superseded UMAP/HDBSCAN labels**, not the mixture-model biotypes. The whole-brain contrast under the current labels **has never been run**. If asked, say both halves of that.

---

## Slide 26 — Genetics: where we are

**On the slide**
- UKB GWAS, REGENIE two-step, 15 phenotypes (case/control, biotype, brain PCs, 11 regions)
- **λ_GC = 1.00–1.08 — no inflation**
- **0 genome-wide-significant hits**; 25 suggestive (p < 10⁻⁵)
- Recurring themes at suggestive loci: **HPA-axis (CRHR1)**, **calcium signalling**, **synaptic vesicle machinery (NSF)**
- Hypothesis-generating only — and on a ~180K-SNP pruned panel

**Figure** — `/Volumes/Extreme_SSD/payam_data/Tinception/GWAS/plots/tinnitus_subtype.manhattan.png` with `tinnitus_vs_controls.qq.png` inset. Equivalent to **v2 slides 51–52**.

**Say** — "Clean lambda, no hits. I'm showing this in one slide rather than five because that's proportionate to the evidence. The interesting design point is that a *neuroanatomically defined* phenotype is a better-powered question than the case/control one — it's more homogeneous than a yes/no symptom report."

**Watch out — two things.**

1. **This GWAS did not use the biotypes you showed on slides 19–25.** The phenotype comes from `biotypes/ukb.csv` (n = 1,512 — the "n = 1,497" in the summary), which holds the **superseded UMAP/HDBSCAN** labels. Nothing in the repo has been re-pointed at `biotypes_gmm/`. Say "a data-driven subtype GWAS" and don't imply continuity with the mixture-model result. The same applies to `plots/correlations/putamen_vs_THI.pdf`, which `12_run_correlations.py` builds off the same old 132/153 split.
2. If asked about PRS: one of 555 tests survived correction (bipolar PRS × right CA1-body, p_adj = 0.047). Say plainly that this is what you'd expect by chance.

---

## Slide 27 — What's next

**On the slide**
- **New datasets incoming** → full pipeline re-run
- **Test the axis properly, not just describe it**: derive **high-frequency PTA** (the audiograms exist at 9 of 10 sites, up to 16 kHz) and ask whether the axis is really an extended-high-frequency axis
- Demonstrate heterogeneity **without clustering** — a variance-ratio test of patient vs control deviation spread makes the point with no cluster solution to defend
- Formal discreteness tests the pipeline has never had: **Hartigan dip test**, simulated-unimodal null for the ARI
- Merge the clinical archive — **TFI, TQ, HADS, BDI, hyperacusis, LDL, duration, laterality, psychoacoustics** all exist per-site and were never merged
- In progress: structural covariance gradients + Allen transcriptomics · mode-based (eigenmode) morphometry
- Open question from this analysis: **does hearing asymmetry explain any of it?**

**Figure** — `plots/ssa/brain_diff_grad_3.pdf` (gradient difference map) as a "coming soon" visual; optionally `MBM/lh_mbm_result.png`.

**Say** — "The honest summary of the extra analyses is that they're built and not yet interpretable. The eigenmode analysis in particular currently returns chance-level coverage, so I'm not going to show you a map of it. What I would rather spend the new data on is the hearing question — because if this axis turns out to be an extended-high-frequency axis, that is a much more useful thing to know than a subtype label."

**Watch out** — Present MBM and gradients as *in progress*, not as results. The MBM threshold map is uncorrected p≤0.05 with 5.6%/6.0% coverage — the null rate.

---

## Slide 28 — Conclusions

**On the slide**
1. Structural differences in tinnitus are **real, replicable, and small** (d ≈ 0.10–0.35)
2. They run **TI > CO** — larger volumes, not atrophy — across both cohorts and every method
3. Convergent circuit: **putamen → VA thalamus → medial/basal amygdala → MDm**
4. **Rigor changes the count**: 44 regions by nominal criteria, 3 under a Bayesian model, HDIs still crossing zero
5. **Heterogeneity is the tractable path forward** — but as a **continuous, hearing-organised axis of deviation**, not two discrete subtypes. The same axis appears independently in resting-state EEG.

**Figure** — the circuit schematic from slide 4 or a simple putamen→VA→amygdala/MDm diagram. `NEEDS MAKING`.

**Say** — "If you take one thing away: the reason this literature doesn't converge probably isn't that everyone measured badly. It's that we've all been averaging over a heterogeneous population and reporting the residue."

---

## Slide 29 — Acknowledgements

**On the slide**
- 13 recruiting sites and their PIs — Delano · Neff · Meyer · Müller · Kleinjung · Schoisswohl · Sedley · Cederroth
- Every participant across both cohorts
- UK Biobank (application number)
- Funders
- Contact + a note that the analysis code is available

**Figure** — logo strip.

**Say** — "This is a consortium result. Ten sites contributed matched controls, which is genuinely hard to organise, and it's the reason the study has the power it does."

---

# Deep dive — what the subtypes are, and what they are not

*This is the material to have in your head for slides 19–25 and for the questions that will follow them. It goes further than the slides do, deliberately.*

## 1. The split is a difference of magnitude, not of pattern

This is the single most important thing to understand about the result, and it is not what the earlier version of this deck said.

Biotype 0 sits **above** the normative mean on 10 of 11 clustering features; Biotype 1 sits marginally below. There is no region that is elevated in one group and reduced in the other. There is no second, alternative profile. What the mixture model found is a group of patients whose left limbic and thalamic nuclei deviate further from age-, sex-, hearing- and TIV-adjusted expectation, and a larger group whose nuclei deviate less — in the same direction, on essentially every feature.

That matters for three reasons. It explains why every whole-brain contrast between the groups is null: a uniform offset of d ≈ 0.4 concentrated in small subcortical nuclei is not something a voxelwise TFCE contrast at n = 58 vs 212 will find. It explains why the silhouette is 0.12: a shifted distribution and its parent overlap heavily by construction. And it means the honest description is a **dimension**, because "differs in degree along a common axis" is what a dimension is.

The EEG project reached the identical conclusion from the other direction, and its wording is worth borrowing: the missed patients deviate "in the same direction… but the magnitude is smaller", and "there is no feature on which FN subjects show a distinctively different (third) pattern."

## 2. "Biotype" is the wrong word, and switching it strengthens the talk

A biotype claim is a claim about kinds — that patients fall into categories with different underlying biology. Nothing here establishes that. What is established is that patients vary continuously in deviation magnitude and that this variation is structured.

Say **"axis of deviation"** or **"audiological axis"**. Use k = 2 as an operational convenience for describing the ends of it, and say that is what you are doing.

This is not a retreat. A dimensional claim is *harder to refute and easier to use*: it survives the silhouette question, the "is n = 58 a cluster or a tail?" question, and the replication question, and it makes a sharper prediction — that deviation magnitude should correlate continuously with hearing, which is testable in any cohort with an audiogram.

## 3. The dominant axis of heterogeneity in tinnitus neuroimaging is audiological, not symptomatic

Both projects find this, independently, and it is the most interesting scientific claim in the talk.

- **MRI**: the group that deviates more has PTA 27.7 vs 18.1 dB, p_FDR = 0.0001. And the relationship is graded, not binary — the proportion of patients in the deviant group rises **8% → 16% → 28% → 34%** across PTA quartiles, monotonically. (Computed while building the slide 24 figure; it was not in any previous document.)
- **EEG**: the classifier's confidence tracks HF-PTA (r = 0.249, p < 0.0001) and PTA4 (r = 0.198) but is flat against THI (r = −0.015) and TFI (r = 0.017). The proportion detected rises 63% → 85% across HF-PTA quartiles.

That the dose-response appears in both, on the same axis, in cohorts that share no participants, is the strongest single piece of evidence in this section — and it is a *gradient* in both, which is the other half of the argument.

The implication is uncomfortable and worth saying out loud: **the neural signal we can measure tracks the audiological substrate, not the thing patients actually complain about.** Two patients with identical THI can sit at opposite ends of this axis. If you are looking for a biomarker of tinnitus *distress*, neither of these modalities is currently delivering one.

That is a genuine contribution to the TRI conversation, because the field's implicit model is that a stronger percept produces a bigger brain signal. Neither dataset supports that.

## 4. Why controlling for PTA is not enough — and why that is a solvable problem, not a mystery

The obvious objection to slide 22 is: you removed hearing twice and it came back, so either your removal failed or your clusters are hearing.

There is a better explanation, and the EEG work points straight at it. **PTA is a 0.5–4 kHz average. It cannot see above 4 kHz at all.** The EEG effect localises specifically to **HF-PTA**, and is essentially absent for PTA4 as a discriminator once HF-PTA is in the model. If the mechanism lives in the extended high frequencies — which is where noise-induced and age-related cochlear damage appears first, and where tinnitus pitch usually matches — then residualising PTA4 removes the part of hearing that does not matter and leaves the part that does.

This reframes "hearing keeps coming back" from an embarrassment into a **hypothesis with a specific test**: derive HF-PTA (4 + 8 kHz for ~257 patients, extended frequencies for a smaller subset) and ask whether the axis is really an extended-high-frequency axis. If it is, that is a better result than the one currently on the slide.

## 5. Age and hearing are collinear, and the 8.5-year gap is doing work

Biotype 0's median age is 68; Biotype 1's is 59.5. Age is in the normative model, so the Z-scores are age-adjusted — yet an 8.5-year difference survives. Either the age–volume relationship is not fully captured by a B-spline mean function, or age is carrying something the model does not have (cumulative noise exposure, comorbidity, duration).

Two consequences. First, "older, more hearing loss, more deviation" is a coherent single story, and you cannot currently separate its three components — say so. Second, the THI difference is exposed: older patients with longer-standing tinnitus habituate, and report lower handicap. Age alone could produce the THI gap without any biology.

What would break the tie: model age and hearing jointly against the mixture posterior; include tinnitus duration once it is merged; and check whether deviation magnitude still tracks hearing *within* an age stratum.

## 6. Where the two projects disagree — and why the disagreement is informative

They agree on the hearing axis and on its continuity. They disagree on distress.

- **EEG**: the two groups are indistinguishable on THI and TFI (d = 0.10 and 0.01).
- **MRI**: Biotype 0 has median THI 14, Biotype 1 has 28 (p_FDR = 0.023).

The likeliest explanation is not biological. Tinception's THI is **massively confounded with site** — means run from 17.2 at cogtail to 58.9 at triple — and **cogtail supplies 43% of Biotype 0**. The MRI "distress" result may simply be the cogtail effect. TIDE collected THI/TFI uniformly across all seven sites, so its null is the cleaner measurement.

The honest reading is that **the EEG result is the better-measured one, and the MRI distress finding should be treated as unproven until it is redone with site in the model.** Saying that at TRI costs you nothing and buys a lot of credibility, especially with anyone who has tried to pool THI across sites themselves.

## 7. Why the whole-brain contrasts are null, and why that was predictable

Expect the question and do not be defensive about it. Three reasons, in order of importance:

1. The axis was **defined on 11 subcortical nuclei**. There is no reason a voxelwise gray-matter contrast should recover it.
2. **Power**: 58 vs 212 with d ≈ 0.4 in structures of a few hundred voxels, under FWE correction across the whole brain. This is not a test that was ever going to succeed.
3. VBM at 2 mm cannot resolve the amygdalar and thalamic nuclei the axis is built from — that is precisely why the nuclei-level segmentation exists.

The correct statement is: **the whole-brain contrast is uninformative here, not evidence of absence.** And note the provenance issue — the null maps on disk were computed from the *old* labels; the current-label contrast has never been run.

## 8. What would make this a real subtype claim

Worth having as a slide box and as an answer. In rough order of how convincing each would be:

- **Membership predicts something the model never saw** — treatment response, EEG signature, longitudinal trajectory. This is the only one that would settle it.
- **The split survives in UK Biobank tinnitus cases but not in controls.** Right now it does neither: controls split 23.8/76.2 and cases 26.0/74.0.
- **A Hartigan dip test rejects unimodality** on the mixture posterior or on PC1.
- **Bootstrap ARI holds up when the PCA and scaler are refitted inside each resample** — the current procedure fits them once on the full sample, which inflates stability.
- **The axis reproduces on an independent feature set** — cortical thickness, or the EEG deviation profile in the same people.
- **It survives hearing being modelled properly**, i.e. HF-PTA rather than PTA4, and it survives PTA-matched pairs (the EEG project's SA2 design).

Until at least two of those hold, "we find a continuous axis of deviation organised by hearing" is the claim the data support, and it is a perfectly good one.

---

# Discussion — prepared talking points

**1. Small effects are a finding, not a failure.**
Cohen's d of 0.1–0.35 means the distributions overlap almost completely. No structural measure here is remotely diagnostic at the individual level. That argues against a focal structural lesion and for tinnitus as a functional/network disorder with modest structural correlates — or for genuine heterogeneity diluting a larger effect present in a subgroup. The VA thalamus result (9.1% vs 2.1% extreme deviation) is the one piece of direct evidence for the second reading.

**2. Heterogeneity as the explanation for a non-converging literature.**
If patients vary continuously in how far their brains sit from normative expectation, a study of N=30 samples that distribution unpredictably — and if the variation is organised by hearing, then a study's result depends on the audiometric profile of whoever happened to walk into that clinic. That predicts exactly the pattern the field shows: many findings, low replication, no consistent direction, and — as Makani's meta-analysis found — sign reversals between normal-hearing and hearing-loss strata. This is a testable version of the heterogeneity claim rather than a rhetorical one, and the axis framing makes it sharper than a subtype framing would: it predicts a *continuous* relationship between cohort hearing profile and effect size, not a mixture of two study types.

**3. Why larger volumes?**
Be honest that you don't know. Four candidates worth naming: activity-dependent plasticity from chronic aberrant input; neuroinflammatory/glial changes (which increase apparent volume on T1); selection effects in who seeks treatment and enrols; and reverse causation — pre-existing structure as a risk factor. Cross-sectional data cannot separate these, and saying so is stronger than picking one.

**4. What each cohort buys.**
The clinical cohort has audiograms, THI, and diagnostic certainty; UKB has power and independence. The findings that appear in both (putamen, cerebellum, direction of effect) are the trustworthy ones. The findings that appear only in the clinical cohort (the STS cluster, VA thalamus, THI relationships) are either clinical-severity-dependent — plausible, since UKB captures milder self-reported tinnitus — or false positives. You cannot currently distinguish these, and the new datasets are the way to.

**5. On the rigor point.**
The 44-vs-3 gap is generalisable beyond tinnitus. Treating a nominal p<0.05 in a 285-vs-285 sample as equal evidence to p<0.05 in a 2,180-vs-2,176 sample is a category error, and a hierarchical model makes that visible. This is a methods contribution worth claiming.

**6. What the new data must do.**
Not simply add N to the same group contrast — that only shrinks confidence intervals around effects you already know are small. It should be aimed at the axis: enough patients with **extended-high-frequency audiometry** to test whether this is really an EHF axis, enough with uniformly-collected distress instruments to test the THI question without site confounding it, and ideally enough with a second modality in the same people to ask whether structural deviation and EEG deviation are the same axis or two.

**7. A methods point worth making explicitly.**
Two of our projects independently produced a beautiful-looking cluster solution that turned out to be a single site — silhouette 0.52 in the MRI cohort, 0.415 in the EEG cohort. In both cases the artefactual solution had *better* internal validity metrics than the real one. Clustering quality indices are not a defence against site effects in multi-site data; they reward them. Leave-one-site-out should be routine in this literature and currently is not.

---

# Anticipated Q&A

**Q: Isn't this all just hearing loss?**
Groups were matched on PTA at recruitment, not just covaried post hoc. The putamen effect survives adding high-frequency PTA as an extra covariate; the cerebellar vermis effect does not — which shows the covariate is doing real work. And the biotypes were derived with PTA partialled out before clustering, yet still separate on PTA afterwards.

**Q: Why the putamen? Tinnitus is an auditory disorder.**
Two reasons it's less surprising than it looks. The putamen sits in a limbic–striatal gating circuit repeatedly implicated in tinnitus models — the "noise-cancellation" framework in particular. And our normative analysis independently flags the VA thalamic nucleus, the primary thalamic output target of the putaminal/pallidal circuit. Two methods, one circuit. That said, we did not find primary auditory cortex, and I won't pretend otherwise.

**Q: Why doesn't the cortical thickness finding replicate in UKB?**
Three candidates, and I can't rule out the third. Phenotype: UKB is self-report, which attenuates small effects. Acquisition: thickness estimation in sulcal regions like STS is sensitive to sequence differences. Or it's a false positive — CWP = 0.047 in a single dataset without replication has a non-trivial false-positive probability, and that's why I labelled it hypothesis-generating.

**Q: Your subtypes — did you just recover the audiogram?**
Partly, and I think that's the interesting part rather than the objection. Hearing was a covariate inside the normative model and was regressed out again before clustering, and the groups still differ on PTA at p = 0.0001. My reading is that PTA is the wrong hearing variable: it's a 0.5-to-4-kilohertz average and it can't see above 4 kHz at all. In our EEG cohort the same axis localises specifically to high-frequency PTA, and there we could match patients exactly on high-frequency hearing and the neural difference got *larger*, not smaller. So the honest answer is: this is an audiological axis, we haven't yet measured hearing well enough to say which part of it, and that's the next analysis.

**Q: But you used to say the less-hearing-loss group was the more distressed one.**
I did, and I don't think that claim survives scrutiny. THI in our cohort is heavily confounded with site — site means run from 17 to 59 — and one site supplies 43% of the smaller group and has the lowest THI of any site. The test was unadjusted, on 43 versus 124 patients, with an 8.5-year age difference on top. Until it's redone with site in the model I'd treat it as unproven. Our EEG cohort, which collected distress instruments uniformly across all seven sites, finds no distress difference at all — and I'd trust that measurement more.

**Q: Two clusters with a silhouette of 0.12 — is that real structure?**
I'd say what we have is an axis, not two clusters. The silhouette is 0.12, BIC prefers two components over one by only about 34 units, a third of resamples give a substantially different partition, and we never ran a unimodality test. What the data support is that patients differ continuously in deviation magnitude and that the variation is organised by hearing. Two components is where BIC lands once you forbid tiny clusters — I use it to describe the ends of the axis, not to claim two kinds of patient. And I'd note that the cleanest-looking solution we found, silhouette 0.52, was a single scanner.

**Q: Is n = 58 out of 270 a cluster, or just the tail of one distribution?**
On current evidence I can't distinguish those, and I'd lean toward the tail. The two groups don't have different profiles — one sits further from normative expectation than the other, in the same direction, on ten of eleven features. That's what a tail looks like. To call it a cluster I'd want a dip test rejecting unimodality, or membership predicting something the model never saw.

**Q: Do the subtypes differ in brain structure at the whole-brain level?**
The whole-brain contrasts on disk are null, but I'd call them uninformative rather than negative. The axis was defined on eleven subcortical nuclei; 58 versus 212 at d around 0.4, under whole-brain FWE, in structures VBM can't resolve at 2 mm — that test was never going to succeed. I should also be precise: those maps were computed from an earlier clustering solution we discarded. The contrast under the current labels hasn't been run yet.

**Q: Does the UK Biobank replicate the subtypes?**
Not as currently computed, and I'd rather tell you that than not. When you apply the model to UK Biobank properly, controls split 24/76 and cases split 26/74 — essentially the same. Two caveats: the analysis had a bug that pooled cases and controls, and the Biobank deviation scores were never residualised for hearing the way the clinical ones were, so the model is being applied across a preprocessing mismatch. Both are fixable and both are at the top of my list. But right now I can't tell you this axis is tinnitus-specific rather than ordinary variation in subcortical morphometry.

**Q: How does this relate to your EEG work?**
Directly, and that's the reason I'm reasonably confident in the axis. Separate consortium, separate people, resting-state EEG instead of structural MRI, a supervised classifier instead of unsupervised clustering. A third of tinnitus patients there are confidently classified as controls — not borderline, median predicted probability 0.17 — and those patients have better high-frequency hearing and identical distress scores. Same axis, same direction, same continuity: their clustering silhouette is 0.07, ours is 0.12. Two cohorts and two modalities agreeing that the dominant axis of variation inside the tinnitus group is hearing rather than severity is a lot more persuasive than either result alone.

**Q: How do you handle 13 different sites and two field strengths?**
Site is in every design matrix as a fixed effect; UKB derivatives were ComBat-harmonised on controls. But covariate adjustment isn't sufficient on its own — the biotyping analysis found a site that survived adjustment and still drove the clustering. That's why we run leave-one-site-out sensitivity checks, and why I'd argue everyone doing multi-site clustering should.

**Q: Only three regions with HDIs that cross zero — is there anything here at all?**
That's the honest state of the evidence for individual regions. What's more robust than any single region is the *pattern*: consistent direction across 79 of 97 regions, consistent across two cohorts and four analysis methods, anatomically concentrated in one circuit. I'd defend the pattern before I'd defend any one region.

**Q: Multiple comparisons across the whole project?**
Each analysis is corrected within itself — TFCE, FDR, Bonferroni as appropriate. There is no correction across the number of *analyses*, which is a standard limitation of multi-modal discovery work and one I'll state in the paper. It's also the reason I'm presenting the Bayesian cross-dataset result as the headline rather than the longer nominal list.

**Q: Cause or consequence?**
Cannot be determined from cross-sectional data. Both directions are plausible for the putamen and for the CRHR1/stress-axis genetic signal. Longitudinal imaging around tinnitus onset is the study that would settle it, and it isn't this one.

**Q: What about the genetics — is it worth continuing with no hits?**
Lambda is clean, so the null is informative rather than technical: at this N, tinnitus has no common variant of large effect. The lever is phenotype refinement, not sample size alone — which is why the data-driven subtype arm is the more interesting one, and why it needs FUMA annotation next. One caveat I'd state: that GWAS used an earlier clustering solution, not the one I showed you, so it needs rerunning before I'd read anything into it.

**Q: What will the new datasets change?**
Mainly they let us test the axis rather than describe it: whether it's really an extended-high-frequency hearing axis once we measure hearing properly, whether it replicates within tinnitus cases in the Biobank, and whether position on it predicts anything clinically. Adding N to the same group contrast would mostly just tighten intervals around effects we already know are small.

**Q: If this is all hearing, is there anything tinnitus-specific here at all?**
That is the right question and I don't have a clean answer. What I can say is that in the EEG cohort, matching patients exactly on high-frequency hearing made the neural separation larger rather than smaller, and hearing contributed nothing once group membership was known — so within tinnitus patients, the axis isn't reducible to the audiogram. What we haven't yet shown is that the same axis doesn't exist in people without tinnitus, and the Biobank result currently suggests it might. That's a specificity test we owe the field, and it's the analysis I'd run first.

---

# Pre-flight checklist

## Figures that must be made (none of these exist)

| Figure | For slide | Note |
|---|---|---|
| Pipeline diagram | 8 | Replaces the ten repeated "Outlook" navigator slides in v2 |
| Heterogeneity schematic — **one distribution with a shifted tail**, not two separate bumps | 4, reused 28 | The conceptual spine of the talk. Draw it as a continuum: the old two-bump version now contradicts slides 21 and 23 |
| **Forest plot: the 11 subtype Cohen's d values with CIs** | 21 | Data confirmed available; the table has never been written to disk — only PDFs exist |
| "Evidence for discreteness" summary panel | 23 | Optional — `01_bic_selection.pdf` + `02_bootstrap_ari.pdf` will do |
| Circuit schematic (putamen → VA → amygdala/MDm) | 28 | Ties the whole argument together |
| Consortium logo strip | 1, 29 | |

Also absent from the results tree, if you ever want them: **no PRS figure exists at all**; **UKB SBM was never rendered** (only `SBM/UKB/sig_{lh,rh}.mgh`); `plots/rois_norm/ukb/` is **empty**; the Bayesian ROI analysis has **no figure of any kind**.

## Numbers to re-verify

- [ ] **Fix the UKB transfer — two defects, not one.** Add the group filter to `load_ukb_features()` (`src/main/13_biotype_gmm.py:214`), *and* create `ukb_covars.csv` so UKB Z-scores are residualised for hearing and site the way the clinical ones are. Then rerun and compare cases against controls as a specificity test. Slide 25 reports the preliminary answer (cases 26.0% vs controls 23.8%).
- [ ] **Redo the subtype–THI test with site in the model** before making any distress claim (warning item 10), and fix the merge bug that drops antinomics THI from `subcortical_roi/master.csv` (167 → 205 patients).
- [ ] **Cite San-Martín et al. 2025 and Zainul Abidin et al. 2025** — see the box below. Both are 2025 papers that directly touch your headline claims; a TRI audience will know them.
- [ ] Optionally extend `src/figures/slide02_literature.csv` with the 9 commented-out rows at its foot (direction confirmed, per-group N still unverified — mostly paywalled)
- [ ] Confirm the sex coding (0/1) in `subcortical_roi/master.csv` before showing any sex breakdown
- [ ] Decide one N per module and footnote the rest: 569 (ROI/SBM/MBM) / 572 (VBM) / 560 (SSA)
- [ ] Regenerate and archive the VBM cluster table — `vbm_results/` stores no stats text, so the numbers on slides 9–11 were re-derived with FSL for this plan
- [ ] `BIOTYPING_RESULTS.md` says 13 sites; the data have 10, and the **primary biotyping run used 9** (270 patients, london_2 excluded)
- [ ] Fix the mis-attributed 86.7% — it is 13/15 of **london_2's patients**, not 86.7% of the minority cluster (which was 13/27 = 48.1%)
- [ ] Report bootstrap ARI as **median 0.747, IQR 0.25–0.84, 67% ≥ 0.60**, not the median alone
- [ ] State that clustering runs on **6 principal components**, not the 11 raw Z-scores
- [ ] Disclose that **site remains significant** between subtypes after excluding london_2 (χ² = 22.09, p_FDR = 0.0079)
- [ ] Check whether any TIDE participant also appears in Tinception before calling the two cohorts independent on slide 24

## Analyses to run next, in value order

Not needed for the talk — this is what turns the axis from a description into a result. Several already exist as working code in the EEG project and can be ported rather than rewritten: `~/Tinnorm/src/35_fn_sensitivity_analysis.py` (the SA1–SA5 battery), `33_fn_subgroup_analysis.py` (clinical effect sizes, clustering checks), `27_deviation_counts.py` (extreme-deviation indices).

1. **Derive HF-PTA and re-run the whole analysis against it.** Per-frequency thresholds exist in the per-site `master_files/*_master.xlsx` for 9 of 10 cohort sites — 517 of 569 subjects have 0.5–8 kHz, so a 4+8 kHz HF-PTA is available for ~257 patients. Extended frequencies go to 16 kHz at tinmeg, 14 kHz at tinspect, 12.5 kHz at antinomics, 12 kHz at neuropren; tinspect's `Audiogramme.xlsx` already carries a precomputed `Mittelwert Hochton`. **This is the single highest-value test in the list** and the direct bridge to the EEG result. Note `triple` (52 subjects) has PTA only, and cogtail stores one binaural value mislabelled with an `RH-` prefix.
2. **Fix the THI merge bug, then redo the distress test with site in the model.** antinomics THI (n = 38) is in `master_files/master.csv` but lost in `subcortical_roi/master.csv` and `VBM_design/covars.csv`. That takes THI from 167 to 205 patients. Until this is done, do not claim a distress difference.
3. **Demonstrate heterogeneity without clustering.** Test whether the *spread* of normative Z is larger in patients than controls (Levene / Brown-Forsythe, or a variance ratio per feature). The EEG project reports ~3× on its top feature (control SD 0.42 vs patient SD 1.38). This makes the heterogeneity argument with no cluster solution to defend, and it sidesteps the entire "your clusters aren't real" line of attack.
4. **Formal discreteness tests, which the pipeline has never had.** Hartigan dip test on the mixture posterior and on PC1; a simulated-unimodal null for the ARI; and refit the PCA and scaler *inside* each bootstrap resample — currently they are fitted once on the full sample, which inflates stability.
5. **Port the EEG sensitivity battery to the MRI cohort.** SA2: PTA-matched subtype pairs — does the Z separation survive audiometric equating? SA3: subtype membership rate by PTA quartile — is there a dose-response? SA5: `posterior ~ subtype + PTA + age` — is the axis independent of hearing once membership is known?
6. **Fix and re-run the UKB transfer** (group filter + `ukb_covars.csv`), then use controls as the specificity test.
7. **Run the whole-brain contrast under the current labels.** It has never been done; the maps on disk are from the discarded solution.
8. **Re-point the downstream analyses at `biotypes_gmm/`** — the GWAS and `12_run_correlations.py` still consume the old UMAP/HDBSCAN labels — and commit `13_biotype_gmm.py`, which is untracked.
9. **Write the per-region subtype Z table to disk.** It exists only as PDFs, which is why slide 21 did not exist until now.
10. **Merge the clinical archive.** `/Volumes/Extreme_SSD/payam_data/Tinception_repo/metadata/` holds TFI, TQ (erlangen), MDI (triple, inditms), HADS and PHQ-15 (london, new-castle), BDI/BAI (tinspect), hyperacusis (GÜF/HQ), LDL/UCL, tinnitus duration, laterality, pitch and loudness matching, MML, residual inhibition, a full cognitive battery (cogtail, n = 144), education, handedness and hearing-aid use. **None of it was ever merged into the analysis cohort.** Most instruments span only 1–2 sites, but after harmonisation **laterality (~110 patients), duration (~90–110), sound tolerance (~75) and hearing-aid use (~90)** are within reach of a real test. This is what would let the axis be characterised clinically rather than only audiologically — and it is the cheapest large win available.
11. **A free external check**: cogtail carries its own a-priori audiological label (`CG / TI-HL / TI-HA / NT-HL / NT-HA`, n = 144). Does mixture-model membership align with it? Agreement would be meaningful; cogtail supplies 43% of Biotype 0.

## Files not to grab by accident

- `biotypes_gmm/biotype_assignments_ukb.csv` (top level) → **stale**, contradicts `ukb_transfer/`
- `biotypes_gmm/biotype_assignments_main.csv` (top level) → the **diagnostic** run (258/27), not the primary (58/212)
- `biotypes_gmm/biotype_stats_ukb.csv` and `ukb_transfer/biotype_stats_ukb.csv` → both **1-byte empty**
- `VBM_design/thi.csv` and `biotypes/` → the **rejected UMAP/HDBSCAN** labels. **Note these are still what the biotype GWAS and `12_run_correlations.py` actually consume** — nothing in the repo reads `biotypes_gmm/`, and `13_biotype_gmm.py` is untracked in git
- `vbm_results/Bio_TFCE/` → null, **and computed from the rejected labels** (dated a month before the mixture model ran). The current-label whole-brain contrast does not exist
- `biotypes_gmm/sensitivity_excl_london_2/` → byte-identical duplicate of `primary/`; harmless, but don't cite it as an independent sensitivity analysis
- `vbm_results/cluster_masks/cluster_uvula-rh_mask.nii.gz` → **0 voxels, broken**
- `phenotypes/phenocode-*.tsv` → **chromosome 1 only**, not GWAS results for those traits
- `vbm_norm` / `sbm_norm` `statistics_test.csv` → R²/SMSE/EXPV numerically exploded; Rho is fine

## ⚠️ Two 2025 papers you must position yourself against

The literature search for the slide 2 figure surfaced two recent papers that bear directly on your headline claims. Neither is in the current deck. A TRI audience will know both.

**1. San-Martín et al. 2025** — *Sci Rep*, PMID 41271804 · 51 tinnitus vs 40 controls, **PTA-matched**
Reports **increased putamen, caudate, pallidum, accumbens and a striatum composite** in tinnitus (putamen t = −3.6, p = 0.0016), with **no difference** in auditory areas, ACC or DLPFC.

- This is **independent corroboration of your headline putamen finding, in the same direction, from a hearing-matched sample.** It substantially strengthens slide 9.
- It is the **ANDES cohort** — Paul Delano, Chile — which is one of your own recruiting sites. Check whether this overlaps your Andes subjects (48 CO / 58 TI). If it does, it is not independent replication and must be described as a partner analysis of shared data. **Find this out before the talk**; getting it wrong in either direction is bad.
- Either way: cite it on slide 9 and say the direction agrees.

**2. Zainul Abidin et al. 2025** — *Brain Commun*, PMID 40496672 · **UK Biobank, 3,121 tinnitus vs 3,121 matched controls**
Reports **increased bilateral thalamus** (strongest effect), **hippocampus, left pallidum, precuneus thickness, right pars triangularis**. Amygdala, hippocampus, putamen and accumbens were significant before covariate adjustment and **lost significance after** it. Heschl's gyrus changes were attributed to **hearing difficulty, not tinnitus**.

- This is the **largest tinnitus structural study in existence** and it uses **the same resource as your validation arm**. You must position your UKB analysis relative to it on slides 7 and 11 — otherwise the first question will be "how does this differ from Zainul Abidin?"
- Points of contact with your own results: their thalamus increase aligns with your MDm/thalamic findings; their "significant until adjusted" subcortical set is a close analogue of your 44-vs-3 problem; their hearing-difficulty dissociation supports your PTA-control design.
- Points of difference to state plainly: their tinnitus phenotype is self-report only, they use FreeSurfer IDPs rather than your nuclei-level segmentation, and they have no clinical cohort and no biotyping. **Your subtype analysis and your fine-grained nuclei segmentation are what is genuinely new** — lead with that when you contrast.

## Before the deck leaves your laptop

- [ ] **Delete slides 66–67 of `wtinception_v2.pptx`** — slide 66 carries plaintext lab-equipment credentials. Rotate them if still active.
- [ ] Fix the raw-Rho bug in `11_compute_ssa_4.py` before quoting any raw correlation from `ssa_results/genes.csv` (spin-test columns are unaffected)
- [ ] Back up `/Volumes/Extreme_SSD/payam_data/Tinception/plots/` — **no figures are tracked in git**; they exist only on the SSD
- [ ] Embed fonts; the PDFs are vector and will scale, but check they render on the conference machine

---

## Figure quick-reference

All paths relative to `/Volumes/Extreme_SSD/payam_data/Tinception/`.

| Slide | Figure |
|---|---|
| 2 | `plots/talk_tri/slide02a_literature_landscape.pdf`, `plots/talk_tri/slide02b_power_curves.pdf` |
| 3 | `plots/audiometry/audio.pdf`, `plots/audiometry/sites.pdf` |
| 6 | `plots/demographics/group_distribution.pdf`, `plots/demographics/sex_distribution.pdf` |
| 7 | `plots/demographics/ukb_distributions.pdf` |
| 9 | `plots/vbm_group/ti_gt_co_{x,y,z}.pdf`, `plots/vbm_group/boxplots_ti_gt_co.pdf` |
| 10 | `plots/vbm_group/ti_gt_co_pta_hf.pdf`, `plots/vbm_group/co_gt_ti_pta_hf.pdf` |
| 11 | `plots/vbm_group/ukb/ti_gt_co_x_thr_0.005.pdf`, `.../ti_gt_co_z_thr_0.05.pdf` |
| 12 | `plots/sbm/brain.pdf`, `plots/sbm/boxplots.pdf` |
| 13 | `plots/sbm/brain_doss.pdf` |
| 14 | `plots/rois/cohen.pdf`, `plots/rois/heatmap.pdf`, `plots/rois/thalamic_nuclei_manhattan_ukb.pdf` |
| 15 | `plots/talk_tri/slide15_bayesian_forest.pdf` (+ `plots/tier_heatmap_v2.pdf` as backup) |
| 16 | `plots/rois/boxplots_Amygdala.pdf`, `plots/rois/boxplots_Thalamic-nuclei.pdf` |
| 17 | `plots/rois_norm/main/Thalamic-nuclei_extreme_{count,stats}.pdf` |
| 18 | `plots/VBM_Normative_Metrics_Final.pdf`, `plots/SBM_Normative_Metrics.pdf` |
| 19 | `plots/biotypes_gmm/primary/01_bic_selection.pdf`, `.../02_bootstrap_ari.pdf` (+ `00_pca_variance.pdf`) |
| 20 | `plots/biotypes_gmm/primary/09_site_composition.pdf` (+ `diagnostic_all_sites/02_bootstrap_ari.pdf` for the side-by-side) |
| 21 | `plots/biotypes_gmm/primary/05_zprofile_partialed.pdf`, `.../07_radar_primary_(excl._london_2).pdf` (+ forest plot to be made) |
| 22 | `plots/biotypes_gmm/primary/06_clinical_clinical.pdf` |
| 23 | reuse slide 19's two figures, or a summary panel to be made |
| 24 | `plots/talk_tri/slide24_mri_eeg_convergence.pdf` |
| 25 | `plots/biotypes_gmm/ukb_transfer/08_ukb_transfer_summary.pdf`, `.../04_umap_ukb.pdf` — **regenerate after the fix** |
| 26 | `GWAS/plots/tinnitus_subtype.manhattan.png`, `GWAS/plots/tinnitus_vs_controls.qq.png` |
| 27 | `plots/ssa/brain_diff_grad_3.pdf`, `MBM/lh_mbm_result.png` |
| backup | `plots/correlations/PTA_THI_stats.png`, `plots/QC.pdf`, `plots/biotypes_gmm/primary/05_zprofile_raw.pdf`, `.../03_dendrogram.pdf`. **Not** `putamen_vs_THI.pdf` — it uses the discarded labels |

All slide 19–25 paths above were confirmed present on the SSD while writing this plan.

---

## Appendix — the two cohorts side by side

Keep this as a backup slide for slide 24. Everything in the EEG column comes from `~/Tinnorm/{FN_SUBGROUP_ANALYSIS,DISCUSSION,ANALYSIS_RESULTS}.md`.

| | Tinception (structural MRI) | TIDE (resting-state EEG) |
|---|---|---|
| Cohort | 10 sites, 569 subjects; **270 patients** in the subtype analysis | 7 sites, 544 subjects; **268 patients** |
| Method | unsupervised — normative Z → PCA → Gaussian mixture | supervised — LOSO LightGBM on source-space PLI, AUC 0.82 |
| Normative model | PCNToolkit BLR, B-spline (deg 3, 5 knots), SinhArcSinh warp, controls only | **identical configuration** |
| Deviates from norms | Biotype 0, n = 58 (21.5%) | detected, n = 183 (68.5%) |
| Looks normative | Biotype 1, n = 212 (78.5%) | missed, n = 85 (31.5%) |
| Hearing, deviant group | PTA 27.7 dB | HF-PTA 48.8 dB |
| Hearing, normative-looking group | PTA 18.1 dB | HF-PTA 38.6 dB |
| Distress difference | THI 14 vs 28 — **probably a site effect** | THI d = 0.10, TFI d = 0.01 — **none** |
| Hyperacusis | not measured in the merged cohort | HQ d = 0.03 — none |
| Nature of the split | magnitude offset, same direction, 10/11 regions | magnitude offset, same direction, all top features |
| Discreteness | silhouette 0.121 | silhouette 0.072 |
| Cleanest solution was… | a single site (london_2), silhouette 0.52 | a single site (Austin), silhouette 0.415 |
| Direction of the group effect | larger volumes (TI > CO) | hyperconnectivity (18:1 positive deviation) |
| Top region | putamen; auditory cortex absent | right pars orbitalis; auditory cortex ranks 8th |

---

*Built from the repository source, the results archive on the external SSD, and `wtinception_v2.pptx`. Every statistic was re-derived from the raw result files rather than carried over from earlier documentation; where the two disagreed, the raw files won and the discrepancy is flagged above.*
