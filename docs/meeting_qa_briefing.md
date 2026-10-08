# Meeting briefing: anticipated questions, answers and background reading

Prepared 2026-10-08 to support the single-cell statistics and T2D literature talks.

## 0. How to use this pack, and what it can and cannot support

- Part 1 gives you the ten sentences to say. Part 2 is a numbers cheat sheet. Part 3 explains the concepts well enough to answer follow-up questions. Part 4 is a bank of about 60 questions with short answers, deeper answers, and what you must not claim. Parts 5 to 7 cover what you do not know, a reading list by question type, and a glossary.
- **Evidence labels.** **[P]** = from one of the 12 core papers, as I extracted it. **[G]** = from a GEO sample record I read directly. **[S]** = from my own toy simulation (Part 3A). **[B]** = general background from standard statistics or biology, not from the 12 papers; verify before quoting. **[I]** = my own judgement.
- **Limits.** I read the 12 papers through an extraction tool that returns passages and summaries, not every page. Check any number against the paper before you quote it in a meeting. Items I could not extract are listed in Part 5.
- If you do not know an answer, the best response is shown in Part 5. Admitting what you did not check is stronger than guessing.

### The 12 core papers (keys used below)
LT19 Luecken and Theis 2019 (tutorial). HEU23 Heumos et al. 2023 (best practices). SQU21 Squair et al. 2021. JUN22 Junttila et al. 2022. SCH21 Schmid et al. 2021 (scPower). LUE22 Luecken et al. 2022 (integration benchmark). ANT25 Antonsson and Melsted 2025. DAN22 Dann et al. 2022 (Milo). BUT21 Büttner et al. 2021 (scCODA). ELG23 Elgamal et al. 2023. LAW17 Lawlor et al. 2017. SEG16 Segerstolpe et al. 2016.

---

## 1. The ten sentences to say

1. The unit of evidence is the donor, not the cell. Cells from one person are not independent replicates. [P: SQU21, JUN22, LT19]
2. In a real null test, naive tests gave up to 40% false positives; pseudobulk and mixed models stayed near the expected 5%. [P: JUN22] (40% is the maximum reported.)
3. All six top methods in the Squair benchmark were pseudobulk. [P: SQU21]
4. Pseudobulk means summing raw counts per donor and cell type, then using a standard count model (edgeR, DESeq2, limma-voom) with covariates and FDR control. [P: HEU23, ELG23]
5. Batch should be in the model. Do not test on batch-corrected values. [P: LT19, ANT25]
6. The two integration benchmarks measure different things, so they do not contradict each other. [P: LUE22, ANT25]
7. Cell-type proportions sum to one, so they need compositional methods (scCODA, Milo) at donor level. No independent benchmark exists for them. [P: BUT21, DAN22, HEU23]
8. In our main islet dataset, chemistry, islet centre and BMI all track with disease, so they must be modelled. [G: GSE221156]
9. Studies with more donors report fewer T2D beta-cell genes (248 from 8 donors; 79 from 46). That is a reason to check how each study handled replicates, not proof of an error. [P: LAW17, ELG23]
10. The supplied papers support islet-focused gaps. They do not cover adipose, muscle or liver, so adipose gaps need more literature. [I]

---

## 2. Numbers cheat sheet (check against the source before quoting)

### Statistics papers
| Fact | Value | Source |
|---|---|---|
| Squair benchmark size | 18 gold-standard datasets, 14 methods (7 single-cell, 6 pseudobulk, 1 mixed) | SQU21 |
| Top performers | All 6 top methods were pseudobulk | SQU21 |
| In vivo check | 5 of 6 pseudobulk calls validated vs 3 of 13 for a Wilcoxon-based list | SQU21 |
| Variance underestimation | Shuffling replicate labels reduced variance estimates in 98.2% of genes in one dataset | SQU21 |
| Junttila size | 18 methods; 1,280 simulated datasets (6 to 40 samples, 100 to 2,000 cells per sample) | JUN22 |
| Junttila null test | B cells from 14 healthy donors, split at random into 2 groups, 30 times | JUN22 |
| Junttila result | Naive and latent-variable methods up to 40% false positives; pseudobulk and mixed models near 0.05; summing beat averaging | JUN22 |
| scPower example | About 74% power with 10 vs 10 donors and 3,000 cells per cell type, for large effects (84 genes, median absolute log fold change 2.8) | SCH21 |
| Integration benchmark | 16 methods, 68 method and preprocessing combinations (590 runs), 13 tasks, 85 batches, about 1.2M cells, 14 metrics | LUE22 |
| Scaling effect | Raised batch removal in 79% of comparisons, lowered biology conservation in 72% | LUE22 |
| Calibration study | 8 methods, 25 random pseudobatch splits; Harmony changed data least; MNN, scVI, LIGER altered most | ANT25 |
| Seurat and MNN artifacts | Over 800 DE genes where the original had 179 and 154 | ANT25 |
| scCODA | Average Matthews correlation 0.64 in low-sample simulations | BUT21 |

### Islet T2D studies
| Fact | Value | Source |
|---|---|---|
| Lawlor 2017 | 8 donors (5 ND, 3 T2D); 1,050 cells; Fluidigm C1; 248 beta, 138 alpha, 24 delta T2D genes | LAW17 |
| Segerstolpe 2016 | 10 donors (6 healthy, 4 T2D); 3,386 sequenced, 2,209 kept; Smart-seq2; INS and FXYD2 down in T2D beta cells | SEG16 |
| Elgamal 2023 | 67 downloaded (29 ND, 17 T2D, 10 T1D, 9 Aab+), 65 in map; 192,203 cells; T2D 84 genes (79 in beta); T1D 1,808 genes | ELG23 |
| Elgamal model | Pseudobulk DESeq2 Wald; covariates sex, age, BMI, chemistry, tissue source; BH FDR 0.10 | ELG23 |
| Elgamal biology | TSHR, SLC4A4, TNFRSF11B up; mitochondrial and oxidative-phosphorylation genes down | ELG23 |

### Dataset records (GEO, read directly)
| Accession | Fact | Source |
|---|---|---|
| GSE221156 (Bandesh, islets) | 17 T2D, 17 ND, 14 pre-diabetic donors; 54 libraries from 48 donors; 42 single-donor plus 12 pooled libraries from 6 donor pairs | G |
| GSE221156 chemistry | Single-donor libraries: T2D 7 V2 and 7 V3; ND 4 V2 and 13 V3; PD 1 V2 and 10 V3 | G |
| GSE221156 centre | UPenn: 4 libraries, 3 T2D, none ND | G |
| GSE221156 BMI | Mean 28.4 (ND), 30.8 (PD), 34.0 (T2D) | G |
| GSE342773 (visceral adipose) | 21 Korean adults: Lean 8, Obesity 5, Obesity+T2D 8; snRNA-seq 11 channels | G |
| GSE268904 (subcutaneous adipose) | 14 snRNA-seq samples: 4 lean, 7 obese, 3 obese diabetic | G |
| GSE278526 | scRNA-seq of adipose stromal cells; 58 recruited, 10 profiled; 4 pooled libraries | G |
| GSE249089 (METSIM) | 84 participants, snRNA-seq, 4 donors pooled per run; raw files withheld; T2D labels unverified | G |
| GSE154126 (Avrahami) | 22 donors; adult controls 4, T2D 10; only 11 to 54 cells per adult donor | G |

---

## 3. Concept explainers (so you can answer follow-ups)

### 3A. Pseudoreplication, with a number you can quote [S]
**The idea.** Each donor has their own baseline for each gene (genetics, age, how the islets were isolated). All cells from that donor share it. A cell-level test treats 2,000 cells from 10 donors as 2,000 independent measurements, so it mistakes donor baselines for disease effects.

**Toy simulation (my own; `scripts/pseudoreplication_sim.py`).** 10 donors per group, 200 cells per donor, negative-binomial counts, a donor-level random effect on the mean, **no true condition effect**. I tested 1,500 genes at p < 0.05:

| Donor-to-donor variation (log SD) | Donors per group | Cells per donor | Cell-level t-test | Cell-level Wilcoxon | Pseudobulk (sum, then t-test) |
|---|---|---|---|---|---|
| 0.3 | 10 | 200 | 57% | 55% | 5% |
| 0.3 | 10 | 1,000 | 81% | 79% | 5% |
| 0.1 | 10 | 200 | 16% | 16% | 5% |
| 0.3 | 3 | 200 | 58% | 57% | 3% |

How to explain it: with no real difference, the expected false-positive rate is 5%. The cell-level tests give far more, and **more cells make it worse**, because the test becomes more certain about a difference that is only donor noise. The pseudobulk rate stays near 5%. The size of the problem depends on how much donors differ (compare 0.1 and 0.3). **Caution:** this is a toy model with settings I chose. It illustrates the mechanism; it is not a measurement of any real dataset and is not from the 12 papers. The real numbers are SQU21 and JUN22.

**Follow-up: "Does more sequencing help?"** More cells per donor improves each donor's estimate but not the number of independent donors, so it does not fix pseudoreplication (see the 1,000-cell row).

### 3B. Why the negative binomial
Counts are non-negative integers. A Poisson model says variance equals the mean. Real gene counts have variance greater than the mean (overdispersion), partly because true expression differs between cells. The negative binomial adds a dispersion parameter, so variance = mean + dispersion × mean². edgeR and DESeq2 fit NB generalised linear models. **[B]** For droplet UMI data, Svensson (2020) argued that zero-inflation is usually unnecessary because most zeros are sampling zeros at low counts; check that paper before relying on it.

### 3C. How pseudobulk works, step by step
1. Pick a cell type (for example beta cells).
2. For each donor, sum raw counts over that donor's beta cells, one value per gene.
3. Filter genes with too few counts (ELG23: at least 5 counts in at least half of samples per group).
4. Build a design matrix: condition plus covariates (sex, age, BMI, chemistry, centre). **It must be full rank**: if a covariate is perfectly aligned with condition (for example all T2D donors from one centre), the model cannot separate them.
5. Fit an NB GLM (edgeR, DESeq2) or a weighted linear model (limma-voom). These tools borrow strength across genes by shrinking dispersion estimates (empirical Bayes) so that few-donor designs are stable.
6. Test with a Wald, likelihood-ratio or quasi-likelihood F-test, then correct for multiple testing within the cell type.
- **Why sum, not average?** JUN22 found summing counts beat averaging. Summing keeps the count scale that the NB model needs and weights donors by how many cells they contribute. **[I]** It does not model uncertainty from very different cell numbers per donor; filtering donors with few cells is a common safeguard.
- **Limit:** needs enough cells per donor and cell type. Avrahami donors have 11 to 54 cells each [G], so cell-type pseudobulk would be very sparse.

### 3D. False discovery rate, in plain words
A p-value of 0.05 means 5% of true nulls will look significant. Test 15,000 genes and about 750 false hits are expected. Benjamini-Hochberg controls the expected **proportion** of false discoveries among the genes you call. FDR 0.10 means about 10% of the called genes are expected to be false. ELG23 used 0.10, which is looser than the common 0.05. FDR is controlled **within each cell type** in the template. A significant gene with a tiny fold change can be real but unimportant, so report effect sizes too.

### 3E. Batch, donor and disease: why design beats software
Batch effects are technical differences between runs (chemistry, centre, day). If disease status is mixed across batches, a model can separate the two. If batch and disease are aligned, **nothing can**. GSE221156 shows partial alignment [G]: T2D single-donor libraries are 7 V2 and 7 V3, non-diabetic 4 V2 and 13 V3, so V2 is overrepresented in T2D. Islet centre is also uneven: UPenn has 4 libraries, 3 T2D and none non-diabetic. The model approach is to include chemistry and centre as covariates and accept reduced power.

**Why integrate at all?** Embeddings such as Harmony help clustering and visualisation so that the same cell type groups together across batches. They are not meant for testing.

**Why the two benchmarks look contradictory.** LUE22 asks: does the method remove real batch effects and keep biology? ANT25 asks: when there is no batch, does the method leave the data alone? scVI-type methods can do well on the first and poorly on the second.

### 3F. Compositional data
Proportions across cell types sum to 1. If beta cells fall from 45% to 29%, the other types' shares must rise even if their absolute numbers did not change. Testing each type separately would call false increases. scCODA models all types jointly (Dirichlet-multinomial) and states effects relative to a reference type. Milo avoids pre-set clusters by testing overlapping neighbourhoods on a kNN graph with an NB GLM. Both need donor-level counts. A chi-square test on cell counts is pseudoreplication. **Dissociation bias:** SEG16 found cell dissociation and sorting shifted proportions (more alpha, fewer beta than in tissue), so proportions measured in dissociated cells are biased before any statistics.

### 3G. Power, in one paragraph
Power is the chance of detecting a real effect. It depends on number of donors, cells per donor, depth and effect size. SCH21's main result: for a fixed budget, many shallow cells beat few deep ones. Its worked example needs about 10 donors per group and 3,000 cells per cell type for about 74% power, but the effect sizes were very large (median absolute log fold change 2.8) and came from a leukaemia study. Real T2D effects are much smaller (ELG23 found only 84 genes with 46 donors), so power will be lower, and donor number sets the ceiling.

### 3H. Normalisation, HVGs, PCA, UMAP in one view
Normalisation removes depth differences between cells (shifted log is the common default; scran and SCTransform are alternatives). HVG selection keeps the genes that vary most (often about 2,000; a default, not a rule). PCA compresses them to tens of components. A kNN graph is built in PCA space; clustering runs on the graph; UMAP draws the graph in two dimensions. **UMAP is a picture**: distances and blob sizes are not quantitative. **Raw counts are used for testing.** LT19 warns that log transformation can introduce false DE when size factors differ between groups.

### 3I. Clustering and annotation caution
Louvain and Leiden maximise modularity on the graph; Leiden fixes a flaw where Louvain can produce disconnected communities **[B]**. Resolution sets the number of clusters and has no single correct value. LT19: marker p-values do not validate clusters, because the clusters and the markers come from the same data (circularity). HEU23: use a three-step annotation (automated, manual markers, verification).

### 3J. Pathway analysis
ORA tests whether a thresholded gene list overlaps a gene set more than expected (hypergeometric test). GSEA ranks all genes and tests whether a set clusters at the top or bottom, using a permutation null. GSVA and ssGSEA score each sample per gene set. **No benchmark of these is in the 12 core papers.** **[I]** The safest practice is to rank genes by the donor-level pseudobulk statistic and run enrichment on that ranking.

### 3K. QC
Per sample, look at genes detected, UMI counts and mitochondrial fraction together; set permissive thresholds per sample and never tune them to make a test work (LT19). Remove doublets (Scrublet, scDblFinder) and handle ambient RNA (SoupX, CellBender). In islets, ambient insulin and glucagon can appear in other cell types **[I]**, so a lone marker can mislead.

---

## 4. Question bank

Format: **Q** the question. **Short** the answer to say first. **If pressed** the deeper answer. **Evidence** where it comes from. **Do not claim** the limit.

### A. Study design and data

**A1. Why is the sample size 20 and not 20,000?**
Short: Cells from one person share genetics, age, batch and ambient RNA, so they are not independent replicates. If pressed: See 3A; with no real effect a cell-level test flagged 57% of genes in my toy model, and more cells made it worse. Evidence: SQU21, JUN22, LT19; toy numbers [S]. Do not claim: that the toy numbers apply to a real dataset.

**A2. Does that mean more cells per donor are useless?**
Short: No, they improve each donor's estimate, but they do not add independent donors. If pressed: SCH21 says many shallow cells beat few deep ones for a fixed budget, and the donor number still sets the ceiling. Do not claim: a specific cells-per-donor target for T2D; the example used large effects.

**A3. How many donors do you need?**
Short: It depends on the effect size; I cannot give one number for T2D. If pressed: SCH21's worked example needed about 10 donors per group and 3,000 cells per type for 74% power at very large effects; T2D effects are smaller (ELG23: 84 genes at 46 donors). Do not claim: that 10 donors per group is enough for T2D.

**A4. What is the right control group?**
Short: Non-diabetic donors, ideally matched for age, sex and BMI. If pressed: In GSE221156 BMI rises with disease [G], so T2D partly overlaps with obesity. Adipose dataset GSE342773 has an Obesity versus Obesity+T2D design [G] that holds adiposity closer to constant. Do not claim: that covariate adjustment fully separates BMI from T2D.

**A5. Why not use bulk RNA-seq?**
Short: Bulk gives an average; a change in one cell type can be hidden or confused with a change in composition. If pressed: single-cell resolves types, but costs more and is noisier per cell. Evidence: background [B].

**A6. scRNA-seq or snRNA-seq?**
Short: The 12 papers do not compare them. If pressed: SEG16 reports dissociation shifted proportions. **[B]** Adipocytes are generally too large for standard scRNA-seq, so adipose studies typically use nuclei. The repo's literature search found no direct comparison for islets or adipose. Do not claim: that either is better for T2D.

**A7. Are the datasets independent?**
Short: Not necessarily. If pressed: the large HPAP and PanKbase islet sets overlap, so check donor IDs before calling something a replication. Evidence: [I], noted in the plan slides. Do not claim: that two named cohorts are independent without checking donors.

**A8. How were donors counted in Bandesh?**
Short: 48 donors from 54 libraries. If pressed: 42 single-donor libraries plus 12 libraries from 6 donor pairs; six donors appear only in pooled libraries and need demultiplexing; some donors appear both alone and in pools and must be counted once. Evidence: [G].

### B. Pipeline

**B1. What is the standard pipeline?**
Short: Count matrix, QC, normalise, variable genes, PCA, neighbours, clustering, annotation. If pressed: see slide "pipeline" and Part 3H, 3I, 3K. Evidence: LT19, HEU23.

**B2. How do you choose QC thresholds?**
Short: Per sample, from the distributions, permissively. If pressed: HEU23 uses sample-wise MAD thresholds; LT19 says do not tune thresholds to improve a result. ELG23 also removed 13,036 barcodes manually on top of Scrublet [P]. Do not claim: a numeric cut-off; they are tissue and sample specific.

**B3. What about doublets and ambient RNA?**
Short: Doublet tools simulate artificial doublets and score real cells by their neighbours; ambient tools estimate the background from empty droplets. If pressed: ELG23 used SoupX and Scrublet. Evidence: HEU23, ELG23.

**B4. Which normalisation is best?**
Short: The 12 papers do not settle it. If pressed: LT19 recommends scran then log; HEU23 shifted-log with deviance feature selection; SCTransform is not covered. Do not claim: a winner.

**B5. How many highly variable genes?**
Short: Commonly about 2,000, which is a default. If pressed: LT19 says 1,000 to 5,000. Evidence: LT19.

**B6. Is UMAP good for quantifying anything?**
Short: No. If pressed: distances and cluster sizes are not meaningful; use PCA space and models. Evidence: LT19, HEU23.

**B7. Louvain or Leiden?**
Short: Leiden is now more common (HEU23, ELG23), but no paper here compares them. If pressed: LT19 (2019) recommended Louvain. Do not claim: Leiden is proven superior in this set.

**B8. How do you choose the clustering resolution?**
Short: Try several and check against markers and stability. If pressed: ELG23 used 0.5 for its islet map. Do not claim: 0.5 is generally right.

**B9. How do you know a cluster is a real cell type?**
Short: Markers plus independent support. If pressed: LT19 says marker p-values are circular. HEU23: automated plus manual plus verification. Evidence: LT19, HEU23.

**B10. Why not trust automated annotation alone?**
Short: It depends on the reference matching your tissue and platform. If pressed: no benchmark of CellTypist, SingleR or Azimuth is in the 12 papers. Do not claim: any tool is more accurate than manual annotation.

### C. Statistics

**C1. Why pseudobulk rather than a cell-level test?**
Short: It respects the donor as the replicate. If pressed: see 3A and 3C; SQU21's top six were all pseudobulk; JUN22's null test found up to 40% false positives for naive tests versus about 5%. Evidence: SQU21, JUN22.

**C2. Is Wilcoxon always wrong?**
Short: Wilcoxon on cells is wrong for condition comparisons across donors because it ignores donors. If pressed: Wilcoxon on pseudobulk samples is not benchmarked in these papers. For marker discovery between clusters it is the common tool but its p-values are circular. Do not claim: Wilcoxon on pseudobulk is good or bad; flag it as needing more literature.

**C3. Why not a mixed model?**
Short: It is a valid alternative. If pressed: JUN22 found pseudobulk at least as good, and HEU23 lists both. Mixed models keep cells and add a donor random effect; they can be slower and unstable with few donors **[I]**. Do not claim: that one is proven better.

**C4. What does "sum, not average" mean?**
Short: Add raw counts across a donor's cells for each gene. If pressed: JUN22 found summing beat averaging. Evidence: JUN22.

**C5. DESeq2, edgeR or limma-voom?**
Short: JUN22 recommends limma, DESeq2 or edgeR with summed counts for sensitivity, and ROTS for precision. If pressed: ELG23 used DESeq2. Evidence: JUN22, ELG23.

**C6. Why FDR 0.10 in Elgamal?**
Short: It is their choice; it is looser than the common 0.05. If pressed: it trades more false positives for more discoveries given limited donors. Evidence: ELG23; the interpretation is [I]. Do not claim: that the authors justified it in a particular way.

**C7. Should you correct per cell type or globally?**
Short: The ELG23 template corrects within each cell type. If pressed: global correction across all cell types would be stricter. The papers do not compare the two. Do not claim: one is correct.

**C8. What covariates should be included?**
Short: Sex, age, BMI, chemistry and centre or tissue source (ELG23 template). If pressed: BMI is correlated with T2D so the T2D coefficient is conditional on BMI. Evidence: ELG23; GSE221156 [G].

**C9. What if batch and disease are completely aligned?**
Short: Then nothing can separate them. If pressed: LUE22 and ANT25 imply that corrections cannot rescue a confounded design; fix it with design or treat the result as unconfirmed. Do not claim: any method resolves complete confounding.

**C10. Why not test on batch-corrected data?**
Short: Correction can create artifacts. If pressed: ANT25 found MNN, scVI and LIGER altered null data; Seurat and MNN reported over 800 DE genes where the original had 179 and 154. LT19 says test measured data. Evidence: ANT25, LT19.

**C11. How do you know your pipeline is calibrated?**
Short: Split controls at random into two groups and count significant genes. If pressed: JUN22 used this with 14 healthy donors; ANT25 with pseudobatches. If you permute, permute donor labels, not cell labels. Evidence: JUN22, ANT25.

**C12. What is the difference between statistical and biological significance?**
Short: A small p-value can accompany a tiny effect. If pressed: report effect sizes with FDR; shrinkage (for example apeglm) stabilises fold changes for low counts **[B]**.

**C13. What about the effect sizes in the islet papers?**
Short: I did not extract them. If pressed: check the supplementary tables of ELG23 before quoting fold changes. Do not claim: any fold-change number.

**C14. How reliable is 74% power?**
Short: It applies to the example's large effects. If pressed: median absolute log fold change was 2.8, from leukaemia data; priors are reliable up to about the pilot's sample size. Evidence: SCH21.

**C15. Does 40% mean the false-positive rate is typically 40%?**
Short: No, it is the maximum reported. If pressed: JUN22 says up to 40% for naive and latent-variable methods. Evidence: JUN22.

**C16. Is pseudobulk reliable with only 3 T2D donors?**
Short: Treat results as exploratory. If pressed: LAW17 had 3 T2D donors; GSE268904 has 3 diabetic samples [G]. In my toy model the pseudobulk test stayed at or below 5% false positives with 3 donors per group, but with very low power **[S]**. Do not claim: confident discoveries from 3 donors.

### D. Batch and integration

**D1. Harmony or scVI?**
Short: Depends on purpose. If pressed: LUE22: scVI-type methods lead on complex tasks; ANT25: Harmony changed null data least and is their only recommendation. For testing, put batch in the model. Evidence: LUE22, ANT25.

**D2. What does "poorly calibrated" mean?**
Short: The method changes data or finds differences when no real batch exists. If pressed: ANT25 split identical data into fake batches and counted differences. Evidence: ANT25.

**D3. Why did scaling matter?**
Short: Scaling raised batch removal in 79% of comparisons but lowered biology conservation in 72%. Evidence: LUE22.

**D4. How do you evaluate integration?**
Short: Batch mixing plus biology conservation plus null calibration. If pressed: LUE22 used 14 metrics. **[B]** Common metrics are kBET, LISI and silhouette scores. Do not claim: that good mixing alone means success.

**D5. Is integration needed for pseudobulk?**
Short: Not for the test: include batch as a covariate. If pressed: integration can still help clustering. Evidence: ELG23 used Harmony for the map but covariates in the DE model.

### E. Composition

**E1. Why not just compare percentages?**
Short: Percentages sum to 100, so one change shifts the rest. If pressed: see 3F; in the illustration, only beta cells fall yet alpha and other shares rise. Evidence: BUT21; the illustration uses made-up numbers.

**E2. scCODA or Milo?**
Short: scCODA models predefined types with a reference; Milo tests graph neighbourhoods. If pressed: scCODA does not model donor-to-donor response variability; Milo needs replicates and avoids complete batch-condition confounding. Evidence: BUT21, DAN22.

**E3. Are these validated for T2D?**
Short: No. If pressed: HEU23 says no independent compositional benchmarks exist; BUT21 and DAN22 used other systems. Do not claim: validation in T2D.

**E4. Can dissociation affect proportions?**
Short: Yes. If pressed: SEG16 found more alpha and fewer beta cells than in tissue. Evidence: SEG16.

### F. T2D studies and biology

**F1. What did the islet studies find?**
Short: ELG23: 84 T2D genes (79 in beta cells); TSHR, SLC4A4, TNFRSF11B up and mitochondrial and oxidative-phosphorylation genes down. LAW17: 248 beta-cell genes. SEG16: INS and FXYD2 down in beta cells. Evidence: ELG23, LAW17, SEG16.

**F2. Why do LAW17 and ELG23 differ so much?**
Short: Possible pseudoreplication, but also different platforms, thresholds and tissue handling. If pressed: I could not confirm how LAW17 handled donors. Do not claim: that LAW17 treated cells as independent.

**F3. Do the studies agree on specific genes?**
Short: I did not extract full gene lists, so I cannot say. If pressed: overlap between LAW17, SEG16 and ELG23 lists needs additional work. Do not claim overlap or conflict.

**F4. What is the biology of T2D beta cells?**
Short: **[B]** Beta-cell dysfunction and insulin resistance are central to T2D; ELG23's mitochondrial and oxidative-phosphorylation decrease is consistent with metabolic stress. If pressed: do not over-interpret a gene list without replication. Do not claim: a mechanism.

**F5. What about the data for adipose tissue?**
Short: No adipose paper is among the 12; I only read the GEO records. If pressed: GSE342773 has 8 Obesity+T2D donors, GSE268904 has 3, GSE278526 profiled 10, GSE249089 (METSIM, n = 84) has unverified T2D labels [G]. Do not claim: adipose findings.

**F6. Why does BMI matter?**
Short: BMI rises with disease in GSE221156 (28.4, 30.8, 34.0) [G], so T2D effects overlap with obesity effects. Evidence: [G]; SEG16 also correlated BMI with expression.

**F7. Are there sex differences?**
Short: Not analysed in the extracted text beyond sex as a covariate (ELG23) or blocking factor (LAW17). If pressed: GSE221156 has 5/17 female ND and 6/14 female T2D [G]. Do not claim: any sex result.

**F8. What about treatment effects (metformin, GLP-1)?**
Short: Not reported in the supplied papers. Do not claim anything.

### G. Gaps and future work

**G1. What is the strongest gap?**
Short: Independent donor-level replication of islet T2D signatures. If pressed: plus separating obesity from T2D, and T2D-specific calibration. Evidence: synthesis.

**G2. Why is adipose not listed as a confirmed gap?**
Short: No supplied adipose paper, so I cannot call it a biological gap. If pressed: only dataset records show small T2D arms. Evidence: [I], [G].

**G3. What would your first study be?**
Short: Donor-level replication of the islet T2D beta-cell signature using GSE221156 and HPAP after checking donor overlap. If pressed: it uses existing data and methods supported by the papers. Evidence: synthesis.

**G4. What is novel about the proposed framework?**
Short: Built-in null-split calibration and explicit confounding handling in a T2D cohort. Do not claim: that no one has done this; I did not systematically search for prior T2D calibration studies.

### H. Challenging questions

**H1. "Are your numbers verified?"**
Short: They come from extracted passages and GEO records; check before citing. If pressed: I did not read every page of each paper. Say which items are unverified (Part 5).

**H2. "Your simulation is made up."**
Short: Yes, it is a toy model to illustrate the mechanism; the real evidence is SQU21 and JUN22. If pressed: the script is `scripts/pseudoreplication_sim.py` and its settings can be changed.

**H3. "Isn't a mixed model the right thing?"**
Short: Both are valid; JUN22 found pseudobulk at least as good. See C3.

**H4. "Did you validate the claims about method X?"**
Short: Only claims from the 12 papers are supported by them. Anything else is flagged as background or needing more literature.

**H5. "Why not use the newest tool?"**
Short: ANT25 found newer integration tools less well calibrated than Harmony. Newer is not automatically better.

**H6. "Can you replicate the Elgamal result?"**
Short: Not yet. The data downloads are blocked in this environment, and HPAP needs access terms. Say so plainly.

**H7. "So is T2D a cell-type-specific disease?"**
Short: In ELG23, 79 of the 84 T2D genes were in beta cells, which suggests the largest changes are in beta cells; other cell types had few calls. Do not claim: other cell types are unaffected; power differs by cell type.

---

## 5. What I do not know, and how to say it

| Not established | What to say |
|---|---|
| How LAW17 and SEG16 handled donors | "I could not confirm that from the methods text I extracted; I would check the supplement." |
| SEG16's DE test | "Reported adjusted p of 0.01 or lower; I did not extract the test." |
| Effect sizes in ELG23, LAW17, SEG16 | "I did not extract fold changes; check the supplementary tables." |
| Software used in each study (Seurat vs Scanpy) | "Not extracted." |
| Overlap of gene lists across studies | "Not assessed." |
| Comparisons of Louvain/Leiden, SCTransform, annotation tools, enrichment methods, scRNA vs snRNA | "These are not covered by the 12 papers; they need additional literature." |
| Adipose T2D biology | "No adipose paper was in the set; I only have GEO records." |
| T2D labels in METSIM (GSE249089) | "Unverified." |
| E-MTAB-5061, E-MTAB-15009, HPAP/PANC-DB, PanKbase, CELLxGENE, HRA002549 | "Hosts were blocked; not verified." |
| Data download | "Blocked by the environment's network settings." |

---

## 6. Reading list by question type (start with the papers named)

| If the question is about | Read first | Then |
|---|---|---|
| Replicates, false positives | JUN22 (Methods: null test; Results: false-positive rates), SQU21 (gold standard design) | LT19 section on DE |
| How to run DE | ELG23 (Methods), HEU23 (DE section) | DESeq2, edgeR, limma-voom papers |
| Batch | ANT25 (design, results), LUE22 | LT19 batch section |
| Composition | BUT21 (model, benchmark), DAN22 | HEU23 abundance section |
| Power | SCH21 | JUN22 simulation design |
| Pipeline steps | LT19 then HEU23 | Scanpy and OSCA tutorials |
| T2D islets | ELG23, then LAW17 and SEG16 | Bandesh (GSE221156; bioRxiv doi:10.1101/2025.01.17.633590) |
| Adipose | None in the set; start with the repo reading lists (`docs/more_sources.md`, `docs/literature_gaps_filled.md`) | Dataset GEO records |
| Pathways | decoupleR, GSEA and GSVA original papers **[B]** | Enrichment on ranked donor-level statistics |

Companion documents in this repo: `docs/paper_findings_and_roadmap.md` (paper-by-paper outline), `docs/literature_review_synthesis.md` (critical synthesis), `docs/geo_verified_metadata.md` (dataset facts), `docs/mentor_plan_part2.md` (glossary and plan).

---

## 7. Glossary (one line each)

- **Count matrix:** genes by cells, integer UMI counts. **UMI:** molecule tag that lets you count molecules. **Sparse:** mostly zeros.
- **Overdispersion:** variance greater than the mean. **Negative binomial:** count model with a dispersion parameter.
- **Library size / size factor:** per-cell depth scaling. **Normalisation:** remove depth effects.
- **HVG:** highly variable gene. **PCA:** linear dimensionality reduction. **UMAP:** nonlinear plot, not a measurement.
- **kNN graph:** cells linked to nearest neighbours. **Leiden/Louvain:** community detection. **Resolution:** controls cluster number.
- **Doublet:** two cells read as one. **Ambient RNA:** background RNA in every droplet.
- **Batch effect:** technical variation between runs. **Confounding:** batch or covariate aligned with disease.
- **Pseudoreplication:** treating correlated observations as independent. **Pseudobulk:** summing a donor's cells per cell type.
- **Design matrix:** the table of condition and covariates a model uses. **Full rank:** every column adds independent information.
- **FDR / Benjamini-Hochberg:** expected proportion of false discoveries among calls. **Effect size:** size of the change, for example log fold change.
- **Compositional data:** proportions summing to 1. **Reference cell type:** baseline used by scCODA.
- **Power:** chance of detecting a real effect. **Null split / A/A test:** split identical groups to check for false positives.
- **Calibration:** whether stated error rates match actual ones. **Ground truth:** known correct answer in a benchmark.
