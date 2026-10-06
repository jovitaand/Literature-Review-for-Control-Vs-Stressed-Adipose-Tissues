# Literature review synthesis: single-cell transcriptomics for Control versus Type 2 Diabetes (T2D)

Prepared 2026-10-06 for a PhD-level literature review or methods background section.

## 0. Read this first: scope and honesty statement

1. **No paper files were attached to the request that produced this document.** I based it on the 12 core papers I read earlier in this project through a full-text extraction tool (summaries and quoted passages, not every page), recorded in `paper_findings_and_roadmap.md`, plus dataset records read directly from GEO (`geo_verified_metadata.md`, `data/metadata/`). Anything I did not extract is written "Not reported" or "Not extracted". "Not extracted" means I did not obtain it from the passages I read. It does not mean the paper omits it.
2. **The supplied literature is strong on statistics and islets, and thin on everything else.** It contains no adipose, muscle or liver T2D paper, no pathway-analysis benchmark, no cell-cell communication or trajectory paper, no Seurat versus Scanpy comparison, and no scRNA-seq versus snRNA-seq comparison. Where a requested comparison cannot be settled from these papers I say so and mark it **[needs additional literature]**. Where I add general expert background it is marked **[background, not from the supplied papers]**.
3. **Labels used for limitations.** **[A]** = reported by the authors, as I extracted it (verify against the paper). **[I]** = my independent assessment.
4. **Numbers.** Every number below was reported in the extracted text or in GEO records. Check each against the paper before you cite it.
5. If you want true PDF-level extraction, send the PDFs (or place them in the repo) and I will fill every "Not extracted" cell.

### Citation key (the 12 core papers)

| Key | Paper |
|---|---|
| LT19 | Luecken and Theis, Mol Syst Biol 2019 (tutorial) |
| HEU23 | Heumos et al., Nat Rev Genet 2023 (best practices across modalities) |
| SQU21 | Squair et al., Nat Commun 2021 (false discoveries in single-cell DE) |
| JUN22 | Junttila, Smolander, Elo, Brief Bioinform 2022 (benchmark of differential state methods) |
| SCH21 | Schmid et al., Nat Commun 2021 (scPower) |
| LUE22 | Luecken et al., Nat Methods 2022 (integration benchmark) |
| ANT25 | Antonsson and Melsted, Genome Res 2025 (calibration of batch correction) |
| DAN22 | Dann et al., Nat Biotechnol 2022 (Milo) |
| BUT21 | Büttner et al., Nat Commun 2021 (scCODA) |
| ELG23 | Elgamal et al., Diabetes 2023 (HPAP integrated islet map) |
| LAW17 | Lawlor et al., Genome Res 2017 (islets, T2D) |
| SEG16 | Segerstolpe et al., Cell Metab 2016 (islets, T2D) |

---

## 1. Paper-by-paper extraction

### 1A. The three primary-data T2D islet studies

| Field | ELG23 | LAW17 | SEG16 |
|---|---|---|---|
| Research question | Build an integrated single-cell map of human islets from HPAP donors and find disease-associated gene changes (T1D, T2D) | Define human islet cell-type signatures and cell-type-specific expression changes in T2D | Transcriptomes of human islet cells in health and T2D, including rare cell types |
| Year | 2023 | 2017 | 2016 |
| Biological system / tissue | Human pancreatic islets | Human pancreatic islets | Human pancreatic islets |
| Organism | Human | Human | Human |
| Condition | Non-diabetic, T2D, T1D, autoantibody-positive | Non-diabetic, T2D | Healthy, T2D |
| Biological samples (donors) | 67 downloaded (29 ND, 17 T2D, 10 T1D, 9 Aab+); 65 in final map | 8 (5 ND, 3 T2D) | 10 (6 healthy, 4 T2D) |
| Cells / nuclei | 192,203 cells | 1,050 cells (622 ND, 428 T2D) | 3,386 sequenced, 2,209 kept |
| scRNA or snRNA | scRNA-seq | scRNA-seq | scRNA-seq |
| Platform | HPAP data (raw FASTQ from PANC-DB); chemistry varies and is modelled; 10x versus other not confirmed | Fluidigm C1, about 3M reads per cell | Smart-seq2, about 750k reads per cell |
| Accession | PANC-DB portal (access terms apply) | SRA SRP075970 / PRJNA323853; processed GEO GSE86473 (SuperSeries) | ArrayExpress E-MTAB-5061 (unverified); portal sandberg.cmb.ki.se/pancreas |
| Main workflow | SoupX; Scrublet plus manual removal; Harmony (donor, chemistry, tissue source); Leiden at resolution 0.5; 10 cell types | Not extracted beyond differential expression below | Not extracted beyond the findings below |
| Statistics | Pseudobulk per donor and cell type; DESeq2 Wald; covariates sex, scaled age, scaled BMI, chemistry, tissue source; gene filter (at least half of samples per group with 5 counts); Benjamini-Hochberg FDR 0.10 | edgeR with sex as a blocking factor; FDR below 5%; genes filtered by number of cells. Donor handling Not extracted | Adjusted p of 0.01 or lower reported; test Not extracted |
| Main findings | T1D 1,808 genes (1,305 in beta); T2D 84 genes (79 in beta). In T2D beta cells TSHR, SLC4A4, TNFRSF11B up; mitochondrial and oxidative-phosphorylation genes down | Signatures for alpha, beta, delta, PP. T2D versus control: 248 genes in beta, 138 alpha, 24 delta (plus 74 acinar, 35 ductal, 28 stellate). 536 islet eQTL target genes checked, 263 cell-type-specific | Programs for rare delta, gamma, epsilon, stellate cells; subpopulations of alpha, beta, acinar. In T2D: INS and FXYD2 down in beta, WFS1 down in alpha. BMI correlated with obesity and diabetes gene expression |
| Main methodological contribution | A donor-level, covariate-adjusted template on the largest islet map | Early cell-type-resolved T2D signatures and eQTL link | Early rare-cell and T2D map; showed dissociation and FACS shift proportions |

### 1B. Methods, benchmark and review papers (no new T2D data)

| Paper | Year | Research question | Design and scale | Main contribution |
|---|---|---|---|---|
| LT19 | 2019 | Which tools and steps for scRNA-seq? | Tutorial; no accession | Per-step best practice: QC with genes, counts and mitochondrial fraction together and permissive; scran then log; Louvain on kNN graph; do not validate clusters with marker p-values; do not test on corrected data |
| HEU23 | 2023 | Updated best practice across modalities | Review using independent benchmarks | Sample-wise QC (MAD), SoupX/CellBender, scDblFinder, shifted-log, deviance HVG, Leiden, three-step annotation, pseudobulk or mixed model for DE, scCODA/Milo for abundance. Flags unresolved issues: no independent compositional benchmarks, low DE tool consensus |
| SQU21 | 2021 | Why do single-cell DE tests disagree, and which are right? | 18 gold-standard datasets (matched bulk and single-cell from the same purified cells); 14 methods (7 single-cell, 6 pseudobulk, 1 mixed) | All six top methods were pseudobulk. Single-cell tests found hundreds of DE genes with no perturbation, biased to high expression. In one in vivo check 5 of 6 pseudobulk calls validated versus 3 of 13 for a Wilcoxon list |
| JUN22 | 2022 | Which methods control false positives across subjects? | 18 methods on 1,280 simulated datasets (6 to 40 samples, 100 to 2,000 cells per sample); real null test: B cells from 14 healthy donors split at random into 2 groups, 30 times | Naive and latent-variable methods up to 40% false positives; pseudobulk and mixed models near 0.05. Summing beat averaging. Pseudobulk at least as good as mixed models |
| SCH21 | 2021 | How many donors, cells and reads are needed? | Power model (scPower); worked example with large effects (84 genes, median absolute log fold change 2.8) | Shallow sequencing of many cells beats deep sequencing of fewer for a fixed budget. About 74% power needed 3,000 cells per cell type per person and 10 donors per group in the example |
| LUE22 | 2022 | Which integration methods work? | 16 methods, 68 method and preprocessing combinations (590 runs), 13 tasks, 85 batches, about 1.2M cells, 14 metrics | scANVI, Scanorama, scVI, scGen good on complex tasks; Harmony and Seurat v3 on simple ones. Trade-off between batch removal and biological conservation |
| ANT25 | 2025 | Do correction methods leave data alone when no batch effect exists? | Public datasets randomly split into pseudobatches 25 times; 8 methods | MNN, scVI, LIGER often altered data substantially; ComBat, ComBat-seq, BBKNN, Seurat introduced artifacts; Harmony changed it least and is the only method they recommend |
| DAN22 | 2022 | Abundance testing that does not depend on cluster boundaries | Milo, simulations; ageing mouse thymus; cirrhotic human liver | kNN-graph neighbourhoods, negative-binomial GLM on per-sample counts, spatial FDR |
| BUT21 | 2021 | Abundance testing for compositional data | scCODA, simulations; supercentenarian B-cell decrease re-found | Hierarchical Dirichlet-multinomial with reference cell type; of compared methods only scCODA, ALDEx2 and ALR-transformed tests controlled FDR in all scenarios |

Public accessions: not applicable to most methods papers (simulations or reused public data); for these "Not extracted" individually.

### 1C. Dataset records read directly from GEO (papers themselves not read)

These are relevant because the supplied papers contain **no** adipose or Bandesh data, so the records are my only evidence about those cohorts.

| Accession | Tissue and design (GEO record) | T2D-relevant size |
|---|---|---|
| GSE221156 (Bandesh) | Islets, 10x scRNA-seq; 17 T2D, 17 ND, 14 pre-diabetic donors; 54 libraries from 48 donors; V2 and V3 chemistry | 17 T2D donors |
| GSE342773 | Visceral adipose, 21 Korean adults: Lean 8, Obesity 5, Obesity+T2D 8; snRNA-seq on 11 channels | 8 Ob+T2D |
| GSE268904 | Deep subcutaneous adipose snRNA-seq, 14 samples: 4 lean, 7 obese, 3 obese diabetic | 3 diabetic |
| GSE278526 | Adipose stromal cells, scRNA-seq, obese versus obese+T2D; 58 recruited, 10 profiled; 4 pooled libraries | Few |
| GSE249089 | METSIM subcutaneous adipose snRNA-seq, 84 men, 4 donors pooled per run (demuxlet); raw files withheld; T2D labels unverified | Unknown |
| GSE154126 (Avrahami) | Islets, 22 donors (adult controls 4, T2D 10); only 11 to 54 cells per adult donor | 10 T2D |

---

## 2. Methodological comparison matrix

Only the three data papers and the two best-practice papers report concrete choices. SQU21, JUN22, SCH21, LUE22, ANT25, DAN22 and BUT21 are benchmarks or methods and are discussed under the relevant step.

| Step | LT19 (2019) | HEU23 (2023) | ELG23 (2023) | LAW17 (2017) | SEG16 (2016) |
|---|---|---|---|---|---|
| scRNA vs snRNA | Not addressed | Not extracted | scRNA | scRNA | scRNA; reports FACS/dissociation shifts proportions [A] |
| Platform | n/a | n/a | HPAP, chemistry varies (modelled) | Fluidigm C1 | Smart-seq2 |
| Fresh vs frozen | Not addressed | Not extracted | Not extracted | Not extracted | Not extracted |
| Replication | n/a | Replicates are donors | 29 ND, 17 T2D | 5 ND, 3 T2D | 6 healthy, 4 T2D |
| QC metrics | Genes, counts, mito fraction jointly; permissive; per sample; no data peeking | Sample-wise MAD thresholds | Not extracted; 4,382 barcodes by Scrublet plus 13,036 manually removed | Not extracted | 3,386 to 2,209 cells; thresholds Not extracted |
| Doublets / ambient | Mentioned | scDblFinder; SoupX, CellBender | Scrublet; SoupX | Not extracted | Not extracted |
| Normalization | scran then log(x+1); avoid unit-variance scaling | Shifted log | Not extracted | Not extracted | Not extracted |
| HVG | 1,000 to 5,000 | Deviance-based | Not extracted | Not extracted | Not extracted |
| Batch handling | Do not test corrected data | Harmony simple; scVI, scANVI, scGen, Scanorama complex | Harmony (donor, chemistry, tissue source) for the map; covariates in DE | Not extracted | Not extracted |
| Dim. reduction | PCA to summarise; UMAP to look | Same | Not extracted | Not extracted | Not extracted |
| Clustering | Louvain on kNN (k 5 to 100); subclustering can create false structure | Leiden, several resolutions | Leiden, resolution 0.5 | Not extracted | Not extracted |
| Annotation | Markers; no marker p-value validation | Automated, manual, verify | 10 cell types; method Not extracted | Signatures for alpha, beta, delta, PP | Rare-cell programs |
| Differential expression | Test measured, not corrected or imputed data | Pseudobulk (limma, edgeR, DESeq2) or mixed model; correct for multiple testing | Pseudobulk DESeq2 Wald, BH FDR 0.10 | edgeR, sex blocking, FDR 5% | Adj. p 0.01 or lower; test Not extracted |
| Pathway analysis | Not addressed | Not extracted | Reports OXPHOS and mitochondrial genes down; enrichment method Not extracted | eQTL target overlap | Not extracted |
| Cell-cell comm., trajectory, velocity, multiomics, spatial, ML | Not addressed | Mentioned as modalities; detail Not extracted | Not reported | Not reported | Not reported |

**Not comparable from the supplied set:** Seurat versus Scanpy (software was not extracted for any data paper), QC numeric thresholds, and ribosomal-gene handling.

---

## 3. Advantages and disadvantages of major methods

For each method: the evidence base is stated. Where the supplied papers do not decide the question I say so and do not declare a winner.

### 3.1 Cell-level DE versus pseudobulk DE (strong evidence)
- **Problem solved:** cells from one donor are not independent; testing them as independent inflates significance. [SQU21, JUN22]
- **Why used:** best performer in SQU21 (all six top methods) and well calibrated in JUN22; recommended by HEU23; used by ELG23.
- **Advantages:** respects the donor as replicate; simple, standard count models (edgeR, DESeq2, limma); summed counts beat averaged ones [JUN22].
- **Disadvantages:** needs enough cells per donor and cell type (Avrahami donors contribute 11 to 54 cells, so pseudobulk would be very sparse [I]); loses within-donor cell-level heterogeneity; needs enough donors.
- **Works well:** many donors, hundreds of cells per donor per type, a well-defined cell type.
- **Misleading when:** few donors; rare cell types; batch confounded with condition; clusters defined from the same data.
- **Control versus T2D:** appropriate and the default. T2D donor numbers (3 to 17 in the supplied data) are the binding constraint.

### 3.2 Pseudobulk versus mixed models (moderate evidence)
JUN22 found pseudobulk at least as good as mixed models; HEU23 lists both. A mixed model keeps cells and models a donor random effect. It can be slower and less stable with few donors **[I]**. Neither paper gives a clear winner, so state this as an open choice, not a settled one.

### 3.3 Wilcoxon versus negative-binomial models (partial evidence)
SQU21 reports a Wilcoxon-based list validating only 3 of 13 versus 5 of 6 for pseudobulk, and single-cell tests biased to highly expressed genes. The key lever is the **unit of analysis**, not Wilcoxon itself; a Wilcoxon test on pseudobulk samples is a different procedure that these papers do not benchmark here. Negative-binomial models (DESeq2, edgeR) are the count-appropriate choice used by ELG23. **[needs additional literature]** for Wilcoxon on pseudobulk with small n.

### 3.4 Harmony versus Seurat integration versus scVI (two benchmarks that answer different questions)
- LUE22: scANVI, Scanorama, scVI, scGen lead on complex tasks; Harmony and Seurat v3 on simple ones; batch removal trades off against biology conservation.
- ANT25: with no true batch effect, MNN, scVI and LIGER altered data most, Seurat and others added artifacts, Harmony least. Authors recommend only Harmony.
- **Reading together:** LUE22 tests whether real batch effects are removed while keeping biology; ANT25 tests whether the method leaves null data alone. A method can be strong on the first and poorly calibrated on the second. For T2D, integration is for clustering and visualisation. For testing, model batch in the regression [LT19, HEU23].
- **Misleading when:** batch and disease are confounded. No method can separate them; design must.

### 3.5 Louvain versus Leiden (weak direct evidence)
LT19 (2019) recommends Louvain; HEU23 and ELG23 use Leiden. No supplied paper compares them directly, and the original Leiden paper is not in the set **[needs additional literature]**. The practical point in both reviews is the same: resolution is a choice, so test several and verify with markers.

### 3.6 Manual versus automated annotation (partial evidence)
HEU23 recommends a three-step workflow (automated, manual markers, verification). LT19 warns that marker p-values do not validate clusters because clusters and markers come from the same data. No supplied paper benchmarks CellTypist, SingleR or Azimuth **[needs additional literature]**. In islets, ambient INS and GCG can mimic hormone expression in other cells **[I]**, which makes marker-only annotation fragile.

### 3.7 Normalization: SCTransform versus conventional log (not decided by the supplied papers)
LT19 suggests scran then log; HEU23 suggests shifted-log with deviance-based feature selection. Neither extracted passage discusses SCTransform **[needs additional literature]**. LT19 does warn that log transformation can introduce false DE when size factors differ between groups; for pseudobulk this concern is reduced because counts, not logs, are modelled **[I]**. **[background, not from the supplied papers]** SCTransform is a regularised negative-binomial approach; no recommendation can be made here.

### 3.8 PCA + UMAP versus alternatives
Both reviews treat PCA as a summary and UMAP as a picture, not a measurement. Alternatives (t-SNE, diffusion maps, scVI latent space) are not evaluated in the set **[needs additional literature]**. Never test or quantify on the 2D embedding.

### 3.9 GSEA versus ORA versus GSVA (no supplied evidence)
No supplied benchmark exists. The repo's earlier search found no single-cell-specific GSEA benchmark and only listed decoupleR as an unread lead (`literature_gaps_filled.md`). Defensible default **[I]**: rank genes by the donor-level pseudobulk statistic and run enrichment on the ranked list, not on a thresholded list from a cell-level test. Flag as **[needs additional literature]**.

### 3.10 scRNA-seq versus snRNA-seq
SEG16 reports that dissociation and FACS shifted cell-type proportions (more alpha, fewer beta than in tissue) [A]. No supplied paper compares scRNA and snRNA directly in any tissue. **[background, not from the supplied papers]** Adipocytes are generally too large and buoyant for standard scRNA-seq, which is why adipose studies commonly use nuclei; the repo's literature search found no single-nucleus versus single-cell comparison for adipose or islets (`literature_gaps_filled.md`). **[needs additional literature]**.

### 3.11 Seurat versus Scanpy
Not addressed by the supplied papers. Both implement the same pipeline; choose by language and ecosystem. No claim of superiority.

### 3.12 Compositional methods: scCODA versus Milo
- scCODA [BUT21]: models proportions jointly with a reference cell type; needs predefined types; does not model donor-to-donor variability in response [A].
- Milo [DAN22]: kNN neighbourhoods, negative-binomial GLM; needs replicates; avoid complete batch-condition confounding; a neighbourhood is not necessarily a distinct subpopulation [A].
- HEU23: no independent benchmarks exist for compositional analysis. Treat both as reasonable, neither as validated for T2D.

---

## 4. Statistical methodology critique

### 4.1 Core principle
The biological replicate is the donor. Cells are nested within donors, share genetics, age, batch and ambient RNA, and are not independent [SQU21, JUN22, LT19]. Treating 20,000 cells from 10 + 10 donors as 20,000 samples inflates significance. The consequence for Control versus T2D is a gene list that reflects donor and batch differences rather than disease.

### 4.2 Per-paper assessment

| Paper | True replicate | Cells treated as independent? | Donor effect handled? | Pseudobulk | Batch / confounding | Multiple testing | Effect sizes |
|---|---|---|---|---|---|---|---|
| ELG23 | Donor | No (donor-level) | Yes, pseudobulk | Yes, DESeq2 Wald | Chemistry and tissue source as covariates; Harmony used for the map; DE on pseudobulk counts | BH FDR 0.10 (looser than the common 0.05) [I] | Not extracted |
| LAW17 | Donor | **Not determinable from extracted methods.** Only sex blocking and a cells-based gene filter were seen; supplement not checked | Not extracted | Not extracted | Not extracted | FDR below 5% | Not extracted |
| SEG16 | Donor | **Not determinable**; test not extracted | Not extracted | Not extracted | Not extracted | Adj. p 0.01 or lower | Not extracted |
| SQU21 | Donor (explicitly) | Shows cell-level tests fail | Yes | Yes | n/a | n/a | Reports validation rates |
| JUN22 | Sample | Shows naive tests inflate (up to 40%) | Yes | Yes | n/a | Nominal error 0.05 | n/a |
| SCH21 | Donor | n/a (power) | Models donors | n/a | n/a | n/a | Example uses a leukaemia-derived large effect size, not T2D |
| LUE22, ANT25 | n/a | n/a | n/a | n/a | They address batch | n/a | n/a |
| DAN22, BUT21 | Donor/sample | Count per-sample | Milo yes; scCODA does not model response variability [A] | n/a | Milo: avoid confounding [A] | FDR (spatial FDR; Bayesian credible effects) | Log fold change per neighbourhood / effect on abundance |

### 4.3 Specific concerns that could change a Control versus T2D conclusion
1. **Possible pseudoreplication in the older islet studies [I].** LAW17 (8 donors) reports 248 beta-cell genes; ELG23 (46 donors in the T2D/ND comparison) reports 79. The shrinkage is consistent with the SQU21 and JUN22 warnings, but the studies also differ in platform, thresholds and tissue handling. It is a reason to check, not proof. For LAW17 I could not establish how donors were handled.
2. **Disease confounded with batch.** In GSE221156 single-donor libraries, T2D has 7 V2 and 7 V3 libraries, ND has 4 V2 and 13 V3, PD 1 and 10; UPenn islets include 3 T2D and no ND donors [GEO record]. ELG23 models chemistry and tissue source; a naive analysis would attribute chemistry to T2D.
3. **BMI tracks with disease.** In GSE221156 mean BMI rises from 28.4 (ND) to 30.8 (PD) to 34.0 (T2D) [GEO record]. T2D effects partly overlap with obesity unless BMI is modelled. ELG23 uses scaled BMI as a covariate, but BMI and diabetes are correlated, so covariate adjustment cannot fully separate them **[I]**.
4. **Very low power in small studies.** With 3 T2D donors (LAW17, GSE268904) donor-level inference is exploratory [SCH21 logic]. SCH21's 74% figure applies to very large effects; T2D effects are smaller, so expect less power.
5. **Cells per donor.** Avrahami adult donors contribute 11 to 54 cells [GEO record], so cell-type pseudobulk would be sparse **[I]**.
6. **Gene-count thresholds in ELG23.** The expression filter (5 counts in at least half of samples per group) is reasonable but can drop genes that are on or off in one group **[I]**.
7. **Proportions.** SEG16 reports dissociation shifts; compositional effects (BUT21) mean a beta-cell loss raises other cell types' shares even if their absolute numbers are unchanged.
8. **Circularity.** Clusters defined from the data and then tested on the same data yield optimistic p-values for markers [LT19].
9. **Significance versus biological significance.** Effect sizes were not extracted for the data papers; gene-count comparisons across studies mix different thresholds.

---

## 5. Strengths and limitations matrix

[A] author-reported (as extracted), [I] my own assessment.

| Paper | Main method | Major strengths | Major limitations | Statistical concerns | Biological limitations | Potential improvement |
|---|---|---|---|---|---|---|
| ELG23 | Harmony map; pseudobulk DESeq2 | Largest islet map (192,203 cells); donor-level model; explicit covariates | [I] FDR 0.10; heavy manual barcode removal (13,036); T2D 17 donors versus 29 ND | [I] BMI and T2D correlated; unbalanced groups | [I] Cadaveric islets; one cohort | Independent replication; sex-stratified analysis; sensitivity to FDR |
| LAW17 | edgeR per cell type | Early cell-type-resolved T2D signatures; eQTL link (536 genes, 263 cell-type-specific) | [I] 3 T2D donors; 1,050 cells | [I] Donor handling not established | [I] Small cohort | Re-analyse with donor-level pseudobulk |
| SEG16 | Smart-seq2 plus differential expression | Full-length data; rare cell types; BMI link | [A] Dissociation and FACS shift proportions | [I] Test not extracted; 4 T2D donors | [A] Proportions biased | Compare with tissue-level data |
| LT19 | Tutorial | Clear step-wise guidance | [A] Dated tool choices (2019) | None for DE; consistent with later work | Not tissue-specific | Update with HEU23 |
| HEU23 | Review | Uses independent benchmarks | [A] No benchmarks for compositional analysis; DE tools low consensus | Recommends pseudobulk | Not tissue-specific | Add disease-specific validation |
| SQU21 | Benchmark | Gold-standard matched bulk | [I] Not a T2D benchmark; transfer to T2D untested | Strong | Transfer to T2D untested | Test in T2D cohorts |
| JUN22 | Simulation plus null split | Real null test (14 donors); many methods | [I] B-cell null test; 100 to 2,000 cells per sample simulated | Strong | Not T2D | Replicate in islets and adipose |
| SCH21 | Power model | Practical design guidance | [A] Needs discrete cell types; priors reliable only up to pilot size | Example effect sizes unlike T2D | Not T2D | Pilot-based T2D power analysis |
| LUE22 | Integration benchmark | 16 methods, 13 tasks, about 1.2M cells | [I] Does not test null calibration | Trade-off between batch removal and biology | Atlas tasks, not T2D | Pair with ANT25 |
| ANT25 | Null calibration | Direct test of over-correction | [I] Tests one property (null calibration); recommends only Harmony | Strong, narrow | Not T2D | Test on confounded designs |
| DAN22 | Milo | Cluster-independent; FDR across batch effects | [A] Needs replicates; not always distinct subpopulation | Avoid confounding | Mouse and liver examples | Apply to islets |
| BUT21 | scCODA | Compositional; FDR control | [A] Predefined types; no donor-level response variability | Reference choice affects interpretation | Not T2D | Apply with donors in T2D |

---

## 6. Methodological evolution of the field

Chronological order of the supplied papers: SEG16 (2016), LAW17 (2017), LT19 (2019), SQU21, SCH21, BUT21 (2021), JUN22, DAN22, LUE22 (2022), HEU23, ELG23 (2023), ANT25 (2025).

| Area | 2016 to 2019 | 2021 to 2023 | 2025 | Status |
|---|---|---|---|---|
| QC | Per-cell metrics; permissive per-sample thresholds [LT19] | MAD-based per-sample thresholds; ambient RNA (SoupX, CellBender) and doublet tools [HEU23, ELG23] | Not addressed | Newer approach addresses ambient RNA, which is particularly relevant to INS and GCG in islets **[I]** |
| Normalization | scran then log [LT19] | Shifted-log; deviance-based feature selection [HEU23] | Not addressed | Debate open; SCTransform not covered |
| Clustering | Louvain [LT19] | Leiden at several resolutions [HEU23, ELG23] | Not addressed | Resolution choice unresolved |
| Annotation | Marker-based; warn against circularity [LT19] | Three-step workflow [HEU23] | Not addressed | No supplied benchmark |
| Batch correction | Test measured data [LT19] | Benchmarked [LUE22]; recommended by task type [HEU23] | Calibration shows over-correction [ANT25] | Evolving toward modelling batch and testing calibration |
| Differential expression | Cell-level (older islet studies, as far as extracted) | Pseudobulk shown superior [SQU21, JUN22]; adopted [HEU23, ELG23] | Not addressed | Clearest methodological shift |
| Abundance | Raw proportions | scCODA, Milo | Not addressed | No independent benchmark [HEU23] |
| Pathway analysis | Not addressed | Not addressed | Not addressed | **[needs additional literature]** |
| Cell-cell comm., trajectory, multimodal, ML | Not covered | scVI-type models appear only as integration methods [LUE22, HEU23] | Not covered | **[needs additional literature]** |

**Older approaches that remain useful:** log-normalised, per-sample QC (LT19); Harmony for simple tasks (LUE22, ANT25). **Newer approaches that fix specific weaknesses:** pseudobulk fixes pseudoreplication; scCODA and Milo address compositionality and cluster dependence; null-split calibration (ANT25, JUN22) exposes over-correction and false positives; Leiden and deviance feature selection are claimed improvements with limited direct evidence in this set. The newest method is not automatically best: ANT25 found newer integration tools less well calibrated than Harmony.

---

## 7. Contradictions and disagreements

| Issue | What the papers say | Why they may differ | Resolution |
|---|---|---|---|
| Louvain vs Leiden | LT19 Louvain; HEU23, ELG23 Leiden | Five years of tool development; no head-to-head in the set | Not settled; test multiple resolutions |
| Normalization | LT19 scran + log; HEU23 shifted-log + deviance | Different benchmark evidence | Not settled |
| Which integration method | LUE22: scVI, scANVI, Scanorama, scGen lead on complex tasks; ANT25: only Harmony recommended | They ask different questions (batch removal and biology conservation versus calibration on null data) | Not contradictory; use both and model batch for testing |
| Integrate first or model batch | HEU23 and ELG23 use Harmony; LT19 and ANT25 argue against testing on corrected values | ELG23 uses Harmony for the map but pseudobulk covariates for DE | Consistent if integration is only for clustering |
| Pseudobulk vs mixed model | JUN22: pseudobulk at least as good; HEU23: both acceptable | Different benchmarks | Open choice |
| Number of T2D genes | LAW17 248 beta genes (8 donors); ELG23 79 beta genes (46 donors) | Possible pseudoreplication, platform, thresholds, tissue handling, donor numbers | Cannot attribute to a single cause |
| Gene-level T2D findings | ELG23: TSHR, SLC4A4, TNFRSF11B up, OXPHOS down; SEG16: INS, FXYD2 down in beta, WFS1 down in alpha | Gene lists from LAW17 and SEG16 were not extracted in full, so overlap cannot be assessed here | **[needs additional literature]**: extract full gene lists and compare |
| QC strictness | LT19 permissive; HEU23 MAD per-sample; ELG23 manually removed 13,036 barcodes in addition to Scrublet | Different goals and datasets | No consensus; thresholds must be documented |

---

## 8. Research gaps

### Biological gaps
- Which T2D changes are cell-type-specific versus shared across cell types, beyond beta cells, is only partly answered in the supplied papers (ELG23: 79 of 84 genes in beta cells).
- Disease progression: pre-diabetic (PD) donors exist in GSE221156 (14) but no supplied paper analyses PD versus ND versus T2D as a continuum.

### Dataset gaps
- No supplied T2D paper covers adipose, muscle or liver. GEO records show small T2D arms (3 in GSE268904, 8 Ob+T2D in GSE342773, 10 profiled in GSE278526) and one larger cohort whose labels are unverified (GSE249089, n = 84).
- Sex and ancestry: GSE221156 has 5/17 female ND and 6/14 female T2D; no supplied paper reports sex-stratified T2D results.

### Statistical gaps
- Few T2D donors; unbalanced groups; effect sizes seldom reported in the extracted text.
- No benchmark for compositional methods (HEU23).
- No power analysis built on real T2D effect sizes (SCH21 example is leukaemia-derived).
- Confounding of batch and disease is not solved by integration (LUE22, ANT25 imply it cannot be).

### Computational gaps
- Pipelines differ between papers (QC, normalization, clustering); no common reference pipeline for T2D.
- Pathway analysis, cell-cell communication, trajectory and RNA velocity are absent from the supplied set **[needs additional literature]**.

### Integration gaps
- Multi-omics (scATAC, proteomics, metabolomics, spatial, clinical) are not covered by the supplied papers. LAW17 links islet eQTLs to cell-type expression, the only genetic integration present. **[needs additional literature]**

### Validation gaps
- ELG23 and LAW17 findings are computational; experimental or independent replication were not extracted. HPAP and PanKbase overlap, so replication needs donor-ID checking.

### Reproducibility gaps
- Code availability and QC thresholds were not extracted for any data paper.
- The GSE249089 raw files are withheld for privacy; GSE221156 mixes pooled and single-donor libraries requiring demultiplexing; metadata completeness varies.

---

## 9. T2D-specific gaps (what the supplied literature does and does not support)

| Topic | Does the supplied literature support this as a gap? |
|---|---|
| Islets, beta cells | Well covered by three papers; **gap** is replication with donor-level methods and independent cohorts (ELG23 vs LAW17/SEG16 differences) |
| Adipose tissue, adipocytes, stromal, endothelial, macrophages | **No supplied paper.** GEO records show small T2D arms and an obesity-versus-T2D design in GSE342773; evidence is about datasets only. **[needs additional literature]** to claim a biological gap |
| Skeletal muscle, liver | No supplied paper. Repo leads exist (E-MTAB-15009: 10 controls, 9 T2D, 38 biopsies, 135,225 nuclei; liver leads unchecked) |
| Immune populations, inflammatory pathways | Not covered |
| Insulin resistance, metabolic and mitochondrial dysfunction | Partly: ELG23 reports mitochondrial and OXPHOS genes down in T2D beta cells; no adipose or muscle support |
| Cellular heterogeneity | SEG16 reports beta-cell subpopulations; no replication in the supplied set |
| Disease progression | PD donors exist in GSE221156; no supplied analysis |
| Sex differences | Sex used as a covariate (ELG23) or blocking factor (LAW17); no stratified result extracted. **Supported as an underexplored area, but only weakly** |
| Obesity versus T2D | **Supported:** SEG16 correlates BMI with expression; GSE221156 BMI rises with disease; GSE342773 has Obesity versus Obesity+T2D arms. No supplied paper separates the two |
| Treatment effects | Not reported in the supplied papers; medication metadata not extracted |

---

## 10. Future research opportunities (realistic and publishable)

1. **Donor-level re-analysis and replication of islet T2D effects**
   - Question: which T2D beta-cell genes replicate across independent cohorts at donor level?
   - Limitation: gene counts fall from hundreds to dozens as donors increase (LAW17 vs ELG23); methods differ.
   - Method: pseudobulk per donor and cell type with covariates (sex, age, BMI, chemistry, center); GSE221156 as discovery and HPAP as replication after checking donor overlap.
   - Data: GSE221156 (17 T2D), HPAP/PANC-DB.
   - Methods: edgeR/DESeq2/limma-voom; BH FDR; effect-size reporting.
   - Contribution: a replicated T2D beta-cell signature with effect sizes.
   - Limitation: access to HPAP; overlap between cohorts; unequal chemistry.
2. **Quantify pseudoreplication on real T2D data**
   - Question: how much do cell-level and donor-level analyses disagree in T2D islets?
   - Existing limitation: the claim is supported by non-T2D benchmarks only (SQU21, JUN22).
   - Method: apply cell-level, pseudobulk and mixed-model analyses to the same T2D datasets and run null splits within controls.
   - Data: LAW17, SEG16, Avrahami, GSE221156.
   - Contribution: T2D-specific evidence for the donor-as-replicate rule.
   - Limitation: low cells per donor in the older sets.
3. **Separate obesity from T2D in adipose snRNA-seq**
   - Question: which adipose cell-type changes are T2D-specific once obesity is held constant?
   - Existing limitation: BMI tracks with disease; supplied papers do not address adipose.
   - Data: GSE342773 (Obesity 5 vs Obesity+T2D 8), GSE249089 (if T2D labels exist), GSE268904.
   - Methods: pseudobulk with BMI; scCODA/Milo for composition; scPower for planning.
   - Contribution: a design that controls adiposity.
   - Limitation: small n; may be exploratory. **[needs additional literature]**
4. **Compositional analysis with donor-level replication across scRNA and snRNA**
   - Question: how robust are T2D cell-type proportion changes to dissociation bias and to compositional closure?
   - Methods: scCODA and Milo with sensitivity to the reference type; compare tissues.
   - Contribution: practical guidance where HEU23 says no benchmark exists.
   - Limitation: dissociation bias cannot be removed by statistics alone.
5. **Calibration of batch handling in confounded T2D designs**
   - Question: does modelling batch or integrating first produce more reliable T2D calls when chemistry and centre are confounded with disease?
   - Method: A/A null splits (ANT25) within controls; compare modelling batch with Harmony and scVI embeddings.
   - Data: GSE221156 (V2/V3 and center structure).
   - Contribution: an empirical test on a real confounded T2D dataset.
   - Limitation: one dataset.
6. **Power planning for T2D with realistic effect sizes**
   - Question: how many donors are needed for T2D effects as small as those observed (for example 84 genes at 46 donors)?
   - Method: scPower calibrated with effect sizes from ELG23 or GSE221156 pseudobulk.
   - Contribution: design guidance for adipose and muscle.
   - Limitation: needs real effect-size priors.

Items needing new literature before they are defensible: pathway/enrichment choice, cell-cell communication, adipose-specific biology, scRNA versus snRNA comparison.

---

## 11. Detailed literature review outline

Items marked **[add lit]** need papers beyond the supplied set.

### 1. Introduction
1.1 Single-cell transcriptomics and cellular heterogeneity. 1.2 Why tissues rather than averages matter in T2D (islet, adipose, muscle, liver). **[add lit]** 1.3 Motivation: donors as the unit of evidence in Control versus T2D; the statistical problem in one paragraph (SQU21, JUN22).

### 2. Technologies
2.1 scRNA-seq platforms: Smart-seq2 (SEG16), Fluidigm C1 (LAW17), droplet (ELG23). 2.2 snRNA-seq and why adipose uses nuclei **[add lit]**. 2.3 Dissociation bias (SEG16). 2.4 Pooling and demultiplexing (GSE221156, GSE249089).

### 3. Computational workflow
3.1 QC: per-sample thresholds, ambient RNA, doublets (LT19, HEU23, ELG23). 3.2 Normalization and features (LT19, HEU23; SCTransform **[add lit]**). 3.3 Dimensionality reduction (PCA, UMAP limits). 3.4 Clustering: Louvain, Leiden, resolution (LT19, HEU23, ELG23). 3.5 Annotation and circularity (LT19, HEU23).

### 4. Statistical analysis
4.1 Biological replication and pseudoreplication (SQU21, JUN22). 4.2 Pseudobulk and mixed models. 4.3 Batch effects and integration (LUE22, ANT25). 4.4 Power (SCH21). 4.5 Multiple testing and effect sizes (ELG23, LAW17). 4.6 Covariates: sex, age, BMI, chemistry, centre.

### 5. Biological interpretation
5.1 Pathway analysis (GSEA, ORA, GSVA, decoupleR) **[add lit]**. 5.2 Cell composition (BUT21, DAN22, HEU23). 5.3 Cell-cell communication **[add lit]**. 5.4 Trajectory and RNA velocity **[add lit]**.

### 6. Advanced and emerging methods
6.1 scVI/scANVI and deep-learning integration (LUE22, ANT25). 6.2 Reference mapping and foundation models **[add lit]**. 6.3 Multimodal, spatial, multiome **[add lit]**.

### 7. Critical comparison of existing T2D studies
Table of donors, platform, test, replicate unit, gene count (LAW17, SEG16, ELG23, plus datasets from GEO). Statistical limitations; confounding; reproducibility.

### 8. Research gaps
Biological, dataset, statistical, computational, integration, validation, reproducibility (Section 8 above).

### 9. Future directions
Proposed analyses (Section 10 above); recommended analytical framework (Section 12).

### 10. Conclusion
State of the field; strongest approaches (donor-level pseudobulk, modelled batch); remaining challenges; research opportunity.

---

## 12. Recommended conceptual framework for your research

**Principle: the donor is the unit of evidence; the cell is a measurement within it.**

1. **Design audit (before any analysis).** Tabulate donors per group, and chemistry, centre, sex, age, BMI, HbA1c, pooled libraries. Flag confounding (for example GSE221156 chemistry and centre). Decide the comparison (T2D vs ND; and Obesity vs Obesity+T2D in adipose).
2. **Power check.** Use scPower with effect sizes from published pseudobulk T2D results to judge feasibility. Declare small cohorts exploratory.
3. **Cell-level processing, per sample.** Ambient RNA and doublet handling; per-sample QC; do not tune thresholds to the result. Keep raw counts.
4. **Embedding for annotation only.** HVG, PCA, optional Harmony for clustering, Leiden at several resolutions; annotate with automated plus marker plus verification. Treat UMAP as a picture.
5. **Donor-level testing.** Sum raw counts per donor and cell type; edgeR/DESeq2/limma-voom (mixed model as a sensitivity analysis) with sex, age, BMI, batch; report effect sizes and FDR.
6. **Composition.** scCODA and Milo at donor level; test sensitivity to reference type; state dissociation bias.
7. **Calibration.** Run null splits within controls (JUN22, ANT25) to estimate false-positive rate in your own data.
8. **Interpretation.** Enrichment on ranked donor-level statistics; label as hypothesis-generating until replicated.
9. **Replication.** Independent cohort with donor-ID overlap checked.

---

## 13. What the Literature Tells Us and What It Still Cannot Tell Us

**Strongest evidence.**
- Donor, not cell, is the replicate. Two independent benchmarks agree: SQU21 (six top methods all pseudobulk) and JUN22 (naive tests up to 40% false positives in a real null split versus about 5% for pseudobulk and mixed models).
- Integration methods differ in what they protect: LUE22 and ANT25 test different properties, and neither supports testing on corrected data.
- Compositional data need compositional methods (BUT21, DAN22), though no independent benchmark exists (HEU23).

**Most reliable analytical approach in this literature.** Per-sample QC with ambient RNA and doublet handling; donor-level pseudobulk with edgeR, DESeq2 or limma; batch and covariates in the model; FDR control; donor-level composition testing; replication in an independent cohort. ELG23 is the closest published template.

**Unresolved methodological debates.** Louvain versus Leiden; scran-log versus shifted-log versus SCTransform; mixed model versus pseudobulk; which integration method; manual versus automated annotation; GSEA versus ORA versus GSVA; scRNA versus snRNA. The supplied papers do not settle any of these directly.

**Important limitations of current T2D studies.** Few T2D donors in the older islet studies (3 and 4), a single large cohort (ELG23), confounding of chemistry, centre and BMI with disease (GSE221156), dissociation bias in proportions (SEG16), limited reporting of effect sizes in what I extracted, and no replication across independent cohorts documented in the supplied papers.

**Strongest research gaps.** (i) Donor-level, covariate-adjusted, independently replicated T2D signatures in islets; (ii) separating obesity from T2D in adipose, where cohorts are small and the supplied literature is silent; (iii) T2D-specific quantification of pseudoreplication and calibration; (iv) compositional analysis with donor-level replication.

**Where a new study can contribute.** A study that (a) uses a donor-level model with explicit batch and BMI handling, (b) quantifies its own false-positive rate by null splitting, (c) reports effect sizes and replication, and (d) extends to adipose or muscle with an obesity-matched design would be both publishable and methodologically defensible. Claims about adipose biology, pathway methods, cell-cell communication and scRNA-versus-snRNA require additional literature that was not supplied.
