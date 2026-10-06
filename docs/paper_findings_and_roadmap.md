# What the 12 papers say, and what to read next

Prepared 2026-10-06. Findings below were read from the papers' full text on PubMed Central (numbers and quoted ideas are theirs).
Where I could not extract something, I say so. Items marked "verify" were not checked in this session.

## Part 1. The 60-second version for a technical researcher

Single-cell RNA-seq turns a tissue into a big sparse table of counts: genes by cells. The biology question here is about people,
not cells: does type 2 diabetes (T2D) change which cells are present, and what each cell type expresses?
The papers fit together into one argument.

1. A standard pipeline exists and is well documented (papers 1, 2): clean the cells, normalize, pick informative genes, reduce dimensions, cluster, label.
2. The unit of evidence is the donor, not the cell (papers 3, 4, 5). Cells from one person are not independent. Tests that treat 20,000 cells as 20,000 samples
   find hundreds of "differences" even when nothing differs. Add up each donor's cells (pseudobulk) and test across donors.
3. Donors also differ by batch (kit, lab, site), so you must correct for batch, but correction can damage the signal (papers 9a, 9b).
4. Cell-type proportions are a second outcome and need their own statistics because they sum to 1 (papers 10a, 10b).
5. Applied to islets, three T2D papers (6, 7, 8) show what is found in practice: the best-powered study (46 donors in the T2D/ND comparison) reports 84 T2D-associated genes,
   mostly in beta cells; older small studies report hundreds.

If you remember one idea: effective sample size = number of donors.

## Part 2. A primer in your language

| Single-cell idea | Analogy you already know |
|---|---|
| Count matrix, genes by cells | A sparse count feature matrix; overdispersed, so Poisson fails and negative binomial is the usual model |
| Donor | The group or cluster in a hierarchical model; observations are nested inside it |
| Pseudoreplication | Treating correlated rows as iid: an inflated n, anti-conservative p-values, like train-test leakage for inference |
| Pseudobulk | Aggregate to the true independent unit, then fit a regression (edgeR, DESeq2, limma-voom) |
| Batch effect | Domain shift between runs or sites |
| Integration | Domain alignment; the risk is aligning away the label you want to study |
| Calibration test (paper 9b) | An A/A test: split identical data randomly in two and check you find nothing |
| Clustering and UMAP | Unsupervised grouping plus a 2D picture; UMAP distances and cluster sizes are not quantitative |
| Compositional data | Data on the simplex: one proportion rising forces others to fall (closure) |
| FDR | Control of the expected fraction of false positives among reported genes, applied per cell type |

## Part 3. Findings outline, paper by paper

### Theme A. The standard pipeline

**1. Luecken and Theis, Mol Syst Biol 2019 (tutorial).**
- Problem: many tools, little guidance.
- QC: look at genes per cell, counts per cell and mitochondrial fraction together; be permissive; set thresholds per sample when quality differs; do not tune thresholds to improve a test ("data peeking").
- Normalization: scran suggested; log(x+1) afterwards; avoid scaling genes to unit variance; log transforms can introduce false differential expression when size factors differ between groups.
- Features and reduction: 1,000 to 5,000 highly variable genes; PCA for summarizing; UMAP for looking; a 2D plot is not a summary of the data.
- Clustering: Louvain community detection on a kNN graph (k typically 5 to 100); subclustering can create patterns from noise.
- Annotation: do not use marker p-values to validate clusters (circular, because clusters and markers come from the same data).
- Differential expression: test on measured data, not on batch-corrected or imputed data; do not rely on tools to fix confounded designs.
- Caveat: written in 2019; later papers update the tool choices.

**2. Heumos et al., Nat Rev Genet 2023 (best practices across modalities).**
- Updates tool choices using independent benchmarks: sample-wise QC with median absolute deviations; ambient RNA removal (SoupX, CellBender); doublets (scDblFinder);
  shifted-log normalization; deviance-based gene selection; Leiden clustering at several resolutions; three-step annotation (automated, manual markers, verification).
- Integration: Harmony for simple tasks; scVI, scANVI, scGen or Scanorama for complex ones.
- Differential expression: pseudobulk with summed counts and limma, edgeR or DESeq2, or a mixed model; correct for multiple testing.
- Abundance: use scCODA or related methods (Milo for graph-based).
- What is unresolved: no independent benchmarks for compositional analysis; differential-expression tools show low consensus and a trade-off between sensitivity and precision.

### Theme B. Donors are the replicates

**3. Squair et al., Nat Commun 2021.**
- Question: why do popular single-cell tests disagree, and which are right?
- Design: 18 "gold standard" datasets with matched bulk and single-cell data from the same purified cells; 14 methods (7 single-cell, 6 pseudobulk, 1 mixed model).
- Result: all six top methods were pseudobulk. Single-cell methods found hundreds of differentially expressed genes with no perturbation present, biased to highly expressed genes.
  In one in vivo check, 5 of 6 pseudobulk calls validated versus 3 of 13 for a Wilcoxon-based list. Shuffling replicate labels reduced variance estimates in 98.2% of genes in one dataset.
- Take-away: inference must be at the level of biological replicates.

**4. Junttila, Smolander, Elo, Brief Bioinform 2022.**
- Design: 18 methods on 1,280 simulated datasets (6 to 40 samples, 100 to 2,000 cells per sample) plus muscat-based simulations; real null test: B cells from 14 healthy donors split at random into 2 groups, 30 times.
- Result: naive and latent-variable methods had false-positive proportions up to 40% in the null test; pseudobulk methods and mixed models stayed near the nominal 0.05 error. Summing counts beat averaging. Pseudobulk was at least as good as mixed models.
- Recommendation: ROTS with summed counts for precision; limma, DESeq2 or edgeR with summed counts for sensitivity; avoid Seurat's latent-variable tests.

**5. Schmid et al., Nat Commun 2021 (scPower; the index lists it under an earlier preprint title, same DOI).**
- Question: how many individuals, cells per person and reads are needed?
- Result: shallow sequencing of many cells beats deep sequencing of fewer cells for a fixed budget. Worked example: for a differential-expression study using large effect sizes (84 genes, median absolute log fold change 2.8),
  about 74% power needs 3,000 cells per cell type per person and 10 donors per group.
- Caveats: effect sizes in real T2D data are much smaller than that example, so expect lower power; needs discrete cell types; priors only reliable up to roughly the pilot's sample size.
- A correction to my earlier shorthand ("donors, not cells, drive power"): the paper's headline is the cells-versus-depth trade-off; the donor-level message comes from papers 3 and 4.

### Theme C. Batch and integration

**9a. Luecken et al., Nat Methods 2022 (integration benchmark).**
- Design: 16 methods, 68 method and preprocessing combinations (590 runs), 13 tasks, 85 batches, about 1.2 million cells, 14 metrics split into batch removal and biological conservation.
- Result: scANVI, Scanorama, scVI and scGen do well on complex tasks; Harmony and Seurat v3 on simpler ones. There is a clear trade-off between removing batch and keeping biology. Scaling raised batch removal in 79% of comparisons and lowered biology conservation in 72%; selecting highly variable genes helped.

**9b. Antonsson and Melsted, Genome Res 2025 (calibration of correction).**
- Design: public datasets randomly split into two pseudobatches 25 times (no true batch effect); 8 methods.
- Result: MNN, scVI and LIGER often altered the data a lot; ComBat, ComBat-seq, BBKNN and Seurat introduced detectable artifacts (for example Seurat and MNN reported over 800 differentially expressed genes where the original had 179 and 154); Harmony changed the data least.
- Their conclusion: Harmony is the only method they recommend.
- Reconcile 9a and 9b: they ask different questions. 9a asks whether a method removes real batch effects and keeps biology; 9b asks whether it leaves data alone when there is no batch effect. Use both. Never run differential expression on corrected values (paper 1); model batch in the regression instead.

### Theme D. Cell-type abundance

**10a. Dann et al., Nat Biotechnol 2022 (Milo).**
- Problem: cluster-based abundance testing depends on how clusters were drawn and fails on continuous states.
- Method: overlapping neighbourhoods on a kNN graph; negative-binomial GLM on per-sample counts; a spatially weighted FDR.
- Findings: sensitive and controls FDR across batch effects in simulations; applied to ageing mouse thymus and cirrhotic human liver.
- Caveats: needs biological replicates; avoid complete confounding of batch and condition; a neighbourhood is not necessarily a distinct subpopulation.

**10b. Büttner et al., Nat Commun 2021 (scCODA).**
- Problem: proportions sum to 1, so one type shrinking makes the rest look larger.
- Method: hierarchical Dirichlet-multinomial model with a spike-and-slab prior; effects are interpreted relative to a reference cell type.
- Findings: average Matthews correlation 0.64 in low-sample simulations; of the compared methods only scCODA, ALDEx2 and ALR-transformed tests controlled FDR in all scenarios; re-found the known B-cell decrease in supercentenarians.
- Caveats: needs pre-defined cell types; does not model donor-to-donor variability in response.

### Theme E. What the T2D islet papers found

**6. Elgamal et al., Diabetes 2023 (HPAP integrated map).**
- Data: 67 donors downloaded (29 non-diabetic, 17 T2D, 10 T1D, 9 autoantibody-positive); 65 in the map; 192,203 cells.
- Pipeline: SoupX for ambient RNA; Scrublet doublets (4,382 barcodes plus 13,036 more removed manually); Harmony using donor, chemistry and tissue source; Leiden at resolution 0.5; 10 annotated cell types.
- Differential expression: pseudobulk counts per donor and cell type, DESeq2 Wald test, covariates sex, scaled age, scaled BMI, chemistry and tissue source, genes tested only if at least half of samples per group had 5 counts, Benjamini-Hochberg FDR 0.10.
- Result: T1D 1,808 genes (1,305 in beta cells); T2D 84 genes (79 in beta cells). In T2D beta cells TSHR, SLC4A4 and TNFRSF11B went up; mitochondrial and oxidative-phosphorylation genes went down.
- Why it matters: this is the modelling template to copy.

**7. Lawlor et al., Genome Res 2017.**
- Data: 8 donors (5 non-diabetic, 3 T2D), 1,050 cells, Fluidigm C1.
- Findings: signatures for alpha, beta, delta and PP cells; T2D versus control 248 genes in beta, 138 in alpha, 24 in delta (plus 74 acinar, 35 ductal, 28 stellate); 536 islet eQTL target genes checked, 263 with cell-type-specific expression.
- Method as stated: edgeR with sex as a blocking factor, FDR below 5%, genes filtered by number of cells. The extracted methods text mentions no donor term or per-donor aggregation. I did not check the supplement, so do not conclude the unit was cells; do ask how replicates were handled before citing the gene counts.

**8. Segerstolpe et al., Cell Metab 2016.**
- Data: 6 healthy and 4 T2D donors; 3,386 cells sequenced, 2,209 kept; Smart-seq2.
- Findings: programs for rare delta, gamma, epsilon and stellate cells; subpopulations of alpha, beta and acinar cells; in T2D, INS and FXYD2 down in beta cells, WFS1 down in alpha cells; BMI correlated with expression of obesity and diabetes genes.
- Lesson: dissociation and FACS shifted cell-type proportions (more alpha, fewer beta than in tissue), so proportions from dissociated cells are biased.
- I could not extract this paper's differential-expression test; adjusted p of 0.01 or lower is what it reports.

**The pattern across 6, 7, 8.** Gene lists shrink as the donor count grows: hundreds of genes from 8 donors, 84 from 46. That is consistent with the warning in papers 3 and 4,
but the studies also differ in tissue handling, platform and thresholds, so treat it as a reason to check, not proof.

## Part 4. How to explain it (talk outline, about 10 minutes)

1. The question: how does T2D change islets, cell type by cell type? (60 s)
2. The data shape: counts by cells; cells nested in donors; 17 versus 17 donors, not 245,000 samples. (90 s)
3. The pipeline in one slide: QC, normalize, variable genes, PCA, neighbours, clusters, labels; one sentence on what each assumes. (90 s)
4. The central statistical problem: pseudoreplication; evidence from papers 3 and 4 (null split: 40% false positives vs about 5%). (2 min)
5. Design: power depends on donors and cells (paper 5); our donors fix the ceiling. (60 s)
6. Batch: needed but dangerous; the two benchmarks answer different questions; model batch, do not test on corrected data. (90 s)
7. Composition: proportions sum to 1; use scCODA or Milo at donor level. (60 s)
8. What the T2D studies found and how their methods differ; our template is Elgamal. (90 s)
9. Our plan and the traps in our dataset (chemistry and center confounded with disease). (60 s)

Questions a sharp listener will ask, with short answers:
- "Why not just use a mixed model?" Paper 4 found pseudobulk at least as good; paper 2 lists both as acceptable.
- "Why not integrate everything first?" Paper 9b: correction can create artifacts; paper 1: do not test on corrected data.
- "Is a cluster a cell type?" Only after markers and replication support it (paper 1).
- "What if we have only 3 T2D donors?" Then the result is exploratory; paper 5 and paper 3 explain why.

## Part 5. Reading roadmap

Time is rough and assumes you read selectively. "Can explain" lines are self-tests.

### Stage 0. Biology you need first (1 to 2 days)
- Read: a short molecular-biology refresher on gene, transcript, exon, intron, UMI; one review of islet biology and T2D pathophysiology (verify: search PubMed for a recent review of T2D beta-cell dysfunction).
- Resource you already have: the glossary and biology list in `mentor_plan_part2.md`.
- Can explain: what a count in the matrix represents; what INS, GCG, SST mark; why T2D is not one disease state (donor BMI, age, HbA1c vary).

### Stage 1. Pipeline literacy (week 1 to 2)
- Papers 1 then 2 (above). Companion docs: Scanpy tutorials and the Bioconductor book Orchestrating Single-Cell Analysis (OSCA) (verify the current links).
- Can explain: each step's assumption; why QC thresholds are per sample.

### Stage 2. QC details that matter in islets (week 2 to 3)
- Ambient RNA and doublets: SoupX, CellBender, Scrublet, scDblFinder (papers named in 1, 2, 6; search by tool name for the original papers).
- Can explain: why ambient INS and GCG create false signals in other cell types.

### Stage 3. Statistics for case-control (week 3 to 5), the core
- Papers 3, 4, 5, then Elgamal's methods section (paper 6).
- Method references to read for the model you will actually run: the DESeq2 paper, the edgeR user guide, limma-voom (verify titles), and the Benjamini-Hochberg FDR idea.
- Concepts to learn: negative binomial model, size factors, design matrices, why a donor-level design matrix must be full rank, confounding of batch and condition.
- Can explain: pseudoreplication with a numeric example; what pseudobulk sums; what FDR 0.10 means in Elgamal.

### Stage 4. Batch and integration (week 5 to 6)
- Papers 9a and 9b; Harmony's original paper (verify); the scVI/scANVI paper already in `literature_reading_list.md`.
- Can explain: why 9a and 9b disagree on scVI and what each measures.

### Stage 5. Abundance (week 6 to 7)
- Papers 10a and 10b; the sccomp preprint listed in `literature_gaps_filled.md`.
- Can explain: closure; what a reference cell type means in scCODA.

### Stage 6. The T2D literature (week 7 to 9)
- Papers 6, 7, 8, then Bandesh (doi:10.1101/2025.01.17.633590; GEO GSE221156) and HPAP (Cell Metab 2022, listed on PANC-DB; verify).
- Task: write a table of each study's donors, platform, test and gene count, and mark which unit was the replicate.

### Stage 7. Interpretation (week 9 to 10)
- GSEA, GO and Reactome, decoupleR (listed in `literature_gaps_filled.md`); rank genes by the pseudobulk statistic.
- Can explain: why enrichment on a thresholded list from a cell-level test is unreliable.

### Stage 8. Later or optional
- Trajectory and RNA velocity; cell-cell communication (CellChat); scATAC and GWAS integration; spatial transcriptomics; adipose-specific papers.

### Concept checklist (tick when you can explain it aloud)
Count matrix; UMI; library size; dropout and sparsity; overdispersion; size factor; normalization; highly variable gene; PCA; kNN graph; Leiden resolution; UMAP limits;
marker gene; annotation circularity; doublet; ambient RNA; batch effect; donor effect; confounding; pseudoreplication; pseudobulk; design matrix; log fold change; p-value; FDR;
compositional data; reference cell type; calibration (A/A test); power; effect size.

## Part 6. Honest limits of this summary
- I read the 12 papers through an extraction tool that returns summaries and quoted passages, not every page. Numbers above were reported in those passages; check any number you will cite against the paper.
- I did not verify Segerstolpe's testing method or whether Lawlor's supplement handled donors; both are flagged.
- Several resources in the roadmap (OSCA, DESeq2, edgeR, limma-voom, Harmony, SoupX and others) are named from memory; search by title before relying on a link.
