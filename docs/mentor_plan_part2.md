# Mentor plan, part 2: remaining deliverables

Verification labels: VERIFIED = read in paper text; ABSTRACT = abstract only; UNVERIFIED = memory, check before use.
Deliverable 3 (12-week roadmap) is in `mentor_plan_part1.md`. Links are in `links.md`.

---
## Deliverable 1: ranked public datasets for Control vs T2D

I cannot honestly give 10 verified datasets. I verified 4 from paper text and have leads for the rest. Do not treat leads as confirmed.

| Rank | Dataset | Tissue | Donors (ND/T2D) | Cells | Status | Role |
|---|---|---|---|---|---|---|
| 1 | HPAP, Elgamal 2023, doi:10.2337/db23-0130 | Islets | 29 / 17 of 67 downloaded (65 in final map) | 192,203 | VERIFIED | Main |
| 2 | Lawlor 2017, GSE86473 | Islets | 5 / 3 | 1,050 | VERIFIED | Learning |
| 3 | Segerstolpe 2016 | Islets | 6 / 4 | 2,209 | VERIFIED counts; accession UNVERIFIED | Learning, replication |
| 4 | Avrahami 2020, GSE154126 | Islets | 4 adult ND + 4 younger groups / 10 T2D | not retrieved | VERIFIED donors | Replication |
| 5 | Bandesh 2025/2026, doi:10.1038/s44318-026-00744-w | Islets | not retrieved | not retrieved | ABSTRACT | Advanced |
| 6 | Xin 2016, GSE81608 | Islets | not retrieved | 1,492 | ABSTRACT; accession UNVERIFIED | Extra replication |
| 7 | Camunas-Soler 2020 patch-seq, doi:10.1016/j.cmet.2020.04.005 | Islets, with electrophysiology | not retrieved | not retrieved | ABSTRACT | Advanced (links function to expression) |
| 8 | Alpha-cell scRNA+snRNA-seq, doi:10.1038/s41467-025-62670-5 | Islets | not retrieved | not retrieved | ABSTRACT | Advanced; snRNA-seq example |
| 9 | Sex differences, 52 donors with and without T2D, doi:10.21203/rs.3.rs-4607352/v1 | Islets, scRNA + snATAC | 52 total | not retrieved | ABSTRACT (preprint) | Validation if data are public |
| 10 | Adipose: no verified T2D-labelled set | Adipose | none | none | Gap | Leads in `literature_reading_list.md` F |

Unsuitable or caution:
- Datasets without diabetes labels (e.g. healthy-only islet atlases, mouse atlases) cannot answer your question.
- Mouse db/db or HFD islet data are not human T2D; use only for hypothesis support.
- Obesity-only adipose snRNA-seq studies do not isolate T2D; any adipose finding would be confounded by BMI.

---
## Deliverable 2: 18 papers and resources (priority order)

Full DOIs and PMC links are in `links.md`. Short list:
1 Luecken & Theis best practices (whole). 2 Best practices across modalities (skim). 3 Confronting false discoveries (whole).
4 Multi-subject DS benchmark (whole). 5 Multi-sample power/design (methods, figures). 6 Elgamal 2023 (whole, copy the pipeline).
7 Lawlor 2017 (results, methods). 8 Segerstolpe 2016 (results). 9 Atlas-level integration benchmark + Genome Res 2025 critique (abstract, figures).
10 Milo + scCODA (methods overview). 11 Avrahami 2020 (dedifferentiation result). 12 Islet identity genesets (annotation markers).
13 decoupleR (tutorial). 14 CellChat Nat Protoc 2024 (protocol; optional). 15 Comparison of CCC methods (Nat Commun 2022; optional).
16 Bandesh 2025/2026 (methods, once you are on HPAP). 17 Camunas-Soler patch-seq (advanced). 18 Human adipose integrated map (adipose context).
Reviews from 2022-2026 on single-cell transcriptomics in diabetes: NOT yet found; search PubMed for "single-cell diabetes review" before citing any.

---
## Deliverable 4: end-to-end workflow

Teaching format for each step: what, why, assumptions, failure modes, validation, expected plot, how to say it to your PI.

### 4.1 Data and metadata
- What: one donor table (ID, diagnosis, sex, age, BMI, batch, isolation site, culture time) and one count matrix per dataset.
- Why: all inference is donor-level; metadata drives the model.
- Fail: batch perfectly aligned with diagnosis (then disease and batch cannot be separated).
- Validate: cross-tab diagnosis x batch x sex.
- Plot: donor-level bars of cell counts per donor.
- PI line: "Here are the donors; this is what we can and cannot adjust for."

### 4.2 QC
- What: remove empty droplets, dying cells, doublets; correct ambient RNA in droplet data.
- Why: low-quality cells form fake clusters; ambient INS/GCG makes alpha cells look INS-positive and creates false DEGs (ambient-RNA abstracts in `literature_reading_list.md` section B).
- Assumptions: thresholds are tissue-specific; do not copy PBMC cutoffs. Check MT fraction per donor and cell type.
- Fail: filtering with one global threshold that removes a whole donor or cell type.
- Validate: per-donor distributions before/after; cells removed per donor.
- Plot: violin of genes/UMIs/MT% per donor.
- PI line: "QC removed X% of cells, evenly across diagnosis."

### 4.3 Normalization
- What: put cells on a comparable scale. Standard: divide by total counts, scale, log1p. Alternatives: Pearson residuals; for Smart-seq2, length-normalized values.
- Why: cells differ in capture efficiency and depth, not biology.
- Assumption: most genes are unchanged and cell types have similar total RNA. Violated when a cell type has much more RNA.
- Fail: normalized values used for pseudobulk testing (use raw counts there).
- Validate: library size no longer drives PC1.
- PI line: "Normalization removes depth, not biology."

### 4.4 HVG, PCA, neighbours, clustering, UMAP
- HVG: keep genes with variance beyond expectation (about 2000), reducing noise.
- PCA: compress to 30-50 components; check an elbow plot and whether PCs track donor rather than biology.
- kNN graph: each cell links to its k nearest neighbours in PC space (try k 15-30).
- Leiden: finds densely connected groups on that graph. Resolution sets granularity: higher means more, smaller clusters.
- UMAP: 2D picture for display only. Distances and cluster sizes in UMAP are not quantitative; never test on UMAP coordinates.
- Fail: choosing resolution until you see a "T2D cluster" (circular).
- Validate: run several resolutions; check cluster stability, marker specificity, and that clusters contain multiple donors.
- Plot: UMAP coloured by cell type, donor, diagnosis side by side.
- PI line: "Clusters are hypotheses until markers and replication support them."

### 4.5 Annotation
- Use canonical markers (INS beta, GCG alpha, SST delta, PPY PP, PRSS1 acinar, KRT19 ductal, plus endothelial/stellate/immune) and the islet identity genesets paper; optionally reference mapping.
- Validate: dotplot; hormone ambient check; compare proportions with the original paper.

### 4.6 Cell-type abundance
- Unit: donor. Compute proportions per donor, then use scCODA, Milo (donor as replicate) or sccomp.
- Why: proportions sum to 1, so one type rising makes others fall.
- Fail: chi-square test on cell counts (pseudoreplication).
- Plot: per-donor proportions as points with group means.

### 4.7 Differential expression: Control vs T2D
- Default: pseudobulk. For each cell type, sum raw counts per donor; filter cell types with too few cells per donor and too few donors per group; model counts with edgeR (quasi-likelihood), DESeq2 or limma-voom; design: diagnosis + sex + age (+ batch/BMI if donors allow).
- Alternative: mixed model with donor random effect (e.g. NEBULA, glmmTMB, dreamlet). The Bioinformatics 2024 paper abstract says pseudobulk with proper offsets matches GLMMs; use mixed models mainly when donors contribute varying cell numbers and you want covariates handled per cell.
- Report: log2 fold change, raw p, FDR (Benjamini-Hochberg) per cell type, and decide how to handle multiplicity across cell types.
- Fail: p-values from per-cell Wilcoxon; treating UMAP separation as evidence; unadjusted covariates.
- Validate: p-value histogram is flat plus a spike near 0; MA/volcano plots; leave-one-donor-out stability; permutation of donor labels.
- Elgamal found 84 T2D genes across cell types with 17 T2D vs 29 ND (VERIFIED), so expect modest effects, not hundreds.
- PI line: "Effective sample size is donors; here is the donor-level result."

### 4.8 Pathways and regulators
- Run GSEA/GO/Reactome/KEGG on the ranked pseudobulk statistics (not on a thresholded list when possible); decoupleR for TF and pathway activity.
- Fail: interpreting enrichment from a ranked list built on cell-level p-values.

### 4.9 Replication
- Repeat the identical pipeline in an independent cohort; report direction concordance and effect-size correlation of DEGs, not just overlap counts.
- Caveat: platforms differ (Smart-seq2/C1 vs droplet), so compare signs and ranks, not magnitudes. Small cohorts have little power, so lack of replication is weak evidence.

### 4.10 Advanced, later
- Trajectory/pseudotime, RNA velocity, integration across datasets, scATAC, spatial. Not needed for the first project.

---
## Deliverable 5: research question and aims

Question: How does type 2 diabetes alter cell-type composition and cell-type-specific transcription in human pancreatic islets, and are these changes reproducible across independent donor cohorts?

- Aim 1, cell-type identification: annotate islet cell types in HPAP using canonical markers; evaluate against published labels. Analysis: clustering stability, marker specificity.
- Aim 2, abundance: test donor-level proportion differences (scCODA/Milo/sccomp), adjusting for sex and age. Statistics: compositional model with donor as replicate; FDR across cell types.
- Aim 3, cell-specific DE: pseudobulk edgeR/limma-voom per cell type with covariates. Statistics: FDR per cell type; leave-one-donor-out sensitivity.
- Aim 4, pathways: GSEA/decoupleR on ranked pseudobulk statistics for beta, alpha, delta cells; check stress, inflammation, mitochondrial themes against literature.
- Aim 5, replication: repeat Aims 3-4 in Lawlor, Segerstolpe and Avrahami; test directional concordance (sign test, rank correlation).

Realism check: HPAP access and FASTQ processing may take weeks. A safe internship scope is Aims 1-3 on one dataset plus Aim 5 on a small one.

---
## Deliverable 6: portfolio structure

```
t2d-islet-sc/
  README.md              question, data, key results, how to run
  environment.yml        pinned versions (python, scanpy, edgeR via rpy2 or R env)
  data/
    README.md            source URL, accession, download date, checksums
    raw/                 never edited (git-ignored)
    metadata/donors.csv  donor table
  notebooks/
    01_qc.ipynb ... 08_replication.ipynb
  src/                   reusable functions (QC, pseudobulk)
  results/figures, results/tables
  reports/               2-page summary, slide deck for PI
  docs/                  literature notes, decisions log (why each threshold)
```
Practices: fixed random seeds, one notebook per step, a decisions log, no large data in git, a short methods paragraph per analysis.

---
## Deliverable 7: common mistakes (statistically invalid ones first)

Invalid:
1. Treating cells as independent replicates (pseudoreplication).
2. Using per-cell Wilcoxon p-values to claim Control vs T2D differences.
3. Permuting cell labels instead of donor labels.
4. Chi-square on cell counts for abundance.
5. Pseudobulk built from normalized/logged values instead of raw counts.
6. Batch or site perfectly confounded with diagnosis, and not stating it.
7. Clustering, then testing the clusters for diagnosis differences on the same data without acknowledging circularity (Cell Syst 2019 paper on post-clustering inference).
8. Integration that removes the disease signal, then claiming no effect.
9. Ignoring multiple testing across genes and cell types.
10. Reporting a cell type with 3 cells in a donor as a donor-level measurement.

Common technical:
11. Copying QC thresholds from other tissues.
12. Ignoring ambient hormone RNA in islets.
13. Reading UMAP distances and cluster sizes quantitatively.
14. Tuning resolution to get the cluster you want.
15. Annotating by one marker gene.
16. Ignoring sex, age, BMI and donor-level covariates.
17. Overinterpreting a p-value without effect size or replication.
18. Skipping donor-level plots.
19. Reporting pathway enrichment from an arbitrary thresholded list.
20. Mixing Smart-seq2 and droplet data without addressing technology effects.

---
## Deliverable 8: glossary (50 terms)

1 Transcriptome: all RNA in a sample. 2 Gene expression: RNA produced from a gene. 3 Transcript: an RNA copy of a gene. 4 Exon: retained coding segment. 5 Intron: removed segment (reads from nuclei often include introns). 6 mRNA. 7 Poly(A): tail used for capture in most protocols. 8 Read: one sequenced fragment. 9 FASTQ: raw read file with quality. 10 Alignment: mapping reads to a genome. 11 UMI: tag that counts original molecules. 12 Cell barcode: tag identifying the cell. 13 Count matrix: genes x cells table of counts. 14 Library size: total counts in a cell. 15 Sequencing depth: reads per cell. 16 Dropout: gene present but not detected. 17 Sparsity: fraction of zeros. 18 scRNA-seq: whole cells. 19 snRNA-seq: nuclei; used for tissues like adipose. 20 Smart-seq2: full-length plate protocol. 21 10x Chromium: droplet protocol. 22 Droplet: oil bubble holding one cell. 23 Empty droplet: no cell, ambient RNA only. 24 Doublet: two cells in one barcode. 25 Ambient RNA: free RNA contaminating droplets. 26 QC metrics: genes, UMIs, MT fraction. 27 Mitochondrial genes: MT- genes; high fraction suggests damage. 28 Ribosomal genes. 29 Normalization. 30 Log1p. 31 HVG: highly variable genes. 32 PCA. 33 kNN graph. 34 Leiden/Louvain: community detection. 35 Resolution. 36 UMAP/t-SNE: 2D visualization. 37 Marker gene. 38 Cell type annotation. 39 Batch effect: technical difference between runs. 40 Donor effect: differences between individuals. 41 Biological replicate: independent donor/animal. 42 Pseudoreplication. 43 Pseudobulk: summed counts per donor per cell type. 44 Differential expression. 45 Log2 fold change. 46 p-value. 47 FDR / adjusted p-value (Benjamini-Hochberg). 48 Mixed-effects model: includes donor random effect. 49 Differential abundance / compositional data. 50 GSEA/gene set enrichment, with TF activity and cell-cell communication as the next-layer interpretations.

Extra for T2D biology (not counted): beta cell, alpha cell, delta cell, PP cell, insulin, glucagon, somatostatin, insulin resistance, beta-cell dedifferentiation, ER stress, oxidative stress, islet, HbA1c.

---
## Open items
- Verify GEO metadata for GSE86473, GSE154126; E-MTAB-5061; HPAP access path.
- Find 2022-2026 diabetes single-cell reviews.
- Find a T2D-labelled adipose or other metabolic dataset.
- Confirm pydeseq2/NEBULA/dreamlet specifics from official documentation before using.
