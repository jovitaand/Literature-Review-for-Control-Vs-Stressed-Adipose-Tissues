# Mentor plan, part 1: datasets, papers, 12-week roadmap, first analysis

Status 2026-10-05. Verification level is marked on every claim: VERIFIED = read in paper text; ABSTRACT = abstract only;
UNVERIFIED = memory, must be checked. GEO/NCBI/EBI were unreachable from my sandbox, so no GEO page has been read directly.

## 1. Five datasets for Control vs T2D (human islets)

| Rank | Dataset | Donors (ND / T2D) | Cells | Tech | Access | Best for |
|---|---|---|---|---|---|---|
| 1 | Elgamal 2023, HPAP (doi:10.2337/db23-0130) | 65 in final map; 29 ND, 17 T2D of 67 downloaded (also 10 T1D, 9 ND Aab+) | 192,203 | HPAP scRNA-seq (chemistry not confirmed) | Raw FASTQ via PANC-DB (data-access terms apply). Authors' processed object: UNVERIFIED | Main project |
| 2 | Lawlor 2017 (doi:10.1101/gr.212720.116) | 8: 5 / 3 | 1,050 (622 ND, 428 T2D) | Fluidigm C1 | SRA SRP075970; processed GEO GSE86473 | Learning |
| 3 | Segerstolpe 2016 (doi:10.1016/j.cmet.2016.08.020) | 10: 6 healthy / 4 T2D | 2,209 | Smart-seq2 | Portal sandberg.cmb.ki.se/pancreas; E-MTAB-5061 UNVERIFIED | Validation, learning |
| 4 | Avrahami 2020 (doi:10.1016/j.molmet.2020.101057) | 22 total: 4 adult ND, 10 adult T2D, 8 younger ND | not retrieved | Fluidigm C1 | GEO GSE154126 | Validation of beta-cell dedifferentiation |
| 5 | Bandesh 2025/2026 (doi:10.1101/2025.01.17.633590; 10.1038/s44318-026-00744-w) | not retrieved | not retrieved | not retrieved | not retrieved | Advanced; read methods first |

Why not the others:
- Xin 2016 (1,492 cells, GSE81608 UNVERIFIED): donor counts not retrieved; small and older. Fine as extra validation once checked.
- No human adipose dataset with T2D vs control labels is verified. Adipose leads in `literature_reading_list.md` section F are obesity/insulin-resistance studies, not T2D case-control. Do not use adipose as the main tissue until a labelled dataset is confirmed.
- Donor reality check: Lawlor (3 T2D) and Segerstolpe (4 T2D) are too small for donor-level inference. Use them to learn the pipeline and for exploratory replication, not as evidence.
- Cell counts per donor are heavily uneven in these plate-based sets; check this before any proportion analysis.

## 2. Ten papers to read first (in order)

| # | Paper | DOI | Read | Learn |
|---|---|---|---|---|
| 1 | Current best practices in scRNA-seq analysis: a tutorial (Mol Syst Biol 2019) | 10.15252/msb.20188746 | Whole | The full pipeline and its decisions. |
| 2 | Best practices for single-cell analysis across modalities (Nat Rev Genet 2023) | 10.1038/s41576-023-00586-w | Skim, then reference | Benchmark-backed defaults. |
| 3 | Confronting false discoveries in single-cell differential expression (Nat Commun 2021) | 10.1038/s41467-021-25960-2 | Whole | Why cell-level tests fail; pseudobulk. |
| 4 | Benchmarking methods for differential states from multi-subject scRNA-seq (Brief Bioinform 2022) | 10.1093/bib/bbac286 | Whole | Abstract: pseudobulk best overall. |
| 5 | Design and power analysis for multi-sample single cell genomics (Nat Commun 2021) | 10.1038/s41467-021-26779-7 | Methods and figures | Donors, not cells, drive power. |
| 6 | Elgamal 2023 integrated islet map (Diabetes) | 10.2337/db23-0130 | Whole; copy the pipeline | A donor-aware T2D workflow. |
| 7 | Lawlor 2017 (Genome Res) | 10.1101/gr.212720.116 | Results and methods | Your first dataset. |
| 8 | Segerstolpe 2016 (Cell Metab) | 10.1016/j.cmet.2016.08.020 | Results | Second dataset; compare findings. |
| 9 | Benchmarking atlas-level data integration (Nat Methods 2022) | 10.1038/s41592-021-01336-8 | Abstract, figures | When integration is needed. Read with the Genome Res 2025 critique (10.1101/gr.279886.124). |
| 10 | Milo (Nat Biotechnol 2022) and scCODA (Nat Commun 2021) | 10.1038/s41587-021-01033-z ; 10.1038/s41467-021-27150-6 | Intro and methods overview | Differential abundance. |

Reviews from 2022-2026 on single-cell transcriptomics in diabetes: NOT yet found. I have not searched for reviews specifically; that is the next literature task.

## 3. Statistical stance (non-negotiable)

10 control + 10 T2D donors with 20,000 cells is n = 20, not 20,000. Cells from one donor share genotype, age, BMI,
isolation batch, and ambient RNA, so they are nested within donors and correlated.
- Defensible for publication: pseudobulk (sum counts per donor per cell type) with edgeR, DESeq2 or limma-voom, with donor as the unit; or mixed models with a donor random effect, where the literature (see reading list section D) argues they perform comparably if specified well. Include covariates (sex, age, BMI, batch) when donor numbers allow.
- Not defensible: Wilcoxon or t-test across cells, reporting cell-level p-values as evidence of disease effect. Fine for marker finding within one dataset, not for Control vs T2D.
- Permutation tests: permute donor labels, never cell labels.
- Abundance: use donor-level proportions with a compositional method (scCODA, Milo, sccomp); proportions sum to 1, so one cell type rising makes others fall.

## 4. 12-week roadmap

| Wk | Concepts | Read | Do | Output |
|---|---|---|---|---|
| 1 | DNA/RNA, transcription, exons/introns, RNA-seq, UMI, count matrix, MT/ribosomal genes | Paper 1 sections on experimental setup | Open a count matrix in pandas, compute library size, genes per cell | One-page glossary in your own words |
| 2 | Islet biology: alpha/beta/delta, insulin, glucose handling, T2D: beta-cell dysfunction, ER and oxidative stress, inflammation | Lawlor intro/discussion; Avrahami abstract | List marker genes per islet cell type (INS, GCG, SST, PPY, PRSS1, KRT19) | Marker table with sources |
| 3 | Cell vs nucleus, Smart-seq2 vs 10x, dropout, sparsity, depth, batch and donor effects | Paper 1; technology comparison papers (reading list part 2, section I) | Load Lawlor into AnnData; explore sparsity | Plots: counts/cell, genes/cell, MT fraction |
| 4 | QC, normalization, HVGs | Paper 1 QC/normalization sections | Run QC filters; justify thresholds per donor | Before/after QC figures, written thresholds |
| 5 | PCA, kNN graph, Leiden resolution, UMAP | Paper 1 clustering section | Cluster at several resolutions; check stability | Clustering notebook; resolution sweep |
| 6 | Marker genes, annotation | Islet identity genesets paper (10.1038/s41467-022-29588-8) | Annotate cell types; check for doublets and ambient hormone signal | Annotated UMAP, dotplot |
| 7 | Pseudoreplication, pseudobulk | Papers 3, 4 | Pseudobulk by donor x cell type; edgeR/limma or pydeseq2; compare with a naive cell-level test | Volcano plots, inflation comparison |
| 8 | Power, abundance | Papers 5, 10 | Donor-level proportions; scCODA or Milo | Abundance plot with donor points |
| 9 | Multiple testing, GSEA/GO/Reactome, decoupleR | decoupleR (10.1093/bioadv/vbac016) | Pathways on pseudobulk results | Enrichment figure, table |
| 10 | Replication | Segerstolpe paper | Repeat steps on Segerstolpe; compare direction of effects | Concordance plot, short report |
| 11 | Batch/integration, HPAP structure | Paper 9 and Genome Res critique | Request PANC-DB access; read Elgamal methods | Access request, data-plan note |
| 12 | Project design | Paper 6 | Write aims and analysis plan | 3-page proposal for your PI |

## 5. Recommended first dataset and analysis

Dataset: Lawlor 2017 (GSE86473 processed; SRA SRP075970 raw). Reason: small, both labelled groups, plate-based so no
droplet ambient-RNA complication, and the accession was read in the paper itself. Weakness: only 3 T2D donors.

First analysis:
1. Load the processed matrix; confirm the donor and diagnosis labels from GEO sample metadata (I could not read them).
2. Per-donor QC, then normalize, HVGs, PCA, kNN, Leiden, UMAP.
3. Annotate alpha, beta, delta, PP, acinar, ductal; color UMAP by donor and diagnosis to see donor effects.
4. Pseudobulk by donor x cell type; test beta cells T2D vs ND with 5 vs 3 donors. Treat as exploratory.
5. Also run the naive cell-level test and compare; the gap is the lesson on pseudoreplication.

Tool choice: start with Scanpy/AnnData (you already know Python), and learn R only for the pseudobulk step if needed
(edgeR/limma-voom or muscat are the field standard; I believe pydeseq2 is a Python alternative but I have not verified it
this session). Learn Seurat later, mainly to read other people's code.

## Open items
- Verify GEO metadata for GSE86473, GSE154126 and E-MTAB-5061, and the HPAP access path.
- Search for 2022-2026 reviews on single-cell transcriptomics in diabetes.
- Find a T2D-labelled adipose or other metabolic-tissue dataset.
- Remaining deliverables: top-10 dataset list, full workflow, aims, portfolio structure, mistakes list, 50-term glossary, step-by-step teaching notes.
