# Reading list: Control vs T2D single-cell transcriptomics

Every entry below was returned by a literature index search (title, DOI, abstract) on 2026-10-05.
I did NOT read full texts except where noted in `datasets_verified_so_far.md`. Takeaways are quoted or
paraphrased from abstracts only. Verify details before citing. "Order" is a suggested reading order.

## A. Foundational: workflow, QC, normalization

| Order | Paper | DOI | Why read it |
|---|---|---|---|
| 1 | Current best practices in single-cell RNA-seq analysis: a tutorial (Mol Syst Biol 2019) | 10.15252/msb.20188746 | Start here. End-to-end workflow and the decisions at each step. |
| 2 | Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data (Nat Protoc 2021) | 10.1038/s41596-020-00409-w | Second pass with a protocol focus. |
| 3 | Bioconductor workflow for single-cell RNA sequencing: normalization, dimensionality reduction, clustering, lineage inference (F1000Res 2017) | 10.12688/f1000research.12122.1 | R/Bioconductor route if you choose R over Python. |
| 4 | Analysis and visualization of single-cell sequencing data with Scanpy and MetaCell: a tutorial (Methods Mol Biol 2024) | 10.1007/978-1-0716-3642-8_17 | Python/Scanpy route, matches your background. |

## B. Quality control: ambient RNA and doublets

| Paper | DOI | Why read it |
|---|---|---|
| Decontamination of ambient RNA in single-cell RNA-seq with DecontX (Genome Biol 2020) | 10.1186/s13059-020-1950-6 | Ambient RNA method. |
| Understanding and mitigating the impact of ambient mRNA contamination (PLoS One 2025) | 10.1371/journal.pone.0332440 | Abstract states ambient transcripts appear among DEGs; compares CellBender and SoupX. Directly relevant to islets (insulin/glucagon spillover). |
| FastCAR: fast correction for ambient RNA to facilitate differential expression (BMC Genomics 2023) | 10.1186/s12864-023-09822-3 | Abstract: ambient RNA levels are highly sample-specific, which matters for case vs control. |
| Souporcell: clustering by genotype without reference genotypes (Nat Methods 2020) | 10.1038/s41592-020-0820-1 | Donor demultiplexing and cross-genotype doublets. |

## C. Integration and batch correction (needed because donors are batches)

| Paper | DOI | Takeaway from abstract |
|---|---|---|
| A benchmark of batch-effect correction methods for scRNA-seq (Genome Biol 2020) | 10.1186/s13059-019-1850-9 | Recommends Harmony, LIGER, Seurat 3. |
| Benchmarking atlas-level data integration in single-cell genomics (Nat Methods 2022) | 10.1038/s41592-021-01336-8 | BBKNN, Scanorama, scVI do well on complex tasks. |
| Batch correction methods ... are often poorly calibrated (Genome Res 2025) | 10.1101/gr.279886.124 | Recommends only Harmony. Read next to the two above; they disagree. |
| Probabilistic harmonization and annotation with deep generative models, scVI/scANVI (Mol Syst Biol 2021) | 10.15252/msb.20209620 | scVI/scANVI model reference. |

## D. Case-control statistics (the most important section for your design)

| Paper | DOI | Takeaway from abstract |
|---|---|---|
| Confronting false discoveries in single-cell differential expression (Nat Commun 2021) | 10.1038/s41467-021-25960-2 | Why cell-level tests inflate false positives. Read first. |
| Differential gene expression analysis for multi-subject scRNA-seq with aggregateBioVar (Bioinformatics 2021) | 10.1093/bioinformatics/btab337 | Naive testing gives many false discoveries; pseudobulk has better FDR control. |
| A balanced measure shows superior performance of pseudobulk methods over mixed models and pseudoreplication (Nat Commun 2022) | 10.1038/s41467-022-35519-4 | Supports pseudobulk. |
| Benchmarking methods for detecting differential states between conditions from multi-subject scRNA-seq (Brief Bioinform 2022) | 10.1093/bib/bbac286 | Abstract: pseudo-bulk methods performed generally best. |
| Pseudobulk with proper offsets has the same statistical properties as GLMMs (Bioinformatics 2024) | 10.1093/bioinformatics/btae498 | Counterpoint on mixed models. |
| Valid post-clustering differential analysis for scRNA-seq (Cell Syst 2019) | 10.1016/j.cels.2019.07.012 | Double-dipping after clustering. |

## E. T2D single-cell papers (islets)

| Paper | DOI | Role |
|---|---|---|
| Xin et al., RNA sequencing of single human islet cells reveals T2D genes (Cell Metab 2016) | 10.1016/j.cmet.2016.08.018 | 1,492 cells, ND vs T2D donors (abstract). |
| Segerstolpe et al., Single-cell transcriptome profiling of human islets in health and T2D (Cell Metab 2016) | 10.1016/j.cmet.2016.08.020 | 6 healthy, 4 T2D, 2,209 cells (read in full). |
| Lawlor et al., Single-cell transcriptomes identify human islet cell signatures ... in T2D (Genome Res 2017) | 10.1101/gr.212720.116 | 5 ND, 3 T2D, 1,050 cells (read in full). |
| Avrahami et al., Single-cell transcriptomics of human islet ontogeny ... beta-cell dedifferentiation in T2D (Mol Metab 2020) | 10.1016/j.molmet.2020.101057 | 22 donors incl. 10 T2D (read in full). |
| Camunas-Soler et al., Pancreas patch-seq links physiologic dysfunction in diabetes to transcriptomic phenotypes (Cell Metab 2020) | 10.1016/j.cmet.2020.04.005 | Links function to expression; advanced. |
| Elgamal et al., An integrated map of cell type-specific gene expression in pancreatic islets (Diabetes 2023) | 10.2337/db23-0130 | 192,203 cells, 17 T2D vs 29 ND (read in full). Best model for your own pipeline. |
| Bandesh et al., Single-cell decoding of human islet cell type-specific alterations in T2D (bioRxiv 2025; EMBO J 2026) | 10.1101/2025.01.17.633590 ; 10.1038/s44318-026-00744-w | Recent, large; read methods. |
| Single-cell mRNA-regulation analysis reveals cell type-specific mechanisms of T2D (Nat Commun 2025) | 10.1038/s41467-025-65060-z | Re-analysis of public islet scRNA-seq. |
| Human pancreatic alpha-cell heterogeneity ... SMOC1 as a beta-cell dedifferentiation gene (Nat Commun 2025) | 10.1038/s41467-025-62670-5 | Combines scRNA-seq and snRNA-seq. |
| Meta-analysis of single-cell expression of human alpha and beta cells in T2D (J Diabetes 2021) | 10.1111/1753-0407.13236 | Example of combining the small public datasets. |

## F. Adipose tissue (relevant to this repo)

| Paper | DOI | Note |
|---|---|---|
| Single-cell analysis of human adipose tissue identifies depot and disease specific cell types (Nat Metab 2020) | 10.1038/s42255-019-0152-6 | scRNA-seq of SVF from obese individuals. |
| An integrated single cell and spatial transcriptomic map of human white adipose tissue (Nat Commun 2023) | 10.1038/s41467-023-36983-2 | Integrates ten studies; >60 subpopulations; a reference for annotation. |
| Human subcutaneous and visceral adipocyte atlases (Nat Genet 2025) | 10.1038/s41588-024-02048-3 | snRNA-seq depot atlases. |
| Spatial mapping reveals human adipocyte subpopulations with distinct sensitivities to insulin (Cell Metab 2021) | 10.1016/j.cmet.2021.07.018 | Insulin sensitivity link. |
| Integrative single-cell analysis of metabolic syndrome in adipose tissue (Diabetol Metab Syndr 2025) | 10.1186/s13098-025-01975-3 | Metabolic-syndrome snRNA-seq reanalysis. |
| Robust snRNA-seq reveals depot-specific dynamics in adipose remodeling during obesity (eLife 2025) | 10.7554/eLife.97981 | Mouse-focused nuclei isolation method; check species before use. |
| Causal effect of obesity on prediabetes and insulin resistance (PLoS Genet 2020) | 10.1371/journal.pgen.1009018 | Bulk plus snRNA-seq, insulin resistance. |

## G. Atlases and portals

| Paper | DOI |
|---|---|
| The Tabula Sapiens (Science 2022) | 10.1126/science.abl4896 |
| Tabula Sapiens 2.0 (bioRxiv 2024) | 10.1101/2024.12.03.626516 |

## Gaps (not yet found or not yet read)
- Experimental design and power for scRNA-seq case-control; sequencing technology comparisons (Smart-seq2 vs 10x vs snRNA-seq).
- Cell type annotation tools; cell-cell communication; pathway/GSEA on pseudobulk; trajectory inference.
- Your remaining items 3 to 5 (message was cut off).
