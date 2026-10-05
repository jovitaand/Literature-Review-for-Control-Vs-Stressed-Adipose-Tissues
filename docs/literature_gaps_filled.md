# Reading list, part 2: filling the gaps

Companion to `literature_reading_list.md`. Same caveat: every entry was returned by a literature-index search
(title, DOI, abstract) on 2026-10-05; I did not read full texts. Takeaways come from abstracts only.

## H. Experimental design and power

| Paper | DOI | Why read it |
|---|---|---|
| Design and power analysis for multi-sample single cell genomics experiments (Nat Commun 2021) | 10.1038/s41467-021-26779-7 | Abstract: models sample size, cells per individual and depth. Core design paper. |
| Hierarchicell: power for differential expression with single-cell data (BMC Genomics 2021) | 10.1186/s12864-021-07635-w | Abstract: existing calculators focus on total cells, not independent units (donors), and overestimate power. |
| Experimental design for single-cell RNA sequencing (Brief Funct Genomics 2018) | 10.1093/bfgp/elx035 | General design considerations. |
| Quantifying the tradeoff between sequencing depth and cell number in scRNA-seq (bioRxiv) | 10.1101/762773 | Budget allocation: cells vs depth. |
| Optimal design of scRNA-seq experiments for cell-type-specific eQTL analysis (Nat Commun 2020) | 10.1038/s41467-020-19365-w | Abstract: different sample/cell/read designs can give similar power. |
| Robustness of scRNA-seq for identifying differentially expressed genes (BMC Genomics 2023) | 10.1186/s12864-023-09487-y | Small cell numbers per cluster. |

## I. Technology comparison (matters for pooling Smart-seq2, C1, 10x, snRNA-seq datasets)

| Paper | DOI | Takeaway from abstract |
|---|---|---|
| Direct comparative analysis of 10X Genomics Chromium and Smart-seq2 (Genomics Proteomics Bioinformatics 2021) | 10.1016/j.gpb.2020.02.005 | Smart-seq2 detected more genes. |
| Comparative analysis of single-cell RNA sequencing methods (Mol Cell 2017) | 10.1016/j.molcel.2017.01.023 | Smart-seq2 most sensitive per cell. |
| Full-length RNA-seq from single cells using Smart-seq2 (Nat Protoc 2014) | 10.1038/nprot.2014.006 | The protocol used by Segerstolpe. |
| Benchmarking single-cell mRNA-seq technologies ... cell types with low RNA content (2022) | 10.7171/3fc1f5fe.dbeabb2a | Sensitivity in low-RNA cell types. |

Gap: I did not find a dedicated single-nucleus vs single-cell comparison for adipose or islets in this search.

## J. Cell type annotation

| Paper | DOI | Why read it |
|---|---|---|
| Automated methods for cell type annotation on scRNA-seq data (Comput Struct Biotechnol J 2021) | 10.1016/j.csbj.2021.01.015 | Overview of approaches. |
| Evaluation of cell type annotation R packages on scRNA-seq data (Genomics Proteomics Bioinformatics 2021) | 10.1016/j.gpb.2020.07.004 | Benchmark of supervised tools. |
| Generation of human islet cell type-specific identity genesets (Nat Commun 2022) | 10.1038/s41467-022-29588-8 | Abstract: meta-analysis of scRNA-seq to define alpha/beta/gamma/delta marker sets. Use for islet annotation. |
| Transcriptomes of the major human pancreatic cell types (Diabetologia 2011) | 10.1007/s00125-011-2283-5 | Older bulk reference for alpha, beta, duct, acinar. |
| Integrated single cell and spatial map of human white adipose tissue (Nat Commun 2023) | 10.1038/s41467-023-36983-2 | Adipose annotation reference (also in section F). |

Skipped: LLM-based annotation benchmark (Brief Bioinform 2025, 10.1093/bib/bbaf622). It exists, but I would not rely on it as a beginner.

## K. Compositional / differential abundance (T2D may change cell-type proportions)

| Paper | DOI | Takeaway from abstract |
|---|---|---|
| scCODA: a Bayesian model for compositional single-cell data analysis (Nat Commun 2021) | 10.1038/s41467-021-27150-6 | Handles compositionality and low sample size. |
| Milo: differential abundance testing using k-NN graphs (Nat Biotechnol 2022) | 10.1038/s41587-021-01033-z | Neighbourhood-level DA without hard clusters. |
| Benchmarking differential abundance methods (Genome Biol 2024) | 10.1186/s13059-023-03143-0 | Compares DA methods. |
| Robust differential composition and variability analysis, sccomp (bioRxiv) | 10.1101/2022.03.04.482758 | Alternative composition model. |
| Best practices for single-cell analysis across modalities (Nat Rev Genet 2023) | 10.1038/s41576-023-00586-w | Abstract: summarises benchmarks into best-practice workflows. Good umbrella paper to read alongside section A. |

## L. Pathways and regulator activity

| Paper | DOI | Why read it |
|---|---|---|
| decoupleR: ensemble of methods to infer biological activities from omics data (Bioinform Adv 2022) | 10.1093/bioadv/vbac016 | One framework for pathway/TF activity; works on pseudobulk output. |
| PROGENy: perturbation-response genes reveal signaling footprints (Nat Commun 2018) | 10.1038/s41467-017-02391-6 | Pathway activity from footprint genes. |
| Transfer of regulatory knowledge from human to mouse (BBA Gene Regul Mech 2019) | 10.1016/j.bbagrm.2019.194431 | DoRothEA and PROGENy benchmark. |

Gap: search for single-cell-specific GSEA tools returned niche packages only. For T2D, running GSEA/decoupleR on pseudobulk DE results (per section D) is the safer route; I have no benchmark paper to cite for that choice.

## M. Cell-cell communication

| Paper | DOI | Takeaway from abstract |
|---|---|---|
| Inference and analysis of cell-cell communication using CellChat (Nat Commun 2021) | 10.1038/s41467-021-21246-9 | The CellChat method. |
| CellChat for systematic analysis of cell-cell communication ... (Nat Protoc 2024) | 10.1038/s41596-024-01045-4 | Protocol, good for a beginner. |
| Comparison of resources and methods to infer cell-cell communication from single-cell RNA data (Nat Commun 2022) | 10.1038/s41467-022-30755-0 | Compares resources and methods (LIANA). Read before choosing a tool. |
| ESICCC: evaluation, selection and integration of CCC inference methods (Genome Res 2023) | 10.1101/gr.278001.123 | Abstract: benchmarks 18 ligand-receptor methods. |
| Comparative analysis of cell-cell communication at single-cell resolution (Nat Biotechnol 2023) | 10.1038/s41587-023-01782-z | Condition-comparison angle, useful for control vs T2D. |

## N. Trajectory / dedifferentiation

| Paper | DOI | Why read it |
|---|---|---|
| A comparison of single-cell trajectory inference methods (Nat Biotechnol 2019) | 10.1038/s41587-019-0071-9 | Abstract: compares 29 methods. |
| Cell-state trajectories in T2D beta cells | see Avrahami 2020 and the SMOC1 alpha-cell paper in `literature_reading_list.md` section E | Applied examples of dedifferentiation analysis. |

Index quirk: one search returned a record with this benchmark's abstract attached to an unrelated BMJ title
(doi 10.1136/bmj.325.7378.1449). Ignore it; it is a metadata error in the index.
Trajectory methods are optional for a first Control vs T2D analysis; treat as advanced.

## O. Genetics link (T2D GWAS to cell types)

| Paper | DOI | Why read it |
|---|---|---|
| Single cell chromatin accessibility reveals pancreatic islet cell type- and state-specific regulatory programs of diabetes risk (Nat Genet 2021) | 10.1038/s41588-021-00823-0 | Abstract: T2D GWAS enrichment in beta-cell states. |
| 3D chromatin maps of the human pancreas reveal lineage-specific regulatory architecture of T2D risk (Cell Metab 2022) | 10.1016/j.cmet.2022.08.014 | Variant-to-gene links. |
| Cell-type-specific cis-eQTLs in pancreatic cell types identify novel risk genes for T2D (Brief Bioinform 2025) | 10.1093/bib/bbaf531 | Abstract: 328 cell-type-specific cis-eQTLs. |
| EPIC: inferring relevant cell types for complex traits from GWAS and scRNA-seq (PLoS Genet 2022) | 10.1371/journal.pgen.1010251 | Method to link GWAS to cell types. |
| Implicating T2D effector genes with promoter-focused Capture-C (Diabetologia 2024) | 10.1007/s00125-024-06261-x | Metabolic cell models. |

## Still open
- Single-nucleus vs single-cell comparisons specific to adipose or islets.
- A benchmark specifically supporting pathway analysis choice on single-cell/pseudobulk data.
- Your original items 3 to 5 (message was cut off before them).
