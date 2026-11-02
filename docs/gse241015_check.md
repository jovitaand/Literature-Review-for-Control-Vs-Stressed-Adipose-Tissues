# GSE241015 (Chavkin et al., obesity-induced EndoMT): criteria check from the two shared files (2 Nov 2026)

Two Drive files were shared and **both were read in full** (small enough for the connector):
- `GSE241015_HAMEC_EndoMT_RNAseq.csv.gz` (968,524 bytes): 15,083 genes x 3 control + 3 EndoMT samples (human gene symbols).
- `GSE241015_Adipose_EndoMT_RNAseq.csv.gz` (456,278 bytes): 8,933 genes x 3 EC + 3 EndoMT samples (mouse gene symbols such as Cstdc2, Defb42).

Both are **bulk RNA-seq result tables** (expression values per sample plus p-value, q-value and log2 fold change). They are not single-cell data, and there are no cells to count.

## What the samples are (GEO sample records, matching the file columns)
| File | Samples | Biology |
|---|---|---|
| HAMEC | 3 control + 3 treated | **cultured human adipose microvascular endothelial cells** (HAMEC cell line), untreated vs TGF-b/TNF-a/IFN-g for 4 days (induced EndoMT). No donors, no tissue depot, no disease status. |
| Adipose | 3 endothelial (CD31+eYFP+) + 3 EndoMT (CD31-eYFP+) | **mouse** visceral adipose tissue, lineage-traced Cdh5-CreERT2;Rosa26-eYFP mice on a high-fat/high-sucrose diet for 12 months |

## Requirement by requirement
| Requirement | Result | Evidence |
|---|---|---|
| Control group (lean or obese without diabetes) | **NOT MET** | HAMEC "control" = untreated cultured cells; the mouse comparison is endothelial vs EndoMT cells, both from high-fat-fed mice |
| Type 2 diabetes group | **NOT MET** | no diabetes in either file |
| Endothelial cells | **Met only as bulk / cultured cells** | HAMEC are endothelial cells in culture; the mouse EC fraction is sorted cells |
| More than 10,000 endothelial cells | **NOT APPLICABLE / NOT MET** | bulk RNA-seq, 3 + 3 samples per file; no cell counts |
| Human adipose tissue | **NOT MET for tissue** | human data are a cultured cell line; the tissue data are mouse |
| Both subcutaneous and visceral | **NOT MET** | mouse visceral fat only (the paper also studied mouse subcutaneous fat by flow cytometry, not in these files) |
| Gene list | Partly checked | see below |

## Gene list in the two files
**HAMEC (human, control vs EndoMT; log2 fold change, q-value):** all 8 genes present.
ANGPT2 -0.53 (q < 0.001); CDH5 -0.35 (q = 0.002); ICAM1 -0.78 (q < 0.001); VCAM1 -1.08 (q < 0.001); ERG +0.08 (q = 0.078); HIF1A +0.04 (q = 0.17); NOS3 +0.30 (q = 0.18, very low expression, about 5 to 8); PECAM1 -0.08 (q = 0.13).
**Adipose (mouse):** 5 of 8 genes present (Cdh5, Hif1a, Icam1, Pecam1, Vcam1); Nos3, Erg and Angpt2 are absent from the table. None is significant (q 0.52 to 1.0) with 3 vs 3 samples.
These are EndoMT (cytokine-induced) vs control comparisons, **not** type 2 diabetes vs control, so they cannot answer the project question.

**Verdict: does not meet the problem statement** (cultured human cells and mouse bulk RNA-seq; no cells to count; no diabetes; no depot comparison; no lean/obese human donors). The paper's human single-cell part reuses public visceral adipose datasets from obese individuals, which are not in these files.
