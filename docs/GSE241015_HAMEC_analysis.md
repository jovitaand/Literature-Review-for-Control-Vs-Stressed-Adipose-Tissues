# GSE241015 HAMEC file: analysis (2 Nov 2026)

File: `GSE241015_HAMEC_EndoMT_RNAseq.csv.gz` (read in full). Bulk RNA-seq of **cultured human adipose microvascular endothelial cells (HAMEC)**: 3 untreated controls vs 3 treated with TGF-b/TNF-a/IFN-g for 4 days ("EndoMT"). 15,083 genes; values are TPM-like (each sample sums to about 0.9 to 1.0 million); the file also holds p-value, q-value, fold change and log2 fold change.
Outputs: `docs/GSE241015_HAMEC_analysis.xlsx` (summary, gene list, marker validation, top 50 up/down, sample correlation) and `docs/GSE241015_HAMEC_analysis.png`.

## Data quality
- No missing values, no duplicated genes. 9 genes have `#DIV/0!` in the fold-change columns (both groups near zero; a spreadsheet artefact, not a result).
- Replicates agree: within-group Spearman correlation 0.990 to 0.993; between-group 0.971 to 0.976. PC1 separates treated from control and explains 67% of variance.
- The file's log2 fold change agrees with one recomputed from the group means (r = 0.93), and its p-values rank like a Welch test on log2(TPM + 1) (Spearman 0.81).
- 2,390 genes have q < 0.05 (1,288 up, 1,093 down); 550 also have |log2 fold change| > 1.

## Your gene list
| Gene | Control TPM | EndoMT TPM | log2 FC | q | Level | Reading |
|---|---|---|---|---|---|---|
| PECAM1 | 2061 | 1947 | -0.08 | 0.13 | high | no clear change |
| CDH5 | 675 | 529 | -0.35 | 0.002 | high | lower |
| NOS3 | 5.6 | 6.9 | +0.30 | 0.18 | very low | no clear change; barely expressed |
| HIF1A | 276 | 284 | +0.04 | 0.17 | moderate | no clear change |
| ERG | 120 | 127 | +0.08 | 0.078 | moderate | no clear change (borderline) |
| VCAM1 | 8.6 | 4.1 | -1.08 | 0.0002 | very low | lower |
| ICAM1 | 231 | 134 | -0.78 | 0.0005 | moderate | lower |
| ANGPT2 | 427 | 296 | -0.53 | 0.0003 | moderate | lower |
Four genes are significantly lower after treatment (CDH5, VCAM1, ICAM1, ANGPT2); four show no clear change (PECAM1, NOS3, HIF1A, ERG). All 8 are detected in this cell line (NOS3 and VCAM1 at very low levels).

## Did the treatment behave like "EndoMT"?
- **Interferon response: strong.** IDO1 (+8.6 log2), CXCL10 (+7.3), CXCL11 (+8.1), GBP1 (+4.2), STAT1 (+2.9), IRF1 (+2.8) and MHC class II genes (HLA-DRB1 +12.7, CD74 +9.2) rise massively.
- **Classic mesenchymal genes: not up.** ACTA2 (-0.08), TAGLN (-0.07), FN1 (-0.02), SNAI1 (+0.16) are unchanged; only SNAI2 rises (+0.9); VIM, SERPINE1, CCN2 and MMP2 fall.
- **Endothelial identity genes: mostly kept.** PECAM1, CLDN5, TEK are unchanged; CDH5 and VWF drop modestly; FLT1, KDR and EGFL7 rise slightly.
- **Early TNF-response genes are low at day 4.** SELE falls 5.5 log2, CCL2 and CXCL8 fall about 1.3 log2, and ICAM1 and VCAM1 fall, although TNF-a is in the treatment. This fits a late, post-peak time point but was not tested here.
So the transcriptional change is dominated by the interferon response, with a weak mesenchymal signature at the RNA level.

## Against the problem statement
This file cannot meet it: cultured cells, 3 vs 3, no donors, no diabetes, no depot, no cell counts. The comparison is **control vs cytokine-induced EndoMT, not control vs type 2 diabetes**, so the gene-list changes above must not be read as diabetes effects.
