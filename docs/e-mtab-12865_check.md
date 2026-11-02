# E-MTAB-12865 (EMBL-EBI BioStudies): criteria check (2 Nov 2026)

Source: the BioStudies study page, its sample table (SDRF) and the file listing, plus the linked paper (Dhumale et al., Sci Rep 2023, "CD31 defines a subpopulation of human adipose-derived regenerative cells with potent angiogenic effects"). The count matrices (about 90 to 230 MB each) were not downloaded.

## The dataset
Single-cell RNA-seq (10x Genomics 3' v3, NovaSeq 6000) of **CD31-enriched adipose-derived regenerative cells** (CD31 magnetic sorting) from **4 donors**: Pt1 male 64 y, Pt2 female 69 y, Pt3 female 69 y, Pt4 female 52 y. Cells came from **subcutaneous lipoaspirates** processed with a Celution device. 8,000 cells were loaded per sample. Experiment design: "cell type comparison" (not a disease comparison). Files per donor: barcodes, features and matrix (.mtx), plus BAM.

## Cell counts (estimated from the barcode files, which are plain text at about 19 bytes per cell)
| Donor | Barcode file | Estimated cells |
|---|---|---|
| Pt1 | 151 KB | about 8,100 |
| Pt2 | 181 KB | about 9,800 |
| Pt3 | 84 KB | about 4,500 |
| Pt4 | 64 KB | about 3,400 |
| **Total** | | **about 25,900** (estimate, sizes are rounded) |
The share of these cells that are endothelial is not stated. CD31 is also present on some immune cells, so the endothelial count is lower than the total.

## Requirement by requirement
| Requirement | Result | Evidence |
|---|---|---|
| Control group (lean or obese without diabetes) | **NOT MET** | no BMI, disease or control/case design; four donors in a single group |
| Type 2 diabetes group | **NOT MET** | diabetes status not stated for any donor |
| Endothelial cells | **Met in principle** | CD31+ enriched; the paper describes 9 endothelial subsets |
| Over 10,000 endothelial cells | **Not established** | about 25,900 cells in total; endothelial-only count not given; would need more than about 40% endothelial |
| Human adipose tissue | **Met** | human subcutaneous adipose tissue (lipoaspirate), 4 donors |
| Both subcutaneous and visceral | **NOT MET** | subcutaneous only |
| Gene list | Not checked | needs the matrices; the paper reports ANGPT2 among proteins enriched in CD31+ cells |

**Verdict: does not meet the problem statement.** No control vs type 2 diabetes design, no BMI or diabetes information, subcutaneous only, and an unconfirmed endothelial count. It is a description of CD31+ cell types in 4 donors, useful for endothelial subtypes but not for a diabetes comparison.
