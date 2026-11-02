# SCP1903 (Single Cell Portal, human adipose single-nucleus RNA-seq): criteria check from the files you uploaded (2 Nov 2026)

Files read in full: `metadata.tsv` (35,488 nuclei, 18 columns), `v1-cluster.tsv` (UMAP coordinates), `raw_barcodes.tsv.gz` (35,488; identical to the metadata cells), `raw_genes.tsv.gz` (27,202 genes) and `genes.tsv.gz` (2,000 variable genes). **No count matrix was included**, so gene expression cannot be measured; cell counts, groups and depots come directly from the metadata. Reproduce with `scripts/scp1903_metadata_check.py`; tables in `docs/SCP1903_criteria_check.xlsx`.

## The dataset
Human adipose tissue, 10x Genomics 3' v2 single-nucleus RNA-seq, **6 donors, 11 samples** (5 donors with both depots, donor 5 visceral only). Nuclei: visceral 28,356, subcutaneous 7,132. Groups (`wtstatus`): lean 30,627 nuclei (3 donors), obese 4,861 nuclei (3 donors). Clusters: ASC1 10,066, ADIPO 6,815, ASC2 5,013, MAC 4,102, **EC 3,113**, ASC3 2,831, TCELL 1,466, **EC-P 888**, PC 840, **LEC 354**.

## Requirement by requirement
| Requirement | Result | Evidence |
|---|---|---|
| Control group (lean or obese without diabetes) | **Met (lean)** | 3 lean donors (and 3 obese); diabetes status is not recorded for anyone |
| Type 2 diabetes group | **NOT MET** | no diabetes field; the `disease` column reads "obesity disorder" for every nucleus, lean included |
| Endothelial cells | **Met** | clusters EC (3,113), EC-P (888) and lymphatic LEC (354) |
| Over 10,000 endothelial cells | **NOT MET** | **4,001 blood endothelial nuclei** (4,355 with lymphatic) of 35,488 |
| Human adipose tissue | **Met** | human, 6 donors |
| Both subcutaneous and visceral | **Met** | 28,356 visceral and 7,132 subcutaneous nuclei |
| Gene list | **Partly checked** | 7 of 8 genes are in the gene list; **PECAM1 is absent** (CDH5, NOS3, HIF1A, ERG, VCAM1, ICAM1, ANGPT2 present). Only ERG, VCAM1 and ICAM1 are among the 2,000 variable genes. Expression levels need the count matrix. |

## Blood endothelial nuclei by depot and group
| Depot | Lean | Obese |
|---|---|---|
| Subcutaneous | 1,252 (2 donors) | 890 (3 donors) |
| Visceral | 1,711 (3 donors) | **148** (3 donors) |
Per sample the endothelial share is 4.9% to 10.4% in visceral and 22.5% to 35.5% in subcutaneous samples. Visceral obese samples hold only 17 to 110 endothelial nuclei each, too few for a donor-level comparison. The dataset is dominated by lean visceral samples (26,369 of 35,488 nuclei).

## Verdict
**Does not meet the problem statement.** It satisfies lean control, human adipose, both depots and endothelial cells present, but has **no type 2 diabetes group** and only **4,001 endothelial nuclei** (below 10,000). It is a lean vs obese comparison across both depots.
To check expression of the gene list, the count matrix (portal download, may need a login) is needed; PECAM1 may be missing because of how genes were filtered or named.
