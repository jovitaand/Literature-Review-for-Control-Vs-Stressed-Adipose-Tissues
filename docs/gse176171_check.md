# GSE176171 (Emont et al., Nature 2022): criteria check, redone with the paper's tables (2 Nov 2026)

**The Drive file `GSE176171.gz` (45,122,572 bytes) could not be read** (connector limit 10 MB; Drive and GEO file hosts blocked; an attempt to read it through another route was refused by the environment and not retried). Everything below comes from the GEO series and sample records, the GEO file listing, the Single Cell Portal page (SCP1376) and the paper's Extended Data tables (PMC9504827). The figures from the tables were read from the web page by an automated extraction, so they were cross-checked by arithmetic (below).

## Verified numbers
- Human data come from **two** experiments: single-nucleus RNA-seq (**GSE176171**, 13 people) and Drop-seq single cells (a different accession, GSE176067, 9 people).
- Human nuclei/cells by depot and BMI (Extended Data Table 2): visceral snRNA-seq 26,092 (BMI < 30) + 7,454 (30 to 40) + 46,539 (> 40) = 80,085; subcutaneous snRNA-seq 24,831 + 12,429 + 20,339 = 57,599; subcutaneous Drop-seq 19,826 + 8,639 = 28,465. Total 166,149 = the SCP1376 human total (check passed). **GSE176171 alone holds the 137,684 nuclei.**
- Endothelial (blood) cells: **15,483** (subcutaneous 5,726 + visceral 9,757 = 15,483, check passed); lymphatic endothelial 2,901; pericytes 1,464; smooth muscle 2,886. Sum 22,734 = SCP1376 human vascular total (check passed).
- The 15,483 includes the Drop-seq cells. **The endothelial count for the GSE176171 nuclei alone is not stated** and is lower.

## Requirement by requirement
| Requirement | Result | Evidence |
|---|---|---|
| Human adipose tissue | **Met** | 13 people (10 paired SAT + VAT, 3 SAT only; 10 women, 3 men), 24 samples in GEO |
| Both subcutaneous and visceral | **Met** | 13 subcutaneous and 11 visceral samples |
| Control group | **Met, by BMI only** | BMI < 30 nuclei: 26,092 (visceral) + 24,831 (subcutaneous); lean (< 25) people in GEO records: 3 (BMI 22.2, 24.0, 24.3); obese (BMI >= 30): 9 |
| Type 2 diabetes group | **NOT MET** | The paper states that subjects with a diagnosis of diabetes, or taking insulin-sensitizing drugs (thiazolidinediones, metformin), were excluded (BIDMC), and that the Pittsburgh cohort excluded type 1 and type 2 diabetes. There is no T2D group by design. (The 43-person insulin sensitive vs resistant cohort is a separate adipocyte experiment, not in GSE176171.) |
| Endothelial cells present | **Met** | blood endothelial 15,483 (both experiments) |
| More than 10,000 endothelial cells | **Partly / not established** | Met only when both depots and both experiments are combined (15,483). Subcutaneous 5,726 and visceral 9,757 are each under 10,000, and the nuclei-only (GSE176171) count is lower and not stated. |
| Gene list (PECAM1, CDH5, NOS3, HIF1A, ERG, VCAM1, ICAM1, ANGPT2) | Not checked | needs the count matrix (738 MB at GEO) |

## Verdict
**Does not meet the problem statement.** The decisive failure is the absence of any type 2 diabetes group (excluded by the study design). It is also borderline on the 10,000-endothelial-cell requirement for the nuclei dataset. It remains usable as a lean vs obese resource with both depots.

## Still unknown
What the shared Drive file contains (it matches no single file in the GEO listing: cell_metadata is 27 MB, the human counts matrix 738 MB). Run `notebooks/GSE176171_criteria_check.ipynb` on it.
