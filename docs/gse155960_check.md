# GSE155960 (Hildreth et al., Nature Immunology 2021): criteria check (2 Nov 2026)

Drive file shared: `GSE155960_RAW.tar`, 165,826,560 bytes (same size as the GEO file). **Its contents could not be read** here (Drive connector limit 10 MB; Drive and GEO file hosts blocked). The check uses the GEO series and 12 sample records and the paper (PMC8102391). Run `notebooks/GSE155960_criteria_check.ipynb` on the tar for the cell counts and gene-list detection.

| Requirement | Result | Evidence |
|---|---|---|
| Human adipose tissue | **Met** | Homo sapiens, 12 libraries (3 lean + 3 obese donors, CD45-negative and CD45-positive sorts), 10x Chromium |
| Control group | **Met** | 3 lean (BMI < 25) and 3 obese (BMI > 30) donors; "healthy" = no history of cardiovascular or liver disease, diabetes or immune disorder |
| Type 2 diabetes group | **NOT MET** | no diabetic donors; diabetes was an exclusion |
| Endothelial cells | **Present, count not stated** | endothelial cells were identified; the paper gives 110,250 cells in total (82,577 in one figure legend); the endothelial cluster size is not given |
| More than 10,000 endothelial cells | **Not established** | only 6 donors; endothelial cells are within the CD45-negative libraries; unlikely but unverified |
| Both subcutaneous and visceral | **NOT MET** | deep subcutaneous abdominal adipose tissue only |
| Gene list | Not checked | needs the sample count files (notebook does this) |

**Verdict: does not meet the problem statement.** No type 2 diabetes group, no visceral depot, and the endothelial count is unknown. It is a lean vs obese subcutaneous resource.
Per-sample files are dense gene x cell tables (`GSM47171xx_CD45N-L1.csv.gz` etc.), so the notebook needs High-RAM.
