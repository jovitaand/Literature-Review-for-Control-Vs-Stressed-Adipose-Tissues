# GSE128889 (Merrick et al., mesenchymal progenitor hierarchy): criteria check (2 Nov 2026)

Drive file shared: `GSE128889_family.soft.gz` (4,271 bytes). It **was read in full** (small enough for the connector). It is a GEO metadata file (series and sample records); it contains **no count data**, so cell and gene-level questions cannot be answered from it. Saved at `data/metadata/GSE128889_family.soft`; checked with `scripts/check_geo_soft.py`.

## What the file contains
12 samples: **1 human** single-cell sample (GSM3717979: subcutaneous abdominal stromal vascular fraction, 31-year-old woman, 10x Genomics) and **11 mouse** samples (2 mouse single-cell SVF samples and 9 bulk RNA-seq samples of sorted CD142+, DPP4+ and ICAM1+ progenitor populations). The series is a SuperSeries of GSE128890 and GSE128891. Count files are in `GSE128889_RAW.tar` (not provided).

| Requirement | Result | Evidence in the file |
|---|---|---|
| Control group (lean or obese without diabetes) | **NOT MET** | no BMI, lean or obese label on the single human sample |
| Type 2 diabetes group | **NOT MET** | no diabetes information at all |
| Endothelial cells | **Not applicable / not found** | the study targets adipose progenitor cells; no endothelial annotation in the records |
| More than 10,000 endothelial cells | **NOT MET (cannot be)** | one human sample; the whole human library is a single donor; cells filtered to > 500 genes and < 5% mitochondrial reads, counts not stated |
| Human adipose tissue | **Partly** | one human sample (11 of 12 samples are mouse) |
| Both subcutaneous and visceral | **NOT MET** | subcutaneous only (human abdominal; mouse inguinal/axillary); no visceral |
| Gene list | Not checked | no counts in this file; `ICAM1` appears only as a sorting marker of mouse progenitor cells |

**Verdict: does not meet the problem statement** on every condition that matters (single donor, no disease or BMI information, subcutaneous only, progenitor-cell focus, mostly mouse and bulk data).
