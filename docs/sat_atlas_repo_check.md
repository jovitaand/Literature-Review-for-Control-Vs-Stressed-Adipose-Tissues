# SAT vascular atlas, GitHub repository Kalucka-lab/SAT_Atlas: criteria check (2 Nov 2026)

Repository cloned and read in full (62 MB, commit cd70f40, Feb 2026): **R/Python code and notebooks only, no data**. The scripts read objects such as `All_Data_Atlas.RDS` and `Vascular_IntegratedObj_annotated.RDS` that are not in the repository. So cell counts and gene values cannot be computed from it; the design, conditions and methods can be read from the code and the paper (AlZaim et al., Nature Metabolism 2026; preprint PMC12485908).

## What the code shows
- README: "R and python scripts ... performed on **subcutaneous adipose tissue** samples."
- Metadata column `Condition` has three values used throughout the code: **Lean, Obese, Diabetic Obese**. Differential expression is run per vascular cell type for Obese vs Lean, Diabetic Obese vs Obese and Diabetic Obese vs Lean.
- Vascular labels: venous endothelial cells 1 to 3, capillary endothelial cells 1 and 2, arterial endothelial cells, sub-endothelial cells, lymphatic endothelial cells, pericytes, VSMCs.
- Method: Seurat `FindMarkers` with the likelihood-ratio test on scaled data, covariates sex and study/chemistry, per cluster. This tests **cells, not donors** (pseudoreplication), so its p-values are too small; a per-donor pseudobulk analysis would be more reliable.
- Integration: Harmony across studies/chemistries. Condition definitions (paper): obesity = BMI > 30; non-obese < 30; type 2 diabetes = previous T2D history plus antidiabetic medication.
- The code uses PECAM1, CDH5 and VWF as endothelial markers (EndMT score), and ERG in a transcription-factor analysis. **NOS3, HIF1A, VCAM1, ICAM1 and ANGPT2 are not analysed anywhere in the repository.**

## Requirement by requirement
| Requirement | Result | Evidence |
|---|---|---|
| Control group (lean or obese without diabetes) | **Met** | Lean and Obese groups. In the in-house 14-donor cohort: 4 non-obese (BMI 23.3 +/- 0.7), 7 obese (BMI 53.2 +/- 8.9). |
| Type 2 diabetes group | **Met** | "Diabetic Obese": 3 donors in the in-house cohort (BMI 51.3). The donors per group in the 65-donor atlas are not stated. |
| Endothelial cells | **Met** | 8 endothelial populations, plus pericytes and VSMCs |
| Over 10,000 endothelial cells | **Likely met, not stated** | about 68,500 vascular cells in total; the endothelial-only count is not given |
| Human adipose tissue | **Met** | human, 65 donors |
| Both subcutaneous and visceral | **NOT MET** | subcutaneous only |
| Gene list | **Cannot be checked from the repository** | needs the Seurat object; `scripts/SAT_atlas_criteria_check.R` does it |

The in-house cohort (4 / 7 / 3 donors) matches the GSE268904 counts verified earlier (lean 4, obese 7, obese diabetic 3).

**Verdict: the best match so far on diabetes (a lean group, an obese group and a diabetic group, large vascular cell numbers) but it fails on "both subcutaneous and visceral" (subcutaneous only), and the repository contains no data.** To run the gene list and obtain endothelial counts per donor, the atlas object is needed (portal states matrices are downloadable; not verified). Notebook: `notebooks/SAT_atlas_criteria_check.ipynb`.
