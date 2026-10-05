# Better public datasets for Control vs T2D (update, 2026-10-05)

This supersedes the dataset ranking in `mentor_plan_part1.md` and `mentor_plan_part2.md`. Labels: VERIFIED = read in paper/portal text;
ABSTRACT = abstract only; UNVERIFIED = not checked. GEO, EBI and CELLxGENE pages themselves could not be opened from the compiling
environment; accessions come from the papers' own data-availability statements.

## What changed
Two larger, freely downloadable islet datasets turned up that beat the small plate-based sets for statistics, and one of them does
not need an access application. The 2017-era datasets (Lawlor, Segerstolpe) drop to learning/replication only.

## Ranked

| Rank | Dataset | Tissue | Donors | Cells | Tech | Where | Verification |
|---|---|---|---|---|---|---|---|
| 1 | Bandesh et al. 2025/2026, Stitzel lab (Jackson Lab) | Human islets | 48: 17 non-diabetic, 14 pre-diabetic, 17 T2D | 245,878 | droplet-based scRNA-seq (10x version not confirmed) | GEO **GSE221156**; BioProject PRJNA913127; CELLxGENE collection; code on Zenodo 14656366 | VERIFIED (data-availability text) |
| 2 | PanKbase integrated map (Vu et al., bioRxiv 2026) | Human islets | 140 donors, 191 samples: 69 controls, 12 autoantibody+, 11 pre-diabetes, 12 T1D, 36 T2D | 448,935 | integrated from HPAP, IIDP, Prodo and other sources | Final annotated object on Zenodo record 15596314; PanKbase portal | VERIFIED (paper text); preprint |
| 3 | HPAP / Elgamal 2023 | Human islets | 29 ND, 17 T2D (of 67) | 192,203 | HPAP scRNA-seq | PANC-DB; registration and data-use terms | VERIFIED |
| 4 | Skeletal muscle snRNA-seq, Hansen et al. 2025 (J Physiol) | Human skeletal muscle (non-islet) | controls n=10; T2D n=9 | 135,225 nuclei from 38 biopsies | snRNA-seq | ArrayExpress **E-MTAB-15009** | VERIFIED (data-availability text) |
| 5 | Lawlor 2017 | Islets | 5 ND, 3 T2D | 1,050 | Fluidigm C1 | GEO GSE86473 | VERIFIED |
| 6 | Segerstolpe 2016 | Islets | 6 healthy, 4 T2D | 2,209 | Smart-seq2 | sandberg.cmb.ki.se/pancreas | VERIFIED counts |
| 7 | Avrahami 2020 | Islets | 10 adult T2D, 4 adult ND, 8 younger ND | not retrieved | Fluidigm C1 | GEO GSE154126 | VERIFIED donors |

Leads, not verified: "Latent plasticity of the human pancreas" (57 donors across development, health, T2D; over 4 million cells and nuclei,
doi:10.1016/j.cmet.2026.07.023, ABSTRACT); the sex-differences islet study (52 donors, scRNA+snATAC, ABSTRACT); the "largest single-cell islet
dataset, 650,000 cells across 121 donors" alpha-to-beta transdifferentiation study (doi:10.1101/2025.02.14.638309, ABSTRACT); PBMC scRNA-seq in T2D
(Front Endocrinol 2024, doi:10.3389/fendo.2024.1397661, ABSTRACT; blood, not a primary disease tissue).

## What each is best for
- **Best for the main project: Bandesh GSE221156.** 17 vs 17 donors in one lab with a single processing pipeline, droplet data, GEO deposit, and a free
  CELLxGENE copy. The authors report 14 cell types present in every donor and 511 DEGs in T2D beta cells (VERIFIED, abstract/results text).
  The CELLxGENE collection is split into four instances (all cells, beta, alpha, delta), so you can start with a smaller beta-cell file.
- **Best for scale and validation: PanKbase.** 36 T2D and 69 control donors, but it is an integration of several sources including HPAP, so it is
  NOT independent of Elgamal/HPAP. Do not use HPAP and PanKbase as two separate replications.
- **Independence check:** donor overlap between Bandesh and PanKbase is unknown. Check donor IDs before claiming replication.
- **Best for learning:** Lawlor (small, plate-based, easy), then Bandesh beta-cell subset.
- **Best for advanced work:** PanKbase object (donor-level modelling across phenotypic groups); the muscle snRNA-seq set if you want a second tissue.
- **Adipose:** still no verified T2D-labelled human adipose set. A Broad Single Cell Portal adipose atlas exists (SCP1376, human and mouse white fat across body weight; obesity-focused, not T2D) and a
  human adipose SVF Drop-seq study (SCP133, 10 subcutaneous samples) exists, but neither is a T2D case-control design.

## Cautions
- **Muscle study design:** 38 biopsies from 19 participants (10 control, 9 T2D) before/after 2 weeks of high-intensity training. That implies two biopsies per person: repeated measures. Use donor as a random effect or analyse the baseline biopsy only. Donor number is 19, not 38.
- **Pre-diabetes group (Bandesh):** you can drop it or model it as a third group; it is not "control".
- **Preprints:** Bandesh (bioRxiv, later EMBO J 2026) and PanKbase (bioRxiv) include different versions; use the published or latest version and note it.
- **Confounders:** age, BMI, sex, HbA1c, and isolation site differ by donor; check them in the supplementary donor tables before modelling.

## How to download (untested)
- **Bandesh:** GEO page https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221156 (supplementary files). Or open the CELLxGENE collection
  https://cellxgene.cziscience.com/collections/58e85c2f-d52e-4c19-8393-b854b84d516e and download each dataset as .h5ad (Download button on the dataset row). Python alternative: the cellxgene-census package (verify the current API in the official docs).
- **PanKbase:** https://zenodo.org/records/15596314 (final annotated dataset); portal https://pankbase.org; analysis files https://data.pankbase.org/analysis-sets/PKBDS1349YHGQ/.
  Zenodo: `wget` the file links listed on the record page; check file sizes first.
- **Muscle:** https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-15009 (Files tab).
- **Code for Bandesh pipeline:** https://zenodo.org/records/14656366 (a ready reference pipeline to read alongside the paper).

## Updated recommendation
1. Start on Lawlor only if you want a tiny warm-up; otherwise go straight to the Bandesh beta-cell subset (CELLxGENE, 17 vs 17 donors).
2. Main analysis: Bandesh, donor-level pseudobulk per cell type, T2D vs ND (and PD as a sensitivity group).
3. Replication: Elgamal/HPAP or the PanKbase object restricted to donors not in Bandesh; directional concordance of DEGs.
4. Secondary: Lawlor/Segerstolpe/Avrahami only as exploratory support.
5. Verify GEO sample labels and donor tables before any inference.

Still open: a T2D-labelled adipose set, and 2022-2026 reviews on single-cell work in diabetes.
