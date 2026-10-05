# GEO metadata read directly (2026-10-05): what is verified now, and what is still blocked

## 1. Download status
- `www.ncbi.nlm.nih.gov` is now reachable, so GEO series and sample records can be read. I read them (below).
- The data **files** are still not downloadable. GEO download links return HTTP 301 and redirect to `ftp.ncbi.nlm.nih.gov`
  (verified: `.../geo/download/?acc=GSE249089&format=file&file=...` redirects to `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249089/suppl/...`).
  That host returns 403 from this environment, so all 27 manifest files failed (0 bytes). figshare, Zenodo, CELLxGENE, EBI also still return 403.
- **To unblock:** add `ftp.ncbi.nlm.nih.gov` (and `figshare.com`, `zenodo.org`, `cellxgene.cziscience.com`, `datasets.cellxgene.cziscience.com`, `www.ebi.ac.uk`)
  to Allowed domains in the environment's Network access, then start a new session. Or run `scripts/download_manifest.py` on your own machine.

## 2. What I extracted instead (committed to `data/metadata/`)
| File | Content |
|---|---|
| `GSE221156_samples.csv` | All 54 Bandesh GEO samples with disease state, sex, age, BMI, HbA1c, race, chemistry, islet center, cause of death |
| `GSE268904_samples.csv` | 14 adipose snRNA-seq samples with sex, age, BMI, condition |
| `GSE342773_samples.csv` | 39 samples: 20 bulk whole-tissue, 8 ex vivo explant, 11 snRNA 10x channels with donor assignment |

## 3. Verified facts per dataset (from the GEO records)
| Accession | What GEO says | Use |
|---|---|---|
| **GSE221156** (Bandesh, PMID 41986506) | Islets from 17 T2D, 17 ND and 14 PD cadaveric donors, 10x scRNA-seq. Donor IDs listed by group | Main islet dataset |
| **GSE342773** | Visceral adipose from 21 Korean adults: Lean n=8 (HbA1c <6.5%, BMI <23), Obesity n=5, Obesity+T2D n=8. Bulk RNA-seq on 20 donors (8/5/7); snRNA-seq on 11 10x channels | Best adipose T2D case-control found: compare Ob vs Ob+T2D to control for BMI |
| **GSE268904** (PMID 41040169) | 14 deep subcutaneous adipose snRNA-seq samples, Göttingen: 4 lean, 7 obese, 3 obese diabetic | Too few T2D (n=3) for donor-level inference |
| **GSE278526** | scRNA-seq (not snRNA-seq) of adipose stromal cells, omental and subcutaneous, obese vs obese+T2D; 58 recruited, 10 profiled by scRNA-seq; 4 pooled libraries | Manifest wrongly calls it snRNA-seq; donor labels live in the barcode metadata file |
| **GSE249089** (PMID 38297378) | 84 METSIM participants, subcutaneous adipose snRNA-seq, 10x v3.1, 4 donors pooled per run and demultiplexed with Demuxlet; stated purpose includes obesity and T2D; raw files withheld for privacy | Large (n=84); T2D labels must come from the metadata file |
| **GSE81608** (PMID 27667665) | Xin 2016, 1,492 cells, 1,600 samples; supplementary RPKM matrix | Accession now confirmed |
| **GSE86473** (PMID 27864352) | A SuperSeries of Lawlor subseries (662 samples) | Confirmed as Lawlor |
| **GSE154126** (PMID 32739450) | Avrahami: 1 newborn, 5 toddlers, 2 adolescents, 4 adult controls, 10 T2D; cells per donor listed (adult controls 22, 22, 37, 54; T2D 24, 13, 11, 30, 15, 23, 11, 19, 34, 24). Supplementary CPM and read-count matrices | Small matrices; very few cells per donor |

## 4. Bandesh (GSE221156) design traps, from the sample table
- 54 libraries from 48 donors: 42 single-donor libraries and 12 libraries from 6 pooled donor pairs (Islet47/48, 57/58, 59/60, 70/71, 84/85, 118/119).
  Six donors (Islet70, 71, 84, 85, 118, 119) exist only in pooled libraries, so their cells need genotype or hashtag demultiplexing; check how the authors assigned them.
- Donors 47, 48, 57, 58, 59, 60 appear both as single libraries and in pools: treat them as one donor, not two.
- **Chemistry is confounded with disease.** Among single-donor libraries: T2D has 7 V2 and 7 V3; ND has 4 V2 and 13 V3; PD has 1 V2 and 10 V3. Include chemistry as a batch term or you will mistake it for disease.
- **Islet center varies:** Scharp-Lacy 22, Wisconsin 7, Southern California 5, Miami 3, UPenn 4 (3 of them T2D, none ND). Another batch factor.
- **Covariates (single-donor libraries, means):** ND age 45.4, BMI 28.4, HbA1c 5.2, 5/17 female. PD age 50.3, BMI 30.8, HbA1c 5.9, 3/11 female. T2D age 48.5, BMI 34.0, HbA1c 7.8, 6/14 female.
  BMI rises with disease state, so T2D effects partly overlap with obesity effects; model BMI and say so.
- One ND donor has no center listed (NA).
- The counts of 17/14/17 donors in the GEO series text are for donors; libraries differ, so always aggregate to donor level.

## 5. Corrections to earlier documents
- Xin accession GSE81608 is now verified (earlier marked unverified).
- Lawlor's GSE86473 is a SuperSeries; the cell data are in its subseries.
- Avrahami's cell counts per donor are now known: adult controls and T2D donors contribute only 11 to 54 cells each (the newborn contributes 84), so cell-type-specific pseudobulk will be very sparse.
- GSE278526 is scRNA-seq of stromal cells, not snRNA-seq as the manifest says.
- Bandesh details I earlier called "10x version not confirmed": both V2 and V3 chemistry are used.

## 6. Still unverified (hosts blocked)
E-MTAB-5061 (Segerstolpe), E-MTAB-15009 (muscle), HPAP/PANC-DB, PanKbase Zenodo record, CELLxGENE collections, HRA002549 and its figshare bundle.
