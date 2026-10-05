# More data sources for Control vs T2D single-cell transcriptomics (2026-10-05)

Labels: VERIFIED = numbers read in paper text or abstract returned by the literature index; LEAD = title/abstract only, donor counts and
accessions not checked; PORTAL = a place to search, listing seen only in search snippets. No GEO or portal page was opened (blocked).
Where a manifest entry you uploaded (data/manifest.csv) overlaps, it is noted.

## 1. Network reality check (tested from the compiling session)
Every candidate host below returned HTTP 403 (egress policy) when tested: zenodo.org, cellxgene.cziscience.com, datasets.cellxgene.cziscience.com,
cells.ucsc.edu, singlecell.broadinstitute.org, data.humancellatlas.org, pankbase.org, data.pankbase.org, hpap.pmacs.upenn.edu, www.ebi.ac.uk,
www.ncbi.nlm.nih.gov, figshare.com. github.com answered but hosts no data here. So these sources are fine for you to download on your own machine,
but the sandbox cannot fetch them until Network access is widened. Add all of the above to Allowed domains if you want me to download.

## 2. Islet datasets (human, T2D vs control) beyond the first list
| Source | What it adds | Status |
|---|---|---|
| HCA Data Portal copies of islet studies: Segerstolpe (project ae71be1d-ddd8-4feb-9bed-24c3ddb6e1ad), Lawlor (c6ad8f9b-d26a-4811-b2ba-93d487978446), "Single cell RNA-seq of human pancreatic endocrine cells from juvenile, adult control and T2D donors" = Avrahami (99101928-d9b1-4aaf-b759-e97958ac7403), "Massively parallel single-cell RNA-seq of 26,677 pancreatic islet cells from healthy and T2D donors" (1c6a960d-52ac-44ea-b728-a59c7ab9dc8e) | A second, standardized place to get the same islet studies, plus one 26,677-cell healthy/T2D study I had not listed | LEAD: project IDs seen in the HCA portal update list; donor counts not checked. URL pattern: https://data.humancellatlas.org/explore/projects/<id> |
| HPAP-derived papers listed on PANC-DB: Weng et al., Nat Genet 2023, "Integration of single-cell multiomic measurements across disease states with genetics identifies mechanisms of beta cell dysfunction in T2D" (PMC10550816); Diabetologia 2026, "Single-cell profiling of pancreatic islets maps subtype-associated molecular alterations in T2D" (PubMed 42678439); Sci Rep 2025, same-donor scRNA-seq and snRNA-seq (PubMed 41102292) | Methods papers and extra analyses on the same HPAP donors | LEAD; same donors as HPAP, so NOT independent of Elgamal/PanKbase |
| Human islet multiome / T2D genetics: Nat Genet 2021 chromatin accessibility (doi:10.1038/s41588-021-00823-0); bioRxiv single-cell multiome of islets and diabetes risk (PubMed 39149326) | scATAC plus RNA for the T2D GWAS link | LEAD |
| "Single-cell mRNA-regulation analysis reveals cell type-specific mechanisms of T2D" (Nat Commun 2025, doi:10.1038/s41467-025-65060-z) | Re-analysis of existing public islet datasets; useful as a cross-check of which public sets others combine | LEAD (abstract: "repurposed datasets"?). Check its methods for the list of GSE accessions |
| scRiskCell, iMeta 2025 (PMC12371254) | T2D islet risk-cell framework built on public data | LEAD |
| Imaging mass cytometry of human pancreas in T2D, Cell Rep 2021 | Protein-level, not RNA; validation of composition changes | LEAD |
| "Single-cell atlas of human pancreatic islet and acinar endothelial cells in health and diabetes", Nat Commun 2025 (doi:10.1038/s41467-024-55415-3) | Endothelial cells in diabetes; 3 pancreases for enrichment | LEAD |

## 3. Adipose tissue (the tissue your repo is named for)
| Source | Details | Status |
|---|---|---|
| HRA002549 (GSA-Human, China) | Visceral adipose SVF scRNA-seq; paper (Chin Med J 2025, doi:10.3760/cma.j.cn112137-20250324-00714) compares macrophages in obesity and obesity with T2D. Your manifest also lists a processed figshare bundle | Abstract confirms the dataset is used for obesity vs obesity+T2D; group sizes not checked |
| METSIM snRNA-seq, GSE249089 (in your manifest) | Paper: Diabetol Metab Syndr 2025, doi:10.1186/s13098-025-01975-3, "snRNA-seq of subcutaneous adipose tissue from 84 individuals with MetS from the METSIM cohort" | VERIFIED from abstract; check whether T2D status is recorded per participant |
| Front Immunol 2026, "A mitochondrial-stress adipocyte-macrophage circuit sustaining metaflammation in human T2D adipose tissue" (doi:10.3389/fimmu.2026.1768845, PMC13181277) | Human T2D adipose; likely has a T2D vs control comparison | LEAD; open the data-availability section |
| Nat Metab 2026, "Defining the vascular niche of human adipose tissue across metabolic states" (doi:10.1038/s42255-026-01475-2) | Adipose endothelial cells across metabolic states | LEAD |
| HCA Data Portal adipose projects: "Single cell full-length transcriptome of human subcutaneous adipose tissue" (f77290ae-0d7b-4239-b0fe-3cf2c9e8858d); integrated single cell and spatial WAT map (57916660-af5a-44d5-a7a9-2e84b65f8a68) | Reference atlases, obesity/lean, not T2D case-control | LEAD |
| Single Cell Portal SCP1376 (human and mouse WAT atlas) and SCP133 (10 human subcutaneous SVF samples) | Reference; obesity-focused | PORTAL |
| Adipose snRNA-seq "Adipocyte subsets with distinct metabolic profiles" (bioRxiv 2025, doi:10.1101/2025.09.14.673351; 65,668 nuclei, SAT vs IAT) | Depot comparison | LEAD |
| Mouse db/db adipose snRNA-seq (bioRxiv 2024, doi:10.1101/2024.02.06.578860) | Mouse T2D model; hypothesis support only | Not human |

## 4. Other metabolically relevant tissues
| Tissue | Source | Status |
|---|---|---|
| Skeletal muscle | Hansen et al., J Physiol 2025: 10 controls, 9 T2D, 38 biopsies, 135,225 nuclei, ArrayExpress E-MTAB-15009 (repeated measures) | VERIFIED |
| Skeletal muscle | "Full-length single nuclei atlas of human skeletal muscle insulin resistance" (bioRxiv 2026, doi:10.64898/2026.02.26.708250; clamp-phenotyped) | LEAD |
| Kidney (diabetic kidney disease, T2D complication) | Nat Genet 2024, doi:10.1038/s41588-024-01802-x: 81 samples across healthy, diabetic and hypertensive kidneys (scRNA, snRNA, Visium, CosMx, snATAC) | VERIFIED from abstract |
| Kidney | snRNA-seq of 5 healthy vs 6 DKD kidneys, Diabetes 2025, doi:10.2337/db25-0272 | VERIFIED from abstract; very small |
| Kidney | Early diabetic nephropathy scRNA landscape, listed on HCA (project 577c946d-6de5-4b55-a854-cd3fde40bff2); spatial DKD atlas (Nature 2026, doi:10.1038/s41586-026-10363-4); DKD snRNA+snATAC bioRxiv 2022 (doi:10.1101/2022.01.28.478204) | LEAD |
| Liver | Int J Biol Sci 2024, doi:10.7150/ijbs.99176, "Reconstruction of the hepatic microenvironment ... underlying type II diabetes through scRNA-seq" | LEAD; species and donor numbers unchecked |
| Liver (reference only) | Human liver sc/sn/spatial atlas (Hepatol Commun 2022, doi:10.1002/hep4.1854) | Not T2D-labelled |
| Blood (PBMC) | Front Endocrinol 2024, doi:10.3389/fendo.2024.1397661; Front Immunol 2024, doi:10.3389/fimmu.2024.1501660 | LEAD; blood is a poor proxy for islet/adipose biology, use for exploratory immune questions |

## 5. Portals to search yourself (filter disease = type 2 diabetes mellitus)
- NCBI GEO DataSets: query like `("type 2 diabetes"[All Fields]) AND ("single cell" OR "single nucleus" OR "snRNA-seq") AND "Homo sapiens"[Organism]`; read each Series Matrix file for diagnosis and donor columns.
- CZ CELLxGENE Discover (cellxgene.cziscience.com): filter Disease = type 2 diabetes mellitus, Organism = Homo sapiens. Downloads are .h5ad with standardized metadata. The Bandesh collection (58e85c2f-d52e-4c19-8393-b854b84d516e) is here.
- Human Cell Atlas Data Portal (data.humancellatlas.org): search "diabetes".
- Single Cell Portal (Broad): search "diabetes" or "adipose".
- PanKbase (pankbase.org) and PANC-DB (hpap.pmacs.upenn.edu): islet-focused, with donor metadata.
- Common Metabolic Diseases Genome Atlas, CMDGA (cmdga.org): appeared in search results; curated metabolic-disease genomics resource. I did not open it.
- ArrayExpress/BioStudies: search "type 2 diabetes single-cell".

## 6. Recent reviews on single-cell work in diabetes (the open item from before)
- "Single-Cell Multi-Omics in Type 2 Diabetes Mellitus: Revealing Cellular Heterogeneity and Mechanistic Insights", Int J Mol Sci 2025, doi:10.3390/ijms262211005 (PMC12652634). Review of single-cell multi-omics in islets. LEAD: read the abstract and skim the dataset tables.
- "Adipose tissue macrophage heterogeneity in the single-cell genomics era", Mol Cells 2024, doi:10.1016/j.mocell.2024.100031 (PMC10960114). Adipose macrophages; relevant to obesity and T2D.
- "Single cell and single nucleus RNA sequencing in liver tissues: applications and prospects", Front Genet 2026, doi:10.3389/fgene.2026.1781941. Liver methods, including model organisms.
- Best-practice reviews already in the reading list remain the core methods references.

## 7. What I recommend
1. Keep Bandesh (GSE221156) as the main dataset and PanKbase/HPAP as validation (see `better_datasets.md`).
2. For adipose, rank by credibility: first METSIM GSE249089 (84 participants; check T2D labels), then HRA002549 (obesity vs obesity + T2D), then the 2026 Front Immunol human T2D adipose paper once you read its data-availability section.
3. For a non-islet T2D tissue with verified numbers, use the skeletal muscle snRNA-seq (E-MTAB-15009), handling repeated measures.
4. Treat kidney, liver and PBMC as exploratory.
5. Before using any source, open its GEO/portal page and confirm diagnosis labels and donor counts. I could not.
