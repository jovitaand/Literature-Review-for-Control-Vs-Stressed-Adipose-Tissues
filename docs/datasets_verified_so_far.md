# Control vs T2D single-cell datasets: verification log (work in progress)

Status as of 2026-10-05. Only facts read from a primary source (paper full text on PMC) are marked VERIFIED.
GEO, NCBI, and EBI were unreachable from the sandbox, so GEO pages were NOT read directly. Anything marked
UNVERIFIED comes from memory or abstracts and must be checked on GEO/ArrayExpress before use.

## Pancreatic islets (human, scRNA-seq)

| Dataset | Paper | Donors | Cells | Tech | Accession | Status |
|---|---|---|---|---|---|---|
| Lawlor 2017 | Genome Res, doi:10.1101/gr.212720.116 | 8 (5 ND, 3 T2D) | 1,050 (622 ND, 428 T2D) | Fluidigm C1, ~3M reads/cell | Raw: SRA SRP075970 / PRJNA323853. Processed: GEO GSE86473 | VERIFIED (paper text) |
| Segerstolpe 2016 | Cell Metab, doi:10.1016/j.cmet.2016.08.020 | 10 (6 healthy, 4 T2D) | 2,209 | Smart-seq2, ~750k reads/cell | Portal: sandberg.cmb.ki.se/pancreas. ArrayExpress E-MTAB-5061 UNVERIFIED | Counts VERIFIED; accession UNVERIFIED |
| Avrahami 2020 | Mol Metab, doi:10.1016/j.molmet.2020.101057 | 22 (1 newborn, 5 toddler, 2 adolescent, 4 ND adult, 10 T2D adult) | not retrieved | Fluidigm C1, SMART-seq, median 1.7M aligned reads | GEO GSE154126 | VERIFIED except cell count |
| Elgamal 2023 (HPAP integrated map) | Diabetes, doi:10.2337/db23-0130 | 67 downloaded (29 ND, 17 T2D, 10 T1D, 9 ND Aab+); 65 in final map | 192,203 | HPAP scRNA-seq; raw FASTQ from PANC-DB | PANC-DB portal (data access terms apply) | VERIFIED (paper text); 10x vs other chemistry not confirmed |
| Xin 2016 | Cell Metab, doi:10.1016/j.cmet.2016.08.018 | ND and T2D organ donors (counts not retrieved) | 1,492 (alpha, beta, delta, PP) | not confirmed | GEO GSE81608 UNVERIFIED | Abstract only |
| Bandesh 2025/2026 | bioRxiv doi:10.1101/2025.01.17.633590; EMBO J doi:10.1038/s44318-026-00744-w | not retrieved | not retrieved | not retrieved | not retrieved | Title/abstract only; read the paper |

Notes for a beginner:
- Lawlor, Segerstolpe, Xin are small (4-10 donors, plate-based). Good for learning and for quick replication, weak for donor-level statistics.
- Elgamal/HPAP has the most T2D donors (17 vs 29 ND) and is the best candidate for the main analysis, but needs PANC-DB access and re-processing from FASTQ.

## Adipose tissue

No human adipose scRNA/snRNA dataset with explicit T2D vs control labels has been verified yet.
Leads (all UNVERIFIED for T2D labels, donor counts, and accessions):
- Single-cell analysis of human adipose tissue, depot and disease specific cell types, Nat Metab 2020, doi:10.1038/s42255-019-0152-6.
- Human adipose snRNA-seq study using dataset HRA002549 (GSA, China), obesity vs obesity+T2D macrophage analysis, doi:10.3760/cma.j.cn112137-20250324-00714.
- Insulin-resistance snRNA-seq in subcutaneous adipose, PLoS Genet 2020, doi:10.1371/journal.pgen.1009018.
- Emont et al. 2022 (Nature) human/mouse adipose snRNA-seq: obesity-focused, not T2D.

## Benchmarks / methods leads (for part 2, not yet read)
- Benchmarking differential state methods for multi-subject scRNA-seq: Nat Commun/Brief Bioinform, doi:10.1093/bib/bbac286 (pseudobulk performed best).
- Pseudobulk with offsets vs GLMM: Bioinformatics, doi:10.1093/bioinformatics/btae498.

## Next steps
1. Confirm GSE/E-MTAB accessions and per-sample labels on GEO/ArrayExpress (needs network access to those hosts).
2. Read Bandesh 2025/2026 and Elgamal methods for chemistry, QC, and donor metadata.
3. Search Single Cell Portal, HCA, and CELLxGENE for T2D-labelled adipose/liver/muscle data.
