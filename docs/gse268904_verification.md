# GSE268904: cross-check against the four criteria

Read 2026-10-28 from the GEO series record, all 14 GEO sample records, and the BioProject page (PRJNA1119437). The paper (PubMed 41040169, Nature Metabolism 2026) and the bioRxiv preprint could not be opened (hosts blocked), and the data files (`ftp.ncbi.nlm.nih.gov`) could not be downloaded, so cell counts could not be checked.

| Criterion | Verdict | Evidence |
|---|---|---|
| Human adipose tissue | **Met** | All 14 samples: Homo sapiens, "Subcutaneous Adipose tissue" (deep subcutaneous, Göttingen) |
| Endothelial cells | **Partly verified** | Series summary: the study focuses on the vascular niche and identifies 8 endothelial populations (AdEC, AdLEC subtypes). But the GEO deposit is raw per-sample nuclei matrices of all cell types, with no cell-type labels, so endothelial nuclei must be annotated by you. The "8 populations" come from the atlas, which also includes 7 public datasets that are not in this series |
| Control vs type 2 diabetes | **Partly met** | Three groups (series design text): lean n = 4, obese n = 7, obese diabetic n = 3. T2D = history of T2D plus antidiabetic medication. All three T2D donors are obese, so the like-for-like contrast is obese (7) vs obese diabetic (3). Lean vs diabetic is confounded with obesity (BMI about 23 vs about 51) |
| More than 10,000 cells | **Not verifiable here** | GEO records give no cell or nucleus counts. The series summary says "nearly 70,000 vascular cells from 65 donors" but that is the merged atlas (includes 7 public datasets). A secondary web summary of the preprint said over 120,000 nuclei from 14 donors, which I could not confirm. Endothelial-only counts are unknown |

## Verified design details (GEO)
- Platform: 10x Genomics Chromium Single Cell 3' v3.1, nuclei (FANS-sorted), DNBSEQ-G400, CellRanger v6.1.2 with introns, GRCh38. Same chemistry for all 14, so no chemistry confounding.
- Files per sample: barcodes, features, matrix.mtx (raw counts). Series supplement about 999 MB. 14 SRA experiments.
- Groups (series text): lean 2F/2M, age 45.75 ± 13.27, BMI 23.27 ± 0.72; obese 3F/4M, age 38 ± 5.78, BMI 53.22 ± 8.94; obese diabetic 2F/1M, age 42 ± 6.16, BMI 51.3 ± 6.16.
- Donors: lean from hemicolectomy, hernia/hiatus and bypass-type surgeries; obese from cosmetic abdominoplasty, sleeve and bypass surgeries. Surgery type differs between lean and obese groups.
- Per-sample table: `data/metadata/GSE268904_samples.csv`; raw series record: `data/metadata/GSE268904_series_record.txt`.

## What this means for analysis
- The T2D comparison is 3 vs 7 donors: exploratory only; use donor-level pseudobulk and report effect sizes, not strong claims.
- Sample names fall into two series ("9S, 10S, 13S, 17S" and "SAT1 to SAT60"). **[I]** Check whether this reflects two processing batches; two of three diabetic donors are in the first series and one is in the second.
- To meet "more than 10,000 endothelial cells" you will likely need the atlas (65 donors) or to count endothelial nuclei yourself after annotation. The atlas data are advertised as an interactive resource in the paper; its location was not verified.
- Next steps: allow `ftp.ncbi.nlm.nih.gov` (or run `scripts/download_manifest.py` locally), count barcodes per sample, annotate endothelial clusters (for example PECAM1, VWF, CLDN5), and count nuclei per group.

## Update 2026-10-29: can endothelial cell counts be found?
Attempted again. The count files (`ftp.ncbi.nlm.nih.gov`), the paper (PMC, Nature) and the preprint (bioRxiv) are still unreachable, so **no endothelial cell or nucleus count could be verified**. The SRA experiment pages are reachable and give sequencing depth only (`data/metadata/GSE268904_sra_spots.csv`):
- 14 libraries, 309M to 371M read pairs each, 4.75 billion in total (lean 1.33B, obese 2.39B, obese diabetic 1.04B). Depth is similar across groups.
- Reads per sample do not give nuclei per sample. **[I]** Under a plausible 30,000 to 100,000 reads per nucleus, the total would be roughly 50,000 to 150,000 nuclei, which is consistent with the secondary report of over 120,000 nuclei (unconfirmed). Whether endothelial nuclei exceed 10,000 depends on their share (they would need to be roughly 7% to 20% of nuclei under that range); that cannot be settled without counting.
- The 14 per-sample matrices contain all cell types; endothelial nuclei must be identified by annotation (for example PECAM1, VWF, CLDN5).

## Update 2026-10-28: real counts from the Colab run
Source: `GSE268904_consolidated.csv` produced by `scripts/GSE268904_extract.R` on the real `GSE268904_RAW.tar` (per-sample table in `data/metadata/GSE268904_real_counts_per_sample.csv`).
- All 14 samples read; 131,641 barcodes in total (one row each, none dropped); barcodes are unique within sample; every barcode has at least 500 UMIs (the deposited files are already filtered). Median genes per nucleus 1,569.
- Nuclei per group: lean 37,418 (4 samples); obese 58,141 (7); obese diabetic 36,082 (3). Five samples individually exceed 10,000 nuclei. **The more-than-10,000-cells criterion is met.**
- Endothelial nuclei: **not established.** Crude marker rules give very different totals: all of PECAM1, VWF and CLDN5 detected 9,380; at least 4 of 6 markers 27,679; at least 3 of 6 38,452. Per-sample endothelial fraction (strict rule) varies from 0.6% to 18%. FLT1 and PTPRB are detected in many non-endothelial nuclei, and ambient RNA inflates the looser rules.
- Obese diabetic arm alone: about 2,100 (strict) to 6,200 (at least 4 of 6) endothelial nuclei from 3 donors.
- Next step: cluster and annotate (Seurat section, or add non-endothelial marker genes through `GENES_EXTRA`) to get a defensible endothelial count.
