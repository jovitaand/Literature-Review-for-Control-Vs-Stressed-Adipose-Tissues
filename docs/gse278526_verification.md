# GSE278526: cross-check against the four criteria

Read 2026-10-28 from the GEO series record, the 4 GEO sample records and the BioProject page (PRJNA1167724). The paper and the data files could not be opened (hosts blocked), so cell counts and the barcode metadata file contents are unverified.

| Criterion | Verdict | Evidence |
|---|---|---|
| Human adipose tissue | **Met** | Omental (visceral) and subcutaneous abdominal adipose stromal cells, Homo sapiens |
| Endothelial cells | **Partly verified** | Series summary: omental endothelial cells show "a lipid-overload, inflammatory, antigen presentation profile in T2D". The deposit is stromal cells sorted as CD45+ and CD45-; the raw matrices have no cell-type labels, so endothelial cells must be annotated. Endothelial counts unknown |
| Control vs type 2 diabetes | **Met, with a major donor-identity risk** | Obese without T2D ("Ob") vs obese with a clinical T2D diagnosis ("Ob+T2D"): 10 patients, 5 per group, paired omental and subcutaneous samples. There is no lean control group. But the 5 patients per group were pooled into one 10x lane, so there are only 4 libraries (omental Ob, omental Ob+T2D, subcutaneous Ob, subcutaneous Ob+T2D) |
| More than 10,000 cells | **Not verifiable here** | No counts in the GEO records. BioProject: 132 Gbases of sequence over 4 SRA experiments; stated average 25,726 reads per cell. A very rough, assumption-based reading is tens of thousands of cells in total **[I]**; do not cite it as a count |

## Details read from GEO
- Platform: 10x Genomics Single Cell 3' (Chromium), Cell Ranger v3.1, GRCh38; instrument fields say NovaSeq 6000 (record) and NextSeq 500 (protocol), which disagree.
- Tissue handling: frozen in CellSTOR, thawed, dissociated, flow sorted (CD45+/- events, CD235a- and live cells); stromal fraction only.
- Libraries: GSM8548235 (BRI-1456, omental Ob), GSM8548236 (BRI-1457, omental Ob+T2D, labelled "metabolic disease"), GSM8548237 (BRI-1458, subcutaneous Ob), GSM8548238 (BRI-1459, subcutaneous Ob+T2D).
- Files: per library barcodes, features, matrix.mtx; series files `GSE278526_RAW.tar` and `GSE278526_cell_barcodes_metadata.tsv.gz` (contents not read). Supplement about 450 MB.
- De-identification: BAM files had SNPs removed (BAMboozle) before upload, so genotype-based demultiplexing of the pooled donors from public BAMs is not possible.

## What this means
- **Donor-level pseudobulk depends on the barcode metadata file.** If it assigns each cell to a donor, you get 5 vs 5 donors per depot. If not, the only replicates are 4 libraries and any gene-level test is pseudoreplicated. Earlier notes in this repo said donor labels are in that file; that is **not verified**.
- Because donors are pooled, batch is perfectly aligned with group and depot at the library level.
- Stromal sorting (CD45+/-) means endothelial proportions are not tissue proportions; compare composition cautiously.

## Comparison with GSE268904
| | GSE268904 | GSE278526 |
|---|---|---|
| Material | Nuclei, subcutaneous, 14 per-donor samples | Cells, omental + subcutaneous stromal, 4 pooled libraries |
| T2D contrast | obese (7) vs obese diabetic (3) donors | Ob (5) vs Ob+T2D (5) donors, if labelled |
| Donor identity | One matrix per donor (clear) | Needs the metadata file (unverified) |
| Endothelial focus | Vascular niche is the paper's focus | Omental endothelial cells highlighted |
| Cell counts | Not in GEO | Not in GEO |
