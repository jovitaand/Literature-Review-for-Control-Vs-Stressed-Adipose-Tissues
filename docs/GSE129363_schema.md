# GSE129363 discovery cohort: merged dataset and schema (3 Nov 2026)

Built from the three files you uploaded (CellAnnotation, barcodes, genes) plus the GEO sample records, by `scripts/merge_gse129363_schema.py`. Outputs in `data/processed/`:

| Table | File | Rows | One row per |
|---|---|---|---|
| cells | `GSE129363_cells_merged.csv` | 38,170 | barcode (= matrix column) |
| genes | `GSE129363_genes_merged.csv` | 33,694 | gene (= matrix row) |
| samples | `GSE129363_samples_merged.csv` | 25 | barcode suffix (= sample) |
Column-by-column dictionary: `docs/GSE129363_data_dictionary.csv`.

## Relationships
```
 genes.tsv.gz (33,694)            barcodes.tsv.gz (38,170)           CellAnnotation (26,350)
   matrix_row  <---- rows ----+     matrix_col  <--- columns ---+        CellID ---- barcode (26,350 of 38,170 match)
   ensembl_id                 |     barcode  ----------------------------+
   symbol (68 rows duplicated)|        |  barcode_suffix (1..25)
                              |        v
                    COUNT MATRIX     samples (25)  -- gsm, person, tissue, condition (annotation vs GEO)
                  (33,694 x 38,170;
                   not provided)
```
- **Keys:** cells.barcode = CellAnnotation.CellID (all 26,350 annotated cells are in the barcode file; 11,820 barcodes were removed by the authors' filtering). cells.barcode_suffix = samples.barcode_suffix. The matrix columns follow barcodes order, its rows follow genes order.
- **Checks passed:** barcodes unique; annotation IDs unique; Tissue in the annotation agrees with the GEO sample record for every annotated cell; each suffix maps to exactly one sample.

## What the merge shows
- 25 samples from **14 people** (11 with both depots, 3 with one). Annotated cells: subcutaneous 13,451, visceral 12,899.
- **Diabetes labels conflict** between the annotation file and GEO for 5 samples (people E07, 04, 09; 5,559 annotated cells). Both labels are kept (`condition_annotation`, `condition_geo`, and the derived `group_*`).
  By annotation label: T2D subcutaneous 7,091, visceral 4,888; control subcutaneous 6,360, visceral 8,011. By GEO label: T2D 8,894 and 8,644; control 4,557 and 4,255.
- All donors are obese: the control group is **obese without diabetes**; there is no lean group.
- **Gene list present:** all 8 genes are in the gene file (matrix rows 1407 VCAM1, 13264 NOS3, 14474 ANGPT2, 23284 HIF1A, 26229 CDH5, 28229 PECAM1, 30601 ICAM1, 33466 ERG).
- 34 gene symbols occur twice; use `symbol_unique` or the Ensembl ID as the feature name.
- **Not in the three files:** the count matrix and any cell-type label. Vascular subtypes (venous, capillary, arterial, sub-AdEC, lymphatic) therefore cannot be formed from these files alone; the R pipeline `scripts/GSE129363_endothelial_subsets.R` does it from the matrix.

## Why 11,820 barcodes are not in the annotation (verified 3 Nov 2026)
The GEO sample records state the authors' filters: cells with fewer than 200 genes, fewer than 200 UMI or more than 20% mitochondrial reads were removed, and doublets were removed by excluding cells with more than 2,500 genes or more than 13,750 UMI. Applying these to the per-barcode counts of the real data reproduces the removal: **all 26,350 kept cells pass every rule, and 11,640 of the 11,820 removed barcodes fail at least one** (overlapping: 11,111 too few genes, 2,729 too few UMI, 1,481 above 20% mitochondrial reads, 10 above the doublet limits). The other 180 have exactly 200 genes, so the authors' cut-off effectively kept only cells with more than 200 genes. Most removed barcodes are empty or poor droplets: the visceral sample Ate-06-VAT lost 88% of its barcodes. Script: reproduce with `data/processed/GSE129363_cells_merged.csv` plus the per-barcode QC table.
