# Dimension reduction settings of the published SAT vascular atlas

Source: AlZaim, Hassan and Kalucka, *Nature Metabolism* (2026), doi:10.1038/s42255-026-01475-2. The Methods text could not be read here (the
article page returned only the abstract, figure legends and references), so the settings below were taken from the authors' code,
which the paper's Code availability statement points to: https://github.com/Kalucka-lab/SAT_Atlas (commit cd70f40, files in `Data_processing/`).
They describe what the code does; they were not checked against the paper's Methods wording.

| Step | Setting in the authors' code | File |
|---|---|---|
| Upper QC limits (all cells) | `nCount_RNA < 60000 & nFeature_RNA < 10000` | 01 |
| Normalisation | `NormalizeData`, LogNormalize, scale factor 10,000 | 01, 02 |
| Variable genes | `FindVariableFeatures`, vst, 2,000 | 01, 02 |
| Scaling | `ScaleData(vars.to.regress = "nCount_RNA")` | 01, 02 |
| PCA | `RunPCA` on the variable genes (default 50 components) | 01, 02 |
| Number of significant PCs | first PC with cumulative variance > 90% and < 5% own variance (co1); last PC after which the drop in variance is > 0.05 percentage points, +1 (co2); use min(co1, co2) | 01, 02 |
| Batch correction | `RunHarmony(group.by.vars = c("samples", "Study_chemistry"), assay.use = "RNA", max.iter.harmony = 100)` | 01, 02 |
| Whole data set: UMAP / neighbours | on the Harmony components chosen by the rule above | 01 |
| Whole data set: clustering | `FindClusters` over resolutions 0.1 to 3.0, chosen with clustree; **0.1** used (`RNA_snn_res.0.1`), cluster 5 then sub-clustered at resolution 1 | 01 |
| Vascular subset (endothelial, lymphatic, mural) | new object from the vascular cells, whole process repeated | 02 |
| Vascular: UMAP and neighbours | on **30** Harmony components (hard-coded `dims = c(1:30)`) | 02 |
| Vascular: clustering | resolutions 0.1 to 1.0, chosen with clustree; **0.7** used (`RNA_snn_res.0.7`) | 02 |
| Vascular markers | `FindAllMarkers`, only.pos, min.pct 0.15, logfc.threshold 0.3 (Wilcoxon) | 02 |
| Sub-endothelial cells | UMAP and neighbours on 15 Harmony components | 03 |
| Cluster naming | by hand from marker genes (reference list `Vascular_markers_summary.xlsx`, not published) | 02 |

How the combined-results notebook uses them: section 2b, option `endo_method = "paper_pipeline"`, on the GSE268904 nuclei. Differences
from the paper: only one study and chemistry, so Harmony corrects for sample only; the 0.1 and 0.7 resolutions are taken over
from the authors' choices on their 329,774-cell atlas and were not re-tuned for GSE268904; the cluster naming uses our own marker
lists (PECAM1, CDH5 and the other tested genes excluded), not the authors' reference list; the mitochondrial filter
(`max_mito_pct`, default 5%) is added because this is a nuclei data set and is not in the authors' code.
