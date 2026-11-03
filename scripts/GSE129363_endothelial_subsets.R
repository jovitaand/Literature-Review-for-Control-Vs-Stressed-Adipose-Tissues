# GSE129363: merge the discovery-cohort files, then form the vascular endothelial subsets
#   venous AdECs, capillary AdECs, arterial AdECs, sub-AdECs and lymphatic AdECs,
# and report the genes of interest: ERG, CD31 (PECAM1), VE-cadherin (CDH5), eNOS (NOS3), VCAM1, ICAM1, HIF-1alpha (HIF1A), ANGPT2.
# For Google Colab (R runtime, High-RAM) or the pipeline. Tested on synthetic data in the real file format; NOT yet run on the real matrix.
#
# IMPORTANT: the subtype marker panels below are curated by the assistant from the literature and from the clusters seen in GSE268904 / GSE129363
# (they are NOT the authors' list, which is not published in the repository). Check the sub-cluster markers table before trusting any label.

# ---- SETTINGS -----------------------------------------------------------
DATA_DIR <- "/content/gse129363_data"      # folder with barcodes, genes, matrix and CellAnnotation files (searched recursively)
OUT_DIR  <- "/content/gse129363_subsets"
USE_GDOWN <- TRUE
drive_ids <- c("1xVFGimI87DMyyUHoJxbBBr_PXS_tU3bt", "1l77NQm3jVmJ8ieoD9KRpWjGFRaMev5hk", "1imnxlaJKfZMK5WUJas-5IXjQbuHNi79y", "190PY9cnmseA9_xXRD6vBNBeCxXY4tBZS")   # barcodes, annotation, genes, matrix
GENES <- c(ERG = "ERG", CD31 = "PECAM1", `VE-cadherin` = "CDH5", eNOS = "NOS3", VCAM1 = "VCAM1", ICAM1 = "ICAM1", `HIF-1alpha` = "HIF1A", ANGPT2 = "ANGPT2")
N_PCS <- 20; CLUSTER_RES <- 0.5; SUB_RES <- 0.8; USE_HARMONY <- TRUE
ENDO_Z_MIN <- 0.5; LYMPH_Z_MIN <- 1.5; LABEL_Z_MIN <- 0.5; MIN_CELLS_PER_SAMPLE <- 10
# -------------------------------------------------------------------------
for (p in c("Matrix", "data.table", "writexl", "Seurat", "ggplot2", "R.utils")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("edgeR", quietly = TRUE)) { install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install("edgeR", ask = FALSE, update = FALSE) }
if (USE_HARMONY && !requireNamespace("harmony", quietly = TRUE)) try(install.packages("harmony", repos = "https://cloud.r-project.org"), silent = TRUE)
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl); library(Seurat); library(ggplot2); library(edgeR)})
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); set.seed(1)

# ---- 1. files
if (USE_GDOWN) { system("pip install -q gdown"); for (id in drive_ids) system(sprintf("cd %s && gdown --id %s", DATA_DIR, id)) }
af <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE); f1 <- function(p) { h <- af[grepl(p, basename(af), ignore.case = TRUE)]; if (!length(h)) stop("No file matching ", p); h[1] }
f_bc <- f1("barcodes"); f_gn <- f1("genes|features"); f_mx <- f1("\\.mtx"); f_an <- f1("annotation")

# ---- 2. MERGE: cells (barcodes + annotation + sample key), genes, matrix  (schema: docs/GSE129363_schema.md)
bc <- fread(f_bc, header = FALSE)[[1]]; ge <- fread(f_gn, header = FALSE, col.names = c("ensembl_id", "symbol")); ge[, matrix_row := .I]
ge[, symbol_unique := make.unique(symbol)]                                   # 34 symbols occur twice
m <- as(readMM(f_mx), "CsparseMatrix"); stopifnot(nrow(m) == nrow(ge), ncol(m) == length(bc)); rownames(m) <- ge$symbol_unique; colnames(m) <- bc
an <- fread(f_an); setnames(an, c("barcode", "sample_name", "condition_annotation", "tissue"))
sample_key <- data.table(sample_name = c("SVF2-SAT", "SVF3-VAT", "SVF5-SAT", "Ate-04-VAT", "Ate-04-SAT", "Ate-05-VAT", "Ate-05-SAT", "Ate-06-VAT", "Ate-06-SAT", "Ate-08-VAT", "Ate-08-SAT", "Ate-09-VAT", "Ate-09-SAT", "Ate-10_VAT", "Ate-10-SAT", "Ate-11-VAT", "Ate-11-SAT", "Ate-12-VAT", "Ate-12-SAT", "Ate-13-VAT", "Ate-13-SAT", "Ate-14-VAT", "Ate-14-SAT", "Ate-15-VAT", "Ate-15-SAT"), person = c("E02", "E03", "E07", "04", "04", "05", "05", "06", "06", "08", "08", "09", "09", "10", "10", "11", "11", "12", "12", "13", "13", "14", "14", "15", "15"), group_geo = c("T2D", "T2D", "T2D", "T2D", "T2D", "T2D", "T2D", "control", "control", "control", "control", "T2D", "T2D", "T2D", "T2D", "control", "control", "control", "control", "T2D", "T2D", "control", "control", "control", "control"))
cells <- merge(data.table(barcode = bc, matrix_col = seq_along(bc)), an, by = "barcode", all.x = TRUE, sort = FALSE)[order(matrix_col)]
cells <- merge(cells, sample_key, by = "sample_name", all.x = TRUE, sort = FALSE)[order(matrix_col)]
cells[, group_annotation := ifelse(condition_annotation == "NonDiabetic", "control", ifelse(condition_annotation == "Diabetic", "T2D", NA))]
cells[, in_annotation := !is.na(sample_name)]
cat("Cells (barcodes):", nrow(cells), "| kept by the authors:", sum(cells$in_annotation), "| genes:", nrow(ge), "\n")
fwrite(cells, file.path(OUT_DIR, "cells_merged.csv")); fwrite(ge, file.path(OUT_DIR, "genes_merged.csv"))
keep <- which(cells$in_annotation); counts <- m[, keep]; rm(m); counts <- counts[Matrix::rowSums(counts > 0) >= 3, ]; cm <- as.data.frame(cells[keep]); rownames(cm) <- cm$barcode
seu <- CreateSeuratObject(counts, meta.data = cm[colnames(counts), c("sample_name", "tissue", "person", "condition_annotation", "group_annotation", "group_geo")]); rm(counts)
G <- GENES[GENES %in% rownames(seu)]; cat("Genes of interest found:", paste(G, collapse = ", "), "| missing:", paste(setdiff(GENES, rownames(seu)), collapse = ", "), "\n")

# ---- 3. all cells: cluster, then keep the VASCULAR ENDOTHELIAL clusters
BROAD <- list(endothelial = c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB","CDH5","KDR","RAMP2"), lymphatic = c("PROX1","MMRN1","PDPN","CCL21","TFF3","LYVE1","FLT4"),
  mural = c("ACTA2","RGS5","MYH11","PDGFRB","NOTCH3","CNN1"), fibroblast = c("PDGFRA","DCN","LUM","COL1A1","COL1A2","CD34"), immune = c("PTPRC","CD14","CD68","LYZ","CD3E","MS4A1","CD74"),
  adipocyte = c("ADIPOQ","PLIN1","LIPE"), mesothelial = c("MSLN","KRT19","KRT8","UPK3B"))
proc <- function(o, res, npcs = N_PCS) {
  o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE); o <- ScaleData(o, features = union(VariableFeatures(o), unname(G)), verbose = FALSE)
  npcs <- min(npcs, ncol(o) - 1, length(VariableFeatures(o))); o <- RunPCA(o, npcs = npcs, verbose = FALSE); red <- "pca"
  if (USE_HARMONY && requireNamespace("harmony", quietly = TRUE) && length(unique(o$sample_name)) > 1) { o <- tryCatch(harmony::RunHarmony(o, group.by.vars = "sample_name", verbose = FALSE), error = function(e) o); if ("harmony" %in% names(o@reductions)) red <- "harmony" }
  o <- FindNeighbors(o, reduction = red, dims = 1:npcs, verbose = FALSE); o <- FindClusters(o, resolution = res, verbose = FALSE); o <- RunUMAP(o, reduction = red, dims = 1:npcs, verbose = FALSE); o$cluster <- paste0("C", as.character(o$seurat_clusters)); o }
clz <- function(o, panels) {                                                   # cluster-level z-scored mean expression of each panel
  lg <- GetAssayData(o, layer = "data"); cl <- sort(unique(o$cluster)); cmn <- sapply(cl, function(c) Matrix::rowMeans(lg[, o$cluster == c, drop = FALSE])); z <- t(scale(t(cmn))); z[is.na(z)] <- 0
  if (ncol(z) < 3) z <- cmn - rowMeans(cmn)                                    # too few clusters to z-score
  sapply(panels, function(g) { g <- g[g %in% rownames(z)]; if (length(g)) colMeans(z[g, , drop = FALSE]) else rep(NA_real_, ncol(z)) }) |> (\(x) { if (is.null(dim(x))) x <- matrix(x, nrow = length(cl), dimnames = list(cl, names(panels))); rownames(x) <- cl; x })() }
seu <- proc(seu, CLUSTER_RES); cat("Clusters (all cells):", length(unique(seu$cluster)), "\n")
zb <- clz(seu, BROAD); vasc_cl <- rownames(zb)[(zb[, "endothelial"] >= ENDO_Z_MIN | zb[, "lymphatic"] >= LYMPH_Z_MIN) & !is.na(zb[, "endothelial"])]
cat("Vascular endothelial clusters:", paste(vasc_cl, collapse = ", "), "\n"); print(round(zb, 2))
seu$broad_class <- ifelse(seu$cluster %in% vasc_cl, "vascular endothelial", "other")
E <- subset(seu, subset = broad_class == "vascular endothelial"); cat("Vascular endothelial cells:", ncol(E), "(", round(100 * ncol(E) / ncol(seu), 1), "% of the kept cells )\n")

# ---- 4. SUBTYPES inside the endothelial subset (panels exclude the 8 genes of interest, to avoid circularity)
SUB <- list(
  venous_AdEC    = c("ACKR1", "NR2F2", "SELP", "SELE", "PLVAP", "IL6", "CCL14"),
  capillary_AdEC = c("CA4", "RGCC", "FABP4", "CD36", "RBP7", "GNG11", "BTNL9", "CD300LG", "APOLD1"),
  arterial_AdEC  = c("GJA5", "GJA4", "SEMA3G", "HEY1", "FBLN5", "SOX17", "EFNB2", "DLL4", "LTBP4"),
  lymphatic_AdEC = c("PROX1", "MMRN1", "PDPN", "CCL21", "TFF3", "LYVE1", "FLT4"),
  sub_AdEC       = c("PLIN1", "PDE3B", "ADH1B", "ACSL1", "ADIPOQ", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PTPRC"))     # mixed endothelial / adipocytic / mesenchymal / immune signature
E <- proc(E, SUB_RES); cat("Sub-clusters:", length(unique(E$cluster)), "\n"); zs <- clz(E, SUB); print(round(zs, 2))
E$subtype <- vapply(E$cluster, function(cl) { s <- zs[cl, ]; if (all(is.na(s)) || max(s, na.rm = TRUE) < LABEL_Z_MIN) "unassigned_EC" else names(s)[which.max(s)] }, "")
# cell-level scores for comparison (noisy on shallow data)
lgE <- GetAssayData(E, layer = "data"); for (nm in names(SUB)) { g <- SUB[[nm]][SUB[[nm]] %in% rownames(E)]; if (length(g) >= 2) E[[paste0("score_", nm)]] <- as.numeric(Matrix::colMeans(lgE[g, , drop = FALSE])) }   # mean log-normalised expression of the panel
cat("\nCells per subtype:\n"); print(sort(table(E$subtype), decreasing = TRUE))
sub_markers <- tryCatch(as.data.table(FindAllMarkers(E, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE))[, head(.SD, 10), by = cluster], error = function(e) data.table())
cluster_label <- data.table(cluster = rownames(zs), nuclei = as.integer(table(E$cluster)[rownames(zs)]), label = vapply(rownames(zs), function(c) E$subtype[E$cluster == c][1], ""), round(as.data.table(zs), 2))

# ---- 5. THE SUBSETS: one object / barcode list per subtype
SUBTYPES <- c("venous_AdEC", "capillary_AdEC", "arterial_AdEC", "sub_AdEC", "lymphatic_AdEC")
subsets <- lapply(SUBTYPES, function(s) { k <- colnames(E)[E$subtype == s]; if (length(k)) subset(E, cells = k) else NULL }); names(subsets) <- SUBTYPES
for (s in SUBTYPES) if (!is.null(subsets[[s]])) { fwrite(data.table(barcode = colnames(subsets[[s]])), file.path(OUT_DIR, paste0("subset_", s, "_barcodes.csv"))); saveRDS(subsets[[s]], file.path(OUT_DIR, paste0("subset_", s, ".rds"))) }
md <- as.data.table(E@meta.data, keep.rownames = "barcode")
n_sub <- rbindlist(lapply(c("group_annotation", "group_geo"), function(g) { x <- copy(md); x[, `:=`(label_source = g, grp = get(g))]; x[, .(cells = .N, samples = uniqueN(sample_name), people = uniqueN(person)), by = .(label_source, subtype, tissue, group = grp)] }))[order(label_source, subtype, tissue, group)]
per_sample <- md[, .N, by = .(subtype, sample_name, person, tissue, group_annotation, group_geo)][order(subtype, -N)]; setnames(per_sample, "N", "cells"); print(n_sub[label_source == "group_annotation"])

# ---- 6. THE GENES OF INTEREST per subtype, tissue and group (detection = at least one count)
cnt <- GetAssayData(E, layer = "counts"); lg <- GetAssayData(E, layer = "data")
gene_tab <- rbindlist(lapply(unique(md$subtype), function(s) rbindlist(lapply(c("SAT", "VAT", "SAT+VAT"), function(tt) rbindlist(lapply(c("control", "T2D", "both groups"), function(gr) {
  sel <- md$subtype == s & (tt == "SAT+VAT" | md$tissue == tt) & (gr == "both groups" | md$group_annotation == gr); k <- which(sel); if (!length(k)) return(NULL)
  data.table(subtype = s, tissue = tt, group = gr, cells = length(k), gene = names(G), symbol = unname(G), cells_detected = as.integer(Matrix::rowSums(cnt[unname(G), k, drop = FALSE] > 0)),
             pct_detected = round(100 * Matrix::rowMeans(cnt[unname(G), k, drop = FALSE] > 0), 1), mean_log_expr = round(Matrix::rowMeans(lg[unname(G), k, drop = FALSE]), 3)) }))))))
print(gene_tab[tissue == "SAT+VAT" & group == "both groups", .(subtype, cells, gene, pct_detected, mean_log_expr)][order(gene, subtype)])

# ---- 7. exploratory pseudobulk T2D vs control per subtype (one value per sample; tissue as covariate). Few samples: EXPLORATORY.
de_one <- function(s) { d <- md[subtype == s]; smp <- d[, .(n = .N, group = group_annotation[1], tissue = tissue[1]), by = sample_name][n >= MIN_CELLS_PER_SAMPLE]
  if (sum(smp$group == "T2D") < 2 || sum(smp$group == "control") < 2) return(data.table(subtype = s, status = sprintf("NOT TESTABLE: %d T2D and %d control samples with >= %d cells", sum(smp$group == "T2D"), sum(smp$group == "control"), MIN_CELLS_PER_SAMPLE)))
  ind <- sparseMatrix(i = seq_len(nrow(d)), j = match(d$sample_name, smp$sample_name), x = 1, dims = c(nrow(d), nrow(smp))); keepc <- d$sample_name %in% smp$sample_name
  pb <- as.matrix(cnt[, d$barcode[keepc], drop = FALSE] %*% ind[keepc, , drop = FALSE]); colnames(pb) <- smp$sample_name; grp <- factor(smp$group, c("control", "T2D"))
  y <- DGEList(pb, group = grp); y <- y[filterByExpr(y) | rownames(y) %in% G, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  des <- if (nlevels(factor(smp$tissue)) > 1) model.matrix(~ factor(smp$tissue) + grp) else model.matrix(~ grp); y <- estimateDisp(y, des); tt <- topTags(glmQLFTest(glmQLFit(y, des), coef = ncol(des)), n = Inf)$table
  r <- data.table(subtype = s, gene = rownames(tt), logFC_T2D_vs_control = round(tt$logFC, 3), p_value = tt$PValue, FDR = tt$FDR, T2D_samples = sum(smp$group == "T2D"), control_samples = sum(smp$group == "control"), status = "EXPLORATORY"); r[gene %in% G] }
de_res <- rbindlist(lapply(SUBTYPES, function(s) tryCatch(de_one(s), error = function(e) data.table(subtype = s, status = paste("ERROR:", conditionMessage(e))))), fill = TRUE); print(de_res)

# ---- 8. pictures and save
sv <- function(f, w, h, p) { png(file.path(OUT_DIR, f), width = w, height = h, res = 110); print(p); dev.off() }
sv("umap_all_cells_vascular.png", 1400, 600, DimPlot(seu, group.by = "cluster", label = TRUE, repel = TRUE, label.box = TRUE) + NoLegend() + DimPlot(seu, group.by = "broad_class"))
sv("umap_endothelial_subtypes.png", 1600, 600, DimPlot(E, group.by = "subtype", label = TRUE, repel = TRUE, label.box = TRUE) + DimPlot(E, group.by = "group_annotation") + DimPlot(E, group.by = "tissue"))
sv("dotplot_genes_by_subtype.png", 1000, 600, DotPlot(E, features = unname(G), group.by = "subtype") + RotatedAxis())
sv("featureplot_genes.png", 1800, 900, FeaturePlot(E, features = unname(G), ncol = 4))
sv("dotplot_subtype_markers.png", 1800, 700, DotPlot(E, features = unique(unlist(lapply(SUB, head, 5))), group.by = "subtype") + RotatedAxis())
write_xlsx(Filter(function(x) !is.null(x) && nrow(x) > 0, list(Cells_per_subtype = n_sub, Per_sample = per_sample, Genes_by_subtype = gene_tab, DE_exploratory = de_res, Subcluster_labels = cluster_label, Subcluster_markers_top10 = sub_markers,
  Broad_cluster_scores = as.data.table(round(zb, 2), keep.rownames = "cluster"), Panels = rbindlist(lapply(names(SUB), function(n) data.table(subtype = n, genes = paste(SUB[[n]], collapse = ", ")))),
  Notes = data.table(note = c("Subtype panels are curated by the assistant, not the authors' list; verify with Subcluster_markers_top10.", "Subtypes are assigned per sub-cluster (shallow data make per-cell calls noisy).", "Control = obese without T2D (annotation label); group_geo gives the GEO label; they conflict for 3 people.", "All donors are obese; there is no lean group.")))), file.path(OUT_DIR, "GSE129363_endothelial_subsets.xlsx"))
saveRDS(E, file.path(OUT_DIR, "endothelial_subtyped.rds")); print(list.files(OUT_DIR))
