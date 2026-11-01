# GSE129363: UMAP, cell clusters, where the 8 genes are expressed, and top-50 differentially expressed genes (obese vs obese + T2D).
# For Google Colab (R runtime, High-RAM recommended). Tested on synthetic data in the real file format; NOT yet run on the real matrix.
#
# Why this notebook: the marker-rule check found only 105 endothelial-like cells among 26,350 (0.4%); the data are shallow
# (median 423 UMI, 245 genes per cell), so most endothelial cells miss 3 of 6 markers. Clustering all cells and naming the clusters
# is a more sensitive way to find endothelial cells.

# ---- SETTINGS -----------------------------------------------------------
DATA_DIR <- "/content/gse129363_data"      # the 5 files (downloaded by ID below, or already there)
OUT_DIR  <- "/content/gse129363_umap"
USE_GDOWN <- TRUE
drive_ids <- c("1xVFGimI87DMyyUHoJxbBBr_PXS_tU3bt", "1l77NQm3jVmJ8ieoD9KRpWjGFRaMev5hk", "1imnxlaJKfZMK5WUJas-5IXjQbuHNi79y", "190PY9cnmseA9_xXRD6vBNBeCxXY4tBZS", "1M2Q1NAnSGLhV3HemzSuHUvIy8RF7mTZf")
RULE_CSV <- ""                              # optional: path to GSE129363_cells_one_row_per_barcode.csv from the earlier run (to show the 105 marker-rule cells on the UMAP)
GENE_LIST <- c(ERG = "ERG", CD31 = "PECAM1", `VE-cadherin` = "CDH5", eNOS = "NOS3", VCAM1 = "VCAM1", ICAM1 = "ICAM1", `HIF-1a` = "HIF1A", ANGPT2 = "ANGPT2")
N_PCS <- 20; CLUSTER_RES <- 0.5; USE_HARMONY <- TRUE   # Harmony removes person-to-person differences before clustering
MIN_CELLS_PER_SAMPLE <- 10; MIN_SAMPLES_PER_GROUP <- 2; RELIABLE_SAMPLES_PER_GROUP <- 4; TOP_N <- 50
# -------------------------------------------------------------------------

for (p in c("Matrix", "data.table", "writexl", "Seurat", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("edgeR", quietly = TRUE)) { install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install("edgeR", ask = FALSE, update = FALSE) }
if (USE_HARMONY && !requireNamespace("harmony", quietly = TRUE)) try(install.packages("harmony", repos = "https://cloud.r-project.org"), silent = TRUE)
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl); library(Seurat); library(ggplot2); library(edgeR)})
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); set.seed(1)

# ---- 1. get the five files and find each by name
if (USE_GDOWN) { system("pip install -q gdown"); for (id in drive_ids) system(sprintf("cd %s && gdown --id %s", DATA_DIR, id)) }
all_files <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
find1 <- function(pattern) { h <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]; if (!length(h)) stop("No file matching ", pattern, " in ", DATA_DIR); h[1] }
f_bc <- find1("barcodes"); f_gn <- find1("genes|features"); f_mx <- find1("matrix.*mtx|\\.mtx"); f_an <- find1("annotation")

# ---- 2. load counts, keep the cells the authors kept, attach sample / diabetes / depot
bc <- fread(f_bc, header = FALSE)[[1]]; gn <- fread(f_gn, header = FALSE)
gene <- make.unique(if (ncol(gn) >= 2) gn[[2]] else gn[[1]])
m <- as(readMM(f_mx), "CsparseMatrix"); stopifnot(nrow(m) == length(gene), ncol(m) == length(bc)); rownames(m) <- gene; colnames(m) <- bc
an <- fread(f_an); setnames(an, 1, "barcode"); an <- an[barcode %in% bc]
sample_key <- data.table(sample = c("SVF2-SAT", "SVF3-VAT", "SVF5-SAT", "Ate-04-VAT", "Ate-04-SAT", "Ate-05-VAT", "Ate-05-SAT", "Ate-06-VAT", "Ate-06-SAT", "Ate-08-VAT", "Ate-08-SAT", "Ate-09-VAT", "Ate-09-SAT", "Ate-10_VAT", "Ate-10-SAT", "Ate-11-VAT", "Ate-11-SAT", "Ate-12-VAT", "Ate-12-SAT", "Ate-13-VAT", "Ate-13-SAT", "Ate-14-VAT", "Ate-14-SAT", "Ate-15-VAT", "Ate-15-SAT"), person = c("E02", "E03", "E07", "04", "04", "05", "05", "06", "06", "08", "08", "09", "09", "10", "10", "11", "11", "12", "12", "13", "13", "14", "14", "15", "15"), geo_label = c("diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic"))
md <- merge(an, sample_key, by.x = "SampleName", by.y = "sample", all.x = TRUE, sort = FALSE); stopifnot(!anyNA(md$person))
md[, group_annotation := ifelse(grepl("^non", Condition, ignore.case = TRUE), "control", "disease")]
md[, group_geo := ifelse(geo_label == "diabetic", "disease", "control")]
md <- as.data.frame(md); rownames(md) <- md$barcode
counts <- m[, md$barcode]; counts <- counts[Matrix::rowSums(counts > 0) >= 3, ]; rm(m); gc(verbose = FALSE)
cat("Cells:", ncol(counts), "| genes kept (detected in >= 3 cells):", nrow(counts), "\n")
seu <- CreateSeuratObject(counts, meta.data = md[colnames(counts), c("SampleName", "Condition", "Tissue", "person", "group_annotation", "group_geo")])
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
Gsym <- GENE_LIST[GENE_LIST %in% rownames(seu)]; cat("Genes of interest found:", paste(Gsym, collapse = ", "), "| missing:", paste(setdiff(GENE_LIST, rownames(seu)), collapse = ", "), "\n")

# ---- 3. normalise, reduce, cluster, UMAP
seu <- NormalizeData(seu, verbose = FALSE); seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = union(VariableFeatures(seu), unname(Gsym)), verbose = FALSE)
seu <- RunPCA(seu, npcs = N_PCS, verbose = FALSE); red <- "pca"
if (USE_HARMONY && requireNamespace("harmony", quietly = TRUE)) {
  seu <- tryCatch(harmony::RunHarmony(seu, group.by.vars = "SampleName", verbose = FALSE), error = function(e) { cat("Harmony failed; plain PCA used:", conditionMessage(e), "\n"); seu })
  if ("harmony" %in% names(seu@reductions)) red <- "harmony"
}
cat("Clustering on:", red, "\n")
seu <- FindNeighbors(seu, reduction = red, dims = 1:N_PCS, verbose = FALSE); seu <- FindClusters(seu, resolution = CLUSTER_RES, verbose = FALSE)
seu <- RunUMAP(seu, reduction = red, dims = 1:N_PCS, verbose = FALSE)
seu$cluster <- paste0("C", as.character(seu$seurat_clusters)); cat("Clusters:", length(unique(seu$cluster)), "\n")
if (nzchar(RULE_CSV) && file.exists(RULE_CSV)) { rc <- fread(RULE_CSV); seu$marker_rule_endothelial <- rc$endo_like[match(colnames(seu), rc$barcode)] } else seu$marker_rule_endothelial <- NA

# ---- 4. name the clusters from marker genes (check this against the marker table!)
PANELS <- list(
  endothelial = c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB","CDH5","KDR","RAMP2"), lymphatic = c("PROX1","LYVE1","FLT4","CCL21"),
  fibroblast_progenitor = c("PDGFRA","DCN","LUM","COL1A1","COL1A2","CD34"), mural_pericyte = c("ACTA2","RGS5","MYH11","PDGFRB","NOTCH3"),
  macrophage_monocyte = c("CD14","CD68","CD163","MRC1","LYZ","CSF1R"), T_NK_cell = c("CD3D","CD3E","IL7R","NKG7","GNLY"),
  B_cell = c("MS4A1","CD79A"), mast_cell = c("TPSAB1","KIT"), adipocyte = c("ADIPOQ","PLIN1","LIPE"))
lg <- GetAssayData(seu, layer = "data")
cl_mean <- sapply(sort(unique(seu$cluster)), function(cl) Matrix::rowMeans(lg[, seu$cluster == cl, drop = FALSE]))
z <- t(scale(t(cl_mean)))                                            # per gene, relative to the other clusters
panel_score <- sapply(PANELS, function(g) { g <- g[g %in% rownames(z)]; if (length(g)) colMeans(z[g, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(z)) })
lab <- apply(panel_score, 1, function(x) if (all(is.na(x)) || max(x, na.rm = TRUE) < 0.3) "unassigned" else colnames(panel_score)[which.max(x)])
seu$cell_type <- unname(lab[seu$cluster])
cluster_tab <- as.data.table(panel_score, keep.rownames = "cluster"); cluster_tab[, label := lab[cluster]]; cluster_tab[, cells := as.integer(table(seu$cluster)[cluster])]
cluster_tab[, (names(PANELS)) := lapply(.SD, round, 2), .SDcols = names(PANELS)]
print(cluster_tab[order(label, -cells)])
cat("\nCells per cell type:\n"); print(sort(table(seu$cell_type), decreasing = TRUE))
markers <- tryCatch(as.data.table(FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE))[, head(.SD, 10), by = cluster], error = function(e) data.table())

# ---- 5. pictures: UMAP and where the genes are
sv <- function(f, w, h, p) { png(file.path(OUT_DIR, f), width = w, height = h, res = 110); print(p); dev.off(); print(p) }
sv("umap_clusters_celltype.png", 1700, 600, DimPlot(seu, group.by = "cluster", label = TRUE) + DimPlot(seu, group.by = "cell_type", label = TRUE))
sv("umap_group_depot_sample.png", 2100, 600, DimPlot(seu, group.by = "group_annotation") + DimPlot(seu, group.by = "Tissue") + DimPlot(seu, group.by = "SampleName") + NoLegend())
if (any(!is.na(seu$marker_rule_endothelial))) sv("umap_marker_rule_cells.png", 800, 600, DimPlot(seu, group.by = "marker_rule_endothelial", cols = c(`FALSE` = "grey85", `TRUE` = "red")) + ggtitle("Marker-rule endothelial cells (red)"))
sv("umap_gene_list.png", 1800, 900, FeaturePlot(seu, features = unname(Gsym), ncol = 4))
sv("dotplot_gene_list_by_cluster.png", 1000, 700, DotPlot(seu, features = unname(Gsym), group.by = "cluster") + RotatedAxis())
sv("dotplot_gene_list_by_celltype.png", 1000, 600, DotPlot(seu, features = unname(Gsym), group.by = "cell_type") + RotatedAxis())

# ---- 6. where does each gene sit? mean expression and % cells, per cluster and per cell type
gl <- lg[unname(Gsym), , drop = FALSE]
by_cl <- function(col) rbindlist(lapply(sort(unique(seu[[col]][, 1])), function(v) { k <- which(seu[[col]][, 1] == v)
  data.table(group_by = col, value = v, cells = length(k), gene = names(Gsym), symbol = unname(Gsym),
             mean_log_expr = round(Matrix::rowMeans(gl[, k, drop = FALSE]), 3), pct_cells_expressing = round(100 * Matrix::rowMeans(gl[, k, drop = FALSE] > 0), 1)) }))
gene_where <- rbind(by_cl("cluster"), by_cl("cell_type"))
gene_top <- gene_where[group_by == "cell_type", .SD[which.max(mean_log_expr)], by = gene][, .(gene, symbol, highest_in_cell_type = value, mean_log_expr, pct_cells_expressing)]
cat("\nCell type where each gene is highest:\n"); print(gene_top)

# ---- 7. endothelial cells: counts per group, depot, person (from clusters, not from the strict marker rule)
ENDO <- "endothelial"
E <- as.data.table(seu@meta.data, keep.rownames = "barcode")[cell_type == ENDO]
cat("\nEndothelial cells (cluster-based):", nrow(E), "| marker-rule cells inside them:", sum(E$marker_rule_endothelial, na.rm = TRUE), "of", sum(seu$marker_rule_endothelial, na.rm = TRUE), "\n")
endo_count <- rbindlist(lapply(c("group_annotation", "group_geo"), function(g) { x <- copy(E); x[, `:=`(label_source = g, grp = get(g))]
  x[, .(endothelial_cells = .N, samples = uniqueN(SampleName), people = uniqueN(person)), by = .(label_source, group = grp, depot = Tissue)] }))
print(endo_count[order(label_source, group, depot)])
per_sample <- E[, .(endothelial_cells = .N), by = .(SampleName, person, depot = Tissue, group_annotation, group_geo)][order(-endothelial_cells)]
per_sample[, enough_for_replicate := endothelial_cells >= MIN_CELLS_PER_SAMPLE]; print(per_sample)
if (nrow(E) < 100) cat("WARNING: very few endothelial cells; check the cluster table above, the cluster naming may have missed them.\n")

# ---- 8. differential expression, obese (control) vs obese + T2D (disease): one value per sample (pseudobulk, edgeR)
cnt <- GetAssayData(seu, layer = "counts")
run_de <- function(src, depot_sel) {
  d <- E[SampleName %in% per_sample[enough_for_replicate == TRUE, SampleName]]; if (depot_sel != "SAT+VAT") d <- d[Tissue == depot_sel]
  smp <- d[, .(n_cells = .N, group = get(src)[1], depot = Tissue[1]), by = SampleName]
  n_c <- sum(smp$group == "control"); n_d <- sum(smp$group == "disease"); label <- paste(src, depot_sel, sep = " | ")
  st <- data.table(analysis = label, control_samples = n_c, disease_samples = n_d, endothelial_cells_used = nrow(d))
  if (min(n_c, n_d) < MIN_SAMPLES_PER_GROUP) { st[, status := sprintf("NOT TESTABLE: %d control, %d disease samples with >= %d endothelial cells", n_c, n_d, MIN_CELLS_PER_SAMPLE)]; return(list(st = st, res = NULL)) }
  ind <- sparseMatrix(i = seq_len(nrow(d)), j = match(d$SampleName, smp$SampleName), x = 1, dims = c(nrow(d), nrow(smp)))
  pb <- as.matrix(cnt[, d$barcode, drop = FALSE] %*% ind); colnames(pb) <- smp$SampleName
  grp <- factor(smp$group, c("control", "disease"))
  y <- DGEList(pb, group = grp); y <- y[rowSums(y$counts) >= 3 | rownames(y) %in% Gsym, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  des <- if (depot_sel == "SAT+VAT" && length(unique(smp$depot)) > 1) model.matrix(~ factor(smp$depot) + grp) else model.matrix(~ grp)
  y <- estimateDisp(y, des); tt <- topTags(glmQLFTest(glmQLFit(y, des), coef = ncol(des)), n = Inf)$table
  res <- data.table(analysis = label, gene = rownames(tt), logFC_disease_vs_control = round(tt$logFC, 3), logCPM = round(tt$logCPM, 2), p_value = tt$PValue, p_adj_BH = tt$FDR, control_samples = n_c, disease_samples = n_d)
  res[, in_gene_list := gene %in% GENE_LIST]
  st[, status := ifelse(min(n_c, n_d) < RELIABLE_SAMPLES_PER_GROUP, sprintf("EXPLORATORY: %d control vs %d disease samples", n_c, n_d), "tested")]
  list(st = st, res = res)
}
runs <- list(); for (src in c("group_annotation", "group_geo")) for (dp in c("SAT", "VAT", "SAT+VAT")) runs[[length(runs) + 1]] <- tryCatch(run_de(src, dp), error = function(e) list(st = data.table(analysis = paste(src, dp), status = paste("ERROR:", conditionMessage(e))), res = NULL))
status_tab <- rbindlist(lapply(runs, `[[`, "st"), fill = TRUE); print(status_tab)
de_all <- rbindlist(lapply(runs, `[[`, "res"), fill = TRUE)
top50 <- if (nrow(de_all)) de_all[, head(.SD[order(p_value)], TOP_N), by = analysis][, rank := seq_len(.N), by = analysis] else data.table()
gl_de <- if (nrow(de_all)) de_all[in_gene_list == TRUE][order(analysis, p_value)] else data.table()
overlap50 <- if (nrow(top50)) top50[in_gene_list == TRUE, .(n_in_top50 = .N, genes = paste(gene, collapse = ", ")), by = analysis] else data.table()
if (nrow(top50)) { cat("\nTop 10 of the", TOP_N, "(first testable analysis):\n"); print(head(top50[analysis == top50$analysis[1], .(rank, gene, logFC_disease_vs_control, p_value, p_adj_BH, in_gene_list)], 10)); print(overlap50) }

# ---- 9. supporting cell-level top 50 (Wilcoxon; every cell treated as independent -> p-values too small, use for ranking only)
cell_de <- rbindlist(lapply(c("group_annotation", "group_geo"), function(src) {
  Es <- subset(seu, subset = cell_type == ENDO); Idents(Es) <- Es[[src]][, 1]
  if (!all(c("disease", "control") %in% levels(Idents(Es))) || min(table(Idents(Es))) < 3) return(NULL)
  x <- as.data.table(FindMarkers(Es, ident.1 = "disease", ident.2 = "control", logfc.threshold = 0, min.pct = 0.05, verbose = FALSE), keep.rownames = "gene")
  x[, `:=`(label_source = src, in_gene_list = gene %in% GENE_LIST)][order(p_val)][, rank := .I][]
}), fill = TRUE)
cell_top50 <- if (nrow(cell_de)) cell_de[rank <= TOP_N] else data.table()

# ---- 10. save
notes <- data.table(note = c("Cell types are assigned to clusters by marker genes (cluster_labels sheet). Check them against Cluster_markers_top10 before trusting any count.",
  "Control = obese without T2D; disease = obese with T2D; no lean group. Diabetes labels differ between the annotation file and GEO for people E07, 04, 09, so both are analysed.",
  "Pseudobulk edgeR: one value per sample. Few samples per group means EXPLORATORY results; NOT TESTABLE means too few replicates.",
  "Cell-level top 50 is ranking only: nuclei/cells from one person are not independent.",
  "Data are shallow (median 423 UMI per cell), so low-expressed genes (e.g. NOS3, ANGPT2, ERG) are often undetected."))
sheets <- list(Cluster_labels = cluster_tab, Cells_per_celltype = as.data.table(sort(table(seu$cell_type), decreasing = TRUE)), Gene_where = gene_where, Gene_top_celltype = gene_top,
  Endothelial_counts = endo_count, Endothelial_per_sample = per_sample, DE_status = status_tab, Top50_pseudobulk = top50, Gene_list_DE = gl_de, Overlap_top50 = overlap50,
  Top50_cell_level = cell_top50, Cluster_markers_top10 = markers, Notes = notes)
write_xlsx(Filter(function(x) !is.null(x) && nrow(x) > 0, sheets), file.path(OUT_DIR, "GSE129363_umap_clusters_DE_summary.xlsx"))
if (nrow(de_all)) fwrite(de_all, file.path(OUT_DIR, "DE_pseudobulk_all_genes.csv")); if (nrow(top50)) fwrite(top50, file.path(OUT_DIR, "top50_DE_genes.csv"))
saveRDS(seu, file.path(OUT_DIR, "gse129363_seurat.rds")); print(list.files(OUT_DIR))
