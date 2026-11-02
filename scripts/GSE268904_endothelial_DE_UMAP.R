# GSE268904 endothelial nuclei: UMAP, clusters, expression of the 8 genes of interest per cluster,
# T2D vs control differential expression (donor-level pseudobulk), top-50 DE genes, overlap with the gene list.
# For Google Colab (R runtime). Tested on synthetic data only; NOT yet run on the real archive.
#
# Genes of interest (handwritten list -> symbol): ERG, CD31=PECAM1, VE-cadherin=CDH5, eNOS=NOS3, VCAM1, ICAM1, HIF-1a=HIF1A, ANGPT2

# ---- SETTINGS -----------------------------------------------------------
TAR_PATH <- ""                    # "" = auto-find in /content, else download by Drive ID
WORK_DIR <- "/content/gse268904_raw"; OUT_DIR <- "/content/gse268904_endo"
MIN_GENES <- 200; MIN_UMI <- 500; ENDO_MIN_MARKERS <- 3; ENDO_OVER_OTHER <- 2      # pre-selection of endothelial nuclei (marker rule)
GENES <- c(ERG = "ERG", CD31 = "PECAM1", `VE-cadherin` = "CDH5", eNOS = "NOS3", VCAM1 = "VCAM1", ICAM1 = "ICAM1", `HIF-1a` = "HIF1A", ANGPT2 = "ANGPT2")
CONTAM_Z <- 1.5                  # a cluster is flagged as non-endothelial when an immune / mural / lymphatic / neural marker score is at least this many SDs above the other clusters
N_PCS <- 20; CLUSTER_RES <- 0.4; TOP_N <- 50; MIN_ENDO_NUCLEI <- 50
USE_HARMONY <- TRUE               # remove person-to-person batch effect before clustering (needs the 'harmony' package; falls back if missing)
# -------------------------------------------------------------------------

for (p in c("Matrix", "data.table", "writexl", "Seurat", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("edgeR", quietly = TRUE)) { install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install("edgeR", ask = FALSE, update = FALSE) }
if (USE_HARMONY && !requireNamespace("harmony", quietly = TRUE)) try(install.packages("harmony", repos = "https://cloud.r-project.org"), silent = TRUE)
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl); library(Seurat); library(ggplot2); library(edgeR)})
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# ---- 1. archive and sample table
if (TAR_PATH == "") { h <- list.files("/content", pattern = "GSE268904.*\\.tar$", recursive = TRUE, full.names = TRUE); if (length(h)) TAR_PATH <- h[1] }
if (TAR_PATH == "" || !file.exists(TAR_PATH)) { system("pip install -q gdown"); system("gdown 1mKaM2c4473gV3zCoCV3F1ZknlCIDlAwq -O /content/GSE268904_RAW.tar"); TAR_PATH <- "/content/GSE268904_RAW.tar" }
stopifnot(file.exists(TAR_PATH)); untar(TAR_PATH, exdir = WORK_DIR)
files <- list.files(WORK_DIR, recursive = TRUE, full.names = TRUE)
meta <- data.table(gsm = sprintf("GSM%d", 8302405:8302418),
  donor = c("9S","10S","13S","17S","SAT22","SAT54","SAT1","SAT18","SAT34","SAT39","SAT4","SAT51","SAT58","SAT60"),
  group = c("obese","obese diabetic","lean","obese diabetic","obese diabetic","lean","obese","obese","obese","obese","obese","obese","lean","lean"),
  bmi = c(39.5,61.5,23,52.6,39.7,22.3,41.7,61.3,62.5,61.2,49.7,57,24.3,23.5),
  sex = c("female","female","female","female","male","female","male","female","male","male","male","female","male","male"))
PANELS <- list(endothelial = c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB","CDH5","KDR"), lymphatic = c("PROX1","LYVE1","FLT4"),
  adipocyte = c("ADIPOQ","PLIN1","LIPE","LEP"), fibroblast = c("PDGFRA","DCN","LUM","COL1A1"), immune = c("PTPRC","CD163","MRC1","CD14"), mural = c("ACTA2","RGS5","MYH11","PDGFRB"))
CORE6 <- c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB")

# ---- 2. per person: keep only endothelial nuclei (full gene counts), build one Seurat object
objs <- list(); n_all <- list()
for (i in seq_len(nrow(meta))) {
  g <- meta$gsm[i]; pick <- function(t) files[grepl(g, files) & grepl(t, files)]
  fb <- pick("barcodes"); ff <- pick("features|genes"); fm <- pick("matrix")
  if (length(fb) != 1 || length(ff) != 1 || length(fm) != 1) { cat("Skipping", g, "\n"); next }
  m <- as(readMM(fm), "CsparseMatrix"); feat <- fread(ff, header = FALSE); gene <- make.unique(if (ncol(feat) >= 2) feat[[2]] else feat[[1]])
  bc <- fread(fb, header = FALSE)[[1]]; umi <- Matrix::colSums(m)
  sc <- sapply(PANELS, function(p) { j <- match(p, gene); j <- j[!is.na(j)]; if (!length(j)) rep(0, ncol(m)) else Matrix::colSums(m[j, , drop = FALSE]) / pmax(umi, 1) * 1e4 })
  jc <- match(CORE6, gene); jc <- jc[!is.na(jc)]; nd <- Matrix::colSums(m[jc, , drop = FALSE] > 0)
  oth <- apply(sc[, colnames(sc) != "endothelial", drop = FALSE], 1, max); top <- colnames(sc)[max.col(sc, ties.method = "first")]
  e <- which(umi >= MIN_UMI & nd >= ENDO_MIN_MARKERS & top == "endothelial" & sc[, "endothelial"] >= ENDO_OVER_OTHER * oth)
  n_all[[i]] <- data.table(donor = meta$donor[i], nuclei = sum(umi >= MIN_UMI), endothelial = length(e))
  cat(sprintf("%-6s %-15s %5d endothelial of %6d nuclei\n", meta$donor[i], meta$group[i], length(e), sum(umi >= MIN_UMI)))
  if (length(e) < 10) next
  x <- m[, e, drop = FALSE]; colnames(x) <- paste0(meta$donor[i], "_", bc[e]); rownames(x) <- gene
  s <- CreateSeuratObject(x, project = meta$donor[i]); s$donor <- meta$donor[i]; s$group <- meta$group[i]; s$bmi <- meta$bmi[i]; s$sex <- meta$sex[i]
  objs[[meta$donor[i]]] <- s; rm(m, x); gc(verbose = FALSE)
}
seu <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
seu <- JoinLayers(seu); rm(objs); gc(verbose = FALSE)
seu$t2d <- ifelse(seu$group == "obese diabetic", "T2D", "non-T2D")
seu <- subset(seu, subset = nFeature_RNA >= MIN_GENES)
cat("Endothelial nuclei in the object:", ncol(seu), "from", length(unique(seu$donor)), "people\n")
missing <- setdiff(GENES, rownames(seu)); cat("Genes of interest missing from the data:", if (length(missing)) paste(missing, collapse = ", ") else "none", "\n")

# ---- 3. normalise, reduce, cluster, UMAP
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = union(VariableFeatures(seu), GENES[GENES %in% rownames(seu)]), verbose = FALSE)
seu <- RunPCA(seu, npcs = N_PCS, verbose = FALSE)
red <- "pca"
if (USE_HARMONY && requireNamespace("harmony", quietly = TRUE) && length(unique(seu$donor)) > 1) {
  seu <- tryCatch(harmony::RunHarmony(seu, group.by.vars = "donor", verbose = FALSE), error = function(e) { cat("Harmony failed, using plain PCA:", conditionMessage(e), "\n"); seu })
  if ("harmony" %in% names(seu@reductions)) red <- "harmony"
}
cat("Clustering on:", red, "\n")
seu <- FindNeighbors(seu, reduction = red, dims = 1:N_PCS, verbose = FALSE)
seu <- FindClusters(seu, resolution = CLUSTER_RES, verbose = FALSE)
seu <- RunUMAP(seu, reduction = red, dims = 1:N_PCS, verbose = FALSE)
seu$cluster <- paste0("C", as.character(seu$seurat_clusters))
cat("Clusters found:", length(unique(seu$cluster)), "\n")
print(table(seu$cluster, seu$group))

G <- GENES[GENES %in% rownames(seu)]
pdf(file.path(OUT_DIR, "umap_overview.pdf"), width = 12, height = 5); print(DimPlot(seu, group.by = "cluster", label = TRUE) + DimPlot(seu, group.by = "group") + DimPlot(seu, group.by = "donor")); dev.off()
png(file.path(OUT_DIR, "umap_clusters_group_donor.png"), width = 1800, height = 560, res = 110); print(DimPlot(seu, group.by = "cluster", label = TRUE) + DimPlot(seu, group.by = "group") + DimPlot(seu, group.by = "donor")); dev.off()
png(file.path(OUT_DIR, "umap_gene_panel.png"), width = 1600, height = 800, res = 110); print(FeaturePlot(seu, features = unname(G), ncol = 4)); dev.off()
png(file.path(OUT_DIR, "dotplot_gene_panel_by_cluster.png"), width = 900, height = 600, res = 110); print(DotPlot(seu, features = unname(G), group.by = "cluster") + RotatedAxis()); dev.off()
print(DimPlot(seu, group.by = "cluster", label = TRUE)); print(DimPlot(seu, group.by = "group"))

# ---- 3b. cluster QC: which clusters are NOT endothelial? (the first real run showed T cells, smooth muscle, lymphatic and neural-like clusters inside the endothelial set)
CONTAM <- list(immune = c("PTPRC","DOCK2","IKZF1","ITK","CD247","IL7R","RUNX1","CD96","ZEB2"), mural_smooth_muscle = c("MYH11","CARMN","PRKG1","RGS6","ACTA2","RGS5","PDGFRB"),
               lymphatic = c("PROX1","MMRN1","PKHD1L1","PROX1-AS1","FLT4"), neural_like = c("CSMD3","CNTNAP2","GPC5","NRXN1","PLP1"))
lg0 <- GetAssayData(seu, layer = "data")
cm <- sapply(sort(unique(seu$cluster)), function(cl) Matrix::rowMeans(lg0[, seu$cluster == cl, drop = FALSE])); zc <- t(scale(t(cm)))
cz <- do.call(cbind, lapply(CONTAM, function(g) { g <- g[g %in% rownames(zc)]; if (length(g)) colMeans(zc[g, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(zc)) })); colnames(cz) <- names(CONTAM); rownames(cz) <- colnames(zc)
cluster_qc <- as.data.table(round(cz, 2), keep.rownames = "cluster"); .tab <- table(seu$cluster); cluster_qc[, nuclei := as.integer(.tab[cluster_qc$cluster])]
cluster_qc[, flagged_as := apply(cz, 1, function(x) if (all(is.na(x)) || max(x, na.rm = TRUE) < CONTAM_Z) "keep (endothelial-like)" else paste0("exclude: ", colnames(cz)[which.max(x)]))]
print(cluster_qc[order(cluster)])
core_clusters <- cluster_qc[flagged_as == "keep (endothelial-like)", cluster]; seu$core <- seu$cluster %in% core_clusters
cat("Core endothelial clusters:", paste(core_clusters, collapse = ", "), "| nuclei kept:", sum(seu$core), "of", ncol(seu), "\n")

# ---- 4. where is each gene expressed? mean expression and % of nuclei expressing, per cluster
lg <- GetAssayData(seu, layer = "data")[unname(G), , drop = FALSE]
gene_by_cluster <- rbindlist(lapply(unique(seu$cluster), function(cl) { k <- which(seu$cluster == cl)
  data.table(cluster = cl, nuclei = length(k), gene = names(G), symbol = unname(G),
             mean_log_expr = round(Matrix::rowMeans(lg[, k, drop = FALSE]), 3), pct_expressing = round(100 * Matrix::rowMeans(lg[, k, drop = FALSE] > 0), 1)) }))
top_cluster <- gene_by_cluster[, .SD[which.max(mean_log_expr)], by = gene][, .(gene, symbol, highest_in_cluster = cluster, mean_log_expr, pct_expressing)]
cluster_markers <- tryCatch({ mk <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE); as.data.table(mk)[, head(.SD, 10), by = cluster] }, error = function(e) data.table())
cluster_comp <- as.data.table(as.data.frame.matrix(table(seu$cluster, seu$group)), keep.rownames = "cluster")
print(top_cluster)

# ---- 5. differential expression T2D vs control: donor-level pseudobulk (edgeR), each person counts once, sex as a covariate
cnt_all <- GetAssayData(seu, layer = "counts")
drop_genes <- grepl("^MT-|^RP[SL][0-9]|^MALAT1$", rownames(cnt_all))                 # technical / housekeeping genes are not tested
run_de <- function(cells, analysis) {
  md <- as.data.table(seu@meta.data[cells, ], keep.rownames = "cell"); cnt <- cnt_all[!drop_genes, cells, drop = FALSE]
  donors <- unique(md$donor); pb <- sapply(donors, function(d) Matrix::rowSums(cnt[, md$donor == d, drop = FALSE])); colnames(pb) <- donors
  di <- unique(md[, .(donor, group, sex)])[match(donors, unique(md$donor))]; n_by <- md[, .N, by = donor]; ok <- n_by[N >= MIN_ENDO_NUCLEI, donor]
  one <- function(case, ctrl, label) {
    d <- di[donor %in% ok & group %in% c(case, ctrl)]; if (sum(d$group %in% case) < 2 || sum(d$group %in% ctrl) < 2) return(NULL)
    grp <- factor(ifelse(d$group %in% case, "T2D", "control"), c("control", "T2D")); y <- DGEList(pb[, d$donor], group = grp)
    y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
    sx <- factor(d$sex); des <- if (nlevels(sx) > 1) model.matrix(~ sx + grp) else model.matrix(~ grp)              # sex adjusted
    y <- estimateDisp(y, des); tt <- topTags(glmQLFTest(glmQLFit(y, des), coef = ncol(des)), n = Inf)$table
    data.table(analysis = analysis, contrast = label, gene = rownames(tt), logFC = round(tt$logFC, 3), logCPM = round(tt$logCPM, 2), p_value = tt$PValue, FDR = tt$FDR,
               n_T2D = sum(d$group %in% case), n_control = sum(d$group %in% ctrl), sex_adjusted = nlevels(sx) > 1)
  }
  rbindlist(list(one("obese diabetic", "obese", "T2D vs obese (main)"), one("obese diabetic", "lean", "T2D vs lean"), one("obese diabetic", c("obese", "lean"), "T2D vs all non-diabetic")), fill = TRUE)
}
de <- rbindlist(list(run_de(colnames(seu)[seu$core], "core endothelial clusters"), run_de(colnames(seu), "all pre-selected nuclei")), fill = TRUE)
de[, in_gene_list := gene %in% GENES]; de[, symbol_of := ifelse(in_gene_list, names(GENES)[match(gene, GENES)], NA_character_)]
# is there any signal at all? genes with p < 0.05 vs the 5% expected by chance, and genomic inflation (lambda ~ 1 = no inflation)
signal <- de[, .(genes_tested = .N, p_below_0.05 = sum(p_value < 0.05), expected_by_chance = round(0.05 * .N), fdr_below_0.05 = sum(FDR < 0.05), lambda = round(median(qchisq(p_value, 1, lower.tail = FALSE)) / 0.4549, 2)), by = .(analysis, contrast)]
print(signal)
main <- de[analysis == "core endothelial clusters" & contrast == "T2D vs obese (main)"][order(p_value)]
top50 <- head(main, TOP_N); top50[, rank := .I]
cat("\nTop", TOP_N, "genes (core endothelial clusters, T2D vs obese, sex-adjusted), first 15:\n"); print(head(top50[, .(rank, gene, logFC, p_value, FDR, in_gene_list)], 15))
# composition: share of each cluster within each person (a difference in cell-type mix looks like a difference in expression)
comp <- as.data.table(as.data.frame.matrix(prop.table(table(seu$donor, seu$cluster), 1) * 100), keep.rownames = "donor"); comp <- merge(unique(meta[, .(donor, group, sex)]), comp, by = "donor")

# ---- 6. supporting cell-level test (Wilcoxon). Nuclei from one person are NOT independent, so p-values are inflated: use for ranking only. Core clusters only; technical and sex-linked genes marked.
SEXGENES <- c("XIST", "TSIX", "UTY", "DDX3Y", "KDM5D", "EIF1AY", "RPS4Y1", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "TTTY14", "LINC00278")
Es <- subset(seu, subset = core); Idents(Es) <- "group"
cell_de <- tryCatch({ x <- FindMarkers(Es, ident.1 = "obese diabetic", ident.2 = "obese", logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  x <- as.data.table(x, keep.rownames = "gene"); x <- x[!grepl("^MT-|^RP[SL][0-9]|^MALAT1$", gene)]; x[, in_gene_list := gene %in% GENES]; x[, sex_linked := gene %in% SEXGENES]; x[order(p_val)] }, error = function(e) { cat("Cell-level test failed:", conditionMessage(e), "\n"); data.table() })
cell_top50 <- head(cell_de[sex_linked == FALSE], TOP_N)

# ---- 7. overlap with the gene list
panel_de <- de[gene %in% GENES][, symbol_of := names(GENES)[match(gene, GENES)]][order(contrast, p_value)]
overlap <- data.table(
  set = c("Gene list (8)", "Top 50 pseudobulk (T2D vs obese)", "Top 50 cell-level", "Pseudobulk genes with p < 0.05 (T2D vs obese)", "Pseudobulk genes with FDR < 0.05",
          "Gene list AND top 50 pseudobulk", "Gene list AND top 50 cell-level", "Gene list AND pseudobulk p < 0.05", "Gene list AND pseudobulk FDR < 0.05"),
  n = c(length(GENES), nrow(top50), nrow(cell_top50), sum(main$p_value < 0.05), sum(main$FDR < 0.05),
        sum(top50$gene %in% GENES), sum(cell_top50$gene %in% GENES), sum(main$p_value < 0.05 & main$gene %in% GENES), sum(main$FDR < 0.05 & main$gene %in% GENES)),
  genes = c(paste(GENES, collapse = ", "), "", "", "", "", paste(intersect(top50$gene, GENES), collapse = ", "), paste(intersect(cell_top50$gene, GENES), collapse = ", "),
            paste(main[p_value < 0.05 & gene %in% GENES, gene], collapse = ", "), paste(main[FDR < 0.05 & gene %in% GENES, gene], collapse = ", ")))
print(overlap)

vol <- copy(main); vol[, `:=`(mlog10p = -log10(p_value), rk = .I)]; vol[, lab := ifelse(in_gene_list | rk <= 10, gene, "")]
vp <- ggplot(vol, aes(logFC, mlog10p, colour = in_gene_list)) + geom_point(alpha = 0.5, size = 1) + scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#d1495b"), name = "gene list") +
  geom_text(aes(label = lab), size = 3, vjust = -0.6, check_overlap = TRUE, show.legend = FALSE) + labs(title = "T2D vs obese, core endothelial clusters (donor-level pseudobulk, sex-adjusted)", x = "log2 fold change (T2D / obese)", y = "-log10 p")
png(file.path(OUT_DIR, "volcano_T2D_vs_obese.png"), width = 1000, height = 700, res = 110); print(vp); dev.off(); print(vp)

# ---- 8. save
per_donor <- as.data.table(seu@meta.data)[, .(preselected_nuclei = .N, core_endothelial_nuclei = sum(core)), by = .(donor, group, sex)]
write_xlsx(list(Overlap_with_gene_list = overlap, Top50_pseudobulk = top50, Gene_list_DE = panel_de, Gene_by_cluster = gene_by_cluster, Gene_top_cluster = top_cluster,
  Cluster_markers_top10 = cluster_markers, Cluster_by_group = cluster_comp, Top50_cell_level = cell_top50, All_DE_pseudobulk = de, Signal_check = signal, Cluster_QC = cluster_qc, Composition_by_donor = comp, Per_donor = per_donor,
  Notes = data.table(note = c("Endothelial nuclei were pre-selected by a marker rule (not a formal annotation); clusters then show sub-populations or leftover contaminants.",
   "Pseudobulk edgeR: one value per person. Only 3 people have T2D, so power is very low; use ranks and fold changes as hypotheses.",
   "Cell-level Wilcoxon treats every nucleus as independent (pseudoreplication): p-values are far too small; use only as a supporting list.",
   "Control = obese non-diabetic (main) and lean; obesity and diabetes are partly confounded (BMI).",
   "Clusters can be driven by person if Harmony is off or fails; check Cluster_by_group and the donor UMAP.",
   "Version 2: clusters showing immune, smooth-muscle, lymphatic or neural-like signatures are excluded from the core analysis (Cluster_QC); sex is a covariate; MT-, ribosomal and MALAT1 genes are not tested; sex-linked genes are removed from the cell-level top 50.",
   "Signal_check: if the number of genes with p < 0.05 is not above the number expected by chance and no gene has FDR < 0.05, there is no detectable difference.")),
  Gene_names = data.table(handwritten = names(GENES), symbol = unname(GENES))), file.path(OUT_DIR, "endothelial_DE_UMAP_summary.xlsx"))
fwrite(de, file.path(OUT_DIR, "DE_pseudobulk_all_contrasts.csv")); fwrite(top50, file.path(OUT_DIR, "top50_DE_genes.csv")); fwrite(overlap, file.path(OUT_DIR, "overlap_with_gene_list.csv"))
saveRDS(seu, file.path(OUT_DIR, "endothelial_seurat.rds"))
cat("Saved to", OUT_DIR, ":\n"); print(list.files(OUT_DIR))
