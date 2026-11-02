# GSE155960 (Hildreth et al. 2021, Nature Immunology): criteria check for the file you shared (GSE155960_RAW.tar, 166 MB).
# For Google Colab (R runtime, High-RAM) or the pipeline. Tested on a synthetic tar in the real layout; NOT yet run on the real file.
# The tar holds 12 per-sample files (GSM47171xx_CD45N-L1.csv.gz ...): raw gene counts for every gene and cell, 3 lean + 3 obese donors, CD45-negative and CD45-positive sorts.
#
# Problem statement: control (lean or obese without diabetes) vs type 2 diabetes | endothelial cells > 10,000 | human adipose, subcutaneous AND visceral | gene list.

# ---- SETTINGS -----------------------------------------------------------
TAR_PATH <- ""                         # "" = search /content and the working folder, else download by Drive ID
OUT_DIR <- "/content/gse155960_check"; WORK_DIR <- "/content/gse155960_raw"
GENE_LIST <- c(PECAM1 = "PECAM1", CDH5 = "CDH5", NOS3 = "NOS3", HIF1A = "HIF1A", ERG = "ERG", VCAM1 = "VCAM1", ICAM1 = "ICAM1", ANGPT2 = "ANGPT2")
MIN_ENDO <- 10000; MIN_UMI <- 200; ENDO_MIN_MARKERS <- 3; ENDO_OVER_OTHER <- 2
# -------------------------------------------------------------------------
for (p in c("Matrix", "data.table", "writexl", "R.utils")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl)}); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
if (TAR_PATH == "") { h <- list.files(c("/content", "."), pattern = "GSE155960.*\\.tar$", recursive = TRUE, full.names = TRUE); if (length(h)) TAR_PATH <- h[1] }
if (TAR_PATH == "" || !file.exists(TAR_PATH)) { system("pip install -q gdown"); system("gdown 1sGXAcXDWEqAqaH2p6bfgSYnBdmbn0U7y -O /content/GSE155960_RAW.tar"); TAR_PATH <- "/content/GSE155960_RAW.tar" }
stopifnot(file.exists(TAR_PATH)); untar(TAR_PATH, exdir = WORK_DIR)
files <- list.files(WORK_DIR, pattern = "\\.csv(\\.gz)?$", recursive = TRUE, full.names = TRUE); cat("Sample files found:", length(files), "\n"); print(basename(files))

# ---- 1. who is who (GEO sample records: lean/obese, donor, CD45 sort). Depot and diabetes status are NOT in the records.
sample_key <- data.table(gsm = sprintf("GSM%d", 4717152:4717163),
  sort = rep(c("CD45 negative", "CD45 positive"), each = 6), phenotype = rep(rep(c("lean", "obese"), each = 3), 2), donor = rep(1:3, 4))
sample_key[, file_tag := sprintf("CD45%s-%s%d", ifelse(sort == "CD45 negative", "N", "P"), ifelse(phenotype == "lean", "L", "O"), donor)]
match_row <- function(f) { i <- which(sapply(sample_key$gsm, function(g) grepl(g, basename(f))))
  if (!length(i)) i <- which(sapply(sample_key$file_tag, function(g) grepl(g, basename(f), fixed = TRUE))); if (length(i)) i[1] else NA_integer_ }

# ---- 2. read each sample (dense gene x cell table -> sparse), count cells, find endothelial cells, check the gene list
PANELS <- list(endothelial = c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB","CDH5","KDR","RAMP2","ACKR1","AQP1"), lymphatic = c("PROX1","PDPN","CCL21","MMRN1"),
  fibroblast = c("PDGFRA","DCN","LUM","COL1A1","COL1A2"), mural = c("ACTA2","RGS5","MYH11","PDGFRB"), immune = c("PTPRC","CD14","CD68","LYZ","CD3E","MS4A1"), adipocyte = c("ADIPOQ","PLIN1","LIPE"))
CORE6 <- c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB")
rows <- list(); gene_rows <- list()
for (f in files) {
  k <- match_row(f); cat("Reading", basename(f), "...\n")
  dt <- fread(f, data.table = FALSE, check.names = FALSE)
  first <- dt[[1]]; looks_gene <- mean(grepl("^[A-Za-z][A-Za-z0-9._-]+$", head(first, 50))) > 0.8 && !all(grepl("^[ACGT]{8,}", head(first, 50)))
  mat <- as.matrix(dt[, -1, drop = FALSE]); rownames(mat) <- first; rm(dt)
  if (!looks_gene) { mat <- t(mat) }                        # cells x genes -> genes x cells
  m <- as(Matrix(mat, sparse = TRUE), "CsparseMatrix"); rm(mat); gene <- rownames(m); umi <- Matrix::colSums(m); keep <- umi >= MIN_UMI; m <- m[, keep, drop = FALSE]; umi <- umi[keep]
  if (ncol(m) < 2) { s0 <- if (!is.na(k)) sample_key[k] else sample_key[0]; rows[[length(rows) + 1]] <- data.table(file = basename(f), gsm = s0$gsm, sort = s0$sort, phenotype = s0$phenotype, donor = s0$donor, genes_in_file = length(gene), cells = ncol(m), endothelial_marker_rule = 0L); next }   # no cell passes the UMI floor
  sc <- sapply(PANELS, function(p) { j <- match(p, gene); j <- j[!is.na(j)]; if (!length(j)) rep(0, ncol(m)) else Matrix::colSums(m[j, , drop = FALSE]) / pmax(umi, 1) * 1e4 })
  jc <- match(CORE6, gene); jc <- jc[!is.na(jc)]; nd <- Matrix::colSums(m[jc, , drop = FALSE] > 0)
  oth <- apply(sc[, colnames(sc) != "endothelial", drop = FALSE], 1, max); top <- colnames(sc)[max.col(sc, ties.method = "first")]
  endo <- nd >= ENDO_MIN_MARKERS & top == "endothelial" & sc[, "endothelial"] >= ENDO_OVER_OTHER * oth
  gj <- match(GENE_LIST, gene)
  det_all <- if (any(!is.na(gj))) Matrix::rowSums(m[gj[!is.na(gj)], , drop = FALSE] > 0) else numeric(0)
  det_endo <- if (any(!is.na(gj)) && sum(endo)) Matrix::rowSums(m[gj[!is.na(gj)], endo, drop = FALSE] > 0) else rep(0, sum(!is.na(gj)))
  s <- if (!is.na(k)) sample_key[k] else data.table(gsm = NA, sort = NA, phenotype = NA, donor = NA, file_tag = NA)
  rows[[length(rows) + 1]] <- data.table(file = basename(f), gsm = s$gsm, sort = s$sort, phenotype = s$phenotype, donor = s$donor, genes_in_file = length(gene), cells = ncol(m), median_umi = median(umi), median_genes = median(Matrix::colSums(m > 0)),
                                         endothelial_marker_rule = sum(endo), cells_ge2_core_markers = sum(nd >= 2))
  if (length(det_all)) gene_rows[[length(gene_rows) + 1]] <- data.table(file = basename(f), phenotype = s$phenotype, sort = s$sort, gene = names(GENE_LIST)[!is.na(gj)], cells_detected_all = as.integer(det_all), cells_detected_in_endothelial = as.integer(det_endo), endothelial_cells = sum(endo), cells = ncol(m))
  rm(m); gc(verbose = FALSE)
}
per_sample <- rbindlist(rows, fill = TRUE); print(per_sample)
genes_tab <- rbindlist(gene_rows, fill = TRUE)

# ---- 3. summaries and the verdict
n_cells <- sum(per_sample$cells); n_endo <- sum(per_sample$endothelial_marker_rule)
by_group <- per_sample[, .(samples = .N, donors = uniqueN(donor), cells = sum(cells), endothelial_marker_rule = sum(endothelial_marker_rule)), by = .(phenotype, sort)][order(sort, phenotype)]; print(by_group)
gene_summary <- if (nrow(genes_tab)) genes_tab[, .(cells_detected_all = sum(cells_detected_all), total_cells = sum(cells), cells_detected_in_endothelial = sum(cells_detected_in_endothelial), endothelial_cells = sum(endothelial_cells)), by = .(phenotype, gene)][, `:=`(pct_of_endothelial = round(100 * cells_detected_in_endothelial / pmax(endothelial_cells, 1), 1))][order(gene, phenotype)] else data.table()
print(gene_summary)
verdict <- data.table(requirement = c("Control group (lean or obese without diabetes)", "Type 2 diabetes group", "Endothelial cells present", "More than 10,000 endothelial cells", "Human adipose tissue", "Both subcutaneous and visceral", "Gene list detectable"),
  result = c(sprintf("MET: %d lean and %d obese donors (BMI < 25 vs > 30); donors had no diabetes (paper)", uniqueN(per_sample[phenotype == "lean", donor]), uniqueN(per_sample[phenotype == "obese", donor])),
             "NOT MET: no type 2 diabetes group; the paper's 'healthy' donors had no history of diabetes",
             ifelse(n_endo > 0, sprintf("MET (marker rule): %d cells", n_endo), "NOT FOUND by the marker rule"),
             sprintf("%s: %d by the marker rule (threshold %d); total cells %d", ifelse(n_endo > MIN_ENDO, "MET", "NOT MET"), n_endo, MIN_ENDO, n_cells),
             "MET (human)", "NOT MET: deep subcutaneous abdominal fat only (paper); no visceral samples",
             if (nrow(gene_summary)) sprintf("%d of 8 genes found in the files", uniqueN(gene_summary$gene)) else "genes not found"))
print(verdict)
write_xlsx(list(Verdict = verdict, Per_sample = per_sample, By_group = by_group, Gene_list_detection = gene_summary, Gene_by_sample = genes_tab, GEO_sample_records = sample_key,
  Notes = data.table(note = c("Depot (subcutaneous only) and diabetes status (none) come from the paper, not from the files.", "Endothelial cells here use a strict marker rule on shallow data, so counts are lower bounds; CD45-negative samples hold the endothelial cells.",
   "The paper reports 110,250 cells (82,577 in one figure legend); compare with the cell total above."))), file.path(OUT_DIR, "GSE155960_criteria_check.xlsx")); fwrite(verdict, file.path(OUT_DIR, "verdict.csv")); cat("Saved to", OUT_DIR, "\n")
