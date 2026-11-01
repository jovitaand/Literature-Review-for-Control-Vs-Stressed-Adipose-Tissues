# Differential expression of the 8 genes of interest in ENDOTHELIAL nuclei: type 2 diabetes vs control, GSE268904.
# For Google Colab (R runtime). Tested on synthetic data only; NOT yet run on the real archive.
#
# Genes (names from the handwritten list -> official gene symbols used in the data):
#   ERG -> ERG | CD31 -> PECAM1 | VE-cadherin -> CDH5 | eNOS -> NOS3 | VCAM1 | ICAM1 | HIF-1a -> HIF1A | ANGPT2
# Method in one line: pick endothelial nuclei -> add up each person's counts (pseudobulk) -> compare people, not cells.

# ---- SETTINGS -----------------------------------------------------------
TAR_PATH <- ""                       # leave "" to auto-find in Drive / download by ID
WORK_DIR <- "/content/gse268904_raw"
OUT_DIR  <- "/content/gse268904_genes"
MIN_UMI <- 500; ENDO_MIN_MARKERS <- 3; ENDO_OVER_OTHER <- 2     # same endothelial rule as before
MIN_ENDO_NUCLEI <- 50                # people with fewer endothelial nuclei are flagged and excluded from the tests
GENES <- c(ERG = "ERG", `CD31` = "PECAM1", `VE-cadherin` = "CDH5", eNOS = "NOS3",
           VCAM1 = "VCAM1", ICAM1 = "ICAM1", `HIF-1a` = "HIF1A", ANGPT2 = "ANGPT2")
# -------------------------------------------------------------------------

for (p in c("Matrix", "data.table", "writexl")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl)})
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. get the archive
if (TAR_PATH == "") {
  hits <- list.files("/content", pattern = "GSE268904.*\\.tar$", recursive = TRUE, full.names = TRUE)
  if (length(hits)) TAR_PATH <- hits[1]
}
if (TAR_PATH == "" || !file.exists(TAR_PATH)) {
  system("pip install -q gdown"); system("gdown 1mKaM2c4473gV3zCoCV3F1ZknlCIDlAwq -O /content/GSE268904_RAW.tar")
  TAR_PATH <- "/content/GSE268904_RAW.tar"
}
stopifnot(file.exists(TAR_PATH)); untar(TAR_PATH, exdir = WORK_DIR)
files <- list.files(WORK_DIR, recursive = TRUE, full.names = TRUE)

# ---- 2. who is who (from the GEO sample records)
meta <- data.table(
  gsm = sprintf("GSM%d", 8302405:8302418),
  donor = c("9S","10S","13S","17S","SAT22","SAT54","SAT1","SAT18","SAT34","SAT39","SAT4","SAT51","SAT58","SAT60"),
  group = c("obese","obese diabetic","lean","obese diabetic","obese diabetic","lean","obese","obese","obese","obese","obese","obese","lean","lean"),
  bmi = c(39.5,61.5,23,52.6,39.7,22.3,41.7,61.3,62.5,61.2,49.7,57,24.3,23.5))
meta[, t2d := group == "obese diabetic"]

PANELS <- list(endothelial = c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB","CDH5","KDR"), lymphatic = c("PROX1","LYVE1","FLT4"),
  adipocyte = c("ADIPOQ","PLIN1","LIPE","LEP"), fibroblast = c("PDGFRA","DCN","LUM","COL1A1"),
  immune = c("PTPRC","CD163","MRC1","CD14"), mural = c("ACTA2","RGS5","MYH11","PDGFRB"))
CORE6 <- c("PECAM1","VWF","CLDN5","EGFL7","FLT1","PTPRB")

# ---- 3. per person: find endothelial nuclei, then read the 8 genes in them
donor_rows <- list(); cell_rows <- list(); gene_check <- list()
for (i in seq_len(nrow(meta))) {
  g <- meta$gsm[i]; pick <- function(t) files[grepl(g, files) & grepl(t, files)]
  fb <- pick("barcodes"); ff <- pick("features|genes"); fm <- pick("matrix")
  if (length(fb) != 1 || length(ff) != 1 || length(fm) != 1) { cat("Skipping", g, "(files missing)\n"); next }
  m <- as(readMM(fm), "CsparseMatrix"); feat <- fread(ff, header = FALSE); gene <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
  umi <- Matrix::colSums(m); keep <- umi >= MIN_UMI
  score <- sapply(PANELS, function(p) { j <- match(p, gene); j <- j[!is.na(j)]; if (!length(j)) rep(0, ncol(m)) else Matrix::colSums(m[j, , drop = FALSE]) / pmax(umi, 1) * 1e4 })
  jc <- match(CORE6, gene); jc <- jc[!is.na(jc)]
  ndet <- Matrix::colSums(m[jc, , drop = FALSE] > 0)
  other <- apply(score[, colnames(score) != "endothelial", drop = FALSE], 1, max)
  top <- colnames(score)[max.col(score, ties.method = "first")]
  endo <- keep & ndet >= ENDO_MIN_MARKERS & top == "endothelial" & score[, "endothelial"] >= ENDO_OVER_OTHER * other
  gj <- match(GENES, gene)
  gene_check[[i]] <- data.table(donor = meta$donor[i], gene_symbol = GENES, found_in_gene_list = !is.na(gj))
  e <- which(endo); n_e <- length(e)
  cnt <- if (n_e) as.matrix(m[gj[!is.na(gj)], e, drop = FALSE]) else matrix(0, sum(!is.na(gj)), 0)
  rownames(cnt) <- names(GENES)[!is.na(gj)]
  tot <- sum(umi[e])
  for (nm in names(GENES)) {
    x <- if (nm %in% rownames(cnt)) cnt[nm, ] else rep(NA_real_, n_e)
    donor_rows[[length(donor_rows) + 1]] <- data.table(
      donor = meta$donor[i], group = meta$group[i], t2d = meta$t2d[i], bmi = meta$bmi[i], gene = nm, symbol = GENES[[nm]],
      n_endo_nuclei = n_e, total_counts = sum(x), pseudobulk_cpm = if (n_e) sum(x) / tot * 1e6 else NA_real_,
      pct_nuclei_expressing = if (n_e) 100 * mean(x > 0) else NA_real_,
      mean_log1p_cp10k = if (n_e) mean(log1p(x / umi[e] * 1e4)) else NA_real_)
  }
  cat(sprintf("%-6s %-15s endothelial nuclei: %5d of %6d\n", meta$donor[i], meta$group[i], n_e, sum(keep)))
  rm(m); gc(verbose = FALSE)
}
dl <- rbindlist(donor_rows); gc_tab <- rbindlist(gene_check)
cat("\nGenes missing from the gene list:", if (any(!gc_tab$found_in_gene_list)) paste(unique(gc_tab[found_in_gene_list == FALSE, gene_symbol]), collapse = ", ") else "none", "\n")
dl[, enough_cells := n_endo_nuclei >= MIN_ENDO_NUCLEI]

# ---- 4. compare people (donor level): Mann-Whitney on pseudobulk CPM, plus fold change
compare <- function(a, b, label) {
  rbindlist(lapply(names(GENES), function(nm) {
    d <- dl[gene == nm & enough_cells & !is.na(pseudobulk_cpm)]
    x <- d[eval(a), pseudobulk_cpm]; y <- d[eval(b), pseudobulk_cpm]
    if (length(x) < 2 || length(y) < 2) return(data.table(contrast = label, gene = nm, n_T2D = length(x), n_comparison = length(y)))
    data.table(contrast = label, gene = nm, n_T2D = length(x), n_comparison = length(y),
               mean_cpm_T2D = mean(x), mean_cpm_comparison = mean(y),
               log2_fold_change = log2((mean(x) + 1) / (mean(y) + 1)),
               p_value = suppressWarnings(wilcox.test(x, y, exact = TRUE)$p.value))
  }), fill = TRUE)
}
res <- rbindlist(list(
  compare(quote(group == "obese diabetic"), quote(group == "obese"), "T2D vs obese (main)"),
  compare(quote(group == "obese diabetic"), quote(group == "lean"), "T2D vs lean"),
  compare(quote(group == "obese diabetic"), quote(group != "obese diabetic"), "T2D vs all non-diabetic")), fill = TRUE)
res[, p_adj_BH := p.adjust(p_value, "BH"), by = contrast]
res[, direction := ifelse(is.na(log2_fold_change), NA, ifelse(log2_fold_change > 0, "higher in T2D", "lower in T2D"))]
res[, read_as := ifelse(is.na(p_value), "too few people", ifelse(p_adj_BH < 0.05, "supported (still only 3 T2D people)", "no clear difference"))]
num <- c("mean_cpm_T2D", "mean_cpm_comparison", "log2_fold_change", "p_value", "p_adj_BH"); res[, (num) := lapply(.SD, signif, 3), .SDcols = num]
print(res[contrast == "T2D vs obese (main)"])

# ---- 5. picture: one dot per person
png(file.path(OUT_DIR, "gene_panel_by_group.png"), width = 1400, height = 700, res = 110)
par(mfrow = c(2, 4), mar = c(6, 4, 3, 1)); cols <- c(lean = "#4c78a8", obese = "#f2a541", `obese diabetic` = "#d1495b")
for (nm in names(GENES)) {
  d <- dl[gene == nm & enough_cells]; d[, group := factor(group, c("lean", "obese", "obese diabetic"))]
  stripchart(pseudobulk_cpm ~ group, data = d, vertical = TRUE, method = "jitter", pch = 19, col = cols[levels(d$group)], las = 2, main = nm, ylab = "CPM (endothelial nuclei)", cex = 1.4)
  for (k in 1:3) { v <- d[as.integer(group) == k, pseudobulk_cpm]; if (length(v)) segments(k - .3, mean(v), k + .3, mean(v), lwd = 2) }
}
dev.off()

# ---- 6. save
fwrite(dl, file.path(OUT_DIR, "gene_panel_per_donor.csv")); fwrite(res, file.path(OUT_DIR, "gene_panel_DE_results.csv"))
write_xlsx(list(Results = res, Per_donor = dl, Gene_check = gc_tab, Genes = data.table(handwritten_name = names(GENES), data_symbol = unname(GENES)),
  Notes = data.table(note = c("Endothelial nuclei were picked with a marker rule (not a formal annotation).",
   "Unit of analysis = person (pseudobulk). Cells from one person are not independent.",
   "Only 3 people have type 2 diabetes; the smallest possible two-sided p-value for 3 vs 7 is 0.0167 and 3 vs 4 is 0.057. Treat results as exploratory.",
   "T2D people differ in BMI (39.7 to 61.5) and from lean people (22 to 24); obesity and diabetes are partly confounded.",
   "CPM = gene counts per million endothelial counts. log2 fold change uses mean CPM + 1.",
   "With 8 genes tested and 3 vs 7 people, the smallest possible adjusted p-value is 0.13: nothing can reach significance, so use effect sizes and the raw pattern per person.",
   "Means of low-expressed genes (e.g. NOS3, ANGPT2) are noisy in single-nucleus data."))),
  file.path(OUT_DIR, "gene_panel_DE_summary.xlsx"))
cat("Saved to", OUT_DIR, "\n")
