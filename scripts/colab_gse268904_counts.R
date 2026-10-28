# GSE268904_RAW.tar: nuclei counts, basic QC and a CRUDE endothelial-nuclei estimate, per sample and per group.
# For Google Colab with an R runtime. NOT YET RUN (written without R access): expect small fixes on the first run.
#
# The tar holds, per sample, GSM*_<id>_barcodes.tsv.gz, _features.tsv.gz and _matrix.mtx.gz (10x output, 14 samples).
# Outputs are small CSVs in OUT_DIR; upload only those.

# ---- CONFIG -------------------------------------------------------------
TAR_PATH <- "/content/drive/MyDrive/GSE268904_RAW.tar"   # path to the tar in your mounted Drive
WORK_DIR <- "/content/gse268904"
OUT_DIR  <- "/content/gse268904_qc"
MIN_UMI  <- 500        # crude floor so empty droplets do not count as nuclei; report with and without
# -------------------------------------------------------------------------

for (p in c("Matrix", "data.table")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table)})
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(TAR_PATH)) stop("Tar not found: ", TAR_PATH)
untar(TAR_PATH, exdir = WORK_DIR)
all_files <- list.files(WORK_DIR, recursive = TRUE, full.names = TRUE)
cat("Files extracted:", length(all_files), "\n")

# Sample metadata copied from the GEO records (GSE268904_samples.csv)
meta <- data.table(
  gsm = sprintf("GSM%d", 8302405:8302418),
  sample = c("9S","10S","13S","17S","SAT22","SAT54","SAT1","SAT18","SAT34","SAT39","SAT4","SAT51","SAT58","SAT60"),
  condition = c("obese","obese diabetic","lean","obese diabetic","obese diabetic","lean","obese","obese","obese","obese","obese","obese","lean","lean"))

endo_markers <- c("PECAM1", "VWF", "CLDN5", "EGFL7", "FLT1", "PTPRB")   # crude panel, not an annotation
mito_prefix  <- "^MT-"

res <- list()
for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]; smp <- meta$sample[i]
  fb <- all_files[grepl(paste0(gsm, ".*barcodes"), all_files)]
  ff <- all_files[grepl(paste0(gsm, ".*features"), all_files)]
  fm <- all_files[grepl(paste0(gsm, ".*matrix"), all_files)]
  if (length(fb) != 1 || length(ff) != 1 || length(fm) != 1) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "missing files"); next }
  m <- tryCatch(readMM(fm), error = function(e) NULL)
  if (is.null(m)) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "matrix unreadable"); next }
  bc <- fread(fb, header = FALSE)[[1]]; feat <- fread(ff, header = FALSE)
  gene <- feat[[2]]
  if (ncol(m) != length(bc) || nrow(m) != length(gene)) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "dimension mismatch"); next }
  m <- as(m, "CsparseMatrix")
  umi <- Matrix::colSums(m); ngene <- Matrix::colSums(m > 0)
  mt <- Matrix::colSums(m[grepl(mito_prefix, gene), , drop = FALSE]) / pmax(umi, 1)
  present <- endo_markers[endo_markers %in% gene]
  det <- if (length(present)) Matrix::colSums(m[match(present, gene), , drop = FALSE] > 0) else rep(0, ncol(m))
  keep <- umi >= MIN_UMI
  res[[i]] <- data.table(gsm = gsm, sample = smp, condition = meta$condition[i], status = "ok",
    n_barcodes_file = length(bc), n_genes_file = length(gene),
    n_umi_ge_min = sum(keep), median_umi = median(umi[keep]), median_genes = median(ngene[keep]),
    median_mito_frac = round(median(mt[keep]), 4),
    markers_present = paste(present, collapse = "|"),
    endo_ge2_markers = sum(keep & det >= 2), endo_ge3_markers = sum(keep & det >= 3), endo_ge4_markers = sum(keep & det >= 4))
  cat(smp, ": ", length(bc), "barcodes,", sum(keep), "with UMI >=", MIN_UMI, "\n")
  rm(m); gc(verbose = FALSE)
}
per_sample <- rbindlist(res, fill = TRUE)
fwrite(per_sample, file.path(OUT_DIR, "per_sample_counts.csv"))

ok <- per_sample[status == "ok"]
group <- ok[, .(n_samples = .N, nuclei_total = sum(n_umi_ge_min),
  endo_ge2 = sum(endo_ge2_markers), endo_ge3 = sum(endo_ge3_markers), endo_ge4 = sum(endo_ge4_markers)), by = condition]
fwrite(group, file.path(OUT_DIR, "group_summary.csv"))
fwrite(rbindlist(list(data.table(item = "all samples", nuclei_total = sum(ok$n_umi_ge_min), endo_ge3 = sum(ok$endo_ge3_markers))), fill = TRUE),
       file.path(OUT_DIR, "overall_summary.csv"))
writeLines(c("# Notes", "",
  "- Nuclei = barcodes with total UMI >= MIN_UMI. If the barcodes file is already CellRanger-filtered, n_barcodes_file is the author's count.",
  "- Endothelial counts are a CRUDE marker-detection estimate (>= k of six endothelial genes detected per nucleus), not a cell-type annotation.",
  "- Ambient RNA can inflate the >=2 estimate; use the >=3 and >=4 columns as the more conservative bounds.",
  "- Proper endothelial counts need clustering and annotation (for example Seurat) and would likely differ.",
  capture.output(sessionInfo())), file.path(OUT_DIR, "notes.md"))
print(per_sample); print(group)
cat("Upload the folder", OUT_DIR, "(small CSVs) back to the chat.\n")
