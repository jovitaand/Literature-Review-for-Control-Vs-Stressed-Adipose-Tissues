# GSE268904_RAW.tar -> one consolidated CSV (nuclei counts, QC, crude endothelial estimate)
# Single-file version of notebooks/GSE268904_colab_qc.ipynb. Paste into one Colab cell (R runtime) or run with Rscript.
# Edit TAR_PATH below if the tar is not found automatically. Tested on a synthetic tar, not yet on the real data.

# Edit TAR_PATH after you find the file in the next cell
TAR_PATH <- ""                                   # e.g. "/content/drive/MyDrive/GSE268904_RAW.tar"
WORK_DIR <- "/content/gse268904_raw"             # extracted files (treated as read-only)
OUT_DIR  <- "/content/gse268904_qc"              # small outputs to download or copy to Drive
MIN_UMI  <- 500                                  # barcodes below this total UMI are not counted as nuclei
ENDO_MARKERS <- c("PECAM1", "VWF", "CLDN5", "EGFL7", "FLT1", "PTPRB")  # crude endothelial panel
RUN_SEURAT <- FALSE                              # TRUE = also run the optional clustering section (slow)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)

hits <- list.files("/content/drive", pattern = "GSE268904", recursive = TRUE, full.names = TRUE)
print(hits)
# Pick the .tar path and set TAR_PATH in the previous cell, or do it here:
if (TAR_PATH == "") { tars <- hits[grepl("\\.tar$", hits)]; if (length(tars) >= 1) TAR_PATH <- tars[1] }
cat("Using:", TAR_PATH, "\n")


if (TAR_PATH == "" || !file.exists(TAR_PATH)) {
  system("pip install -q gdown")
  system("gdown 1mKaM2c4473gV3zCoCV3F1ZknlCIDlAwq -O /content/GSE268904_RAW.tar")
  TAR_PATH <- "/content/GSE268904_RAW.tar"
}
stopifnot(file.exists(TAR_PATH)); cat("Using:", TAR_PATH, "\n")


need <- c("Matrix", "data.table", "digest")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(digest)})

untar(TAR_PATH, exdir = WORK_DIR)
files <- list.files(WORK_DIR, recursive = TRUE, full.names = TRUE)
files <- files[!file.info(files)$isdir]
inv <- data.table(
  relative_path = sub(paste0("^", WORK_DIR, "/?"), "", files),
  filename = basename(files),
  size_bytes = file.info(files)$size,
  sha256 = vapply(files, function(f) digest(f, algo = "sha256", file = TRUE), character(1)),
  file_type = ifelse(grepl("barcodes", files), "barcodes",
              ifelse(grepl("features|genes", files), "features",
              ifelse(grepl("matrix", files), "matrix", "other")))
)
inv[, gsm := sub("^(GSM[0-9]+).*", "\\1", filename)]
inv[, duplicate_of_n := .N, by = sha256]
cat("Files extracted:", nrow(inv), " (expected 42 = 14 samples x 3 files)\n")
print(inv[, .N, by = file_type])
print(inv[duplicate_of_n > 1, .(files_with_identical_content = .N), by = file_type])
cat("Identical features.tsv.gz files across samples are expected (same gene list); identical barcodes or matrix files would be suspicious.\n")


meta <- data.table(
  gsm = sprintf("GSM%d", 8302405:8302418),
  sample = c("9S","10S","13S","17S","SAT22","SAT54","SAT1","SAT18","SAT34","SAT39","SAT4","SAT51","SAT58","SAT60"),
  condition = c("obese","obese diabetic","lean","obese diabetic","obese diabetic","lean","obese","obese","obese","obese","obese","obese","lean","lean"),
  sex = c("female","female","female","female","male","female","male","female","male","male","male","female","male","male"),
  age = c(37,41,42,34,50,61,43,29,40,41,30,44,54,26),
  bmi = c(39.5,61.5,23,52.6,39.7,22.3,41.7,61.3,62.5,61.2,49.7,57,24.3,23.5))
# Referential check: every metadata GSM must have 3 files
chk <- inv[, .(n_files = .N, has_barcodes = any(file_type == "barcodes"), has_features = any(file_type == "features"), has_matrix = any(file_type == "matrix")), by = gsm]
print(merge(meta[, .(gsm, sample)], chk, by = "gsm", all = TRUE))

mito_prefix <- "^MT-"
res <- list()
for (i in seq_len(nrow(meta))) {
  gsm <- meta$gsm[i]; smp <- meta$sample[i]
  pick <- function(type) files[grepl(gsm, files) & grepl(type, files)]
  fb <- pick("barcodes"); ff <- pick("features|genes"); fm <- pick("matrix")
  if (length(fb) != 1 || length(ff) != 1 || length(fm) != 1) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "missing files"); next }
  m <- tryCatch(readMM(fm), error = function(e) NULL)
  if (is.null(m)) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "matrix unreadable"); next }
  bc <- fread(fb, header = FALSE)[[1]]; feat <- fread(ff, header = FALSE); gene <- feat[[2]]
  if (ncol(m) != length(bc) || nrow(m) != length(gene)) { res[[i]] <- data.table(gsm = gsm, sample = smp, status = "dimension mismatch"); next }
  m <- as(m, "CsparseMatrix")
  umi <- Matrix::colSums(m); ngene <- Matrix::colSums(m > 0)
  mt <- Matrix::colSums(m[grepl(mito_prefix, gene), , drop = FALSE]) / pmax(umi, 1)
  present <- ENDO_MARKERS[ENDO_MARKERS %in% gene]
  det <- if (length(present)) Matrix::colSums(m[match(present, gene), , drop = FALSE] > 0) else rep(0, ncol(m))
  keep <- umi >= MIN_UMI
  res[[i]] <- data.table(gsm = gsm, sample = smp, condition = meta$condition[i], status = "ok",
    n_barcodes_file = length(bc), n_genes_file = length(gene), n_nuclei_umi_ge_min = sum(keep),
    median_umi = median(umi[keep]), median_genes = median(ngene[keep]), median_mito_frac = round(median(mt[keep]), 4),
    raw_counts_integer = all(m@x == round(m@x)),
    markers_present = paste(present, collapse = "|"),
    endo_ge2 = sum(keep & det >= 2), endo_ge3 = sum(keep & det >= 3), endo_ge4 = sum(keep & det >= 4))
  cat(smp, ":", length(bc), "barcodes;", sum(keep), "with UMI >=", MIN_UMI, "\n")
  rm(m); gc(verbose = FALSE)
}
per_sample <- rbindlist(res, fill = TRUE)
print(per_sample)

ok <- per_sample[status == "ok"]
group <- ok[, .(n_samples = .N, nuclei = sum(n_nuclei_umi_ge_min), endo_ge2 = sum(endo_ge2), endo_ge3 = sum(endo_ge3), endo_ge4 = sum(endo_ge4)), by = condition]
overall <- data.table(condition = "ALL", n_samples = nrow(ok), nuclei = sum(ok$n_nuclei_umi_ge_min), endo_ge2 = sum(ok$endo_ge2), endo_ge3 = sum(ok$endo_ge3), endo_ge4 = sum(ok$endo_ge4))
group <- rbind(group, overall)
print(group)
cat("\nCRUDE estimate: nuclei where >= k of the endothelial markers are detected. Not a cell-type annotation.\n")

if (RUN_SEURAT) {
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat", repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(Seurat))
  seu_res <- list()
  for (i in seq_len(nrow(meta))) {
    gsm <- meta$gsm[i]; smp <- meta$sample[i]
    pick <- function(type) files[grepl(gsm, files) & grepl(type, files)]
    m <- readMM(pick("matrix")); bc <- fread(pick("barcodes"), header = FALSE)[[1]]; gene <- fread(pick("features|genes"), header = FALSE)[[2]]
    rownames(m) <- make.unique(gene); colnames(m) <- bc
    obj <- CreateSeuratObject(m, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset = nCount_RNA >= MIN_UMI & percent.mt < 10)
    obj <- NormalizeData(obj, verbose = FALSE); obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE); obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE); obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
    genes_ok <- intersect(ENDO_MARKERS, rownames(obj))
    obj <- AddModuleScore(obj, features = list(genes_ok), name = "endo", verbose = FALSE)
    key <- intersect(c("PECAM1", "VWF", "CLDN5"), rownames(obj))
    pct <- sapply(levels(obj$seurat_clusters), function(cl) mean(colSums(GetAssayData(obj, layer = "counts")[key, obj$seurat_clusters == cl, drop = FALSE] > 0) >= 1))
    sc  <- tapply(obj$endo1, obj$seurat_clusters, median)
    endo_clusters <- names(sc)[sc > 0.3 & pct[names(sc)] > 0.3]
    n_endo <- sum(obj$seurat_clusters %in% endo_clusters)
    seu_res[[i]] <- data.table(sample = smp, condition = meta$condition[i], nuclei_after_qc = ncol(obj), n_clusters = length(levels(obj$seurat_clusters)),
                               endo_clusters = paste(endo_clusters, collapse = "|"), endothelial_nuclei = n_endo)
    cat(smp, ": ", ncol(obj), "nuclei; endothelial clusters:", paste(endo_clusters, collapse = ","), "->", n_endo, "\n")
    rm(obj, m); gc(verbose = FALSE)
  }
  seu <- rbindlist(seu_res); print(seu)
  print(seu[, .(nuclei = sum(nuclei_after_qc), endothelial = sum(endothelial_nuclei)), by = condition])
} else cat("Skipped (RUN_SEURAT = FALSE)\n")

# ---- wide file info per sample (from the inventory) ----
fi <- dcast(inv, gsm ~ file_type, value.var = c("filename", "size_bytes", "sha256"))
# ---- per-sample table ----
final <- merge(meta, fi, by = "gsm", all.x = TRUE)
final <- merge(final, per_sample[, !c("sample", "condition"), with = FALSE], by = "gsm", all.x = TRUE)
if (exists("seu")) {
  s2 <- copy(seu)[, .(sample, seurat_nuclei_after_qc = nuclei_after_qc, seurat_n_clusters = n_clusters, seurat_endothelial_clusters = endo_clusters, seurat_endothelial_nuclei = endothelial_nuclei)]
  final <- merge(final, s2, by = "sample", all.x = TRUE)
}
final[, row_type := "sample"]
# ---- group and overall rows ----
num <- c("n_barcodes_file", "n_nuclei_umi_ge_min", "endo_ge2", "endo_ge3", "endo_ge4", "seurat_nuclei_after_qc", "seurat_endothelial_nuclei")
num <- intersect(num, names(final))
grp <- final[, c(list(n_samples = .N), lapply(.SD, sum, na.rm = TRUE)), by = condition, .SDcols = num][, row_type := "group"]
ovr <- final[, c(list(n_samples = .N), lapply(.SD, sum, na.rm = TRUE)), .SDcols = num][, `:=`(row_type = "overall", condition = "ALL")]
ovr[, notes := paste("Nuclei = barcodes with total UMI >=", MIN_UMI, "; endo_geK = nuclei with >= K of 6 endothelial genes detected (CRUDE estimate, not an annotation); ambient RNA can inflate endo_ge2;",
                     "identical features files across samples are expected; no data were removed; processed", format(Sys.time(), "%Y-%m-%d"))]
final <- rbindlist(list(final, grp, ovr), fill = TRUE)
lead <- c("row_type", "gsm", "sample", "condition", "sex", "age", "bmi", "n_samples", "status")
setcolorder(final, c(intersect(lead, names(final)), setdiff(names(final), lead)))
out_csv <- file.path(OUT_DIR, "GSE268904_consolidated.csv")
fwrite(final, out_csv)
cat("Wrote", out_csv, "|", nrow(final), "rows x", ncol(final), "columns |", round(file.info(out_csv)$size / 1e3, 1), "KB\n")
print(final[row_type != "sample", c("row_type", "condition", "n_samples", "n_nuclei_umi_ge_min", "endo_ge3"), with = FALSE])
# Copy to Drive
dest <- "/content/drive/MyDrive/gse268904_qc"
dir.create(dest, showWarnings = FALSE); file.copy(out_csv, dest, overwrite = TRUE)
cat("Copied to", dest, "- upload GSE268904_consolidated.csv back to the chat\n")
