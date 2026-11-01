# GSE129363 (Vijay et al.): preliminary analysis of the Discovery Cohort files, for Google Colab with an R runtime.
# Reads: Discovery_Cohort barcodes / genes / matrix / CellAnnotation, and (optionally) SVF_Normalised_Data.
# Writes small readable files to OUT_DIR. Tested on synthetic data only; NOT yet run on the real files.
#
# Question answered: does this dataset meet the four criteria?
#   human adipose tissue | endothelial cells | control vs type 2 diabetes | more than 10,000 cells

# ---- CONFIG -------------------------------------------------------------
DATA_DIR <- "/content/drive/MyDrive"      # folder holding the five files (Drive mounted) - searched recursively
OUT_DIR  <- "/content/gse129363_out"
READ_NORMALISED <- TRUE                   # the 81 MB SVF_Normalised_Data table: only its header + gene column are read
CELL_THRESHOLD <- 10000
MIN_UMI <- 200                            # the authors' own floor (GEO: cells with <200 UMI removed)
ENDO_MIN_MARKERS <- 3; ENDO_OVER_OTHER <- 2
# -------------------------------------------------------------------------

for (p in c("Matrix", "data.table", "writexl")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl)})
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. locate the files (names are matched loosely, so "Discovery Cohort barcodes (1).tsv" also works)
all_files <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
find1 <- function(pattern) {
  hit <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]
  if (!length(hit)) stop("Not found (pattern ", pattern, ") under ", DATA_DIR, ". Check DATA_DIR or copy the file there.")
  if (length(hit) > 1) { cat("Several matches for", pattern, "- using the newest:\n"); print(basename(hit)); hit <- hit[which.max(file.info(hit)$mtime)] }
  hit
}
f_bc  <- find1("barcodes")
f_gn  <- find1("genes|features")
f_mx  <- find1("matrix.*mtx|\\.mtx")
f_an  <- find1("annotation")
f_nrm <- if (READ_NORMALISED) tryCatch(find1("normali[sz]ed"), error = function(e) NA_character_) else NA_character_
cat("Files used:\n"); print(basename(c(f_bc, f_gn, f_mx, f_an, f_nrm)))

# ---- 2. read the 10x-style triplet
bc <- fread(f_bc, header = FALSE)[[1]]
gn <- fread(f_gn, header = FALSE)
gene <- if (ncol(gn) >= 2) gn[[2]] else gn[[1]]
m <- as(readMM(f_mx), "CsparseMatrix")                     # genes x cells
dims_ok <- nrow(m) == length(gene) && ncol(m) == length(bc)
cat(sprintf("Matrix %d genes x %d cells | barcodes file %d | genes file %d | dimensions consistent: %s\n",
            nrow(m), ncol(m), length(bc), length(gene), dims_ok))
if (!dims_ok) stop("Matrix dimensions do not match barcodes/genes. If genes x cells is transposed, use t(m).")

umi <- Matrix::colSums(m); ngene <- Matrix::colSums(m > 0)
mt <- Matrix::colSums(m[grepl("^MT-", gene), , drop = FALSE]) / pmax(umi, 1)

# ---- 3. marker panels (counts per 10,000 UMI), same definitions as the GSE268904 pipeline
PANELS <- list(
  endothelial = c("PECAM1", "VWF", "CLDN5", "EGFL7", "FLT1", "PTPRB"),
  lymphatic   = c("PROX1", "LYVE1", "FLT4"),
  adipocyte   = c("ADIPOQ", "PLIN1", "LIPE", "LEP"),
  fibroblast  = c("PDGFRA", "DCN", "LUM", "COL1A1"),
  immune      = c("PTPRC", "CD163", "MRC1", "CD14"),
  mural       = c("ACTA2", "RGS5", "MYH11", "PDGFRB"))
panel_score <- sapply(PANELS, function(g) {
  idx <- match(g, gene); idx <- idx[!is.na(idx)]
  if (!length(idx)) return(rep(0, ncol(m)))
  Matrix::colSums(m[idx, , drop = FALSE]) / pmax(umi, 1) * 1e4
})
idx_e <- match(PANELS$endothelial, gene); idx_e <- idx_e[!is.na(idx_e)]
endo_det <- if (length(idx_e)) Matrix::colSums(m[idx_e, , drop = FALSE] > 0) else rep(0, ncol(m))
best_other <- apply(panel_score[, setdiff(colnames(panel_score), "endothelial"), drop = FALSE], 1, max)
top_panel <- colnames(panel_score)[max.col(panel_score, ties.method = "first")]
top_panel[rowSums(panel_score) == 0] <- "none"
endo_specific <- endo_det >= ENDO_MIN_MARKERS & top_panel == "endothelial" & panel_score[, "endothelial"] >= ENDO_OVER_OTHER * best_other

sc <- as.data.table(round(panel_score, 2)); setnames(sc, paste0("score_", names(PANELS)))
cells <- cbind(data.table(barcode = bc, umi = umi, genes_detected = ngene, mito_frac = round(mt, 4),
                          barcode_suffix = sub("^.*-", "", bc), endo_markers_detected = endo_det,
                          top_panel = top_panel, endo_specific = endo_specific), sc)

# ---- 4. the authors' cell annotation: find the barcode column, the cell-type column and any sample/disease column
an <- fread(f_an, data.table = TRUE)
cat("CellAnnotation:", nrow(an), "rows; columns:", paste(names(an), collapse = " | "), "\n")
bcol <- names(an)[which.max(sapply(an, function(x) mean(as.character(x) %in% bc)))]
if (mean(as.character(an[[bcol]]) %in% bc) < 0.5) {              # barcodes may be the row names (first, unnamed column)
  bcol <- names(an)[1]; cat("Barcode column assumed to be the first column:", bcol, "\n")
}
cat(sprintf("Barcode column: %s | %.1f%% of annotation rows match the matrix barcodes\n", bcol, 100 * mean(as.character(an[[bcol]]) %in% bc)))
setnames(an, bcol, "barcode"); an[, barcode := as.character(barcode)]
other_cols <- setdiff(names(an), "barcode")
cat_cols <- other_cols[sapply(an[, ..other_cols], function(x) !is.numeric(x) && uniqueN(x) > 1 && uniqueN(x) < 60)]
cat("Candidate label columns (text, 2-59 distinct values):", paste(cat_cols, collapse = ", "), "\n")
guess <- function(pat) { h <- cat_cols[grepl(pat, cat_cols, ignore.case = TRUE)]; if (length(h)) h[1] else NA_character_ }
celltype_col <- guess("cell.?type|singler|label|annot|ident|cluster"); if (is.na(celltype_col) && length(cat_cols)) celltype_col <- cat_cols[1]
sample_col   <- guess("sample|donor|subject|patient|library|orig")
diab_col     <- guess("diab|status|condition|disease|group")
depot_col    <- guess("depot|tissue|fat|vat|sat")
cat("Chosen columns -> cell type:", celltype_col, "| sample:", sample_col, "| diabetes:", diab_col, "| depot:", depot_col, "\n")
cells <- merge(cells, an, by = "barcode", all.x = TRUE, sort = FALSE)
cells[, in_annotation := barcode %in% an$barcode]

# ---- 5. summaries
n_cells <- nrow(cells); n_pass <- sum(cells$umi >= MIN_UMI)
ct <- if (!is.na(celltype_col)) cells[, .(cells = .N, share_pct = round(100 * .N / n_cells, 1)), by = c(celltype_col)][order(-cells)] else data.table()
if (nrow(ct)) setnames(ct, 1, "cell_type_author_label")
endo_lab <- if (!is.na(celltype_col)) grepl("endothel", cells[[celltype_col]], ignore.case = TRUE) else rep(FALSE, n_cells)
by_group <- function(col) if (!is.na(col)) cells[, .(cells = .N, endothelial_by_author_label = sum(endo_lab[.I]), endothelial_by_markers_specific = sum(endo_specific)), by = c(col)][order(-cells)] else data.table()
by_sample <- by_group(sample_col); by_diab <- by_group(diab_col); by_depot <- by_group(depot_col)
by_suffix <- cells[, .(cells = .N, endothelial_by_markers_specific = sum(endo_specific), median_umi = median(umi)), by = barcode_suffix][order(as.integer(barcode_suffix))]

# cross-check the normalised table (header + gene names only)
nrm_info <- data.table(item = "normalised table", value = "not read")
if (!is.na(f_nrm)) {
  hdr <- names(fread(f_nrm, nrows = 0)); g1 <- fread(f_nrm, select = 1L)[[1]]
  nrm_info <- data.table(item = c("file", "columns in header (incl. gene-name column)", "rows (genes)", "header cells found in barcodes file"),
                         value = c(basename(f_nrm), length(hdr), length(g1), sum(hdr %in% bc)))
}

# ---- 6. the four criteria
sym <- function(x) ifelse(x, "MET", "NOT MET")
n_endo_lab <- sum(endo_lab); n_endo_mk <- sum(endo_specific)
crit <- data.table(
  criterion = c("Human adipose tissue", "Endothelial cells present", "Control vs type 2 diabetes", "More than 10,000 cells"),
  result = c("MET (GEO: Homo sapiens; subcutaneous and visceral adipose SVF)",
             if (n_endo_lab > 0 || n_endo_mk > 0) sprintf("MET: %d cells labelled endothelial by the authors; %d by the marker rule", n_endo_lab, n_endo_mk) else "NOT FOUND in this file",
             "PARTLY: all donors are obese; the contrast is obese diabetic vs obese NON-diabetic, with no lean control (from the GEO sample records)",
             sprintf("%s: %d cells in the matrix (%d with UMI >= %d); threshold %d", sym(n_cells > CELL_THRESHOLD), n_cells, n_pass, MIN_UMI, CELL_THRESHOLD)),
  source = c("GEO series/sample records", "this run", "GEO sample records (GSE129363_samples.csv)", "this run"))

# ---- 7. write
fwrite(cells, file.path(OUT_DIR, "GSE129363_cells_one_row_per_barcode.csv"))
fwrite(crit, file.path(OUT_DIR, "criteria_check.csv"))
sheets <- list(Criteria = crit, Cell_types = ct, By_sample = by_sample, By_diabetes = by_diab, By_depot = by_depot,
               By_barcode_suffix = by_suffix, Normalised_table = nrm_info,
               Overview = data.table(item = c("Cells in matrix", "Genes in matrix", "Cells with UMI >= MIN_UMI", "Median UMI per cell", "Median genes per cell",
                                              "Median mito fraction", "Endothelial (author label)", "Endothelial (marker rule)", "Cells matched to annotation"),
                                     value = c(n_cells, nrow(m), n_pass, median(umi), median(ngene), round(median(mt), 4), n_endo_lab, n_endo_mk, sum(cells$in_annotation))),
               Notes = data.table(note = c("Cell-type labels are the authors' (SingleR); the marker rule is an independent crude check.",
                 "Barcode suffix (-1, -2, ...) is the sample index in the Cell Ranger aggregate; mapping it to GSM numbers needs the aggregate CSV, which is not in the five files.",
                 "If the annotation has a sample/diabetes column it is used above; otherwise use By_barcode_suffix and ask the authors or read the paper for the order.",
                 "Reads are raw Cell Ranger counts; the normalised table is only cross-checked for size.")))
sheets <- sheets[sapply(sheets, nrow) > 0]
write_xlsx(sheets, file.path(OUT_DIR, "GSE129363_summary_for_readers.xlsx"))
print(crit); print(head(ct, 15)); print(by_suffix)
cat("Done. Small files to share are in", OUT_DIR, "(skip the per-barcode CSV if large).\n")
