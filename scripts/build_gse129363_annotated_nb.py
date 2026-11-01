import json
C=[]
def md(t): C.append({"cell_type":"markdown","metadata":{},"source":t.strip().splitlines(True)})
def code(t): C.append({"cell_type":"code","metadata":{},"execution_count":None,"outputs":[],"source":t.strip().splitlines(True)})

md("""
# GSE129363: is this dataset suitable? (annotated, plain-language notebook)

**Question:** does this dataset meet our four criteria?
1. human **adipose (fat) tissue**
2. **endothelial cells** (cells lining blood vessels)
3. **control vs type 2 diabetes (T2D)**
4. **more than 10,000 cells**

**How to use:** Runtime > Change runtime type > **R**, then run each cell top to bottom (Shift+Enter).
Each cell has a short "What this does" note above it, and comments (lines starting with `#`) inside.

> Status: tested on synthetic data only, not yet on the real files. Check the numbers it prints against the paper (Vijay et al. 2020, PubMed 32066997).

**Words used here**
- **Cell** = one cell (one barcode = one cell's unique tag).
- **Gene** = a row in the data; the number tells how active that gene is in a cell.
- **UMI** = one counted RNA molecule. More UMIs = better-quality cell.
- **SVF** = stromal vascular fraction: the non-fat-storing cells of fat tissue (vessel cells, fibroblasts, immune cells).
""")

md("""
## Step 1. Get the five files into Colab

**What this does:** downloads the files from Google Drive into Colab's temporary storage.
The five Drive file IDs are below (from your links). Downloading works only if each file is shared as **"Anyone with the link"**.
If it fails, use Option B: click the **folder icon** on the left, then the **Drive icon** to mount your Drive, and set `DATA_DIR` to the folder with the files.
""")
code("""
# ---- SETTINGS (the only things you may need to change) ----
DATA_DIR <- "/content/gse129363_data"    # where the files will be (or already are)
OUT_DIR  <- "/content/gse129363_out"     # where results are saved
USE_GDOWN <- TRUE                        # TRUE = download by Drive ID (Option A); FALSE = files already in DATA_DIR / mounted Drive (Option B)

CELL_THRESHOLD <- 10000   # criterion 4: more than this many cells
MIN_UMI <- 200            # the authors removed cells with fewer than 200 UMIs

# Drive file IDs of the five GSE129363 files (barcodes, CellAnnotation, genes, matrix, SVF_Normalised_Data)
drive_ids <- c(
  "1xVFGimI87DMyyUHoJxbBBr_PXS_tU3bt",
  "1l77NQm3jVmJ8ieoD9KRpWjGFRaMev5hk",
  "1imnxlaJKfZMK5WUJas-5IXjQbuHNi79y",
  "190PY9cnmseA9_xXRD6vBNBeCxXY4tBZS",
  "1M2Q1NAnSGLhV3HemzSuHUvIy8RF7mTZf")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)

if (USE_GDOWN) {
  system("pip install -q gdown")                       # small download helper
  for (id in drive_ids) {
    # gdown keeps the original file name from Drive
    system(sprintf("cd %s && gdown --id %s", DATA_DIR, id))
  }
}
print(list.files(DATA_DIR, recursive = TRUE))           # you should see 5 files
""")

md("""
## Step 2. Load the R tools and find each file by name

**What this does:** installs three helper packages (reading big tables, sparse matrices, writing Excel), then identifies which file is which by its name:
- **barcodes** = list of cell tags (one per cell)
- **genes** = list of gene names (one per row of the matrix)
- **matrix** = the big table of counts (genes x cells), stored compactly
- **CellAnnotation** = the authors' label for each cell (e.g. "endothelial cell")
- **SVF_Normalised_Data** = the same data, normalised (only used for a size check)
""")
code("""
for (p in c("Matrix", "data.table", "writexl"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(Matrix); library(data.table); library(writexl)})

all_files <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)

# helper: find the file whose NAME matches a pattern
find1 <- function(pattern, required = TRUE) {
  hit <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]
  if (!length(hit)) { if (required) stop("Cannot find a file matching '", pattern, "' in ", DATA_DIR) else return(NA_character_) }
  hit[1]
}
f_bc  <- find1("barcodes")
f_gn  <- find1("genes|features")
f_mx  <- find1("matrix.*mtx|\\\\.mtx")
f_an  <- find1("annotation")
f_nrm <- find1("normali[sz]ed", required = FALSE)

files_used <- data.frame(role = c("barcodes", "genes", "matrix", "annotation", "normalised"),
                         file = basename(c(f_bc, f_gn, f_mx, f_an, f_nrm)),
                         size_MB = round(file.size(c(f_bc, f_gn, f_mx, f_an, f_nrm)) / 1e6, 2))
print(files_used)
""")

md("""
## Step 3. Load the cells, genes and counts, and check they agree

**What this does:** reads the three core files and checks that
*number of rows in the matrix = number of genes* and *number of columns = number of cells*.
**The number of columns is your answer to criterion 4** (before any filtering).
""")
code("""
bc <- fread(f_bc, header = FALSE)[[1]]                 # cell tags
gn <- fread(f_gn, header = FALSE)                      # gene list
gene <- if (ncol(gn) >= 2) gn[[2]] else gn[[1]]        # use gene SYMBOLS (2nd column) when present
m <- as(readMM(f_mx), "CsparseMatrix")                 # counts: rows = genes, columns = cells

cat("Matrix:", nrow(m), "genes x", ncol(m), "cells\\n")
cat("Barcodes file:", length(bc), "| genes file:", length(gene), "\\n")
stopifnot(nrow(m) == length(gene), ncol(m) == length(bc))   # stops with an error if the files do not agree
cat("OK: files are consistent.\\n")
""")

md("""
## Step 4. Quality of each cell

**What this does:** for every cell, counts
- `umi` = total RNA molecules counted (cell size / quality)
- `genes_detected` = how many different genes were seen
- `mito_frac` = share of counts from mitochondrial genes (high = stressed or dying cell)
""")
code("""
umi   <- Matrix::colSums(m)
ngene <- Matrix::colSums(m > 0)
mito  <- Matrix::colSums(m[grepl("^MT-", gene), , drop = FALSE]) / pmax(umi, 1)

cat("Median UMI per cell:", median(umi), "| median genes per cell:", median(ngene),
    "| median mito share:", round(median(mito), 3), "\\n")
cat("Cells with at least", MIN_UMI, "UMI:", sum(umi >= MIN_UMI), "of", length(umi), "\\n")
""")

md("""
## Step 5. Our own simple endothelial-cell check (marker genes)

**What this does:** endothelial cells switch on a known set of genes (PECAM1, VWF, CLDN5, EGFL7, FLT1, PTPRB).
We call a cell **endothelial-like** if
1. at least **3 of the 6** endothelial genes are detected, **and**
2. its strongest gene group is the endothelial one, **and**
3. the endothelial score is **at least 2x** the best non-endothelial score (fibroblast, immune, fat, vessel-wall, lymphatic).

This is a rough independent check on the authors' labels, **not** a proper annotation. Scores are counts per 10,000 UMI.
""")
code("""
PANELS <- list(
  endothelial = c("PECAM1", "VWF", "CLDN5", "EGFL7", "FLT1", "PTPRB"),
  lymphatic   = c("PROX1", "LYVE1", "FLT4"),
  adipocyte   = c("ADIPOQ", "PLIN1", "LIPE", "LEP"),
  fibroblast  = c("PDGFRA", "DCN", "LUM", "COL1A1"),
  immune      = c("PTPRC", "CD163", "MRC1", "CD14"),
  mural       = c("ACTA2", "RGS5", "MYH11", "PDGFRB"))   # mural = pericytes / vessel wall cells

# score = marker counts per 10,000 UMI, for each group and each cell
score <- sapply(PANELS, function(g) {
  i <- match(g, gene); i <- i[!is.na(i)]
  if (!length(i)) return(rep(0, ncol(m)))
  Matrix::colSums(m[i, , drop = FALSE]) / pmax(umi, 1) * 1e4
})

# how many of the 6 endothelial genes are seen in each cell
ie <- match(PANELS$endothelial, gene); ie <- ie[!is.na(ie)]
cat("Endothelial markers found in this gene list:", paste(gene[ie], collapse = ", "), "\\n")
endo_n <- Matrix::colSums(m[ie, , drop = FALSE] > 0)

best_other <- apply(score[, colnames(score) != "endothelial", drop = FALSE], 1, max)
top_group  <- colnames(score)[max.col(score, ties.method = "first")]
endo_like  <- endo_n >= 3 & top_group == "endothelial" & score[, "endothelial"] >= 2 * best_other

cat("Endothelial-like cells (marker rule):", sum(endo_like), "\\n")
""")

md("""
## Step 6. The authors' own labels

**What this does:** reads the CellAnnotation file (the authors labelled every cell with a type using SingleR), matches it to our cells by barcode, and guesses which column holds the cell type, and which (if any) hold sample or diabetes status. The printed list of columns lets you check the guess.
""")
code("""
an <- fread(f_an)
cat("Annotation file:", nrow(an), "rows. Columns:", paste(names(an), collapse = " | "), "\\n")

# find the column that holds barcodes (the one matching our barcodes best); else assume the first
share <- sapply(an, function(x) mean(as.character(x) %in% bc))
bcol <- if (max(share) >= 0.5) names(an)[which.max(share)] else names(an)[1]
cat("Barcode column:", bcol, "| matches matrix:", round(100 * max(share), 1), "%\\n")
setnames(an, bcol, "barcode"); an[, barcode := as.character(barcode)]

# columns that look like labels = text with 2-59 distinct values
others <- setdiff(names(an), "barcode")
labcols <- others[sapply(an[, ..others], function(x) !is.numeric(x) && uniqueN(x) > 1 && uniqueN(x) < 60)]
pick <- function(pat) { h <- labcols[grepl(pat, labcols, ignore.case = TRUE)]; if (length(h)) h[1] else NA_character_ }
celltype_col <- pick("cell.?type|singler|label|annot|ident|cluster")   # the real file has NO cell-type column, so this stays NA
sample_col <- pick("sample|donor|subject|patient|library|orig")
diab_col   <- pick("diab|status|condition|disease|group")
depot_col  <- pick("depot|tissue|fat|vat|sat")
cat("Cell type column:", celltype_col, "| sample:", sample_col, "| diabetes:", diab_col, "| depot:", depot_col, "\\n")

# one table, one row per cell: our measures + the authors' labels
cells <- data.table(barcode = bc, barcode_suffix = sub("^.*-", "", bc), umi = umi, genes_detected = ngene,
                    mito_frac = round(mito, 4), endo_markers_detected = endo_n, endo_like = endo_like)
cells <- merge(cells, an, by = "barcode", all.x = TRUE, sort = FALSE)
cells[, in_annotation := barcode %in% an$barcode]
cat("Barcodes in the matrix:", nrow(cells), "| cells kept by the authors (in annotation):", sum(cells$in_annotation), "\\n")
""")

md("""
## Step 7. What cell types are there? (and how many are endothelial)

**What this does:** counts cells per author-label, and draws a bar chart. The endothelial bar is highlighted.
""")
code("""
if (is.na(celltype_col)) {
  cat("The annotation file has no cell-type column (only sample, condition, tissue), so endothelial cells come from our marker rule.\\n")
  is_endo_label <- rep(FALSE, nrow(cells)); ct <- data.table(cell_type_author_label = "no cell-type labels in file", cells = nrow(cells), share_percent = 100)
  print(cells[in_annotation == TRUE, .(cells = .N, endothelial_like = sum(endo_like), pct_endothelial_like = round(100 * mean(endo_like), 1))])
} else {
is_endo_label <- grepl("endothel", cells[[celltype_col]], ignore.case = TRUE)

ct <- cells[, .N, by = c(celltype_col)][order(-N)]
setnames(ct, c("cell_type_author_label", "cells"))
ct[, share_percent := round(100 * cells / sum(cells), 1)]
print(ct)

png(file.path(OUT_DIR, "cell_types.png"), width = 900, height = 500)
par(mar = c(5, 14, 3, 2))
barplot(rev(ct$cells), names.arg = rev(ct$cell_type_author_label), horiz = TRUE, las = 1,
        col = rev(ifelse(grepl("endothel", ct$cell_type_author_label, ignore.case = TRUE), "#d1495b", "#4c78a8")),
        xlab = "Number of cells", main = "GSE129363: cells per author-assigned type (endothelial in red)")
dev.off()
barplot(rev(ct$cells), names.arg = rev(ct$cell_type_author_label), horiz = TRUE, las = 1,
        col = rev(ifelse(grepl("endothel", ct$cell_type_author_label, ignore.case = TRUE), "#d1495b", "#4c78a8")),
        xlab = "Number of cells", main = "Cells per type (endothelial in red)")
}
""")

md("""
## Step 8. Breakdown by sample, diabetes status and fat depot

**What this does:** if the annotation file has sample / diabetes / depot columns, counts cells (and endothelial cells) in each group.
If not, it falls back to the **barcode suffix** (-1, -2, ...), which is the sample's position in the authors' combined run.
> Why this matters: diabetes status belongs to the **person**, so a fair comparison needs enough people per group, not just many cells.
""")
code("""
A <- cells[in_annotation == TRUE]            # analyse the cells the authors kept
A[, endo_lab := if (is.na(celltype_col)) FALSE else grepl("endothel", get(celltype_col), ignore.case = TRUE)]
by_group <- function(cols) {
  if (any(is.na(cols))) return(NULL)
  A[, .(cells = .N, endothelial_author_label = sum(endo_lab), endothelial_marker_rule = sum(endo_like),
        pct_endothelial_marker_rule = round(100 * mean(endo_like), 1)), by = cols][order(-cells)]
}
by_sample <- by_group(sample_col); by_diab <- by_group(diab_col); by_depot <- by_group(depot_col)
by_diab_depot <- by_group(c(diab_col, depot_col))
by_suffix <- cells[, .(barcodes = .N, cells_kept_by_authors = sum(in_annotation), endothelial_marker_rule = sum(endo_like & in_annotation)), by = barcode_suffix][order(as.integer(barcode_suffix))]
for (nm in c("by_diab_depot", "by_diab", "by_depot", "by_sample")) { cat("\\n---", nm, "---\\n"); print(get(nm)) }
cat("\\n--- by barcode suffix (= sample position) ---\\n"); print(by_suffix)
""")

md("""
## Step 9. Final verdict on the four criteria

**What this does:** fills the checklist using the numbers above. Rows 1 and 3 come from the GEO records (not this data):
- all 32 samples are human (*Homo sapiens*), from subcutaneous or visceral fat SVF;
- every donor is obese; the comparison is **obese diabetic vs obese non-diabetic**, with **no lean control**.
""")
code("""
n_raw <- nrow(cells); n_cells <- sum(cells$in_annotation); n_endo_lab <- sum(A$endo_lab); n_endo_mk <- sum(A$endo_like)

criteria <- data.table(
  criterion = c("1. Human adipose tissue", "2. Endothelial cells", "3. Control vs type 2 diabetes", "4. More than 10,000 cells"),
  verdict = c("MET",
              ifelse(n_endo_lab > 0 | n_endo_mk > 0, "MET", "NOT FOUND"),
              "PARTLY (obese diabetic vs obese non-diabetic; no lean control)",
              ifelse(n_cells > CELL_THRESHOLD, "MET", "NOT MET")),
  evidence = c("GEO: Homo sapiens; subcutaneous and visceral adipose SVF",
               if (is.na(celltype_col)) sprintf("no author cell-type labels in the file; %d cells by our marker rule", n_endo_mk) else sprintf("%d cells labelled endothelial by authors; %d by our marker rule", n_endo_lab, n_endo_mk),
               "All donors obese. Condition per cell is in the annotation file; it disagrees with the GEO sample records for 3 people (see Step 8 note)",
               sprintf("%d cells kept by the authors (%d barcodes in the matrix); threshold %d", n_cells, n_raw, CELL_THRESHOLD)))
print(criteria)
""")

md("""
## Step 10. Save results (download from the folder panel on the left)

**What this does:** saves a small Excel workbook for non-technical readers, the criteria CSV and the per-cell CSV.
Share the Excel and `criteria_check.csv` first (small); the per-cell CSV is large.
""")
code("""
sheets <- list(Criteria = criteria, Cell_types = ct, By_diabetes_and_depot = by_diab_depot, By_diabetes = by_diab, By_depot = by_depot,
               By_sample = by_sample, By_barcode_suffix = by_suffix,
               Overview = data.table(item = c("Cells in matrix", "Genes", "Cells with UMI >= MIN_UMI", "Median UMI", "Median genes per cell",
                                              "Endothelial (author label)", "Endothelial (marker rule)"),
                                     value = c(n_cells, nrow(m), sum(umi >= MIN_UMI), median(umi), median(ngene), n_endo_lab, n_endo_mk)))
sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)
write_xlsx(sheets, file.path(OUT_DIR, "GSE129363_summary_for_readers.xlsx"))
fwrite(criteria, file.path(OUT_DIR, "criteria_check.csv"))
fwrite(cells, file.path(OUT_DIR, "GSE129363_cells_one_row_per_barcode.csv"))
print(list.files(OUT_DIR))
""")

nb={"cells":C,"metadata":{"kernelspec":{"display_name":"R","language":"R","name":"ir"},"language_info":{"name":"R"},"colab":{"provenance":[]}},"nbformat":4,"nbformat_minor":5}
json.dump(nb,open("notebooks/GSE129363_colab_annotated.ipynb","w"),indent=1)
print(len(C),"cells")
