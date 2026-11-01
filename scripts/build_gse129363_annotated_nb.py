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


# ======================= PART 2: differential expression in endothelial cells =======================
md(r"""
# PART 2. Differential expression in endothelial cells: obese (control) vs obese with type 2 diabetes (disease)

Both groups are **obese**. The *control* group is obese **without** type 2 diabetes (T2D); the *disease* group is obese **with** T2D. There is no lean or healthy group in this dataset.

**Read this first**
- The files contain **no author cell-type labels**, so endothelial cells are found with the marker rule from Step 5. Step 11 verifies how many there are (it does not assume 50, 105 or any number) and shows how the count changes with stricter rules.
- The unit of replication is the **person / sample**, not the cell: counts are added up per sample (pseudobulk) and compared between samples. Cells from one person are not independent.
- **Two diabetes labels exist** (the cell annotation file and the GEO sample records disagree for 3 people: E07, 04, 09). The analysis is run with both and both are reported.
- Subcutaneous (SAT) and visceral (VAT) fat are analysed **separately**, and also combined (with depot as a covariate) as an exploratory check.
- If a group has too few samples or cells, the notebook says so and marks results exploratory.
""")

md(r"""
## Step 11. How many endothelial cells, in which group, from which person?

**What this does:** builds the endothelial set (cells kept by the authors that pass the marker rule), links every cell to its person, depot and diabetes group under **both** labelling sources, and reports counts per group and per person. It also shows the count under stricter rules so you can see how fragile the number is.
""")
code(r"""
# ---- SETTINGS for Part 2 ----
if (!requireNamespace("edgeR", quietly = TRUE)) { install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install("edgeR", ask = FALSE, update = FALSE) }
suppressPackageStartupMessages(library(edgeR))
GENE_LIST <- c(ERG = "ERG", CD31 = "PECAM1", `VE-cadherin` = "CDH5", eNOS = "NOS3", VCAM1 = "VCAM1", ICAM1 = "ICAM1", `HIF-1a` = "HIF1A", ANGPT2 = "ANGPT2")
MIN_CELLS_PER_SAMPLE <- 10   # a sample needs at least this many endothelial cells to count as a replicate
MIN_SAMPLES_PER_GROUP <- 2   # fewer samples than this in a group = not testable
RELIABLE_SAMPLES_PER_GROUP <- 4   # below this, results are labelled EXPLORATORY

# who is who: sample name -> person, and the diabetes label from the GEO sample records
sample_key <- data.table(
  sample = c("SVF2-SAT", "SVF3-VAT", "SVF5-SAT", "Ate-04-VAT", "Ate-04-SAT", "Ate-05-VAT", "Ate-05-SAT", "Ate-06-VAT", "Ate-06-SAT", "Ate-08-VAT", "Ate-08-SAT", "Ate-09-VAT", "Ate-09-SAT", "Ate-10_VAT", "Ate-10-SAT", "Ate-11-VAT", "Ate-11-SAT", "Ate-12-VAT", "Ate-12-SAT", "Ate-13-VAT", "Ate-13-SAT", "Ate-14-VAT", "Ate-14-SAT", "Ate-15-VAT", "Ate-15-SAT"),
  person = c("E02", "E03", "E07", "04", "04", "05", "05", "06", "06", "08", "08", "09", "09", "10", "10", "11", "11", "12", "12", "13", "13", "14", "14", "15", "15"),
  geo_label = c("diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "diabetic", "diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "diabetic", "diabetic", "non-diabetic", "non-diabetic", "non-diabetic", "non-diabetic"))

# the matrix column of every row of `cells`
cells[, mcol := match(barcode, bc)]
A <- cells[in_annotation == TRUE]
A <- merge(A, sample_key, by.x = sample_col, by.y = "sample", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(A$person))                                   # every sample must be recognised
A[, depot := get(depot_col)]
A[, annotation_label := ifelse(grepl("^non", get(diab_col), ignore.case = TRUE), "control (obese, no T2D)", "disease (obese + T2D)")]
A[, geo_label_grp := ifelse(geo_label == "diabetic", "disease (obese + T2D)", "control (obese, no T2D)")]
E <- A[endo_like == TRUE]
cat("Cells kept by the authors:", nrow(A), "| endothelial by marker rule:", nrow(E), "
")

# 1) how fragile is the endothelial count? (more markers required = stricter)
rule_tab <- data.table(rule = paste0(">= ", 2:6, " of 6 endothelial markers detected"), cells = sapply(2:6, function(k) sum(A$endo_markers_detected >= k)))
rule_tab <- rbind(rule_tab, data.table(rule = "specific rule used here (>= 3 markers + endothelial top family + 2x others)", cells = nrow(E)))
print(rule_tab)

# 2) per group and per person
for (src in c("annotation_label", "geo_label_grp")) {
  cat("
=== Endothelial cells by group and depot, labels from:", ifelse(src == "annotation_label", "CELL ANNOTATION FILE", "GEO SAMPLE RECORDS"), "===
")
  print(E[, .(endothelial_cells = .N, samples = uniqueN(get(sample_col)), people = uniqueN(person)), by = c(src, "depot")][order(get(src), depot)])
}
per_person <- E[, .(endothelial_cells = .N), by = .(person, sample = get(sample_col), depot, annotation_label, geo_label_grp)][order(person, depot)]
per_person[, enough_for_replicate := endothelial_cells >= MIN_CELLS_PER_SAMPLE]
all_samples <- A[, .(cells_kept = .N), by = .(person, sample = get(sample_col), depot)]
per_person <- merge(all_samples, per_person, by = c("person", "sample", "depot"), all.x = TRUE); per_person[is.na(endothelial_cells), endothelial_cells := 0L]
per_person[, enough_for_replicate := endothelial_cells >= MIN_CELLS_PER_SAMPLE]
print(per_person)
cat("
Samples with at least", MIN_CELLS_PER_SAMPLE, "endothelial cells:", sum(per_person$enough_for_replicate), "of", nrow(per_person), "
")
""")

md(r"""
## Step 12. Is each gene detected in the endothelial cells? (per group)

**What this does:** for each gene in your list, counts the endothelial cells in which it is detected (at least one RNA molecule) and the percentage, separately for control and disease, for each depot and both depots together.
**Genes detected in both groups** and **in only one group** are listed separately (Step 14). A gene detected in only one group is still tested in Step 13.
""")
code(r"""
Gsym <- GENE_LIST[GENE_LIST %in% gene]
cat("Genes from the list found in the gene file:", paste(Gsym, collapse = ", "), "| missing:", paste(setdiff(GENE_LIST, gene), collapse = ", "), "
")
ge <- m[match(Gsym, gene), E$mcol, drop = FALSE]                  # genes x endothelial cells (raw counts)
det_one <- function(src) {
  rbindlist(lapply(c("SAT", "VAT", "SAT+VAT"), function(dp) {
    sel <- if (dp == "SAT+VAT") rep(TRUE, nrow(E)) else E$depot == dp
    rbindlist(lapply(unique(E[[src]]), function(gr) {
      k <- which(sel & E[[src]] == gr)
      data.table(label_source = ifelse(src == "annotation_label", "annotation file", "GEO records"), depot = dp, group = gr, endothelial_cells = length(k),
                 gene = names(Gsym), symbol = unname(Gsym),
                 cells_detected = if (length(k)) as.integer(Matrix::rowSums(ge[, k, drop = FALSE] > 0)) else 0L,
                 mean_counts_per_cell = if (length(k)) round(Matrix::rowMeans(ge[, k, drop = FALSE]), 3) else 0)
    }))
  }))
}
detect <- rbind(det_one("annotation_label"), det_one("geo_label_grp"))
detect[, pct_detected := ifelse(endothelial_cells > 0, round(100 * cells_detected / endothelial_cells, 1), NA_real_)]
print(detect[label_source == "annotation file" & depot == "SAT+VAT", .(group, endothelial_cells, gene, cells_detected, pct_detected)][order(gene, group)])
""")

md(r"""
## Step 13. Differential expression: disease vs control, one value per person (pseudobulk, edgeR)

**What this does:** for each depot separately (SAT, VAT) and then both together (depot as a covariate):
1. keeps samples with enough endothelial cells and adds up each gene's counts per sample;
2. tests **disease vs control** with edgeR's quasi-likelihood test (works with very few replicates; genes detected in only one group are still tested);
3. reports log2 fold change (positive = higher in T2D), p-value and BH-adjusted p-value (across all genes tested, and across your gene list alone);
4. lists the **50 most differentially expressed genes** and marks which are in your gene list.
Anything with fewer than the reliable number of samples per group is labelled **EXPLORATORY**; fewer than the minimum is **NOT TESTABLE**.
""")
code(r"""
run_de <- function(src, depot_sel) {
  label <- paste(ifelse(src == "annotation_label", "annotation", "GEO"), depot_sel, sep = " | ")
  d <- E[get(sample_col) %in% per_person[enough_for_replicate == TRUE, sample]]
  if (depot_sel != "SAT+VAT") d <- d[depot == depot_sel]
  smp <- d[, .(n_cells = .N, group = get(src)[1], depot = depot[1], person = person[1]), by = .(sample = get(sample_col))]
  n_c <- sum(smp$group == "control (obese, no T2D)"); n_d <- sum(smp$group == "disease (obese + T2D)")
  status <- data.table(analysis = label, label_source = src, depot = depot_sel, control_samples = n_c, disease_samples = n_d, endothelial_cells_used = nrow(d))
  if (min(n_c, n_d) < MIN_SAMPLES_PER_GROUP) { status[, status := sprintf("NOT TESTABLE: %d control and %d disease samples with >= %d endothelial cells", n_c, n_d, MIN_CELLS_PER_SAMPLE)]; return(list(status = status, res = NULL)) }
  cols <- match(d$mcol, seq_len(ncol(m)))
  ind <- sparseMatrix(i = seq_len(nrow(d)), j = match(d[[sample_col]], smp$sample), x = 1, dims = c(nrow(d), nrow(smp)))
  pb <- as.matrix(m[, d$mcol, drop = FALSE] %*% ind); colnames(pb) <- smp$sample; rownames(pb) <- gene
  pb <- pb[!duplicated(rownames(pb)), , drop = FALSE]
  grp <- factor(ifelse(smp$group == "disease (obese + T2D)", "disease", "control"), c("control", "disease"))
  y <- DGEList(pb, group = grp); y <- y[rowSums(y$counts) >= 3 | rownames(y) %in% Gsym, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  des <- if (depot_sel == "SAT+VAT" && length(unique(smp$depot)) > 1) model.matrix(~ factor(smp$depot) + grp) else model.matrix(~ grp)
  coefn <- ncol(des)
  y <- estimateDisp(y, des); fit <- glmQLFit(y, des); tt <- topTags(glmQLFTest(fit, coef = coefn), n = Inf)$table
  res <- data.table(analysis = label, label_source = src, depot = depot_sel, gene = rownames(tt), logFC_disease_vs_control = round(tt$logFC, 3), logCPM = round(tt$logCPM, 2),
                    p_value = tt$PValue, p_adj_BH_all_genes = tt$FDR, control_samples = n_c, disease_samples = n_d)
  res[, in_gene_list := gene %in% GENE_LIST]
  res[in_gene_list == TRUE, p_adj_BH_gene_list_only := p.adjust(p_value, "BH")]
  status[, status := ifelse(min(n_c, n_d) < RELIABLE_SAMPLES_PER_GROUP, sprintf("EXPLORATORY: only %d control vs %d disease samples", n_c, n_d), "tested")]
  list(status = status, res = res)
}
runs <- list()
for (src in c("annotation_label", "geo_label_grp")) for (dp in c("SAT", "VAT", "SAT+VAT")) runs[[length(runs) + 1]] <- tryCatch(run_de(src, dp), error = function(e) list(status = data.table(analysis = paste(src, dp), status = paste("ERROR:", conditionMessage(e))), res = NULL))
status_tab <- rbindlist(lapply(runs, `[[`, "status"), fill = TRUE); print(status_tab)
de_all <- rbindlist(lapply(runs, `[[`, "res"), fill = TRUE)
if (nrow(de_all)) {
  gl_res <- de_all[in_gene_list == TRUE][order(analysis, p_value)]; print(gl_res[, .(analysis, gene, logFC_disease_vs_control, p_value, p_adj_BH_gene_list_only)])
  top50 <- de_all[, head(.SD[order(p_value)], 50), by = analysis]; top50[, rank := seq_len(.N), by = analysis]
  overlap50 <- top50[in_gene_list == TRUE, .(genes_of_your_list_in_top50 = paste(gene, collapse = ", "), n = .N), by = analysis]; print(overlap50)
} else { gl_res <- top50 <- overlap50 <- data.table(); cat("No analysis was testable with the current thresholds.
") }
""")

md(r"""
## Step 14. Genes detected in both groups, in one group only, or in neither

**What this does:** from Step 12, classifies each gene of your list for every depot and labelling source. This is only about **detection** (at least one RNA molecule in at least one endothelial cell of the group); every gene is still tested in Step 13.
""")
code(r"""
cls <- dcast(detect[, .(label_source, depot, gene, group, cells_detected)], label_source + depot + gene ~ group, value.var = "cells_detected", fill = 0L)
setnames(cls, c("control (obese, no T2D)", "disease (obese + T2D)"), c("cells_detected_control", "cells_detected_disease"), skip_absent = TRUE)
cls[, detected_in := ifelse(cells_detected_control > 0 & cells_detected_disease > 0, "both groups", ifelse(cells_detected_control > 0, "control only", ifelse(cells_detected_disease > 0, "disease only", "neither")))]
print(cls[depot == "SAT+VAT"][order(label_source, gene)])
both_list <- cls[detected_in == "both groups"][order(label_source, depot, gene)]
""")

md(r"""
## Step 15. Supporting cell-level check for your 8 genes (EXPLORATORY ONLY)

**What this does:** when there are too few samples for the person-level test, this compares the 8 genes between all endothelial cells of the two groups (rank test). **It treats every cell as independent, so p-values are far too small.** Use it to see the direction of change only, never as evidence.
""")
code(r"""
lib <- Matrix::colSums(m[, E$mcol, drop = FALSE])
lx <- log1p(sweep(as.matrix(ge), 2, lib / 1e4, "/"))               # log1p counts per 10,000
cell_level <- rbindlist(lapply(c("annotation_label", "geo_label_grp"), function(src) rbindlist(lapply(c("SAT", "VAT", "SAT+VAT"), function(dp) {
  sel <- if (dp == "SAT+VAT") rep(TRUE, nrow(E)) else E$depot == dp
  a <- sel & E[[src]] == "disease (obese + T2D)"; b <- sel & E[[src]] == "control (obese, no T2D)"
  rbindlist(lapply(seq_along(Gsym), function(i) data.table(label_source = src, depot = dp, gene = names(Gsym)[i], cells_disease = sum(a), cells_control = sum(b),
    mean_log_disease = if (sum(a)) round(mean(lx[i, a]), 3) else NA_real_, mean_log_control = if (sum(b)) round(mean(lx[i, b]), 3) else NA_real_,
    p_value_cell_level = if (sum(a) > 1 && sum(b) > 1) suppressWarnings(wilcox.test(lx[i, a], lx[i, b])$p.value) else NA_real_)))
}))))
cell_level[, direction := ifelse(is.na(mean_log_disease) | is.na(mean_log_control), NA, ifelse(mean_log_disease > mean_log_control, "higher in T2D", "lower in T2D"))]
print(cell_level[label_source == "annotation_label" & depot == "SAT+VAT"])
""")

md(r"""
## Step 16. Save the Part 2 results
""")
code(r"""
OUT2 <- file.path(OUT_DIR, "part2_DE"); dir.create(OUT2, showWarnings = FALSE)
notes <- data.table(note = c(
  "Endothelial cells = cells kept by the authors that pass the marker rule (no author cell-type labels exist in the files).",
  "Control = obese without T2D; disease = obese with T2D. No lean group exists.",
  "Diabetes labels differ between the annotation file and GEO for people E07, 04 and 09; both versions are analysed.",
  "Replicate = sample (one person per depot). Pseudobulk edgeR quasi-likelihood. logFC positive = higher in T2D.",
  "Results with fewer than the reliable number of samples per group are EXPLORATORY; below the minimum they are NOT TESTABLE.",
  "SAT and VAT come from the same people in most cases; the combined analysis ignores this pairing (exploratory).",
  "Cell-level p-values (Step 15) treat cells as independent and must not be used as evidence."))
sheets <- list(Status = status_tab, Endothelial_rule_counts = rule_tab, Endothelial_per_sample = per_person, Gene_detection = detect, Detected_both_groups = both_list,
               Detection_class_all = cls, Gene_list_DE = gl_res, Top50_by_analysis = top50, Overlap_top50_gene_list = overlap50, Cell_level_exploratory = cell_level, Notes = notes)
sheets <- Filter(function(x) !is.null(x) && nrow(x) > 0, sheets)
write_xlsx(sheets, file.path(OUT2, "GSE129363_endothelial_DE_summary.xlsx"))
if (nrow(de_all)) fwrite(de_all, file.path(OUT2, "DE_all_genes_all_analyses.csv"))
fwrite(detect, file.path(OUT2, "gene_detection_by_group.csv")); fwrite(per_person, file.path(OUT2, "endothelial_cells_per_sample.csv"))
print(list.files(OUT2))
""")

nb={"cells":C,"metadata":{"kernelspec":{"display_name":"R","language":"R","name":"ir"},"language_info":{"name":"R"},"colab":{"provenance":[]}},"nbformat":4,"nbformat_minor":5}
json.dump(nb,open("notebooks/GSE129363_colab_annotated.ipynb","w"),indent=1)
print(len(C),"cells")
