# GSE176171 (Emont et al. 2022, human + mouse white adipose atlas): criteria check for the file you shared (GSE176171.gz, 45 MB).
# For Google Colab (R runtime) or the pipeline. Tested on synthetic tables only; NOT yet run on the real file.
#
# Problem statement checked:
#   1. control group = obese without diabetes OR lean      2. type 2 diabetes (T2D) group
#   3. endothelial cells / nuclei, more than 10,000        4. human adipose tissue, BOTH subcutaneous and visceral
#   5. the endothelial gene list is detectable             (gene names: PECAM1, CDH5, NOS3, HIF1A, ERG, VCAM1, ICAM1, ANGPT2)

# ---- SETTINGS -----------------------------------------------------------
FILE <- ""                         # path to GSE176171.gz; "" = search /content and the working folder
OUT_DIR <- "/content/gse176171_check"
GENE_LIST <- c(PECAM1 = "PECAM1", CDH5 = "CDH5", NOS3 = "NOS3", HIF1A = "HIF1A", ERG = "ERG", VCAM1 = "VCAM1", ICAM1 = "ICAM1", ANGPT2 = "ANGPT2")
MIN_ENDO <- 10000
# -------------------------------------------------------------------------
for (p in c("data.table", "writexl", "R.utils")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(data.table); library(writexl)}); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (FILE == "") { h <- list.files(c("/content", "."), pattern = "GSE176171.*\\.gz$", recursive = TRUE, full.names = TRUE); if (length(h)) FILE <- h[1] }
stopifnot(file.exists(FILE)); cat("File:", FILE, round(file.size(FILE) / 1e6, 1), "MB\n")

# ---- 1. what kind of file is it?
con <- gzfile(FILE, "rt"); head3 <- readLines(con, n = 3); close(con)
cat("First lines:\n"); print(substr(head3, 1, 200))
kind <- if (grepl("^%%MatrixMarket", head3[1])) "matrix" else "table"
cat("Detected:", kind, "\n")
if (kind == "matrix") stop("This is a count matrix (genes x cells) without cell labels. Also provide GSE176171_cell_metadata.tsv.gz (and the barcodes/features files) so cells can be labelled by cell type, depot and person.")

# ---- 2. read the table and find the useful columns
dt <- fread(FILE, data.table = TRUE); cat("Rows:", nrow(dt), "| columns:", ncol(dt), "\n"); print(names(dt))
pick <- function(pat, exclude = NULL) { h <- names(dt)[grepl(pat, names(dt), ignore.case = TRUE)]; if (!is.null(exclude)) h <- h[!grepl(exclude, h, ignore.case = TRUE)]; if (length(h)) h[1] else NA_character_ }
col <- list(celltype = pick("cell.?type|celltype|annotation|cell_class|^type$|cluster_name|ident"), sample = pick("^sample$|sample_id|sample.?name|orig.ident|library|^id$|donor|patient|subject"),
            depot = pick("depot|region|tissue|fat|site|location"), species = pick("species|organism|genome"), bmi = pick("bmi"), diab = pick("diab|t2d|disease|condition|status"))
cat("Columns chosen:\n"); print(unlist(col))
profile <- rbindlist(lapply(names(dt), function(n) data.table(column = n, class = class(dt[[n]])[1], distinct = uniqueN(dt[[n]]), example = paste(head(unique(as.character(dt[[n]])), 4), collapse = " | "))))
if (!is.na(col$species)) { cat("Species values:\n"); print(table(dt[[col$species]])); hs <- grepl("human|sapiens|^hs|hg", dt[[col$species]], ignore.case = TRUE); if (any(hs)) dt <- dt[hs]; cat("Human rows kept:", nrow(dt), "\n") }

# ---- 3. embed the GEO sample records (BMI, sex, age, depot of the 24 human samples) and try to link them to the cells
sample_key <- data.table(gsm = c("GSM5359325", "GSM5359326", "GSM5359327", "GSM5359328", "GSM5359329", "GSM5359330", "GSM5359331", "GSM5359332", "GSM5359333", "GSM5359334", "GSM5359335", "GSM5359336", "GSM5359337", "GSM5820679", "GSM5820680", "GSM5820681", "GSM5820682", "GSM5820683", "GSM5820684", "GSM5820685", "GSM5820686", "GSM5820687", "GSM5820688", "GSM5820689"), title = c("Human visceral white adipose 01-1", "Human visceral white adipose 01-2", "Human visceral white adipose 253-1", "Human visceral white adipose 254-1", "Human visceral white adipose 255-1", "Human visceral white adipose 256-1", "Human subcutaneous white adipose 01-1", "Human subcutaneous white adipose 02-1", "Human subcutaneous white adipose 04-1", "Human subcutaneous white adipose 253-1", "Human subcutaneous white adipose 254-1", "Human subcutaneous white adipose 255-1", "Human subcutaneous white adipose 256-1", "Human visceral white adipose 09-1", "Human visceral white adipose 10-1", "Human visceral white adipose 12-1", "Human visceral white adipose 13-1", "Human visceral white adipose 266-1", "Human subcutaneous white adipose 09-1", "Human subcutaneous white adipose 10-1", "Human subcutaneous white adipose 11-1", "Human subcutaneous white adipose 12-1", "Human subcutaneous white adipose 13-1", "Human subcutaneous white adipose 266-1"), depot = c("visceral", "visceral", "visceral", "visceral", "visceral", "visceral", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "visceral", "visceral", "visceral", "visceral", "visceral", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous", "subcutaneous"),
  bmi = c(49.3, 49.3, 49.3, 33.1, 25.4, 30.04, 30.04, 23.96, 23.96, 24.27, 24.27, 34.53, 34.53, 45.7, 43.1, 48.7, 43.2, 22.15, 45.7, 43.1, 42.8, 48.7, 43.2, 22.15), sex = c("female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "female", "male", "male", "male", "female", "female", "female", "male", "male", "male"), age = c("29y", "29y", "29y", "57y", "51y", "53y", "53y", "41y", "41y", "73y", "73y", "41y", "41y", "41y", "35y", "36y", "24y", "68y", "41y", "35y", "58y", "36y", "24y", "68y"), surgery = c("VSG", "VSG", "VSG", "Hernia", "CCY", "TAH BSO", "TAH BSO", "TAH BSO", "TAH BSO", "TAH BSO", "TAH BSO", "CCY", "CCY", "VSG", "RYGB", "VSG", "VSG", "Colon polyp", "VSG", "RYGB", "VSG", "VSG", "VSG", "Colon polyp"))
sample_key[, bmi_group := cut(bmi, c(0, 25, 30, 100), c("lean (BMI < 25)", "overweight (25-30)", "obese (BMI >= 30)"), right = FALSE)]
sample_key[, person := paste(bmi, age, sex)]          # samples from the same person share BMI, age and sex
dt[, cell_type := if (!is.na(col$celltype)) as.character(get(col$celltype)) else NA_character_]
dt[, sample_id := if (!is.na(col$sample)) as.character(get(col$sample)) else NA_character_]
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
link <- function(s) { k <- match(norm(s), norm(sample_key$title)); if (anyNA(k)) k[is.na(k)] <- match(norm(s[is.na(k)]), norm(sample_key$gsm)); k }
us <- unique(dt$sample_id); lk <- link(us); cat("Cell-table samples:", length(us), "| linked to GEO sample records:", sum(!is.na(lk)), "\n")
dt[, geo_row := lk[match(sample_id, us)]]
dt[, depot_std := if (!is.na(col$depot)) ifelse(grepl("visc|omen", get(col$depot), ignore.case = TRUE), "visceral", ifelse(grepl("subq|subcut|sc", get(col$depot), ignore.case = TRUE), "subcutaneous", as.character(get(col$depot)))) else sample_key$depot[geo_row]]
dt[is.na(depot_std) & !is.na(geo_row), depot_std := sample_key$depot[geo_row]]
dt[, bmi_v := if (!is.na(col$bmi)) suppressWarnings(as.numeric(get(col$bmi))) else sample_key$bmi[geo_row]]
dt[, bmi_group := cut(bmi_v, c(0, 25, 30, 100), c("lean (BMI < 25)", "overweight (25-30)", "obese (BMI >= 30)"), right = FALSE)]
dt[, diabetes := if (!is.na(col$diab)) as.character(get(col$diab)) else NA_character_]
cat("Diabetes information in the cell table:", if (!is.na(col$diab)) paste(unique(dt$diabetes), collapse = ", ") else "NONE (no diabetes column; groups can only be defined by BMI)", "\n")

# ---- 4. endothelial cells
dt[, is_endo := grepl("endothel", cell_type, ignore.case = TRUE) & !grepl("lymph|LEC", cell_type, ignore.case = TRUE)]
dt[, is_lec := grepl("lymph|^LEC$", cell_type, ignore.case = TRUE)]
cat("Cell types:\n"); print(sort(table(dt$cell_type), decreasing = TRUE))
n_endo <- sum(dt$is_endo); cat("Blood endothelial cells:", n_endo, "| lymphatic endothelial:", sum(dt$is_lec), "\n")
by_depot_group <- dt[is_endo == TRUE, .(endothelial_cells = .N, samples = uniqueN(sample_id)), by = .(depot_std, bmi_group, diabetes)][order(depot_std, bmi_group)]
by_sample <- dt[, .(cells = .N, endothelial = sum(is_endo), depot = depot_std[1], bmi = bmi_v[1], bmi_group = bmi_group[1], diabetes = diabetes[1]), by = sample_id][order(-endothelial)]
print(by_depot_group); print(by_sample)

# ---- 5. gene list: present in the table (as columns) ?
genes_in_table <- intersect(GENE_LIST, names(dt)); cat("Gene-list genes found as columns in this table:", if (length(genes_in_table)) paste(genes_in_table, collapse = ", ") else "none (gene counts are in the matrix file, not in this table)", "\n")

# ---- 6. verdict
has <- function(x) !is.na(x) && length(x) > 0
ctrl_lean <- by_sample[bmi_group == "lean (BMI < 25)" & is.na(diabetes) | grepl("lean|control|normal|no", diabetes, ignore.case = TRUE), .N]
verdict <- data.table(
  requirement = c("Control group (lean or obese without diabetes)", "Type 2 diabetes group", "Endothelial cells present", "More than 10,000 endothelial cells", "Human adipose tissue", "Both subcutaneous and visceral", "Gene list detectable"),
  result = c(sprintf("%s: %d samples are lean (BMI < 25); %d obese (BMI >= 30)", ifelse(sum(by_sample$bmi_group == "lean (BMI < 25)", na.rm = TRUE) > 0, "MET (by BMI)", "NOT MET"), sum(by_sample$bmi_group == "lean (BMI < 25)", na.rm = TRUE), sum(by_sample$bmi_group == "obese (BMI >= 30)", na.rm = TRUE)),
             if (!is.na(col$diab) && any(grepl("t2d|type 2|diabet", dt$diabetes, ignore.case = TRUE) & !grepl("non|no |control", dt$diabetes, ignore.case = TRUE))) "MET (diabetes label present)" else "NOT MET: no type 2 diabetes label",
             ifelse(n_endo > 0, sprintf("MET: %d blood endothelial cells", n_endo), "NOT FOUND"),
             sprintf("%s: %d (threshold %d)", ifelse(n_endo > MIN_ENDO, "MET", "NOT MET"), n_endo, MIN_ENDO),
             ifelse(nrow(dt) > 0, "MET (human rows kept)", "NOT MET"),
             sprintf("%s: depots found = %s", ifelse(all(c("subcutaneous", "visceral") %in% dt[is_endo == TRUE, depot_std]), "MET", "NOT MET / unknown"), paste(unique(dt[is_endo == TRUE, depot_std]), collapse = ", ")),
             if (length(genes_in_table)) "partly: genes found as columns" else "needs the count matrix"))
print(verdict)
write_xlsx(list(Verdict = verdict, Endothelial_by_depot_group = by_depot_group, By_sample = by_sample, Columns = profile, Cell_types = as.data.table(sort(table(dt$cell_type), decreasing = TRUE)), GEO_sample_records = sample_key), file.path(OUT_DIR, "GSE176171_criteria_check.xlsx"))
fwrite(verdict, file.path(OUT_DIR, "verdict.csv")); cat("Saved to", OUT_DIR, "\n")
