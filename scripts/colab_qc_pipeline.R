# Archive inventory and first-pass QC, written to run in Google Colab with an R runtime.
# NOT YET RUN: written without access to R, so expect small fixes on the first run.
#
# Usage in Colab (R runtime):
#   1. Mount Drive:  system("ls /content/drive") after using the Drive button, or copy the ZIP to /content.
#   2. Edit the three paths in CONFIG below.
#   3. source("colab_qc_pipeline.R")
# Only the small files in OUT_DIR are meant to be shared; the raw data never leave Colab.

# ---- CONFIG -------------------------------------------------------------
ARCHIVE   <- "/content/drive/MyDrive/dataset.zip"   # your 1.06 GB archive
RAW_DIR   <- "/content/raw"                         # extracted, treated as read-only afterwards
OUT_DIR   <- "/content/qc_out"                      # small reports (zip this and upload it)
SAMPLE_N  <- 50000                                  # rows read per table for profiling
MAX_DEPTH <- 6                                      # nested-archive depth
# -------------------------------------------------------------------------

pkgs <- c("data.table", "digest", "readxl", "R.utils")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(data.table); library(digest); library(readxl)})

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
errors <- list()
log_err <- function(path, err, action) errors[[length(errors) + 1]] <<- data.table(path = path, error = as.character(err), next_action = action)

# ---- Phase 1: recursive extraction --------------------------------------
is_archive <- function(f) grepl("\\.(zip|tar|tar\\.gz|tgz|gz|7z)$", f, ignore.case = TRUE) && !grepl("\\.(tsv|csv|txt|mtx)\\.gz$", f, ignore.case = TRUE)
extract_one <- function(f, dest) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    if (grepl("\\.zip$", f, ignore.case = TRUE)) utils::unzip(f, exdir = dest)
    else if (grepl("\\.(tar|tar\\.gz|tgz)$", f, ignore.case = TRUE)) utils::untar(f, exdir = dest)
    else if (grepl("\\.7z$", f, ignore.case = TRUE)) { if (Sys.which("7z") == "") stop("7z not installed"); system2("7z", c("x", "-y", paste0("-o", dest), shQuote(f))) }
    else if (grepl("\\.gz$", f, ignore.case = TRUE)) R.utils::gunzip(f, destname = file.path(dest, sub("\\.gz$", "", basename(f))), remove = FALSE, overwrite = TRUE)
    TRUE
  }, error = function(e) { log_err(f, conditionMessage(e), "check password / format / install 7z"); FALSE })
}
extracted <- character(0)
queue <- ARCHIVE; depth <- 0
if (!file.exists(ARCHIVE)) stop("Archive not found: ", ARCHIVE)
ok <- extract_one(ARCHIVE, RAW_DIR); extracted <- c(extracted, ARCHIVE)
while (depth < MAX_DEPTH) {
  all <- list.files(RAW_DIR, recursive = TRUE, full.names = TRUE)
  nested <- setdiff(all[vapply(all, is_archive, logical(1))], extracted)
  if (!length(nested)) break
  for (f in nested) { extract_one(f, paste0(f, "_extracted")); extracted <- c(extracted, f) }
  depth <- depth + 1
}
cat("Extraction passes:", depth, "| nested archives extracted:", length(extracted) - 1, "\n")

# ---- Phase 2: inventory with SHA-256 ------------------------------------
files <- list.files(RAW_DIR, recursive = TRUE, full.names = TRUE)
files <- files[!file.info(files)$isdir]
rel   <- sub(paste0("^", RAW_DIR, "/?"), "", files)
inv <- data.table(
  relative_path = rel, filename = basename(files),
  extension = tolower(tools::file_ext(sub("\\.gz$", "", files))),
  size_bytes = file.info(files)$size,
  sha256 = vapply(files, function(f) tryCatch(digest(f, algo = "sha256", file = TRUE), error = function(e) NA_character_), character(1)),
  is_nested_archive = vapply(files, is_archive, logical(1))
)
inv[, duplicate_group := if (.N > 1) .GRP else NA_integer_, by = sha256]
inv[, duplicate_of_n := .N, by = sha256]

# ---- Phase 3: structural profiling of each readable table ---------------
table_ext <- "(csv|tsv|txt|tab|mtx|xlsx|xls)$"
profile_file <- function(path, rel_path) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  out <- list(readable = FALSE, n_rows = NA_integer_, n_cols = NA_integer_, delimiter = NA, col_names = NA, col_types = NA, note = "")
  tryCatch({
    if (ext %in% c("xlsx", "xls")) {
      sheets <- excel_sheets(path); dt <- as.data.table(read_excel(path, sheet = 1, n_max = SAMPLE_N))
      out$note <- paste("sheets:", paste(sheets, collapse = "|"))
    } else if (ext == "mtx") {
      hdr <- readLines(gzfile(path), n = 3); out$note <- paste(hdr, collapse = " // "); dt <- NULL
      out$readable <- TRUE
    } else if (ext %in% c("csv", "tsv", "txt", "tab")) {
      dt <- fread(path, nrows = SAMPLE_N, data.table = TRUE, showProgress = FALSE)
      out$delimiter <- if (ext == "csv") "," else if (ext == "tsv") "\t" else "auto"
    } else return(out)
    if (!is.null(dt)) {
      out$readable <- TRUE; out$n_cols <- ncol(dt); out$n_rows <- nrow(dt)
      out$col_names <- paste(names(dt), collapse = "|"); out$col_types <- paste(vapply(dt, function(x) class(x)[1], ""), collapse = "|")
    }
  }, error = function(e) { out$note <<- conditionMessage(e); log_err(path, conditionMessage(e), "inspect manually") })
  out
}
prof_rows <- list()
for (i in seq_len(nrow(inv))) {
  if (grepl(table_ext, sub("\\.gz$", "", inv$filename[i]), ignore.case = TRUE)) {
    p <- profile_file(files[i], inv$relative_path[i])
    prof_rows[[length(prof_rows) + 1]] <- data.table(relative_path = inv$relative_path[i], readable = p$readable, sampled_rows = p$n_rows, n_cols = p$n_cols,
      delimiter = p$delimiter, column_names = p$col_names, column_types = p$col_types, note = p$note)
  }
}
profile <- if (length(prof_rows)) rbindlist(prof_rows) else data.table()
inv <- merge(inv, profile, by = "relative_path", all.x = TRUE)

# ---- Phase 5: first-pass QC on each readable table (flag, never remove) --
qc <- list(); dd <- list()
add_qc <- function(file, check, value, severity, action) qc[[length(qc) + 1]] <<- data.table(file = file, check = check, value = as.character(value), severity = severity, recommended_action = action)
for (i in seq_len(nrow(inv))) {
  if (!isTRUE(inv$readable[i])) next
  path <- files[i]; ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  if (!ext %in% c("csv", "tsv", "txt", "tab", "xlsx", "xls")) next
  dt <- tryCatch(if (ext %in% c("xlsx", "xls")) as.data.table(read_excel(path, n_max = SAMPLE_N)) else fread(path, nrows = SAMPLE_N, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) { add_qc(inv$relative_path[i], "empty_table", 0, "warning", "check file"); next }
  nm <- names(dt)
  if (anyDuplicated(nm)) add_qc(inv$relative_path[i], "duplicate_column_names", paste(nm[duplicated(nm)], collapse = "|"), "warning", "rename, keep mapping")
  miss <- vapply(dt, function(x) mean(is.na(x) | (is.character(x) & !is.na(x) & x == "")), numeric(1))
  for (j in which(miss > 0.2)) add_qc(inv$relative_path[i], paste0("missing_over_20pct:", nm[j]), round(miss[j], 3), "info", "review column")
  const <- vapply(dt, function(x) uniqueN(x) <= 1, logical(1))
  for (j in which(const)) add_qc(inv$relative_path[i], paste0("constant_column:", nm[j]), "1 unique value", "info", "likely drop for modelling; keep in raw")
  if (anyDuplicated(dt)) add_qc(inv$relative_path[i], "duplicate_rows_in_sample", sum(duplicated(dt)), "warning", "check whether expected")
  num <- names(dt)[vapply(dt, is.numeric, logical(1))]
  if (length(num)) for (cn in num) { x <- dt[[cn]]; if (any(x < 0, na.rm = TRUE)) add_qc(inv$relative_path[i], paste0("negative_values:", cn), sum(x < 0, na.rm = TRUE), "info", "confirm expected (e.g. log values)") }
  for (j in seq_along(nm)) dd[[length(dd) + 1]] <- data.table(file = inv$relative_path[i], original_column = nm[j], standardized_column = tolower(gsub("[^A-Za-z0-9]+", "_", trimws(nm[j]))),
    original_type = class(dt[[j]])[1], units = "unknown", meaning = "unknown", transformation_history = "none")
}
qc <- if (length(qc)) rbindlist(qc) else data.table(file = character(), check = character(), value = character(), severity = character(), recommended_action = character())
dd <- if (length(dd)) rbindlist(dd) else data.table()

# ---- Write small outputs ------------------------------------------------
fwrite(inv, file.path(OUT_DIR, "file_inventory.csv"))
fwrite(inv[duplicate_of_n > 1], file.path(OUT_DIR, "duplicate_files.csv"))
fwrite(dd, file.path(OUT_DIR, "data_dictionary.csv"))
fwrite(qc, file.path(OUT_DIR, "qc_summary.csv"))
fwrite(if (length(errors)) rbindlist(errors) else data.table(path = character(), error = character(), next_action = character()), file.path(OUT_DIR, "extraction_errors.csv"))
fwrite(data.table(step = "none applied", rationale = "this script only inventories and flags; no data modified"), file.path(OUT_DIR, "transformation_log.csv"))

report <- c("# Processing report", "",
  paste("Files discovered:", nrow(inv)), paste("Nested archives extracted:", length(extracted) - 1),
  paste("Tables profiled (readable):", sum(inv$readable, na.rm = TRUE)), paste("Files with errors:", length(errors)),
  paste("Duplicate-content files:", sum(inv$duplicate_of_n > 1, na.rm = TRUE)), paste("QC flags:", nrow(qc)), "",
  "Row counts are from the first SAMPLE_N rows only, not the whole file. Re-run with a larger SAMPLE_N or count lines for full counts.",
  "No rows, columns or files were removed or changed. Review qc_summary.csv and extraction_errors.csv before standardizing.")
writeLines(report, file.path(OUT_DIR, "processing_report.md"))
cat(paste(report, collapse = "\n"), "\n")
cat("Outputs in", OUT_DIR, "\nTotal size (MB):", round(sum(file.info(list.files(OUT_DIR, full.names = TRUE))$size) / 1e6, 2), "\n")
