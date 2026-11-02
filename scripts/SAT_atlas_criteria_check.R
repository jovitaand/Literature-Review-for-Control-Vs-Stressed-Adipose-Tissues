# Criteria check for the SAT vascular atlas Seurat object (Kalucka lab; code at github.com/Kalucka-lab/SAT_Atlas).
# The repository holds code only; the objects it reads are e.g. Vascular_IntegratedObj_annotated.RDS (metadata columns used in the repo:
# Vascular_labels, Condition = Lean / Obese / Diabetic Obese, gender, Study_chemistry, orig.ident / samples).
# For Google Colab (R runtime, High-RAM) or the pipeline. Tested on a synthetic object with these columns; NOT yet run on the real atlas.

# ---- SETTINGS -----------------------------------------------------------
RDS <- ""                          # path to the atlas .RDS; "" = search /content and the working folder
OUT_DIR <- "/content/sat_atlas_check"
GENE_LIST <- c(PECAM1 = "PECAM1", CDH5 = "CDH5", NOS3 = "NOS3", HIF1A = "HIF1A", ERG = "ERG", VCAM1 = "VCAM1", ICAM1 = "ICAM1", ANGPT2 = "ANGPT2")
MIN_EC <- 10000; MIN_EC_PER_DONOR <- 20
# -------------------------------------------------------------------------
for (p in c("Seurat", "data.table", "writexl", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
if (!requireNamespace("edgeR", quietly = TRUE)) { install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install("edgeR", ask = FALSE, update = FALSE) }
suppressPackageStartupMessages({library(Seurat); library(data.table); library(writexl); library(Matrix); library(edgeR)}); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (RDS == "") { h <- list.files(c("/content", "."), pattern = "\\.RDS$|\\.rds$", recursive = TRUE, full.names = TRUE); h <- h[grepl("Vascular|Atlas|SAT", h, ignore.case = TRUE)]; if (length(h)) RDS <- h[1] }
stopifnot(file.exists(RDS)); cat("Reading", RDS, round(file.size(RDS) / 1e6), "MB\n"); obj <- readRDS(RDS); md <- as.data.table(obj@meta.data, keep.rownames = "cell"); cat("Cells:", nrow(md), "\n")

# ---- 1. which metadata columns?
pick <- function(pat) { h <- names(md)[grepl(pat, names(md), ignore.case = TRUE)]; if (length(h)) h[1] else NA_character_ }
col <- list(label = pick("^Vascular_labels$|label|celltype|cell_type|annotation"), cond = pick("^Condition$|condition|diab|group"), sample = pick("^samples$|^orig.ident$|sample|donor"), sex = pick("gender|sex"),
            study = pick("Study_chemistry|^Study$"), depot = pick("depot|tissue|region|site|fat"))
print(unlist(col)); stopifnot(!is.na(col$label), !is.na(col$cond), !is.na(col$sample))
md[, label := as.character(get(col$label))]; md[, cond := as.character(get(col$cond))]; md[, sample := as.character(get(col$sample))]
if (!is.na(col$depot)) { cat("Depot values:\n"); print(table(md[[col$depot]])) } else cat("No depot column: the repository states the atlas is subcutaneous (SAT) only.\n")
cat("Condition values:\n"); print(table(md$cond))

# ---- 2. cell classes
md[, class := ifelse(grepl("lymph", label, ignore.case = TRUE), "lymphatic endothelial", ifelse(grepl("endothelial", label, ignore.case = TRUE), "blood endothelial (incl. sub-endothelial)", ifelse(grepl("pericyte|vsmc|smooth", label, ignore.case = TRUE), "mural (pericyte / smooth muscle)", "other")))]
by_label <- dcast(md[, .N, by = .(label, cond)], label ~ cond, value.var = "N", fill = 0L); by_class <- dcast(md[, .N, by = .(class, cond)], class ~ cond, value.var = "N", fill = 0L)
by_class[, total := rowSums(.SD), .SDcols = setdiff(names(by_class), "class")]; print(by_class); print(by_label)
ec <- md[class %in% c("blood endothelial (incl. sub-endothelial)", "lymphatic endothelial")]; n_ec <- nrow(ec); n_blood <- sum(md$class == "blood endothelial (incl. sub-endothelial)")
donors <- md[, .(donors = uniqueN(sample), cells = .N), by = cond]; ec_donor <- ec[, .(ec_cells = .N), by = .(cond, sample)][order(cond, -ec_cells)]; ec_donor[, enough := ec_cells >= MIN_EC_PER_DONOR]
print(donors); print(ec_donor[, .(donors = .N, donors_with_enough = sum(enough), ec_cells = sum(ec_cells)), by = cond])

# ---- 3. gene list in endothelial cells (raw counts, RNA assay)
DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cnt <- tryCatch(GetAssayData(obj, layer = "counts"), error = function(e) GetAssayData(obj, slot = "counts")); G <- GENE_LIST[GENE_LIST %in% rownames(cnt)]; cat("Genes found:", paste(G, collapse = ", "), "| missing:", paste(setdiff(GENE_LIST, rownames(cnt)), collapse = ", "), "\n")
ecc <- ec$cell; x <- cnt[unname(G), ecc, drop = FALSE]
det <- rbindlist(lapply(unique(ec$cond), function(cn) { k <- which(ec$cond == cn); data.table(condition = cn, ec_cells = length(k), gene = names(G), cells_detected = as.integer(Matrix::rowSums(x[, k, drop = FALSE] > 0)), pct_detected = round(100 * Matrix::rowMeans(x[, k, drop = FALSE] > 0), 1)) }))
print(det[order(gene, condition)])

# ---- 4. exploratory pseudobulk DE of the gene list per donor (T2D vs obese, T2D vs lean, obese vs lean), blood endothelial cells only
be <- md[class == "blood endothelial (incl. sub-endothelial)"]; xb <- cnt[, be$cell, drop = FALSE]; sm <- unique(be[, .(sample, cond)])
pb <- sapply(sm$sample, function(s) Matrix::rowSums(xb[, be$sample == s, drop = FALSE])); ok <- be[, .N, by = sample][N >= MIN_EC_PER_DONOR, sample]
one <- function(a, b) { d <- sm[sample %in% ok & cond %in% c(a, b)]; if (sum(d$cond == a) < 2 || sum(d$cond == b) < 2) return(NULL)
  y <- DGEList(pb[, d$sample], group = factor(d$cond, c(b, a))); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y); des <- model.matrix(~ y$samples$group); y <- estimateDisp(y, des)
  t <- topTags(glmQLFTest(glmQLFit(y, des), coef = 2), n = Inf)$table; data.table(contrast = paste(a, "vs", b), gene = rownames(t), logFC = round(t$logFC, 3), p_value = t$PValue, FDR = t$FDR, n_a = sum(d$cond == a), n_b = sum(d$cond == b)) }
conds <- unique(sm$cond); lean <- grep("lean", conds, value = TRUE, ignore.case = TRUE)[1]; dob <- grep("diab", conds, value = TRUE, ignore.case = TRUE)[1]; ob <- setdiff(grep("obes", conds, value = TRUE, ignore.case = TRUE), dob)[1]
de <- rbindlist(list(one(dob, ob), one(dob, lean), one(ob, lean)), fill = TRUE); gl_de <- if (nrow(de)) de[gene %in% GENE_LIST][order(contrast, p_value)] else data.table(); print(gl_de)

# ---- 5. verdict
has <- function(pat) any(grepl(pat, md$cond, ignore.case = TRUE))
verdict <- data.table(requirement = c("Control group (lean or obese without diabetes)", "Type 2 diabetes group", "Endothelial cells present", "More than 10,000 endothelial cells", "Human adipose tissue", "Both subcutaneous and visceral", "Gene list detectable"),
 result = c(sprintf("%s: conditions = %s", ifelse(has("lean") || has("^obese$"), "MET", "NOT MET"), paste(unique(md$cond), collapse = ", ")), ifelse(has("diab|t2d"), "MET (a diabetic group exists)", "NOT MET"), ifelse(n_ec > 0, sprintf("MET: %d endothelial cells (%d blood, %d lymphatic)", n_ec, n_blood, n_ec - n_blood), "NOT FOUND"),
            sprintf("%s: %d endothelial cells (threshold %d)", ifelse(n_ec > MIN_EC, "MET", "NOT MET"), n_ec, MIN_EC), "MET (human; per paper)",
            if (!is.na(col$depot) && all(c("subcutaneous", "visceral") %in% tolower(md[[col$depot]]))) "MET" else "NOT MET / not recorded: the atlas is subcutaneous only (repository README)", sprintf("%d of 8 genes present", length(G))))
print(verdict)
write_xlsx(list(Verdict = verdict, Cells_by_class = by_class, Cells_by_label = by_label, Donors_by_condition = donors, EC_per_donor = ec_donor, Gene_detection = det, Gene_list_DE_exploratory = gl_de, Columns_used = data.table(item = names(col), column = unlist(col))), file.path(OUT_DIR, "SAT_atlas_criteria_check.xlsx"))
fwrite(verdict, file.path(OUT_DIR, "verdict.csv")); cat("Saved to", OUT_DIR, "\n")
