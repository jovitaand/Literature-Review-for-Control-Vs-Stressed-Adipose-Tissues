#!/usr/bin/env Rscript
# Combine the per-dataset pseudobulk results of GSE129363 and GSE278526 into one T2D-vs-Control result.
#
# Usage:
#   Rscript scripts/combine_ERG_datasets.R <GSE129363 panel csv> <GSE278526 panel csv> [out_dir]
# Inputs (written by the two notebooks):
#   GSE129363: results_v3/DE/Panel_DE_T2D_vs_Control_all_subsets.csv
#   GSE278526: erg_panel_results_v2/panel_DE_T2D_vs_Control.csv
#
# Design (why it is built this way)
#  * Each dataset keeps its own normalisation, QC and donor-level model. Cells are never merged, because the datasets
#    differ in tissue prep (SVF vs sorted stromal cells), depth and library structure.
#  * The two per-dataset log2 fold changes (T2D relative to Control) are then pooled per depot x EC subset x gene with an
#    inverse-variance meta-analysis (fixed-effect and DerSimonian-Laird random-effects).
#  * A dataset only contributes to a pooled estimate if its test is interpretable: tested, enough donors per group and
#    enough reads for that gene. Excluded rows are still listed, with the reason.
#  * Neither notebook exports a standard error, so SE is recovered from logFC and the two-sided p-value
#    (SE = |logFC| / z, z = qnorm(1 - p/2)). This is an approximation and is stated in the output.
#  * With only two datasets the between-dataset variance (tau^2, I^2) is poorly estimated; both models are reported.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript combine_ERG_datasets.R <GSE129363 csv> <GSE278526 csv> [out_dir]")
f_129363 <- args[1]; f_278526 <- args[2]
out_dir <- if (length(args) >= 3) args[3] else "combined_results"

MIN_COUNTS           <- 50    # minimum total reads of the gene in the tested pseudobulks
MIN_DONORS_PER_GROUP <- 3
FDR_CUTOFF           <- 0.05
PANEL <- c("ERG", "PECAM1", "CDH5", "NOS3", "VCAM1", "ICAM1", "HIF1A", "ANGPT2")

suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. read and harmonise ----------------------------------------------------------------------
harmonise <- function(path, name) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- sub("^n_control$", "n_Control", names(d))
  need <- c("depot", "subset", "gene", "n_Control", "n_T2D", "total_counts", "status", "logFC", "PValue")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop(name, ": missing columns: ", paste(miss, collapse = ", "))
  d %>% transmute(dataset = name, depot, subset, gene,
                  n_Control = as.integer(n_Control), n_T2D = as.integer(n_T2D),
                  total_counts = as.numeric(total_counts), status = as.character(status),
                  logFC = as.numeric(logFC), PValue = as.numeric(PValue)) %>%
    filter(gene %in% PANEL)
}
long <- bind_rows(harmonise(f_129363, "GSE129363"), harmonise(f_278526, "GSE278526"))

# ---- 2. SE from p, eligibility ------------------------------------------------------------------
p_safe <- pmin(pmax(long$PValue, 1e-300), 1 - 1e-12)
long$z  <- qnorm(p_safe / 2, lower.tail = FALSE)
long$se <- abs(long$logFC) / long$z
long$se[!is.finite(long$se) | long$se <= 0] <- NA
long$reason <- with(long, ifelse(status != "tested" | !is.finite(logFC) | !is.finite(PValue), "not tested",
                          ifelse(is.na(se), "no usable standard error",
                          ifelse(pmin(n_Control, n_T2D) < MIN_DONORS_PER_GROUP,
                                 paste0("fewer than ", MIN_DONORS_PER_GROUP, " donors in a group"),
                          ifelse(is.na(total_counts) | total_counts < MIN_COUNTS,
                                 paste0("fewer than ", MIN_COUNTS, " reads (not interpretable)"), "eligible")))))
long$eligible <- long$reason == "eligible"
write.csv(long %>% arrange(gene != "ERG", gene, depot, subset, dataset) %>% select(-z),
          file.path(out_dir, "per_dataset_results_with_eligibility.csv"), row.names = FALSE)

# ---- 3. pool per depot x subset x gene ----------------------------------------------------------
pool_one <- function(d) {
  e <- d[d$eligible, ]
  k <- nrow(e)
  base <- data.frame(k_datasets = k, datasets_used = paste(e$dataset, collapse = " + "),
                     donors_Control = sum(e$n_Control), donors_T2D = sum(e$n_T2D),
                     total_reads = sum(e$total_counts))
  if (k == 0) return(cbind(base, data.frame(method = "none", logFC = NA, se = NA, CI_low = NA, CI_high = NA, p = NA,
                                            Q = NA, I2 = NA, tau2 = NA, fixed_logFC = NA, fixed_p = NA)))
  w <- 1 / e$se^2
  fixed <- sum(w * e$logFC) / sum(w); se_f <- 1 / sqrt(sum(w))
  if (k == 1) {
    return(cbind(base, data.frame(method = "single dataset", logFC = e$logFC, se = e$se,
                                  CI_low = e$logFC - 1.96 * e$se, CI_high = e$logFC + 1.96 * e$se,
                                  p = e$PValue, Q = NA, I2 = NA, tau2 = NA, fixed_logFC = fixed,
                                  fixed_p = 2 * pnorm(-abs(fixed / se_f)))))
  }
  Q <- sum(w * (e$logFC - fixed)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (e$se^2 + tau2)
  re <- sum(wr * e$logFC) / sum(wr); se_r <- 1 / sqrt(sum(wr))
  cbind(base, data.frame(method = "random-effects (DL)", logFC = re, se = se_r,
                         CI_low = re - 1.96 * se_r, CI_high = re + 1.96 * se_r,
                         p = 2 * pnorm(-abs(re / se_r)), Q = Q, I2 = max(0, (Q - (k - 1)) / Q) * 100, tau2 = tau2,
                         fixed_logFC = fixed, fixed_p = 2 * pnorm(-abs(fixed / se_f))))
}
keys <- long %>% distinct(depot, subset, gene)
pooled <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  d <- long[long$depot == keys$depot[i] & long$subset == keys$subset[i] & long$gene == keys$gene[i], ]
  cbind(keys[i, ], pool_one(d))
}))
pooled$fold_change <- 2^pooled$logFC
# Multiple testing: all pooled tests with k >= 2 are corrected together; single-dataset results are not "combined"
ok2 <- pooled$k_datasets >= 2 & is.finite(pooled$p)
pooled$combined_FDR <- NA_real_
pooled$combined_FDR[ok2] <- p.adjust(pooled$p[ok2], method = "BH")
pooled$evidence <- with(pooled, ifelse(k_datasets >= 2, "pooled across both datasets",
                                ifelse(k_datasets == 1, paste("only", datasets_used, "interpretable"),
                                       "not interpretable in either dataset")))
pooled$call <- with(pooled, ifelse(k_datasets == 0, "not testable",
                            ifelse(k_datasets == 1, "single dataset: no replication (exploratory)",
                            ifelse(combined_FDR < FDR_CUTOFF, "significant (combined FDR < 0.05)", "not significant"))))
pooled <- pooled %>% arrange(gene != "ERG", match(gene, PANEL), depot, subset)
write.csv(pooled, file.path(out_dir, "combined_panel_results.csv"), row.names = FALSE)
write.csv(pooled %>% filter(gene == "ERG"), file.path(out_dir, "combined_ERG_results.csv"), row.names = FALSE)

# ---- 4. forest plots ----------------------------------------------------------------------------
forest <- function(gene_name, subsets = NULL, title) {
  d <- long[long$gene == gene_name & is.finite(long$logFC) & !is.na(long$se), ]
  if (!is.null(subsets)) d <- d[d$subset %in% subsets, ]
  if (!nrow(d)) return(NULL)
  d$row <- paste(d$depot, gsub("_", " ", d$subset))
  d$used <- ifelse(d$eligible, "used in pooled estimate", "excluded (see reason)")
  d$lo <- d$logFC - 1.96 * d$se; d$hi <- d$logFC + 1.96 * d$se
  p <- pooled[pooled$gene == gene_name & pooled$k_datasets >= 2 & is.finite(pooled$logFC), ]
  if (!is.null(subsets)) p <- p[p$subset %in% subsets, ]
  p$row <- paste(p$depot, gsub("_", " ", p$subset))
  g <- ggplot(d, aes(logFC, row)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = dataset, alpha = used), orientation = "y", width = 0.25,
                  position = position_dodge(0.6)) +
    geom_point(aes(colour = dataset, alpha = used, shape = used), size = 2.8, position = position_dodge(0.6)) +
    scale_alpha_manual(values = c("used in pooled estimate" = 1, "excluded (see reason)" = 0.35), name = NULL) +
    scale_shape_manual(values = c("used in pooled estimate" = 16, "excluded (see reason)" = 1), name = NULL) +
    scale_colour_manual(values = c(GSE129363 = "#2A7F8E", GSE278526 = "#C8622D"), name = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.box = "vertical") +
    labs(x = paste0(gene_name, " log2 fold change, T2D relative to Control (95% CI)"), y = NULL, title = title,
         subtitle = paste0("Faded hollow points: excluded (< ", MIN_COUNTS, " reads, < ", MIN_DONORS_PER_GROUP,
                           " donors per group).\nBlack diamond: pooled estimate (needs >= 2 interpretable datasets)."))
  if (nrow(p)) g <- g + geom_errorbar(data = p, aes(xmin = CI_low, xmax = CI_high, y = row), orientation = "y",
                                      width = 0.18, colour = "black", inherit.aes = FALSE) +
    geom_point(data = p, aes(logFC, row), shape = 18, size = 4.5, colour = "black", inherit.aes = FALSE)
  g
}
save_forest <- function(g, name, w, h) if (!is.null(g)) {
  ggsave(file.path(out_dir, paste0(name, ".png")), g, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".svg")), g, width = w, height = h, bg = "white")
}
save_forest(forest("ERG", NULL, "ERG: T2D vs Control by depot and EC subset, both datasets"), "ERG_forest_combined", 10, 7)
panel_plots <- lapply(PANEL, function(g) forest(g, "All_blood_ECs", paste(g, "- all blood ECs")))
for (g in seq_along(PANEL)) save_forest(panel_plots[[g]], paste0("forest_", PANEL[g], "_all_blood_ECs"), 9, 3.8)

# ---- 5. plain-language summary ------------------------------------------------------------------
fmt <- function(x, n = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = n))
lines <- c("COMBINED T2D vs Control RESULT (positive logFC = higher in T2D)",
           sprintf("Inclusion rules: tested; >= %d donors per group; >= %d reads for the gene. SE recovered from logFC and p.",
                   MIN_DONORS_PER_GROUP, MIN_COUNTS), "")
e <- pooled[pooled$gene == "ERG", ]
for (i in seq_len(nrow(e))) {
  r <- e[i, ]; where <- paste(r$depot, gsub("_", " ", r$subset))
  lines <- c(lines, if (r$k_datasets == 0) sprintf("ERG %-22s not interpretable in either dataset", where)
             else sprintf("ERG %-22s %s | logFC %s (95%% CI %s to %s), p = %s | %s",
                          where, r$evidence, fmt(r$logFC), fmt(r$CI_low), fmt(r$CI_high), fmt(r$p, 3), r$call))
}
ex <- long[long$gene == "ERG" & !long$eligible, ]
if (nrow(ex)) lines <- c(lines, "", "ERG rows excluded from pooling:",
                         sprintf("  %-10s %-22s %s", ex$dataset, paste(ex$depot, gsub("_", " ", ex$subset)), ex$reason))
sig <- pooled[pooled$call == "significant (combined FDR < 0.05)", ]
lines <- c(lines, "", sprintf("Panel genes significant after combined FDR (pooled tests only): %d", nrow(sig)))
if (nrow(sig)) lines <- c(lines, sprintf("  %s %s: logFC %s, FDR %s", sig$gene, paste(sig$depot, sig$subset), fmt(sig$logFC), fmt(sig$combined_FDR, 3)))
lines <- c(lines, "", "Read with care: with two datasets, heterogeneity (I2, tau2) is unreliable; disease and library are confounded in",
           "GSE278526 (one 10x library per depot x group); diabetes labels in GSE129363 differ between the annotation file and GEO.")
writeLines(lines, file.path(out_dir, "summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
cat("\nWrote results to", normalizePath(out_dir), "\n")
