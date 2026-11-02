suppressPackageStartupMessages(library(limma))
h <- read.csv("hamec.csv", check.names = FALSE, fileEncoding = "UTF-8-BOM"); rownames(h) <- h$target_id
S <- c("CT_1","CT_2","CT_3","EndoMT_1","EndoMT_2","EndoMT_3"); X <- as.matrix(h[, S]); L <- log2(X + 1)
keep <- rowMeans(X) >= 1; cat("genes tested (mean TPM >= 1):", sum(keep), "of", nrow(X), "\n")
des <- model.matrix(~ factor(rep(c("CT","EndoMT"), each = 3), levels = c("CT","EndoMT")))
fit <- eBayes(lmFit(L[keep, ], des), trend = TRUE, robust = TRUE); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt); res <- merge(tt, h[, c("target_id","pval","qval","log2FC")], by.x = "gene", by.y = "target_id")
res$log2FC <- suppressWarnings(as.numeric(res$log2FC)); res$qval <- as.numeric(res$qval)
cat("limma genes adj.P < 0.05:", sum(res$adj.P.Val < 0.05), "| file q < 0.05 among same genes:", sum(res$qval < 0.05, na.rm = TRUE), "\n")
cat("overlap (both significant):", sum(res$adj.P.Val < 0.05 & res$qval < 0.05, na.rm = TRUE), "\n")
cat("Spearman of p-values (limma vs file):", round(cor(res$P.Value, res$pval, method = "spearman", use = "complete.obs"), 3), "| Pearson of log2FC:", round(cor(res$logFC, res$log2FC, use = "complete.obs"), 3), "\n")
G <- c("PECAM1","CDH5","NOS3","HIF1A","ERG","VCAM1","ICAM1","ANGPT2"); g <- res[match(G, res$gene), c("gene","logFC","P.Value","adj.P.Val","log2FC","pval","qval")]
g$mean_CT <- rowMeans(X[G, 1:3]); g$mean_EndoMT <- rowMeans(X[G, 4:6]); g[, 2:7] <- round(g[, 2:7], 4); print(g, row.names = FALSE)
write.csv(g, "gene_list_limma_vs_file.csv", row.names = FALSE)
