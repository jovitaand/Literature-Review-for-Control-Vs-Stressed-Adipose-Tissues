#!/usr/bin/env bash
# Install the R/Python packages the notebooks need, from the OS package manager (works when CRAN is blocked).
set -e
PKGS="r-base-core r-cran-matrix r-cran-data.table r-cran-writexl r-cran-ggplot2 r-cran-seurat r-cran-uwot r-cran-igraph r-cran-r.utils r-bioc-edger python3-openpyxl python3-pil"
apt-get install -y -qq $PKGS >/dev/null 2>&1 || { apt-get update -qq && apt-get install -y -qq $PKGS; }
Rscript -e 'for (p in c("Matrix","data.table","writexl","ggplot2","Seurat","edgeR","R.utils","uwot")) cat(sprintf("%-12s %s\n", p, requireNamespace(p, quietly = TRUE)))'
echo "harmony is not packaged for apt; the notebooks fall back to plain PCA if it is missing (set USE_HARMONY=FALSE to silence the warning)."
