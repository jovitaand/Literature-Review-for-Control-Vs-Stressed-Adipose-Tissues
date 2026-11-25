"""Export a Tabula Sapiens .h5ad file to plain files that R can read.

Usage:  python export_TS_fat_h5ad.py  input.h5ad  output_folder  [preferred_raw_layer]

Writes into output_folder:
  raw_counts.mtx.gz      genes x cells, un-normalised UMI counts (MatrixMarket)
  decontX_counts.mtx.gz  genes x cells, ambient-corrected counts, only if the file has such a layer
  cells.csv              one row per cell: every column of adata.obs
  genes.csv              one row per gene: every column of adata.var, plus the gene symbol
  layers.txt             the layer names found in the file (for the audit)
Nothing is normalised or changed. The original .h5ad is only read.
"""
import sys, os, gzip
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import anndata as ad

h5, out = sys.argv[1], sys.argv[2]
pref = sys.argv[3] if len(sys.argv) > 3 else "raw_counts"
os.makedirs(out, exist_ok=True)
a = ad.read_h5ad(h5)                      # needs RAM: roughly 2-3x the file size
print("cells x genes:", a.shape)
print("layers:", list(a.layers.keys()))
open(os.path.join(out, "layers.txt"), "w").write("\n".join(["X"] + list(a.layers.keys())))

def looks_like_counts(m):
    """True if the first rows hold only whole numbers >= 0 (un-normalised counts)."""
    x = m[:200]
    x = x.data if sp.issparse(x) else np.asarray(x).ravel()
    return x.size > 0 and np.all(x >= 0) and np.allclose(x, np.round(x))

def pick_raw():
    for name in [pref, "raw_counts", "counts"]:
        if name in a.layers and looks_like_counts(a.layers[name]):
            return name, a.layers[name]
    if a.raw is not None and looks_like_counts(a.raw.X):
        return "adata.raw.X", a.raw.X
    if looks_like_counts(a.X):
        return "X", a.X
    sys.exit("No layer with whole-number counts was found. Stop: the raw counts are needed for the pseudobulk test.")

def write_mtx(m, name):
    m = sp.csc_matrix(m.T) if sp.issparse(m) else sp.csc_matrix(np.asarray(m).T)   # genes x cells
    with gzip.open(os.path.join(out, name), "wb") as f:
        mmwrite(f, m, field="integer" if np.allclose(m.data, np.round(m.data)) else "real")

nm, m = pick_raw()
print("raw counts taken from:", nm)
write_mtx(m, "raw_counts.mtx.gz")
open(os.path.join(out, "raw_counts_source.txt"), "w").write(nm)
if "decontXcounts" in a.layers:
    write_mtx(a.layers["decontXcounts"], "decontX_counts.mtx.gz")
    print("decontXcounts layer exported")

obs = a.obs.copy(); obs.insert(0, "cell_id", obs.index.astype(str))
obs.to_csv(os.path.join(out, "cells.csv"), index=False)
var = a.var.copy(); var.insert(0, "feature_id", var.index.astype(str))
for c in ["gene_symbol", "feature_name", "gene_name", "symbol"]:
    if c in var.columns:
        var.insert(1, "symbol", var[c].astype(str)); break
else:
    var.insert(1, "symbol", var.index.astype(str))
var.to_csv(os.path.join(out, "genes.csv"), index=False)
print("done ->", out)
