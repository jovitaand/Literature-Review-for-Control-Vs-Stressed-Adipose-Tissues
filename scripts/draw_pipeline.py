"""Draw the agreed 15-step single-cell analysis outline as a flow diagram.
Usage: python3 scripts/draw_pipeline.py   -> docs/figures/analysis_pipeline.png (300 dpi) and .svg"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

NAVY, TEAL, ORANGE, GREY = "#1F3A5F", "#2A7F8E", "#C8622D", "#F4F6F8"
fig, ax = plt.subplots(figsize=(17, 12)); ax.set_xlim(0, 100); ax.set_ylim(0, 72); ax.axis("off")

def box(x, y, w, h, num, text, color=NAVY, fs=9.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8", fc="white", ec=color, lw=1.8))
    if num:
        ax.add_patch(FancyBboxPatch((x, y), 3.6, h, boxstyle="round,pad=0.25,rounding_size=0.8", fc=color, ec=color))
        ax.text(x + 1.8, y + h / 2, num, color="white", ha="center", va="center", fontsize=12, fontweight="bold")
        ax.text(x + 4.6, y + h / 2, text, ha="left", va="center", fontsize=fs, color="#1b1b1b", linespacing=1.35)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color="#1b1b1b", linespacing=1.35)

def arrow(x1, y1, x2, y2, color="#555", style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=14, lw=1.6, color=color,
                                 connectionstyle=f"arc3,rad={rad}"))

def elbow(pts, color):
    for (x1, y1), (x2, y2) in zip(pts[:-2], pts[1:-1]):
        ax.plot([x1, x2], [y1, y2], color=color, lw=1.6, solid_capstyle="round")
    arrow(*pts[-2], *pts[-1], color=color)

def lane(x, w, title, color):
    ax.add_patch(FancyBboxPatch((x, 9), w, 60.5, boxstyle="round,pad=0.3,rounding_size=1.2", fc=GREY, ec="none"))
    ax.text(x + w / 2, 67.6, title, ha="center", va="center", fontsize=12.5, fontweight="bold", color=color)

ax.text(50, 71, "Adipose-tissue endothelial single-cell analysis: Control vs Type 2 Diabetes (T2D)", ha="center", fontsize=16, fontweight="bold", color=NAVY)

# ---- lane A: whole population (1-7)
lane(1, 31, "A. Whole cell population (steps 1–7)", NAVY)
A = [("1", "Compare complete Control vs T2D\npopulations (cells per group / sample)"),
     ("2", "QC + normalisation\nRemove: genes < 200 · MT% > 20%\nSensitivity: MT 15% · snRNA-seq: 3–5%"),
     ("3", "Elbow plot + PCA\nchoose PCs and clustering resolution"),
     ("4", "UMAP of all cells\nhighlight endothelial population"),
     ("5", "DE across clusters\n(FindAllMarkers, Wilcoxon)"),
     ("6", "Annotate every cluster\nfrom its DEGs"),
     ("7", "Confirm endothelial cluster(s)\nPECAM1 · CDH5 · VWF · CLDN5 …\n(several markers, never PECAM1 alone)")]
ys = [58, 50.2, 43.4, 36.9, 30.4, 23.9, 14]
hs = [5.6, 6.6, 4.9, 4.9, 4.9, 4.9, 6.6]
for (n, t), y, h in zip(A, ys, hs):
    box(2.2, y, 28.4, h, n, t)
for i in range(len(A) - 1):
    arrow(16.4, ys[i], 16.4, ys[i + 1] + hs[i + 1] + 0.1)

# ---- lane B: endothelial subset (8-10)
lane(34, 31, "B. Endothelial subset (steps 8–10)", TEAL)
box(35.2, 58, 28.4, 6.0, "8", "Subset endothelial cluster(s)\nrecompute elbow, PCA, clustering", TEAL)
box(35.2, 49, 28.4, 6.0, "9", "UMAP of endothelial cells\nannotate subclusters from DEGs", TEAL)
box(35.2, 39.5, 28.4, 6.2, "10", "Split into two Seurat objects\nSAT   |   VAT\n(tissues analysed separately)", TEAL)
arrow(49.4, 58, 49.4, 55.1); arrow(49.4, 49, 49.4, 45.8, )
elbow([(30.6, 17.3), (33.4, 17.3), (33.4, 61.0), (35.0, 61.0)], "#555")  # step 7 -> step 8
# QC comparison table
ax.add_patch(FancyBboxPatch((35.2, 10), 28.4, 25, boxstyle="round,pad=0.25,rounding_size=0.8", fc="white", ec=ORANGE, lw=1.8, ls="--"))
ax.text(49.4, 33.4, "Step 2 QC comparison (slide section)", ha="center", fontsize=10.4, fontweight="bold", color=ORANGE)
rows = [("MT cutoff", "Retained", "Removed", "MT loss %"), ("20%", "–", "–", "–"), ("15%", "–", "–", "–"), ("3–5% (nuclei)", "–", "–", "–")]
for r, row in enumerate(rows):
    for c, cell in enumerate(row):
        ax.text([37.6, 47.2, 52.6, 59.4][c], 31.0 - r * 2.2, cell, fontsize=9.2, ha="center" if c else "left",
                fontweight="bold" if r == 0 else "normal", color="#1b1b1b")
ax.text(36.3, 23.4, "Also report per cutoff:\n• effect on clustering\n• effect on endothelial identification\n  (cells and clusters retained)", fontsize=8.9, va="top", linespacing=1.4)
ax.text(36.3, 11.0, "Cutoff to use: 3–5% for snRNA-seq (Dr. Thanu's note);\n15%/20% apply to whole-cell data (e.g. GSE129363).", fontsize=8.5, color=ORANGE, va="bottom", style="italic", linespacing=1.35)
elbow([(30.6, 53.5), (32.4, 53.5), (32.4, 28.0), (35.0, 28.0)], ORANGE)  # step 2 -> QC comparison

# ---- lane C: gene panel + DE (11-15)
lane(67, 32, "C. Gene panel and differential expression (11–15)", ORANGE)
C = [("11", "Endothelial gene panel\nERG · PECAM1 (CD31) · CDH5 (VE-cad) ·\nNOS3 (eNOS) · VCAM1 · ICAM1 · HIF1A · ANGPT2", 58, 7.2),
     ("12", "ERG, T2D vs Control (cell level)\nlog2FC · adjusted p-value · significance", 49.2, 6.0),
     ("13", "Pseudobulk from the original object\nsum counts per sample → edgeR / limma\n(sample = replicate)", 39.4, 7.2),
     ("14", "Pseudobulk separately for VAT and SAT\n(distinct depot profiles)", 31.2, 5.6),
     ("15", "ERG + remaining panel genes:\ndepot-specific and disease-associated\nexpression patterns", 20.8, 7.4)]
for n, t, y, h in C:
    box(68.2, y, 29.6, h, n, t, ORANGE, fs=9.2)
for i in range(len(C) - 1):
    arrow(83, C[i][2], 83, C[i + 1][2] + C[i + 1][3] + 0.1)
elbow([(63.6, 42.6), (66.0, 42.6), (66.0, 61.6), (68.0, 61.6)], TEAL)  # step 10 -> step 11
ax.text(83, 16.4, "Report: effect size (logFC, T2D relative to Control),\nFDR < 0.05, and sample counts per group", ha="center", fontsize=8.8, style="italic", color="#444", linespacing=1.35)

# ---- footer: presentation requirements
ax.add_patch(FancyBboxPatch((1, 0.8), 98, 6.2, boxstyle="round,pad=0.25,rounding_size=0.8", fc="white", ec=NAVY, lw=1.6))
ax.text(50, 5.5, "Presentation requirements", ha="center", fontsize=10.8, fontweight="bold", color=NAVY)
ax.text(50, 3.6, "High-resolution, publication-quality figures  ·  a figure legend for every figure  ·  due end of day", ha="center", fontsize=9.2, color="#1b1b1b")
ax.text(50, 1.9, "QC section: compare MT cutoffs (cells retained / removed, % MT-associated loss, effect on clustering and endothelial identification)", ha="center", fontsize=9.2, color="#1b1b1b")

fig.savefig("docs/figures/analysis_pipeline.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("docs/figures/analysis_pipeline.svg", bbox_inches="tight", facecolor="white")
print("ok")
