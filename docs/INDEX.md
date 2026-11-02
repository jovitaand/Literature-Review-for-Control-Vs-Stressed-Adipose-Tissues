# Project index (consolidated 2 November 2026)

**Start here:** `docs/Dataset_selection_report.docx` is the single consolidated report for the line manager (status at a glance, all datasets, results, review, limitations, next steps).

## Reports and written material
| File | Content |
|---|---|
| docs/Dataset_selection_report.docx | Consolidated manager report |
| docs/gse268904_verification.md | GSE268904 criteria cross-check with real counts |
| docs/gse278526_verification.md | GSE278526 / PRJNA1167724 criteria cross-check |
| docs/GSE268904_one_page_summary.pdf | One-page chart summary for non-technical readers |
| docs/literature_review_synthesis.md, docs/meeting_qa_briefing.md | Earlier literature synthesis and meeting briefing |

## Notebooks (R, Google Colab or the pipeline)
| Notebook | Purpose | Status |
|---|---|---|
| notebooks/GSE268904_colab_qc.ipynb | Extract, QC, marker-based endothelial estimate, consolidated CSV | Run on real data |
| notebooks/GSE268904_gene_panel_DE.ipynb | Eight-gene endothelial differential expression | Synthetic test only |
| notebooks/GSE268904_endothelial_DE_UMAP.ipynb | UMAP, clusters, DE, top 50 for GSE268904 | Synthetic test only |
| notebooks/GSE278526_colab.ipynb | Criteria check and summary for GSE278526 | Synthetic test only |
| notebooks/GSE129363_colab_annotated.ipynb | Annotated check and endothelial DE (marker rule) | Run on real data (Part 1) |
| notebooks/GSE129363_umap_clusters_DE.ipynb | UMAP, clusters, DE, top 50 (version 2) | v1 run on real data; v2 awaits re-run |

## Data and metadata (data/metadata)
GEO series and sample records, sample maps and cell annotation for GSE268904, GSE278526, GSE129363; `atlas_source_datasets_screen.csv` (screening of the datasets behind the integrated vascular atlas).

## Pipeline (pipeline/)
`run_notebook.py`, `runner.R`, `fetch_data.py`, `setup_env.sh`, `summarize_run.py`; see `pipeline/README.md`. Needs Google Drive and GEO file hosts to be allowed to fetch real data.

## Open items
1. Re-run GSE129363 notebook v2; resolve the 3 conflicting diabetes labels (E07, 04, 09).
2. Run the GSE268904 gene-panel / UMAP notebooks on the real archive.
3. Obtain the integrated vascular atlas matrices.
4. Run the GSE278526 notebook on the real files.
