# Notebook pipeline: run notebooks of any length inside Claude's container

Claude's cloud container has **4 CPUs, 15 GB RAM and about 28 GB free disk**, and R (with Seurat, edgeR, data.table) plus Python.
It can run any notebook in `notebooks/`, however long, as a background job with logs and checkpoints.
**The one hard limit is getting the data in** (see "Data access").

## Quick start (what Claude does)
```bash
bash pipeline/setup_env.sh                                  # packages from apt (CRAN is blocked)
python3 pipeline/fetch_data.py GSE129363                    # data -> data/raw/GSE129363 (needs Drive access, see below)
python3 pipeline/run_notebook.py notebooks/GSE129363_umap_clusters_DE.ipynb \
    --set DATA_DIR=data/raw/GSE129363 --set USE_GDOWN=FALSE --set USE_HARMONY=FALSE --background
python3 pipeline/run_notebook.py --status pipeline/runs/<run-id>     # progress, per-cell timing, peak memory
python3 pipeline/summarize_run.py pipeline/runs/<run-id>             # RESULTS.md with previews of the outputs
python3 pipeline/run_notebook.py --resume pipeline/runs/<run-id>     # continue after a failure or a pause
```

## How it works
- `run_notebook.py` turns code cells into `cells/cell_NNN.R` (or `.py`), applies `--set KEY=VALUE` overrides to top-level settings (`DATA_DIR <- ...`), and starts `runner.R`.
- `runner.R` runs the cells **in one R session**, writes `run.log`, `status.json` (state, per-cell seconds, peak memory) and, after every cell, a **checkpoint** of the whole workspace. A failed cell stops the run; fix the cause and `--resume` continues from the failed cell without redoing the earlier ones.
- `--background` detaches the job, so a notebook can run for hours; Claude polls `--status`. `--checkpoint none` saves disk for huge objects; `--cell-timeout` guards a stuck cell; `--stop-after N` runs part of a notebook.
- `install.packages()` is skipped (CRAN is blocked); packages come from `setup_env.sh`. `harmony` is not available there, so use `--set USE_HARMONY=FALSE` (the notebook falls back to plain PCA).
- Outputs go to `pipeline/runs/<run-id>/out/`; `summarize_run.py` writes a short `RESULTS.md` (file list, CSV heads, Excel sheet names) so results can be read without opening large files. Small results can be committed; `data/raw/` and `pipeline/runs/` are git-ignored.

## Data access: the real bottleneck
Checked from inside the container:

| Host | Status |
|---|---|
| Google Drive (`drive.usercontent.google.com`, `drive.google.com`) | **blocked** |
| GEO file server (`ftp.ncbi.nlm.nih.gov`) | **blocked** |
| CRAN | **blocked** |
| GEO web pages, GitHub, PyPI | reachable |

The Google Drive connector can only hand over files under 10 MB, so the 38 MB matrix and 81 MB table of GSE129363, and the 1 GB GSE268904 tar, cannot come through it. Pick one:
1. **Allow the hosts** in the cloud environment's network settings (environment menu in the session title bar, Edit, Network access): add `drive.usercontent.google.com`, `drive.google.com`, `googleusercontent.com` (for Drive) and `ftp.ncbi.nlm.nih.gov` (for GEO). Then `fetch_data.py` works for every dataset. Files must be shared as "Anyone with the link".
2. **Put the files in the repository** (git accepts files up to 100 MB; the GSE129363 files fit). Place them in `data/raw/GSE129363/` and push; Claude pulls them. Not possible for the 1 GB tar.
3. **Run in Colab** and send back the small result files (what has been done so far).

Run `python3 pipeline/fetch_data.py --check-network` any time to see what is reachable.

## Limits to keep in mind
- 15 GB RAM: a full Seurat run on 130,000 nuclei is tight; the GSE268904 notebook pre-selects endothelial nuclei per sample to stay small.
- The container is ephemeral: commit or download the results you need before the session ends.
- Python notebooks run in one process without checkpoints (resume restarts from the chosen cell with a fresh state).
