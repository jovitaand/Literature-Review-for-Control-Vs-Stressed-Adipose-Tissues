#!/usr/bin/env python3
"""Run a Jupyter notebook (R or Python kernel) from the command line, of any length, with logs, checkpoints and a results summary.

Examples
  python3 pipeline/run_notebook.py notebooks/GSE129363_umap_clusters_DE.ipynb --set DATA_DIR=/data/gse129363 --set USE_GDOWN=FALSE --background
  python3 pipeline/run_notebook.py --resume pipeline/runs/<run-id>            # continue after a failure or pause
  python3 pipeline/run_notebook.py --status pipeline/runs/<run-id>
"""
import argparse, json, os, re, subprocess, sys, time, pathlib

HERE = pathlib.Path(__file__).resolve().parent

def load_cells(nb_path):
    nb = json.load(open(nb_path))
    lang = (nb.get("metadata", {}).get("kernelspec", {}).get("language") or nb.get("metadata", {}).get("language_info", {}).get("name") or "R").lower()
    cells = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
    return lang, cells

def fmt(v):
    if re.fullmatch(r"-?[0-9.]+|TRUE|FALSE|NULL|NA|True|False|None", v): return v
    return '"%s"' % v.strip('"')

def patch(src, sets, lang):
    """Override top-level settings such as `DATA_DIR <- "..."` (R) or `DATA_DIR = "..."` (Python)."""
    op = r"<-" if lang == "r" else r"="
    for k, v in sets.items():
        src = re.sub(r"(?m)(^|;\s*)(%s\s*%s\s*)[^;#\n]*" % (re.escape(k), op), lambda m: f"{m.group(1)}{m.group(2)}{fmt(v)}" , src)
    return src

def status(rd):
    f = pathlib.Path(rd) / "status.json"
    if not f.exists(): print("no status yet"); return None
    s = json.load(open(f)); print(json.dumps({k: v for k, v in s.items() if k != "cells"}), "| cells run:", len(s["cells"]),
          "| seconds:", round(sum(c["seconds"] for c in s["cells"])), "| peak MB:", max([c["peak_mem_mb"] or 0 for c in s["cells"]] or [0])); return s

def launch(rd, start, stop_after, ckpt, tmo, lang, background, skip_install):
    log = open(rd / "run.log", "a")
    if lang == "r":
        env = dict(os.environ)
        if skip_install: env["R_PROFILE_USER"] = str(HERE / "no_install.R")
        cmd = ["Rscript", str(HERE / "runner.R"), str(rd), str(start), str(stop_after), ckpt, str(tmo)]
    else:
        cmd = [sys.executable, str(HERE / "runner.py"), str(rd), str(start), str(stop_after), str(tmo)]; env = dict(os.environ)
    if background:
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        (rd / "pid").write_text(str(p.pid)); print(f"started in background: pid {p.pid}\nlog: {rd/'run.log'}\nstatus: python3 pipeline/run_notebook.py --status {rd}"); return 0
    rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT, env=env); return rc

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("notebook", nargs="?"); ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a top-level setting in the notebook")
    ap.add_argument("--run-dir"); ap.add_argument("--resume", metavar="RUN_DIR"); ap.add_argument("--status", metavar="RUN_DIR")
    ap.add_argument("--from-cell", type=int); ap.add_argument("--stop-after", type=int, default=10**6)
    ap.add_argument("--checkpoint", default="each", help="none | each | N (every N cells)"); ap.add_argument("--cell-timeout", type=float, default=7200)
    ap.add_argument("--background", action="store_true"); ap.add_argument("--no-skip-install", action="store_true", help="allow install.packages() (needs CRAN access)")
    a = ap.parse_args()
    if a.status: status(a.status); return 0
    if a.resume:
        rd = pathlib.Path(a.resume); s = json.load(open(rd / "status.json")); meta = json.load(open(rd / "meta.json"))
        start = a.from_cell or (s["last_completed"] + 1); print(f"resuming {rd} at cell {start}")
        return launch(rd, start, a.stop_after, a.checkpoint, a.cell_timeout, meta["lang"], a.background, not a.no_skip_install)
    if not a.notebook: ap.error("give a notebook or --resume/--status")
    lang, cells = load_cells(a.notebook); sets = dict(s.split("=", 1) for s in a.set)
    rd = pathlib.Path(a.run_dir or HERE / "runs" / f"{pathlib.Path(a.notebook).stem}-{time.strftime('%Y%m%d-%H%M%S')}"); rd.mkdir(parents=True, exist_ok=True)
    (rd / "cells").mkdir(exist_ok=True); (rd / "out").mkdir(exist_ok=True)
    sets.setdefault("OUT_DIR", str(rd / "out"))
    for i, c in enumerate(cells, 1): (rd / "cells" / f"cell_{i:03d}.{'R' if lang == 'r' else 'py'}").write_text(patch(c, sets, lang))
    json.dump({"notebook": a.notebook, "lang": lang, "sets": sets, "cells": len(cells), "started": time.ctime()}, open(rd / "meta.json", "w"), indent=1)
    print(f"{a.notebook}: {len(cells)} code cells ({lang}); run dir {rd}")
    return launch(rd, a.from_cell or 1, a.stop_after, a.checkpoint, a.cell_timeout, lang, a.background, not a.no_skip_install)

if __name__ == "__main__": sys.exit(main())
