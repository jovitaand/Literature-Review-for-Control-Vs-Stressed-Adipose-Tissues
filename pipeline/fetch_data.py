#!/usr/bin/env python3
"""Download the data for a dataset listed in pipeline/datasets.json (resumable). Tells you exactly which host is blocked, if any.
   python3 pipeline/fetch_data.py GSE129363 [--dest DIR]        # tries Google Drive (by ID), then any --url given
   python3 pipeline/fetch_data.py --check-network
Large files go to data/raw/ (git-ignored). A file already present with the right size is skipped."""
import argparse, json, os, pathlib, subprocess, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
HOSTS = {"Google Drive": "https://drive.usercontent.google.com", "Google Drive (old)": "https://drive.google.com", "GEO file server": "https://ftp.ncbi.nlm.nih.gov",
         "GEO web": "https://www.ncbi.nlm.nih.gov", "CRAN": "https://cloud.r-project.org", "GitHub": "https://raw.githubusercontent.com", "PyPI": "https://pypi.org"}
def probe(url):
    r = subprocess.run(["curl", "-s", "-o", "/dev/null", "-m", "15", "-I", "-w", "%{http_code}", url], capture_output=True, text=True); return r.stdout.strip()
def check_network():
    for k, u in HOSTS.items():
        c = probe(u); print(f"{k:20s} {u:45s} {'reachable' if c not in ('000', '') else 'BLOCKED by the environment network policy'} ({c})")
def fetch(name, dest, urls):
    cfg = json.load(open(ROOT / "pipeline" / "datasets.json"))[name]; d = pathlib.Path(dest or ROOT / cfg["dest"]); d.mkdir(parents=True, exist_ok=True); bad = 0
    for f in cfg["files"]:
        p = d / f["name"]
        if p.exists() and (f["bytes"] is None or p.stat().st_size == f["bytes"]): print("have", p.name); continue
        cands = ([f"https://drive.usercontent.google.com/download?id={f['drive_id']}&export=download&confirm=t"] if f.get("drive_id") else []) + urls.get(f["name"], [])
        for u in cands:
            print("fetching", p.name, "from", u.split("?")[0][:60])
            subprocess.run(["curl", "-L", "-C", "-", "--fail", "-m", "3600", "-o", str(p), u], capture_output=True)
            if p.exists() and p.stat().st_size > 0 and (f["bytes"] is None or p.stat().st_size == f["bytes"]) and not open(p, "rb").read(15).lower().startswith((b"<!doctype", b"<html")): print("  ok", p.stat().st_size, "bytes"); break
            if p.exists(): p.unlink()
        else:
            bad += 1; print(f"  FAILED: {p.name}. Allow the hosts in the environment's network settings (see --check-network) or copy the file into {d}")
    return bad
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); ap.add_argument("dataset", nargs="?"); ap.add_argument("--dest")
    ap.add_argument("--url", action="append", default=[], metavar="NAME=URL"); ap.add_argument("--check-network", action="store_true"); a = ap.parse_args()
    if a.check_network or not a.dataset: check_network(); sys.exit(0)
    urls = {}
    for kv in a.url: k, v = kv.split("=", 1); urls.setdefault(k, []).append(v)
    sys.exit(1 if fetch(a.dataset, a.dest, urls) else 0)
