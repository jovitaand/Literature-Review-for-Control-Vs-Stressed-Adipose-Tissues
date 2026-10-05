#!/usr/bin/env python3
"""Download every file in data/manifest.csv into data/raw/<dataset_id>/.

Usage:
    python scripts/download_manifest.py            # download all
    python scripts/download_manifest.py --dry-run  # list what would be fetched
    python scripts/download_manifest.py GSE249089  # only one dataset

Resumes partial files (curl -C -), skips files already at roughly the manifest size,
and writes a log to data/raw/download_log.tsv. Needs curl and outbound access to
www.ncbi.nlm.nih.gov and figshare.com.
"""
import csv, subprocess, sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
manifest = root / "data" / "manifest.csv"
out_root = root / "data" / "raw"
dry = "--dry-run" in sys.argv
only = {a for a in sys.argv[1:] if not a.startswith("--")}

rows = list(csv.DictReader(manifest.open(encoding="utf-8")))
rows = [r for r in rows if r["dataset_id"] and (not only or r["dataset_id"] in only)]
total = sum(float(r["approx_size_mb"]) for r in rows)
print(f"{len(rows)} files, about {total:,.0f} MB")

log = []
for r in rows:
    d = out_root / r["dataset_id"]
    dest = d / r["file_name"]
    size_mb = float(r["approx_size_mb"])
    if dest.exists() and dest.stat().st_size >= 0.95 * size_mb * 1e6:
        print(f"skip   {dest.name} (already present)")
        log.append((r["dataset_id"], r["file_name"], "skipped", dest.stat().st_size))
        continue
    print(f"{'would get' if dry else 'get'}   {r['dataset_id']}/{r['file_name']}  ~{size_mb} MB")
    if dry:
        continue
    d.mkdir(parents=True, exist_ok=True)
    rc = subprocess.call(["curl", "-L", "--fail", "-C", "-", "--retry", "3", "-sS",
                          "-o", str(dest), r["download_url"]])
    status = "ok" if rc == 0 else f"curl_exit_{rc}"
    got = dest.stat().st_size if dest.exists() else 0
    print(f"       {status}, {got/1e6:.1f} MB")
    log.append((r["dataset_id"], r["file_name"], status, got))

if not dry:
    out_root.mkdir(parents=True, exist_ok=True)
    with (out_root / "download_log.tsv").open("w") as f:
        f.write("dataset\tfile\tstatus\tbytes\n")
        for row in log:
            f.write("\t".join(map(str, row)) + "\n")
