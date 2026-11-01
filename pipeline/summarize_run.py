#!/usr/bin/env python3
"""Write RESULTS.md for a run: state, timing, files produced, and small previews so results can be read without opening big files."""
import json, pathlib, sys
rd = pathlib.Path(sys.argv[1]); out = rd / "out"; L = []
s = json.load(open(rd / "status.json")); L += [f"# Run {rd.name}", "", f"- state: **{s['state']}** ({s['last_completed']}/{s['total_cells']} cells)", f"- run time: {round(sum(c['seconds'] for c in s['cells']))} s; peak memory {max([c['peak_mem_mb'] or 0 for c in s['cells']] or [0])} MB"]
if s.get("failed_cell"): L.append(f"- FAILED in cell {s['failed_cell']}: see run.log (last lines below)")
L += ["", "## Slowest cells", ""] + [f"- cell {c['cell']}: {c['seconds']} s" for c in sorted(s["cells"], key=lambda c: -c["seconds"])[:5]]
L += ["", "## Files produced", ""]
for f in sorted(out.rglob("*")):
    if f.is_file(): L.append(f"- {f.relative_to(rd)} ({f.stat().st_size/1e6:.2f} MB)")
for f in sorted(out.rglob("*.csv"))[:6]:
    L += ["", f"### {f.name} (first rows)", "```"] + f.read_text(errors="replace").splitlines()[:8] + ["```"]
for f in sorted(out.rglob("*.xlsx"))[:3]:
    try:
        import openpyxl; wb = openpyxl.load_workbook(f, read_only=True); L += ["", f"### {f.name}: sheets", ""]
        for ws in wb: L.append(f"- {ws.title}: {ws.max_row-1 if ws.max_row else 0} rows")
    except Exception as e: L.append(f"(could not read {f.name}: {e})")
if s["state"] != "finished":
    L += ["", "## Last log lines", "```"] + (rd / "run.log").read_text(errors="replace").splitlines()[-25:] + ["```"]
(rd / "RESULTS.md").write_text("\n".join(L)); print("\n".join(L[:12])); print("... full summary:", rd / "RESULTS.md")
