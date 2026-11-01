# Runs Python notebook cells in one shared namespace (no checkpoints). Usage: runner.py RUN_DIR START STOP_AFTER CELL_TIMEOUT
import sys, json, time, glob, traceback, resource
rd, start, stop_after = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
cells = sorted(glob.glob(f"{rd}/cells/cell_*.py")); ns = {"__name__": "__main__"}; st = {"total_cells": len(cells), "last_completed": start - 1, "state": "running", "failed_cell": None, "cells": []}
def save(): json.dump(st, open(f"{rd}/status.json", "w"))
for i in range(start, min(len(cells), stop_after) + 1):
    print(f"\n[runner] ===== cell {i} of {len(cells)} =====", flush=True); t0 = time.time(); ok = True
    try: exec(compile(open(cells[i - 1]).read(), cells[i - 1], "exec"), ns)
    except Exception: traceback.print_exc(); ok = False
    st["cells"].append({"cell": i, "seconds": round(time.time() - t0, 1), "peak_mem_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, "ok": ok})
    if not ok: st["state"], st["failed_cell"] = "failed", i; save(); sys.exit(1)
    st["last_completed"] = i; save()
st["state"] = "finished" if st["last_completed"] >= len(cells) else "paused"; save()
