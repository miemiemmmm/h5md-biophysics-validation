#!/usr/bin/env python
"""
Validation 5 (second half) -- what does carrying the metadata cost to write?

The question is not whether H5MD is faster than XTC. It is:

    how much does embedding the full metadata add to the cost of writing a
    trajectory that would have been written anyway?

Measured by running the *same* .tpr twice, differing only in the output format:

    -o run.h5md   H5MD with the complete biophysics metadata
    -o run.trr    the same trajectory frames, no metadata at all

with the step count held at zero so that the run writes its metadata and a
single frame and stops.  What remains between the two is the metadata.  Each
configuration is repeated and the median taken, because a single mdrun
wall-clock reading is dominated by start-up noise.

    python scripts/validation5_write_overhead.py --repeats 5
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import RESULTS_DIR



def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)



from config.matrix import DATA_ROOT, GMX

ROOT = Path(DATA_ROOT)


def time_run(tpr: Path, out: Path, work: Path, ntomp: int) -> float:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    out.unlink(missing_ok=True)
    t0 = time.perf_counter()
    proc = subprocess.run(
        [GMX, "-quiet", "-nobackup", "mdrun", "-s", str(tpr), "-o", str(out),
         "-nsteps", "0", "-ntmpi", "1", "-ntomp", str(ntomp),
         "-g", str(work / "md.log"), "-e", str(work / "md.edr"),
         "-c", str(work / "out.gro"), "-cpo", str(work / "state.cpt")],
        cwd=work, capture_output=True, text=True, timeout=3600)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"mdrun failed:\n{proc.stderr[-2000:]}")
    return elapsed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--protocol", default="nvt-vrescale",
                    help="compare one protocol across systems")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--ntomp", type=int, default=16)
    ap.add_argument("--out", default=_r("write_overhead.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    targets = {k: v for k, v in manifest.items()
               if "error" not in v and v["protocol"] == args.protocol}

    work = ROOT / "work" / "overhead"
    scratch = ROOT / "overhead"
    scratch.mkdir(parents=True, exist_ok=True)
    rows = []

    for key, entry in sorted(targets.items()):
        tpr = Path(entry["tpr"])
        if not tpr.exists():
            continue
        print(f"--- {key}")
        try:
            h5_times, trr_times = [], []
            h5_out = scratch / "probe.h5md"
            trr_out = scratch / "probe.trr"
            # Interleave the two formats so that any drift in machine load
            # affects both arms equally.
            for i in range(args.repeats):
                h5_times.append(time_run(tpr, h5_out, work, args.ntomp))
                trr_times.append(time_run(tpr, trr_out, work, args.ntomp))
            h5_size = h5_out.stat().st_size
            trr_size = trr_out.stat().st_size
            th5 = statistics.median(h5_times)
            ttrr = statistics.median(trr_times)
            rows.append({
                "key": key, "system": entry["system"],
                "n_particles": _n_particles(h5_out),
                "median_h5md_seconds": th5,
                "median_trr_seconds": ttrr,
                "delta_seconds": th5 - ttrr,
                "relative_write_overhead": (th5 - ttrr) / ttrr if ttrr else None,
                "h5md_bytes": h5_size, "trr_bytes": trr_size,
                "h5md_seconds": h5_times, "trr_seconds": trr_times,
            })
            r = rows[-1]
            print(f"    h5md {th5:.2f}s  trr {ttrr:.2f}s  "
                  f"delta {th5 - ttrr:+.2f}s  "
                  f"({100 * r['relative_write_overhead']:+.1f}%)")
        except Exception as exc:
            print(f"    FAILED: {str(exc)[:300]}")
        finally:
            shutil.rmtree(work, ignore_errors=True)

    shutil.rmtree(scratch, ignore_errors=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"protocol": args.protocol,
                               "repeats": args.repeats, "rows": rows}, indent=2))
    print(f"\nwritten to {out}")
    return 0


def _n_particles(h5_path):
    import h5py
    try:
        with h5py.File(h5_path, "r") as f:
            for sysname in f.get("particles", {}):
                pos = f["particles"][sysname].get("position")
                if pos is not None and "value" in pos:
                    return int(pos["value"].shape[-2])
    except OSError:
        pass
    return None


if __name__ == "__main__":
    raise SystemExit(main())
