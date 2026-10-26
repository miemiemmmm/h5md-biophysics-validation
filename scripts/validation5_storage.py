#!/usr/bin/env python
"""
Validation 5 -- what the metadata costs to store.

Measures, at the HDF5 layer, the stored size of every dataset of each reference
file and attributes it to the metadata groups or to the trajectory. Write time
is not measured here: neither process wall time nor GROMACS's 'Write traj.'
counter separates the one-off metadata write from per-frame output, so it comes
from a timer in the writer (see data/README.md).

    python scripts/validation5_storage.py

Writes <results>/storage.json: a summary, and under ``files`` one entry per
reference file with its system, protocol and byte breakdown.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


import h5py

#: Root groups that carry metadata rather than trajectory frames.
METADATA_GROUPS = (
    "topology",
    "force_field",
    "connectivity",
    "parameters",
    "selections",
    "observables",
)
TRAJECTORY_GROUPS = ("particles",)


def storage_breakdown(path):
    """Bytes of stored data per root group, from the HDF5 layer itself."""
    totals = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            root = name.split("/")[0]
            totals[root] = totals.get(root, 0) + obj.id.get_storage_size()

    with h5py.File(path, "r") as f:
        f.visititems(visit)
        n_particles = None
        for sysname in f.get("particles", {}):
            pos = f["particles"][sysname].get("position")
            if pos is not None and "value" in pos:
                n_particles = pos["value"].shape[-2]
                n_frames = pos["value"].shape[0]
                break
        else:
            n_frames = None

    metadata = sum(totals.get(g, 0) for g in METADATA_GROUPS)
    trajectory = sum(totals.get(g, 0) for g in TRAJECTORY_GROUPS)
    return {
        "per_group_bytes": totals,
        "metadata_bytes": metadata,
        "trajectory_bytes": trajectory,
        "file_bytes": Path(path).stat().st_size,
        "n_particles": n_particles,
        "n_frames": n_frames,
        "metadata_bytes_per_particle": (metadata / n_particles)
        if n_particles
        else None,
        "metadata_over_trajectory": (metadata / trajectory) if trajectory else None,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("storage.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    files = {}
    for key, entry in sorted(manifest.items()):
        h5 = Path(entry.get("h5md", ""))
        if "error" in entry or not h5.exists():
            continue
        record = {"system": entry["system"], "protocol": entry["protocol"]}
        record["file"] = str(h5)
        try:
            record.update(storage_breakdown(h5))
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        files[key] = record

    measured = [v for v in files.values() if "error" not in v]
    summary = {"n_files": len(files), "n_measured": len(measured)}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "files": files}, indent=2))
    print(f"storage measured for {len(measured)} of {len(files)} reference files")
    print(f"written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
