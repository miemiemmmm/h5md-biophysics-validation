#!/usr/bin/env python
"""
Run the conformance validator and the round-trip test over every reference file.

    python scripts/run_validation.py --manifest results/reference_files.json

Writes results/validation.json holding, per reference file: the validator's
findings (Validation 3), the item-by-item round-trip comparison (Validation 2),
and the metadata storage breakdown (Validation 5).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import RESULTS_DIR



def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)



import h5py

from h5md_validator import Validator
from h5md_validator import gmxdump, roundtrip
from config.matrix import DATA_ROOT, GMX

DUMPS = Path(DATA_ROOT) / "dumps"

#: Root groups that carry metadata rather than trajectory frames.  Used for the
#: storage-overhead measurement of Validation 5.
METADATA_GROUPS = ("topology", "force_field", "connectivity", "parameters",
                   "selections", "observables")
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
        "metadata_bytes_per_particle": (metadata / n_particles) if n_particles else None,
        "metadata_over_trajectory": (metadata / trajectory) if trajectory else None,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("validation.json"))
    ap.add_argument("--skip-roundtrip", action="store_true")
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    validator = Validator()
    results = {}

    for key, entry in sorted(manifest.items()):
        if "error" in entry:
            results[key] = {"status": "not_generated", "error": entry["error"],
                            "system": entry["system"], "protocol": entry["protocol"]}
            print(f"  --  {key}: not generated")
            continue
        h5 = Path(entry["h5md"])
        if not h5.exists():
            results[key] = {"status": "missing_file", "system": entry["system"],
                            "protocol": entry["protocol"]}
            continue

        record = {"status": "ok", "system": entry["system"],
                  "protocol": entry["protocol"],
                  "dimensions": entry.get("dimensions", []),
                  "exercises": entry.get("exercises", []),
                  "wall_seconds": entry.get("wall_seconds")}

        report = validator.validate(h5)
        record["validation"] = report.to_dict()

        if not args.skip_roundtrip and Path(entry["tpr"]).exists():
            try:
                dump = gmxdump.from_tpr(entry["tpr"], gmx=GMX, cache_dir=DUMPS)
                rt = roundtrip.roundtrip(dump, h5, entry["tpr"])
                record["roundtrip"] = rt.to_dict()
            except Exception as exc:
                record["roundtrip_error"] = f"{type(exc).__name__}: {exc}"

        try:
            record["storage"] = storage_breakdown(h5)
        except Exception as exc:
            record["storage_error"] = f"{type(exc).__name__}: {exc}"

        results[key] = record
        v = record["validation"]
        r = record.get("roundtrip", {})
        print(f"  {'PASS' if v['conformant'] else 'FAIL'}  {key:44s} "
              f"err={v['counts']['ERROR']:<3d} dev={v['counts']['DEVIATION']:<3d} "
              f"rt: match={r.get('n_match', '-')} mismatch={r.get('n_mismatch', '-')} "
              f"missing={r.get('n_missing_in_file', '-')} "
              f"dropped={r.get('n_dropped_terms', '-')}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))

    ok = [r for r in results.values() if r.get("status") == "ok"]
    conformant = [r for r in ok if r["validation"]["conformant"]]
    faithful = [r for r in ok if r.get("roundtrip", {}).get("faithful")]
    complete = [r for r in ok if r.get("roundtrip", {}).get("all_terms_written")]
    print("\n" + "=" * 72)
    print(f"reference files validated : {len(ok)}")
    print(f"conformant (0 errors)     : {len(conformant)}")
    compared = sum(1 for r in ok if "roundtrip" in r)
    print(f"round-trip faithful       : {len(faithful)} of {compared} compared")
    print(f"every source term written  : {len(complete)} of {compared} compared")
    print(f"written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
