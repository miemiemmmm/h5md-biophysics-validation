#!/usr/bin/env python
"""
Where could force-field provenance come from?

S8 defines thirteen attributes naming the force field a trajectory used, and no
reference file carries any of them (deviations.DEV-FF-PROVENANCE-UNWRITTEN).
Before deciding that the user must supply them by hand, this script measures the
two automatic alternatives.

1. INFER IT FROM THE PARAMETERS the file already holds.  If files built with the
   same force field shared a recognisable parameter signature, provenance would
   need no declaration at all.  Measured as the overlap between the distinct
   bonded parameter tuples of each pair of systems.

2. READ IT FROM THE TOPOLOGY that grompp already parses.  A GROMACS topology
   names its force field in an #include path, which grompp resolves and then
   discards.  Measured as the fraction of systems whose force field and water
   model are recoverable from those paths.

    python scripts/provenance_recoverability.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import re
import sys
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import RESULTS_DIR



def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)



import h5py

from config.matrix import DATA_ROOT

ROOT = Path(DATA_ROOT)

#: the force field each system was actually built with, from systems.py
_FF_FIELD = "force_field"

_INCLUDE = re.compile(r'^#include\s+"([^"]+)"', re.M)
_FF_PATH = re.compile(r"([A-Za-z0-9_.+-]+)\.ff/forcefield\.itp")
_WATER = re.compile(r"\.ff/(tip[0-9a-z]*|spce?|opc3?|tip4pew)\.itp")
_IONS = re.compile(r"\.ff/ions?(?:_[a-z0-9]+)?\.itp")


def parameter_signature(h5_path, kind="bond"):
    """The set of distinct parameter tuples of one interaction type."""
    with h5py.File(h5_path, "r") as f:
        for layer in ("bonded", "nonbonded"):
            ds = f.get(f"force_field/{layer}/{kind}/parameters")
            if ds is not None:
                data = ds[()]
                return {tuple(round(float(data[i][n]), 6) for n in data.dtype.names)
                        for i in range(data.shape[0])}
    return set()


def from_topology(top_path):
    """What the topology's #include paths reveal."""
    text = Path(top_path).read_text(errors="replace")
    joined = " ".join(_INCLUDE.findall(text))
    ff = _FF_PATH.search(joined)
    water = _WATER.search(joined)
    return {"force_field": ff.group(1) if ff else None,
            "water_model": water.group(1) if water else None,
            "ion_parameters": bool(_IONS.search(joined))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=_r("provenance_recoverability.json"))
    args = ap.parse_args()

    from config.systems import SYSTEMS_BY_ID

    # ---- 1. can the parameters identify the force field? -------------------
    signatures, declared = {}, {}
    for sid, system in SYSTEMS_BY_ID.items():
        f = ROOT / "reference_files" / f"{sid}__nvt-vrescale.h5md"
        if not f.exists():
            continue
        try:
            sig = parameter_signature(f)
        except OSError:
            continue
        if sig:
            signatures[sid] = sig
            declared[sid] = getattr(system, _FF_FIELD, "")

    pairs = []
    for a, b in itertools.combinations(sorted(signatures), 2):
        sa, sb = signatures[a], signatures[b]
        overlap = len(sa & sb) / min(len(sa), len(sb))
        pairs.append({"a": a, "b": b, "same_force_field": declared[a] == declared[b],
                      "overlap": overlap})
    same = [p["overlap"] for p in pairs if p["same_force_field"]]
    diff = [p["overlap"] for p in pairs if not p["same_force_field"]]
    worst_same = min(same) if same else None
    best_diff = max(diff) if diff else None
    separable = (worst_same is not None and best_diff is not None
                 and worst_same > best_diff)

    # ---- 2. can the topology reveal it? ------------------------------------
    topo = {}
    for sid, system in SYSTEMS_BY_ID.items():
        # Read the shipped topology directly: the harness work cache is
        # transient, and validation_systems/ is where the systems actually live.
        top = Path(system.location) / system.topology
        if top.is_file():
            topo[sid] = from_topology(top)
    n_ff = sum(1 for v in topo.values() if v["force_field"])
    n_water = sum(1 for v in topo.values() if v["water_model"])

    result = {
        "inference_from_parameters": {
            "n_pairs": len(pairs),
            "worst_same_force_field_overlap": worst_same,
            "best_different_force_field_overlap": best_diff,
            "separable": separable,
            "pairs": sorted(pairs, key=lambda p: -p["overlap"]),
        },
        "recovery_from_topology": {
            "n_systems": len(topo),
            "n_force_field_recoverable": n_ff,
            "n_water_model_recoverable": n_water,
            "per_system": topo,
        },
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))

    print("1. inferring the force field from the parameters in the file")
    print(f"   {len(pairs)} system pairs compared on their bonded parameter sets")
    print(f"   lowest overlap between two systems sharing a force field : "
          f"{worst_same:.1%}")
    print(f"   highest overlap between two systems that do not          : "
          f"{best_diff:.1%}")
    print(f"   separable by this signature: {'yes' if separable else 'NO'}")
    if not separable:
        print("   overlap tracks what is in the box, not which force field "
              "parametrised it, so provenance cannot be inferred from content")
    print()
    print("2. reading it from the topology grompp already parses")
    print(f"   force field named in an #include path : {n_ff}/{len(topo)} systems")
    print(f"   water model named likewise            : {n_water}/{len(topo)} systems")
    missing = [s for s, v in sorted(topo.items()) if not v["force_field"]]
    if missing:
        print(f"   not derivable for: {', '.join(missing)}")
        print("   (externally assembled systems, whose includes name local files)")
    print(f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
