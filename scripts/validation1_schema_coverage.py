#!/usr/bin/env python
"""
Validation 1 -- schema coverage.

Compares the specification model with the corpus of reference files, and
needs no other validation's results. Two questions:

1. ELEMENT COVERAGE.  Of the elements the specification defines, which ones does
   the corpus actually contain?  An element that is specified but never written
   by any run is not evidence that the schema works -- it is an untested branch,
   and coverage.json names it as such.

2. DIMENSION COVERAGE.  Which of the chemical and protocol dimensions of the test
   matrix are exercised at all.

Whether the written items are correct is Validation 2's question.

    python scripts/validation1_schema_coverage.py
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


import h5py
import numpy as np

from h5md_validator import spec


# --------------------------------------------------------------------------
def specified_elements():
    """Every element path the specification defines, with its obligation."""
    out = {}

    def walk_group(node: spec.Group, prefix: str):
        path = f"{prefix}/{node.name}" if prefix else node.name
        out[path] = "required" if node.required else "optional"
        for ds in node.datasets:
            out[f"{path}/{ds.name}"] = "required" if ds.required else "optional"
        for sub in node.groups:
            walk_group(sub, path)

    walk_group(spec.PARAMETERS, "")

    for name in spec.OBSERVABLES:
        out[f"observables/{name}"] = "optional"
    # S6: observables of a subset of particles, in a subgroup per selection group.
    out["observables/<group>"] = "optional"

    for name, itype in spec.INTERACTION_TYPES.items():
        layer = "bonded" if itype.bonded else "nonbonded"
        out[f"force_field/{layer}/{name}/parameters"] = "optional"
        out[f"force_field/{layer}/{name}/parameters_end_state"] = "optional"
        # Bonded lists hold parameter_id (and particle_ids); nonbonded, type_id.
        out[f"topology/<molecule>/interactions/{name}"] = "optional"

    # S8 provenance attributes on the force_field group. These are optional but
    # they are the elements that say *which* force field a trajectory used, so
    # leaving them out of the count would flatter the coverage figure.
    for attr in spec.FORCE_FIELD_ATTRS:
        out[f"force_field/@{attr.name}"] = "optional"
    out["topology/@total_particle_count"] = "required"
    out["topology/@total_residue_count"] = "required"
    out["topology/@system_name"] = "optional"

    for name in (
        "topology/molecule_blocks",
        "topology/<molecule>/particle_data",
        "topology/<molecule>/sequence",
        "topology/<molecule>/exclusions",
        "topology/<molecule>/exclusions/particle_ids",
        "topology/<molecule>/particle_name_table",
        "topology/<molecule>/component_name_table",
        "connectivity/bonds",
        "selections/<group>/ranges",
        "selections/<group>/indices",
        "selections/<group>/selection",
        "selections/<group>/software",
    ):
        out[name] = "optional"
    return out


def observed_values(h5_path):
    """Attribute and scalar values per element path, for the constancy check.

    An element that is present in every file but carries the same value in all
    of them is exercised in form only. A writer that emits a fixed placeholder
    -- 'not declared', say -- would otherwise count as coverage, which would
    flatter the figure without demonstrating any capability.
    """
    values = {}

    def record(path, value):
        if isinstance(value, bytes):
            value = value.decode("utf-8", "replace")
        elif isinstance(value, np.ndarray):
            if value.size > 8:
                return
            value = tuple(value.ravel().tolist())
        elif isinstance(value, np.generic):
            value = value.item()
        values[path] = repr(value)[:120]

    with h5py.File(h5_path, "r") as f:

        def visit(name, obj):
            parts = name.split("/")
            if (
                parts[0] == "topology"
                and len(parts) > 1
                and parts[1] != "molecule_blocks"
            ):
                parts[1] = "<molecule>"
            elif (parts[0] == "selections" and len(parts) > 1) or (
                parts[0] == "observables"
                and len(parts) > 2
                and parts[1] not in spec.OBSERVABLES
            ):
                parts[1] = "<group>"
            elif parts[0] == "particles":
                return
            canonical = "/".join(parts)
            for attr, value in obj.attrs.items():
                record(f"{canonical}/@{attr}", value)
            if (
                isinstance(obj, h5py.Dataset)
                and obj.shape == ()
                and obj.dtype.names is None
            ):
                record(canonical, obj[()])

        f.visititems(visit)
        for attr, value in f.attrs.items():
            record(f"@{attr}", value)
    return values


def observed_elements(h5_path):
    """Element paths actually present in one file, in the same vocabulary."""
    seen = set()
    with h5py.File(h5_path, "r") as f:

        def visit(name, obj):
            for attr in obj.attrs:
                parts = name.split("/")
                if (
                    parts[0] == "topology"
                    and len(parts) > 1
                    and parts[1] != "molecule_blocks"
                ):
                    parts[1] = "<molecule>"
                elif parts[0] == "selections" and len(parts) > 1:
                    parts[1] = "<group>"
                seen.add("/".join(parts) + f"/@{attr}")
            parts = name.split("/")
            if (
                parts[0] == "topology"
                and len(parts) > 1
                and parts[1] != "molecule_blocks"
            ):
                parts[1] = "<molecule>"
            elif (parts[0] == "selections" and len(parts) > 1) or (
                parts[0] == "observables"
                and len(parts) > 2
                and parts[1] not in spec.OBSERVABLES
            ):
                parts[1] = "<group>"
            elif parts[0] == "particles":
                return
            seen.add("/".join(parts))

        f.visititems(visit)
        for root in f:
            seen.add(root)
    return seen


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("coverage.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    ok = {
        k: v
        for k, v in manifest.items()
        if "error" not in v and Path(v.get("h5md", "")).exists()
    }

    # ---- 1. element coverage -------------------------------------------
    specified = specified_elements()
    observed = set()
    per_file_elements = {}
    distinct_values = defaultdict(set)
    for key, rec in ok.items():
        path = rec["h5md"]
        try:
            elems = observed_elements(path)
            vals = observed_values(path)
        except OSError:
            continue
        per_file_elements[key] = elems
        observed |= elems
        for epath, value in vals.items():
            distinct_values[epath].add(value)

    element_rows = []
    for path, obligation in sorted(specified.items()):
        present = path in observed
        n_files = sum(1 for e in per_file_elements.values() if path in e)
        seen_values = distinct_values.get(path)
        row = {
            "element": path,
            "obligation": obligation,
            "observed": present,
            "n_files": n_files,
        }
        if present and seen_values is not None:
            row["n_distinct_values"] = len(seen_values)
            if len(seen_values) == 1 and n_files > 1:
                row["constant_value"] = next(iter(seen_values))
        element_rows.append(row)
    never = [r for r in element_rows if not r["observed"]]
    constant = [r for r in element_rows if r.get("constant_value") is not None]

    # ---- 2. dimension coverage ------------------------------------------
    dims = defaultdict(list)
    for rec in ok.values():
        for d in rec.get("dimensions", []):
            dims[d].append(rec["system"])
    exercises = defaultdict(list)
    for rec in ok.values():
        for e in rec.get("exercises", []):
            exercises[e].append(rec["protocol"])

    not_generated = {k: v["error"] for k, v in manifest.items() if "error" in v}

    summary = {
        "n_reference_files": len(ok),
        "element_coverage": {
            "n_specified": len(specified),
            "n_observed": sum(1 for r in element_rows if r["observed"]),
            "n_never_observed": len(never),
            "never_observed": [r["element"] for r in never],
            "n_constant": len(constant),
            "constant": [
                {
                    "element": r["element"],
                    "value": r["constant_value"],
                    "n_files": r["n_files"],
                }
                for r in constant
            ],
            "rows": element_rows,
        },
        "dimensions": {k: sorted(set(v)) for k, v in sorted(dims.items())},
        "protocol_targets": {k: sorted(set(v)) for k, v in sorted(exercises.items())},
        "not_generated": not_generated,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))

    ec = summary["element_coverage"]
    print(f"reference files              : {summary['n_reference_files']}")
    print(f"specified elements           : {ec['n_specified']}")
    print(f"  observed in the corpus     : {ec['n_observed']}")
    print(f"  never observed             : {ec['n_never_observed']}")
    print(f"chemical dimensions covered  : {len(summary['dimensions'])}")
    if ec["never_observed"]:
        print("\nspecified but never written by any run:")
        for e in ec["never_observed"]:
            print(f"  {e}")
    if ec["constant"]:
        print(
            f"\npresent but carrying one value in every file "
            f"({ec['n_constant']} element(s)) -- exercised in form, not in "
            "substance:"
        )
        for r in ec["constant"]:
            print(f"  {r['element']:52s} = {r['value']}")
    print(f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
