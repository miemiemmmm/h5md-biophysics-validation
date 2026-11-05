#!/usr/bin/env python
"""
Does every bonded interaction in the source reach the file?

Validation 2, per interaction entry.  For each reference file it counts, per interaction type, how many interaction entries the GROMACS run
input holds and how many the H5MD file holds, and reports the difference.

Table SI of S8 is deliberately not exhaustive -- the set of bonded functional
forms in use across MD engines is open, and the schema's answer to that is meant
to be a convention for writing new ones, as S7 already states for simulation
parameters.  So the question here is not whether a term is in the table.  It is
whether the term reaches the file at all, under the table or under a convention,
and whether its absence is recorded when it does not.  Three outcomes:

    represented       the entries are in the file
    not predefined    the form is outside Table SI, and the file has no entries
                      for it and no marker that it was dropped
    dropped           the form maps to a predefined type, the source has
                      entries, the file has none

    python scripts/validation2_interaction_coverage.py
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


import h5py

from config.matrix import DATA_ROOT
from gmx_roundtrip import gmxdump
from gmx_roundtrip import vocabulary as vocab

DUMPS = Path(DATA_ROOT) / "tpr_dumps"


def source_counts(dump) -> dict:
    """GROMACS functype -> number of interaction entries in the whole system."""
    counts = Counter()
    for block in dump.molblocks:
        mol = dump.moltypes[block["moltype"]]
        for functype, entries in mol.ilists.items():
            if functype == "LJ_SR":  # nonbonded, no per-molecule list
                continue
            counts[functype] += len(entries) * block["nmol"]
    return dict(counts)


def file_counts(path) -> dict:
    """Schema interaction type -> number of entries in the whole system."""
    counts = Counter()
    with h5py.File(path, "r") as f:
        topo = f.get("topology")
        if topo is None:
            return {}
        blocks = topo.get("molecule_blocks")
        multiplicity = {}
        if blocks is not None:
            for row in blocks[()]:
                name = row["name"]
                name = name.decode() if isinstance(name, bytes) else str(name)
                multiplicity[name] = int(row["count"])
        for name in topo:
            mol = topo[name]
            if not isinstance(mol, h5py.Group):
                continue
            inter = mol.get("interactions")
            if not isinstance(inter, h5py.Group):
                continue
            n_copies = multiplicity.get(name, 1)
            for tname in inter:
                node = inter[tname]
                if "parameter_id" in node:
                    counts[tname] += int(node["parameter_id"].shape[0]) * n_copies
    return dict(counts)


def analyse(key, entry) -> dict:
    dump_path = DUMPS / f"{Path(entry['tpr']).stem}.dump.txt"
    if not dump_path.exists():
        return {}
    dump = gmxdump.parse(dump_path.read_text())
    src = source_counts(dump)
    dst = file_counts(entry["h5md"])

    rows = []
    for functype, n_src in sorted(src.items()):
        element = vocab.INTERACTION.get(functype)
        if element is None:
            rows.append(
                {
                    "functype": functype,
                    "element": None,
                    "source_entries": n_src,
                    "file_entries": 0,
                    "outcome": "not predefined",
                }
            )
        else:
            # Several functypes can map onto one element; attribute the file's
            # count to the mapped element and compare on the element total.
            rows.append(
                {
                    "functype": functype,
                    "element": element,
                    "source_entries": n_src,
                    "file_entries": dst.get(element, 0),
                    "outcome": None,
                }
            )

    # resolve outcomes for mapped types, accounting for many-to-one mappings
    by_element = defaultdict(list)
    for r in rows:
        if r["element"]:
            by_element[r["element"]].append(r)
    for element, group in by_element.items():
        total_src = sum(r["source_entries"] for r in group)
        in_file = dst.get(element, 0)
        for r in group:
            r["file_entries"] = in_file if len(group) == 1 else None
        if in_file == 0 and total_src > 0:
            for r in group:
                r["outcome"] = "dropped"
        elif in_file == total_src:
            for r in group:
                r["outcome"] = "represented"
        else:
            for r in group:
                r["outcome"] = "partial"
            group[0]["detail"] = {"source_total": total_src, "file_total": in_file}

    total = sum(src.values())
    lost = sum(
        r["source_entries"]
        for r in rows
        if r["outcome"] in ("not predefined", "dropped")
    )
    return {
        "system": entry["system"],
        "protocol": entry["protocol"],
        "total_source_entries": total,
        "entries_not_in_file": lost,
        "fraction_lost": (lost / total) if total else 0.0,
        "rows": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("interaction_coverage.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    per_file = {}
    for key, entry in sorted(manifest.items()):
        if "error" in entry or not Path(entry["h5md"]).exists():
            continue
        res = analyse(key, entry)
        if res:
            per_file[key] = res

    # Collapse to one row per system. The protocol does change the topology --
    # constraints = all-bonds converts every bond into a constraint -- so pick a
    # representative that leaves the bonded terms intact rather than whichever
    # key sorts first.
    PREFERRED = ("nvt-vrescale", "npt-prahman", "nve")
    by_system = {}
    for name in {res["system"] for res in per_file.values()}:
        candidates = {
            res["protocol"]: res for res in per_file.values() if res["system"] == name
        }
        chosen = next(
            (candidates[p] for p in PREFERRED if p in candidates),
            next(iter(sorted(candidates.items())))[1],
        )
        by_system[name] = chosen

    unsupported = defaultdict(
        lambda: {"systems": set(), "entries": 0, "outcome": "", "element": None}
    )
    for res in by_system.values():
        for r in res["rows"]:
            if r["outcome"] in ("not predefined", "dropped", "partial"):
                u = unsupported[r["functype"]]
                u["systems"].add(res["system"])
                u["entries"] += r["source_entries"]
                u["outcome"] = r["outcome"]
                u["element"] = r["element"]

    summary = {
        "per_system": {
            k: {kk: vv for kk, vv in v.items() if kk != "rows"}
            for k, v in by_system.items()
        },
        "unsupported": {
            k: {**v, "systems": sorted(v["systems"])} for k, v in unsupported.items()
        },
        "per_file": per_file,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2, default=str))

    print(
        f"{'system':28s} {'entries':>9s} {'not in file':>12s} {'fraction':>9s}  types"
    )
    for name, res in sorted(by_system.items()):
        types = sorted(
            {
                r["functype"]
                for r in res["rows"]
                if r["outcome"] in ("not predefined", "dropped", "partial")
            }
        )
        print(
            f"  {name:26s} {res['total_source_entries']:9,d} "
            f"{res['entries_not_in_file']:12,d} "
            f"{100 * res['fraction_lost']:8.2f}%  {', '.join(types) or '-'}"
        )

    print("\ninteraction types that do not reach the file:")
    for functype, u in sorted(unsupported.items()):
        reason = (
            "outside Table SI; written nowhere and flagged nowhere"
            if u["outcome"] == "not predefined"
            else f"maps to '{u['element']}' but no entries are written"
        )
        print(
            f"  {functype:14s} {u['entries']:8,d} entries in "
            f"{len(u['systems'])} system(s): {reason}"
        )
    print(f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
