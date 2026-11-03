#!/usr/bin/env python
"""
Validation 2 -- round-trip fidelity.

Compares each reference file with the GROMACS run input it was written from:
every metadata item that ``gmx dump`` reports for the .tpr, and every observable
that ``gmx energy`` reports for the .edr, is looked up in the H5MD file. The two
sides share no code path, so a match means the item survived serialising and
parsing, not that one reader agrees with itself.

Two results, kept separate:

1. PER FILE.  Each compared item with a status of match, mismatch,
   missing_in_file or missing_in_source, and a ``faithful`` verdict.

2. ITEM FIDELITY.  The same comparisons aggregated over all files, per item:
   is it represented in the files, missing, or incorrect?

    python scripts/validation2_roundtrip.py

Writes <results>/roundtrip.json: a summary, ``item_fidelity``, and under
``files`` one entry per reference file.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


from config.matrix import DATA_ROOT, GMX
from h5md_validator import gmxdump, roundtrip

DUMPS = Path(DATA_ROOT) / "tpr_dumps"


def item_fidelity(files: dict) -> dict:
    """Aggregate the per-file comparisons per item, over every file."""
    items = defaultdict(
        lambda: {
            "category": "",
            "files": 0,
            "match": 0,
            "mismatch": 0,
            "missing_in_file": 0,
            "extra_in_file": 0,
        }
    )
    for rec in files.values():
        for comp in rec.get("comparisons", []):
            # collapse the per-molecule-type prefix so items aggregate
            name = comp["key"]
            parts = name.split(".")
            if parts[0] == "topology" and len(parts) > 2:
                parts[1] = "<molecule>"
                name = ".".join(parts)
            e = items[name]
            e["category"] = comp["category"]
            if comp["status"] in ("match", "match_unordered"):
                e["files"] += 1
                e["match"] += 1
            elif comp["status"] == "mismatch":
                e["files"] += 1
                e["mismatch"] += 1
            elif comp["status"] == "missing_in_file":
                e["missing_in_file"] += 1
            elif comp["status"] == "missing_in_source":
                e["extra_in_file"] += 1

    rows = [{"item": k, **v} for k, v in sorted(items.items())]
    represented = [r for r in rows if r["match"] or r["mismatch"]]
    unsupported = [r for r in rows if r["missing_in_file"] and not r["match"]]
    incorrect = [r for r in rows if r["mismatch"]]
    return {
        "n_items": len(rows),
        "n_represented": len(represented),
        "n_unsupported": len(unsupported),
        "n_incorrect": len(incorrect),
        "unsupported": [r["item"] for r in unsupported],
        "incorrect": [r["item"] for r in incorrect],
        "rows": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("roundtrip.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    files = {}
    for key, entry in sorted(manifest.items()):
        h5 = Path(entry.get("h5md", ""))
        if "error" in entry or not h5.exists() or not Path(entry["tpr"]).exists():
            print(f"  --  {key}: not compared")
            continue
        record = {"system": entry["system"], "protocol": entry["protocol"]}
        try:
            dump = gmxdump.from_tpr(entry["tpr"], gmx=GMX, cache_dir=DUMPS)
            # The .edr is the source side of the observables comparison. It is
            # absent from a corpus built before it was kept, in which case
            # observables are simply not compared.
            edr = None
            edr_path = entry.get("edr")
            if edr_path and Path(edr_path).exists():
                edr = gmxdump.from_edr(edr_path, gmx=GMX, cache_dir=DUMPS)
            record.update(
                roundtrip.roundtrip(dump, h5, entry["tpr"], edr=edr).to_dict()
            )
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        files[key] = record
        print(
            f"  {'OK  ' if record.get('faithful') else 'DIFF'}  {key:44s} "
            f"match={record.get('n_match', '-')} "
            f"mismatch={record.get('n_mismatch', '-')} "
            f"missing={record.get('n_missing_in_file', '-')} "
            f"dropped={record.get('n_dropped_terms', '-')}"
        )

    compared = [v for v in files.values() if "error" not in v]
    fidelity = item_fidelity(files)
    summary = {
        "n_compared": len(compared),
        "n_faithful": sum(1 for v in compared if v.get("faithful")),
        "n_all_terms_written": sum(1 for v in compared if v.get("all_terms_written")),
        "n_errors": len(files) - len(compared),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(
            {"summary": summary, "item_fidelity": fidelity, "files": files}, indent=2
        )
    )

    print("\n" + "=" * 72)
    print(
        f"round-trip faithful       : {summary['n_faithful']} of {summary['n_compared']}"
    )
    print(
        f"every source term written : {summary['n_all_terms_written']} of "
        f"{summary['n_compared']}"
    )
    print(f"metadata items compared   : {fidelity['n_items']}")
    print(f"  represented in the files: {fidelity['n_represented']}")
    print(f"  present in source, absent: {fidelity['n_unsupported']}")
    print(f"  represented but incorrect: {fidelity['n_incorrect']}")
    if fidelity["unsupported"]:
        print("\npresent in the GROMACS source but absent from every file:")
        for e in fidelity["unsupported"]:
            print(f"  {e}")
    print(f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
