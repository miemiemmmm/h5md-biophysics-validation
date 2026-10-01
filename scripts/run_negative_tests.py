#!/usr/bin/env python
"""
Validation 4 -- negative / corruption tests.

Generates one deliberately corrupted file per corruption per base file, runs the
validator over each, and checks that the expected finding code was raised.

    python scripts/run_negative_tests.py BASE.h5md [BASE2.h5md ...] \
        --out-dir results/negative --json results/negative_tests.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from matrix import RESULTS_ROOT

#: Results live in ``results/`` by default. RESULTS_DIR redirects every script
#: at once, so an existing set of results can be analysed from anywhere.
RESULTS_DIR = os.environ.get("RESULTS_DIR", RESULTS_ROOT)


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)



from h5md_validator import Validator
from h5md_validator.corrupt import CORRUPTIONS, generate


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("base", nargs="+", help="valid H5MD file(s) to corrupt")
    ap.add_argument("--out-dir", default=_r("negative"))
    ap.add_argument("--json", default=_r("negative_tests.json"))
    ap.add_argument("--keep", action="store_true",
                    help="keep the corrupted files instead of deleting them")
    args = ap.parse_args()

    validator = Validator()
    rows = []

    for base in args.base:
        base = Path(base)
        out_dir = Path(args.out_dir) / base.stem
        print(f"\n=== {base.name} ===")
        for corruption, path, skip_reason in generate(base, out_dir):
            if path is None:
                rows.append({
                    "base": base.name, "corruption": corruption.code,
                    "category": corruption.category, "expected": corruption.expect,
                    "status": "skipped", "reason": skip_reason,
                })
                continue
            report = validator.validate(path)
            codes = {f.code for f in report.findings}
            detected = corruption.expect in codes
            rows.append({
                "base": base.name,
                "corruption": corruption.code,
                "category": corruption.category,
                "description": corruption.description,
                "expected": corruption.expect,
                "status": "detected" if detected else "MISSED",
                "raised": sorted(codes - {f.code for f in report.deviations}),
                "n_errors": len(report.errors),
            })
            mark = "ok " if detected else "MISS"
            print(f"  [{mark}] {corruption.code:32s} expect {corruption.expect}")
            if not args.keep:
                path.unlink(missing_ok=True)
        if not args.keep:
            try:
                out_dir.rmdir()
            except OSError:
                pass

    tested = [r for r in rows if r["status"] != "skipped"]
    detected = [r for r in tested if r["status"] == "detected"]
    missed = [r for r in tested if r["status"] == "MISSED"]
    skipped = [r for r in rows if r["status"] == "skipped"]

    by_category = {}
    for r in tested:
        c = by_category.setdefault(r["category"], {"tested": 0, "detected": 0})
        c["tested"] += 1
        c["detected"] += r["status"] == "detected"

    summary = {
        "n_corruption_kinds": len(CORRUPTIONS),
        "n_files_tested": len(tested),
        "n_detected": len(detected),
        "n_missed": len(missed),
        "n_skipped": len(skipped),
        "detection_rate": (len(detected) / len(tested)) if tested else 0.0,
        "by_category": by_category,
        "missed": [{"base": r["base"], "corruption": r["corruption"],
                    "expected": r["expected"], "raised": r["raised"]} for r in missed],
        "rows": rows,
    }

    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 72)
    print(f"corruption kinds defined : {len(CORRUPTIONS)}")
    print(f"corrupted files tested   : {len(tested)}")
    print(f"detected                 : {len(detected)}")
    print(f"missed                   : {len(missed)}")
    print(f"not applicable (skipped) : {len(skipped)}")
    print(f"detection rate           : {100 * summary['detection_rate']:.1f}%")
    if missed:
        print("\nMISSED:")
        for r in missed:
            print(f"  {r['base']} / {r['corruption']}: expected {r['expected']}, "
                  f"raised {r['raised']}")
    print(f"\nwritten to {out}")
    return 1 if missed else 0


if __name__ == "__main__":
    raise SystemExit(main())
