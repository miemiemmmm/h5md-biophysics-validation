#!/usr/bin/env python
"""
Validation 3 -- internal consistency.

Runs the conformance validator over every reference file: each file is checked
against the specification model on its own, with no other input. A file is
conformant when no ERROR finding is raised; the codes are listed in
docs/check-codes.md.

    python scripts/validation3_conformance.py

Writes <results>/conformance.json: a summary, and under ``files`` one entry per
reference file with its system, protocol, verdict, counts and findings.
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


from h5md_validator import Validator


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out", default=_r("conformance.json"))
    args = ap.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    validator = Validator()
    files = {}
    for key, entry in sorted(manifest.items()):
        h5 = Path(entry.get("h5md", ""))
        if "error" in entry or not h5.exists():
            print(f"  --  {key}: not generated")
            continue
        report = validator.validate(h5)
        files[key] = {
            "system": entry["system"],
            "protocol": entry["protocol"],
            **report.to_dict(),
        }
        c = report.to_dict()["counts"]
        print(
            f"  {'PASS' if report.conformant else 'FAIL'}  {key:44s} "
            f"errors={c['ERROR']:<3d} warnings={c['WARNING']:<3d}"
        )

    summary = {
        "n_files": len(files),
        "n_conformant": sum(1 for v in files.values() if v["conformant"]),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "files": files}, indent=2))

    print("\n" + "=" * 72)
    print(f"reference files validated : {summary['n_files']}")
    print(f"conformant (0 errors)     : {summary['n_conformant']}")
    print(f"written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
