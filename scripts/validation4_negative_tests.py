#!/usr/bin/env python
"""
Validation 4 -- negative / corruption tests.

Generates one deliberately corrupted file per corruption per base file, runs the
validator over each, and checks that the expected finding code was raised.

    python scripts/validation4_negative_tests.py --quick       # config bases
    python scripts/validation4_negative_tests.py --all         # every reference file
    python scripts/validation4_negative_tests.py BASE.h5md [BASE2.h5md ...]

``--quick`` is for routine runs: a handful of bases chosen for differing
composition, enough that most corruptions apply somewhere.
``--all`` sweeps the whole corpus as a corruption's applicability depends on what
elements a file happens to contain. A full sweep shows which kinds are never exercised.

Each corrupted file is checked for having actually changed. A mutator that
reports success but leaves the file semantically identical is reported as
``inert``: that is a broken test, not a validator failure, and counting it as a
miss sends the reader to debug the wrong thing.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import pathlib
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import DATA_ROOT, NEGATIVE_TEST_BASES, RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


from h5md_validator import Validator
from h5md_validator.corrupt import CORRUPTIONS, generate


def resolve_bases(args) -> list:
    """The files to corrupt: explicit arguments, the whole corpus, or the
    configured quick set."""
    if args.base:
        return [Path(b) for b in args.base]
    if args.all:
        manifest = json.loads(Path(args.manifest).read_text())
        return sorted(
            Path(e["h5md"])
            for e in manifest.values()
            if "error" not in e and Path(e["h5md"]).exists()
        )
    ref = Path(DATA_ROOT) / "reference_files"
    return [
        ref / f"{name}.h5md"
        for name in NEGATIVE_TEST_BASES
        if (ref / f"{name}.h5md").is_file()
    ]


def _norm(value):
    value = np.asarray(value).tolist() if not isinstance(value, (bytes, str)) else value
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, list):
        return [_norm(v) for v in value]
    return value


def snapshot(path) -> dict:
    """A value-level view of the file, for deciding whether a mutator bit.

    Records object existence as well as values: a corruption that adds an empty
    group changes nothing a value-only comparison would notice.
    """
    out = {}
    with h5py.File(path, "r") as f:

        def visit(name, obj):
            out[f"{name}#is"] = "group" if isinstance(obj, h5py.Group) else "dataset"
            for key in obj.attrs:
                out[f"{name}@{key}"] = _norm(obj.attrs[key])
            if isinstance(obj, h5py.Dataset):
                out[f"{name}#shape"] = list(obj.shape)
                out[f"{name}#fields"] = list(obj.dtype.names or [])
                # Hash the raw buffer rather than normalising to Python lists:
                # on a 164 MB file the latter costs 16 s, fifty times the
                # validation it is guarding. Attributes are still compared by
                # value, which is where dtype-insensitivity actually matters --
                # a units list rewritten from S4 to S8 holding the same strings
                # must not read as a change.
                try:
                    arr = np.asarray(obj[()])
                    if arr.dtype.hasobject:
                        # Variable-length strings are stored as object arrays:
                        # tobytes() would hash the pointers, which differ
                        # between two copies of the same file. These arrays are
                        # small, so normalise them instead.
                        payload = repr(_norm(arr)).encode()
                    else:
                        payload = np.ascontiguousarray(arr).tobytes()
                    out[f"{name}#data"] = hashlib.sha256(payload).hexdigest()
                except Exception:
                    pass

        f.visititems(visit)
        for key in f.attrs:
            out[f"/@{key}"] = _norm(f.attrs[key])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("base", nargs="*", help="valid H5MD file(s) to corrupt")
    ap.add_argument(
        "--all", action="store_true", help="every reference file in the manifest"
    )
    ap.add_argument(
        "--quick",
        action="store_true",
        help="the bases declared in config/matrix.py (the default)",
    )
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    ap.add_argument("--out-dir", default=str(Path(DATA_ROOT) / "negative"))
    ap.add_argument("--json", default=_r("negative_tests.json"))
    ap.add_argument(
        "--keep",
        action="store_true",
        help="keep the corrupted files instead of deleting them",
    )
    args = ap.parse_args()

    validator = Validator()
    rows = []

    bases = resolve_bases(args)
    if not bases:
        print("no base files; pass paths, or --all, or --quick")
        return 1
    print(f"{len(bases)} base file(s) x {len(CORRUPTIONS)} corruptions")

    for base in bases:
        base = Path(base)
        out_dir = Path(args.out_dir) / base.stem
        pristine = snapshot(base)
        print(f"\n=== {base.name} ===")
        for corruption, path, skip_reason in generate(base, out_dir):
            if path is None:
                rows.append(
                    {
                        "base": base.name,
                        "corruption": corruption.code,
                        "category": corruption.category,
                        "expected": corruption.expect,
                        "status": "skipped",
                        "reason": skip_reason,
                    }
                )
                continue
            # A mutator that reports success but leaves the file semantically
            # unchanged is a broken test, not a missed detection.
            if snapshot(path) == pristine:
                rows.append(
                    {
                        "base": base.name,
                        "corruption": corruption.code,
                        "category": corruption.category,
                        "expected": corruption.expect,
                        "status": "inert",
                        "reason": "the mutator changed nothing in this file",
                    }
                )
                print(f"  [INRT] {corruption.code:32s} changed nothing")
                if not args.keep:
                    path.unlink(missing_ok=True)
                continue
            report = validator.validate(path)
            codes = {f.code for f in report.findings}
            detected = corruption.expect in codes
            rows.append(
                {
                    "base": base.name,
                    "corruption": corruption.code,
                    "category": corruption.category,
                    "description": corruption.description,
                    "expected": corruption.expect,
                    "status": "detected" if detected else "missed",
                    "raised": sorted(codes),
                    "n_errors": len(report.errors),
                }
            )
            mark = "ok " if detected else "MISS"
            print(f"  [{mark}] {corruption.code:32s} expect {corruption.expect}")
            if not args.keep:
                path.unlink(missing_ok=True)
        if not args.keep:
            with contextlib.suppress(OSError):
                out_dir.rmdir()

    # 'inert' is excluded from `tested`: the file was never corrupted, so the
    # validator was never asked anything and a detection rate over those rows
    # would be measuring the wrong population.
    tested = [r for r in rows if r["status"] in ("detected", "missed")]
    detected = [r for r in tested if r["status"] == "detected"]
    missed = [r for r in tested if r["status"] == "missed"]
    skipped = [r for r in rows if r["status"] == "skipped"]
    inert = [r for r in rows if r["status"] == "inert"]

    # A kind that is skipped or inert on every base was never exercised at all.
    exercised = {r["corruption"] for r in detected}
    never = (
        sorted(set(CORRUPTIONS) - exercised)
        if isinstance(CORRUPTIONS, dict)
        else sorted({c.code for c in CORRUPTIONS} - exercised)
    )

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
        "n_inert": len(inert),
        "n_kinds_exercised": len(exercised),
        "kinds_never_exercised": never,
        "detection_rate": (len(detected) / len(tested)) if tested else 0.0,
        "by_category": by_category,
        "missed": [
            {
                "base": r["base"],
                "corruption": r["corruption"],
                "expected": r["expected"],
                "raised": r["raised"],
            }
            for r in missed
        ],
        "inert": [{"base": r["base"], "corruption": r["corruption"]} for r in inert],
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
    print(f"inert (mutator no-op)    : {len(inert)}")
    print(f"kinds exercised          : {len(exercised)} of {len(CORRUPTIONS)}")
    print(f"detection rate           : {100 * summary['detection_rate']:.1f}%")
    if never:
        print(f"\nNEVER EXERCISED on any base ({len(never)}):")
        for code in never:
            print(f"  {code}")
    if inert:
        print("\nINERT -- the mutator changed nothing, so the test is broken:")
        for r in inert:
            print(f"  {r['base']} / {r['corruption']}")
    if missed:
        print("\nMISSED:")
        for r in missed:
            print(
                f"  {r['base']} / {r['corruption']}: expected {r['expected']}, "
                f"raised {r['raised']}"
            )
    print(f"\nwritten to {out}")
    return 1 if missed else 0


if __name__ == "__main__":
    raise SystemExit(main())
