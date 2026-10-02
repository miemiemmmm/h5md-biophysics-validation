#!/usr/bin/env python
"""
Produce the reference H5MD files the validation runs against.

For each (system, protocol) pair: write an .mdp, grompp it into a .tpr, and run
a short mdrun writing H5MD output.  The .tpr and the H5MD file are the two sides
of the round-trip test, so both are kept.

    python scripts/build_reference_files.py --list
    python scripts/build_reference_files.py --system bena_12k --protocol nvt-vrescale
    python scripts/build_reference_files.py --all
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import DATA_ROOT, GMX, PROTOCOLS, PROTOCOLS_BY_ID, RESULTS_DIR
from config.systems import SYSTEMS, SYSTEMS_BY_ID, prepare_system



def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


ROOT = Path(DATA_ROOT)


def run(cmd, cwd=None, check=True, timeout=7200, input=None):
    proc = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True,
                          text=True, timeout=timeout, input=input)
    if check and proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(str(c) for c in cmd)}\n"
            f"--- stdout ---\n{proc.stdout[-3000:]}\n--- stderr ---\n{proc.stderr[-4000:]}"
        )
    return proc


def write_mdp(path: Path, options: dict):
    lines = [f"{k:32s} = {v}" for k, v in options.items()]
    path.write_text("\n".join(lines) + "\n")


def build_tpr(system, protocol, work: Path, force=False) -> Path:
    """grompp one (system, protocol) pair into a .tpr."""
    tpr = ROOT / "tprs" / f"{system.id}__{protocol.id}.tpr"
    if tpr.exists() and not force:
        return tpr
    tpr.parent.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    prep = prepare_system(system)
    mdp = work / "grompp.mdp"
    options = dict(protocol.mdp)
    options.update(prep.get("mdp_overrides", {}))
    write_mdp(mdp, options)

    # Keep grompp's resolved option list: the complete .mdp including every
    # default, which is what scripts/validation1_parameter_surface.py measures.
    # The scratch directory is deleted after each pair, so it lands elsewhere.
    resolved_mdp = ROOT / "resolved_mdp" / f"{system.id}__{protocol.id}.mdp"
    resolved_mdp.parent.mkdir(parents=True, exist_ok=True)

    cmd = [GMX, "-quiet", "-nobackup", "grompp",
           "-f", mdp, "-c", prep["structure"], "-r", prep["structure"],
           "-p", prep["topology"], "-o", tpr, "-po", resolved_mdp, "-maxwarn", "20"]
    if prep.get("index"):
        cmd += ["-n", prep["index"]]
    run(cmd, cwd=work)
    return tpr


def run_mdrun(tpr: Path, out: Path, ntomp=8, force=False) -> dict:
    """Run a short mdrun writing H5MD output; return timing and size."""
    # The .edr is GROMACS's own record of the same observables the H5MD file
    # carries, so it is the source side of the observables comparison. Like the
    # resolved .mdp it has to land outside the scratch directory to survive.
    edr = ROOT / "edr" / f"{out.stem}.edr"
    if out.exists() and edr.exists() and not force:
        return {"cached": True, "size_bytes": out.stat().st_size, "edr": str(edr)}
    out.parent.mkdir(parents=True, exist_ok=True)
    edr.parent.mkdir(parents=True, exist_ok=True)
    work = out.parent / f".run_{out.stem}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    t0 = time.perf_counter()
    try:
        run([GMX, "-quiet", "-nobackup", "mdrun", "-s", tpr, "-o", out,
         "-ntmpi", "1", "-ntomp", str(ntomp),
             "-g", work / "md.log", "-e", edr, "-c", work / "out.gro",
             "-cpo", work / "state.cpt"], cwd=work)
    except Exception:
        # A writer that aborts mid-run leaves a partial, unreadable file behind.
        # Keep it out of the reference corpus; the failure is recorded instead.
        out.unlink(missing_ok=True)
        edr.unlink(missing_ok=True)
        shutil.rmtree(work, ignore_errors=True)
        raise
    elapsed = time.perf_counter() - t0
    shutil.rmtree(work, ignore_errors=True)
    return {"cached": False, "wall_seconds": elapsed,
            "size_bytes": out.stat().st_size, "edr": str(edr)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--system", action="append", default=[])
    ap.add_argument("--protocol", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--ntomp", type=int, default=8)
    ap.add_argument("--manifest", default=_r("reference_files.json"))
    args = ap.parse_args()

    if args.list:
        print("SYSTEMS")
        for s in SYSTEMS:
            print(f"  {s.id:26s} {s.description}")
        print("\nPROTOCOLS")
        for p in PROTOCOLS:
            print(f"  {p.id:26s} {p.description}")
        return 0

    systems = [SYSTEMS_BY_ID[s] for s in args.system] if args.system else SYSTEMS
    protocols = ([PROTOCOLS_BY_ID[p] for p in args.protocol]
                 if args.protocol else PROTOCOLS)
    if not args.all and not args.system and not args.protocol:
        ap.error("pass --all, or select with --system / --protocol")

    manifest_path = Path(args.manifest)
    manifest = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())

    for system in systems:
        allowed = set(system.protocols or [p.id for p in protocols])
        for protocol in protocols:
            if protocol.id not in allowed:
                continue
            # A non-periodic protocol only makes sense for a vacuum system, and
            # a periodic system only under a periodic protocol.
            if (protocol.mdp.get("pbc") == "no") != (not system.periodic):
                continue
            key = f"{system.id}__{protocol.id}"
            work = ROOT / "scratch" / key
            print(f"--- {key}")
            try:
                tpr = build_tpr(system, protocol, work, force=args.force)
                out = ROOT / "reference_files" / f"{key}.h5md"
                info = run_mdrun(tpr, out, ntomp=args.ntomp, force=args.force)
                manifest[key] = {
                    "system": system.id, "protocol": protocol.id,
                    "tpr": str(tpr), "h5md": str(out),
                    "dimensions": list(system.dimensions),
                    "exercises": list(protocol.exercises),
                    **info,
                }
                print(f"    ok  {out.name}  {info.get('size_bytes', 0) / 1e6:.1f} MB")
            except Exception as exc:
                manifest[key] = {"system": system.id, "protocol": protocol.id,
                                 "error": f"{type(exc).__name__}: {exc}"}
                print(f"    FAILED: {str(exc)[:400]}")
            finally:
                shutil.rmtree(work, ignore_errors=True)
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            manifest_path.write_text(json.dumps(manifest, indent=2))

    ok = sum(1 for v in manifest.values() if "error" not in v)
    print(f"\n{ok}/{len(manifest)} reference files available; manifest at {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
