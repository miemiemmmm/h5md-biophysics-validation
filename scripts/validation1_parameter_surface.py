#!/usr/bin/env python
"""
How much of a flat MD parameter file is actually in force?

GROMACS specifies a run as a flat list of key-value pairs. grompp writes the
complete list, defaults included, to mdout.mdp. Many of those options belong to
subsystems that are switched off, and the format marks them in no way at all:
`couple-lambda0 = vdw-q` is written with the same syntax and the same apparent
authority as `dt = 0.002`, even when `free-energy = no` makes it inert.

This measures the ratio, which is the motivation for the conditional subgroups
of S7: a schema that instantiates a subgroup only for the algorithm actually
selected records what was in force, rather than the whole option surface.

    python scripts/validation1_parameter_surface.py [mdout.mdp ...]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import DATA_ROOT, RESULTS_DIR

#: Option families and the master switch that governs them. A family counts as
#: inert when its switch is absent, "no", "none" or zero.


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


FAMILIES = {
    "free energy": (
        "free-energy",
        r"^(free-energy|init-lambda|delta-lambda|fep-lambdas|mass-lambdas|"
        r"coul-lambdas|vdw-lambdas|bonded-lambdas|restraint-lambdas|"
        r"temperature-lambdas|calc-lambda|sc-|couple-|dhdl-|separate-dhdl|"
        r"nstdhdl|init-lambda-state)"),
    "pull code": ("pull", r"^pull"),
    "AWH": ("awh", r"^awh"),
    "enforced rotation": ("rotation", r"^(rotation|rot-)"),
    "QM/MM": ("QMMM", r"^QMMM"),
    "walls": ("nwall", r"^(nwall|wall-)"),
    "ion/water swapping": ("swapcoords", r"^(swapcoords|swap-|split|massw|"
                                         r"solvent-group|coupl-steps|iontype|"
                                         r"threshold|bulk-offset)"),
    "simulated tempering": ("simulated-tempering", r"^(simulated-tempering|sim-temp)"),
    "electric field": ("electric-field", r"^electric-field"),
    "expanded ensemble": ("nstexpanded", r"^(nstexpanded|lmc-|weight-|wl-|"
                                         r"init-lambda-weights|lmc)"),
    "user hooks": (None, r"^(userint|userreal)"),
}

OFF = {"", "no", "none", "0", "false"}


def analyse(path):
    lines = [l for l in Path(path).read_text(errors="replace").splitlines()
             if re.match(r"^[a-zA-Z]", l)]
    settings = {}
    for line in lines:
        if "=" in line:
            k, v = line.split("=", 1)
            settings[k.strip()] = v.strip()

    inert, families = 0, {}
    for name, (switch, pattern) in FAMILIES.items():
        if switch is not None and settings.get(switch, "").lower() not in OFF:
            continue                      # the subsystem is actually in use
        hits = [k for k in settings if re.match(pattern, k)]
        if hits:
            families[name] = len(hits)
            inert += len(hits)
    return {"file": str(path), "n_options": len(settings),
            "n_inert": inert, "inert_families": families,
            "fraction_inert": inert / len(settings) if settings else 0.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mdp", nargs="*", help="mdout.mdp files (default: find one)")
    ap.add_argument("--out", default=_r("parameter_surface.json"))
    args = ap.parse_args()

    paths = args.mdp
    if not paths:
        root = Path(DATA_ROOT)
        paths = [str(p) for p in sorted(root.glob("systems/*/mdout.mdp"))][:1]
    if not paths:
        print("no mdout.mdp found; pass one explicitly", file=sys.stderr)
        return 1

    results = [analyse(p) for p in paths]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2))

    for r in results:
        print(f"{r['file']}")
        print(f"  options written by grompp        : {r['n_options']}")
        for name, n in sorted(r["inert_families"].items(), key=lambda x: -x[1]):
            print(f"      {name:22s} {n:4d}")
        print(f"  belonging to disabled subsystems : {r['n_inert']} "
              f"({100 * r['fraction_inert']:.0f}%)")
        print("  none of which the format marks as inert")
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
