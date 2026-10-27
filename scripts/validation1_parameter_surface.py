#!/usr/bin/env python
"""
How much of a flat MD parameter file is actually in force?

GROMACS specifies a run as a flat list of key-value pairs. grompp writes the
complete list, defaults included, to mdout.mdp. Many of those options belong to
optional features that are switched off, and the format marks them in no way at all:
`couple-lambda0 = vdw-q` is written with the same syntax and the same apparent
authority as `dt = 0.002`, even when `free-energy = no` makes it inert.

Each option is put in one of three classes: represented in the ``parameters``
group, belonging to a feature the run switches off (or GROMACS no longer
supports), or an engine-specific setting the specification leaves out to stay
independent of the MD engine. The inert share is the motivation for the
conditional subgroups of S7: a schema that instantiates a subgroup only for the
algorithm actually selected records what was in force, rather than the whole
option surface.

    python scripts/validation1_parameter_surface.py [mdout.mdp ...]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import DATA_ROOT, RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


#: .mdp options that the ``parameters`` group represents, directly or through the
#: selections and topology it refers to; mirrors ``roundtrip.extract_source``.
SPEC_MAPPED = frozenset(
    {
        "integrator", "nsteps", "dt", "ld-seed", "bd-fric", "emstep", "emtol",
        "comm-mode", "nstcomm", "comm-grps", "mass-repartition-factor",
        "nstxout", "nstvout", "nstfout", "nstenergy", "nstlog", "energygrps",
        "tcoupl", "nsttcouple", "tc-grps", "ref-t", "tau-t", "nh-chain-length",
        "pcoupl", "pcoupltype", "nstpcouple", "ref-p", "compressibility", "tau-p",
        "coulombtype", "rcoulomb", "epsilon-r", "coulomb-modifier", "epsilon-rf",
        "vdwtype", "rvdw", "DispCorr", "vdw-modifier", "rvdw-switch", "pbc",
        "constraints", "constraint-algorithm", "lincs-iter", "lincs-order",
        "shake-tol",
    }
)  # fmt: skip


OFF = {"", "no", "none", "0", "false", "inactive"}


def _off(option):
    """Inactive when the option is absent or every value is "no", "none", 0 ..."""
    return lambda s: all(v.lower() in OFF for v in s.get(option, "").split() or [""])


def _always(s):
    return True


#: Features the run can switch off, each a list of (inactive?, option pattern)
#: parts. Expanded ensemble is a mode of free energy (``free-energy = expanded``);
#: grompp writes its options only in that mode, except the ``init-lambda-*`` and
#: ``init-wl-*`` initial values, which belong to free energy here.
FAMILIES = {
    "free energy": [
        (
            _off("free-energy"),
            (
                r"^(free-energy|init-lambda|delta-lambda|fep-lambdas|mass-lambdas|"
                r"coul-lambdas|vdw-lambdas|bonded-lambdas|restraint-lambdas|"
                r"temperature-lambdas|calc-lambda|sc-|couple-|dhdl-|separate-dhdl|"
                r"nstdhdl|init-wl-|dh_hist)"
            ),
        )
    ],
    "NMR refinement": [
        (_off("disre"), r"^(disre|nstdisreout)"),
        (_off("orire"), r"^(orire|nstorireout)"),
    ],
    "user hooks": [(_always, r"^(userint|userreal|user1-grps|user2-grps)")],
    "non-equilibrium MD": [
        (_off("acc-grps"), r"^(acc-grps|accelerate)$"),
        (_off("freezegrps"), r"^(freezegrps|freezedim)$"),
        (_off("cos-acceleration"), r"^cos-acceleration$"),
        (_off("deform"), r"^deform"),
    ],
    "walls": [(_off("nwall"), r"^(nwall|wall-)")],
    "external interfaces": [
        (_off(opt), rf"^{opt}$")
        for opt in (
            "density-guided-simulation-active",
            "qmmm-cp2k-active",
            "colvars-active",
            "nnpot-active",
            "fmm-backend",
        )
    ],
    "simulated annealing": [(_off("annealing"), r"^annealing")],
    "simulated tempering": [
        (_off("simulated-tempering"), r"^(simulated-tempering|sim-temp)")
    ],
    "electric field": [(_off("electric-field"), r"^electric-field")],
    # Removed from GROMACS, or undocumented, and without effect.
    "removed features": [(_always, r"^(implicit-solvent|adress|Shake-SOR)$")],
    "QM/MM": [(_off("QMMM"), r"^QMMM")],
    # Used only when the topology has shell particles; the test set has none.
    "shell MD": [(_always, r"^(niter|fcstep)$")],
    "pull code": [(_off("pull"), r"^pull")],
    "AWH": [(_off("awh"), r"^awh")],
    "enforced rotation": [(_off("rotation"), r"^(rotation|rot-)")],
    "ion/water swapping": [
        (
            _off("swapcoords"),
            (
                r"^(swapcoords|swap-|split|massw|solvent-group|coupl-steps|iontype|"
                r"threshold|bulk-offset)"
            ),
        )
    ],
    "test particle insertion": [
        (lambda s: s.get("integrator") not in ("tpi", "tpic"), r"^rtpi$")
    ],
    "interactive MD": [(_off("IMD-group"), r"^IMD-group$")],
}


def analyse(path):
    lines = [
        line
        for line in Path(path).read_text(errors="replace").splitlines()
        if re.match(r"^[a-zA-Z]", line)
    ]
    settings = {}
    for line in lines:
        if "=" in line:
            k, v = line.split("=", 1)
            settings[k.strip()] = v.strip()

    mapped = [k for k in settings if k in SPEC_MAPPED]
    inert, families, counted = 0, {}, set(mapped)
    for name, parts in FAMILIES.items():
        hits = []
        for inactive, pattern in parts:
            if not inactive(settings):
                continue  # the feature is actually in use
            new = [k for k in settings if k not in counted and re.match(pattern, k)]
            counted.update(new)
            hits += new
        if hits:
            families[name] = len(hits)
            inert += len(hits)
    engine = sorted(k for k in settings if k not in counted)
    return {
        "file": str(path),
        "n_options": len(settings),
        "n_mapped": len(mapped),
        "n_inert": inert,
        "inert_families": families,
        "fraction_inert": inert / len(settings) if settings else 0.0,
        "n_engine_specific": len(engine),
        "engine_specific": engine,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "mdp",
        nargs="*",
        help="resolved .mdp files (default: all in _data/resolved_mdp)",
    )
    ap.add_argument("--out", default=_r("parameter_surface.json"))
    args = ap.parse_args()

    paths = args.mdp
    if not paths:
        paths = [
            str(p) for p in sorted((Path(DATA_ROOT) / "resolved_mdp").glob("*.mdp"))
        ]
    if not paths:
        print("no resolved .mdp found; pass one explicitly", file=sys.stderr)
        return 1

    results = [analyse(p) for p in paths]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2))

    for r in results:
        print(f"{r['file']}")
        print(f"  options written by grompp        : {r['n_options']}")
        print(f"  represented in parameters        : {r['n_mapped']}")
        print(
            f"  belonging to disabled features   : {r['n_inert']} "
            f"({100 * r['fraction_inert']:.0f}%), none marked as inert"
        )
        for name, n in sorted(r["inert_families"].items(), key=lambda x: -x[1]):
            print(f"      {name:24s} {n:4d}")
        print(f"  engine-specific                  : {r['n_engine_specific']}")
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
