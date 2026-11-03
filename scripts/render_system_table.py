#!/usr/bin/env python
"""
Render the validation system table as camera-ready LaTeX.

    python scripts/render_system_table.py            # to stdout
    python scripts/render_system_table.py --out t.tex

One row per system, grouped by the database it came from. Everything is read
from config/systems.py, so the table cannot disagree with the registry the
validation actually runs on.

Particle counts come from each system's own structure file: the atom count line
of a .gro, or ``gmx dump`` for the systems shipped as a .tpr. Dumping is slow
for the largest systems, so counts are cached in _data/natoms.json and reused.

Used by:
    the paper; not part of any validation stage
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import DATA_ROOT, GMX
from config.systems import SYSTEMS

#: Printed group order, and the caption name of each source.
GROUPS = [
    ("self_generated", "Generated here"),
    ("bena", "BioExcel benchmark suite"),
    ("jcp_benchmark", "GROMACS benchmark set"),
    ("hexabc", "hexABC / MDposit"),
    ("memprotmd", "MemProtMD"),
    ("misato", "MISATO"),
]

#: Force-field identifiers as they should appear in print.
FF_LABEL = {
    "amber99sb": "AMBER99SB",
    "amber99sb-ildn": "AMBER99SB-ILDN",
    "amber99sb-ildn-lipid": "AMBER99SB-ILDN + Berger",
    "amber-hexabc": "AMBER (depositor)",
    "amber-misato": "AMBER ff14SB + GAFF2",
    "charmm27": "CHARMM27",
    "charmm36m": "CHARMM36m",
    "charmm36-mar2014-mod-vs": "CHARMM36 (virtual sites)",
    "gromos54a7": "GROMOS54A7",
    "oplsaa": "OPLS-AA",
    "martini22": "MARTINI 2.2",
}

WATER_LABEL = {
    "tip3p": "TIP3P",
    "tip4p": "TIP4P",
    "spc": "SPC",
    "martini-w": "MARTINI W",
}

#: Systems whose registry description does not reduce to a composition on its
#: own -- the identifier carries the content and the rest restates the force
#: field, which already has a column.
COMPOSITION = {
    "memprotmd_3k1b_cg": "Coarse-grained membrane protein in a DPPC bilayer",
    "memprotmd_5bkf_cg": "Coarse-grained membrane protein in a DPPC bilayer",
    "memprotmd_7e6u_cg": "Coarse-grained membrane protein in a DPPC bilayer",
}


def group_of(system) -> str:
    return str(system.location).split("validation_systems/")[-1].split("/")[0]


def n_particles(system, cache: dict) -> int | None:
    """Atom count from the shipped structure file, cached across runs."""
    if system.id in cache:
        return cache[system.id]
    path = Path(system.location) / system.structure
    if not path.is_file():
        return None
    if path.suffix == ".gro":
        with open(path) as handle:
            handle.readline()
            n = int(handle.readline().strip())
    else:
        # A .tpr is binary; gmx is the only reliable reader of it.
        proc = subprocess.run(
            [GMX, "dump", "-s", str(path)], capture_output=True, text=True, timeout=1800
        )
        n = None
        for line in proc.stdout.splitlines():
            if line.strip().startswith("natoms"):
                n = int(line.split("=")[1])
                break
    cache[system.id] = n
    return n


def tex_escape(text: str) -> str:
    for a, b in (
        ("\\", r"\textbackslash{}"),
        ("_", r"\_"),
        ("%", r"\%"),
        ("&", r"\&"),
        ("#", r"\#"),
    ):
        text = text.replace(a, b)
    # Ions are chemistry, not text: set the charges properly.
    for a, b in (
        ("K+", "K$^+$"),
        ("Na+", "Na$^+$"),
        ("Cl-", "Cl$^-$"),
        ("Mg2+", "Mg$^{2+}$"),
        ("Ca2+", "Ca$^{2+}$"),
    ):
        text = text.replace(a, b)
    for a, b in (
        ("protein-ligand", "protein--ligand"),
        ("protein-DNA", "protein--DNA"),
        ("protein-RNA", "protein--RNA"),
        ("–", "--"),
        ("—", "---"),
    ):
        text = text.replace(a, b)
    return text


def short_description(system) -> str:
    """The composition alone.

    Drops the force-field tail after the first comma, and the ``<SOURCE> <ID>:``
    prefix the deposited systems carry -- the identifier is already the first
    column, so repeating it in prose costs width and says nothing.
    """
    if system.id in COMPOSITION:
        return COMPOSITION[system.id]
    text = system.description.split(",")[0].strip()
    if ": " in text:
        text = text.split(": ", 1)[1].strip()
    # The deposit id often repeats inside the remaining clause as well.
    for part in system.id.split("_"):
        token = part.upper()
        if len(token) == 4 and any(c.isdigit() for c in token):
            text = text.replace(f" {token}", "").replace(f"{token} ", "")
    text = text.strip()
    # Some deposited systems carry a description so terse that stripping the
    # identifier leaves nothing usable; fall back to the clause after it.
    if len(text.split()) < 3:
        tail = system.description.split(",")
        if len(tail) > 1:
            text = f"{text} {tail[1].strip()}".strip()
    return text[0].upper() + text[1:] if text else text


def build(cache: dict) -> str:
    by_group = defaultdict(list)
    for system in SYSTEMS:
        by_group[group_of(system)].append(system)

    header = (
        r"\textbf{System} & \textbf{Composition} & \textbf{Particles} & "
        r"\textbf{Force field} & \textbf{Water} \\"
    )
    rows = [
        r"% Requires \usepackage{booktabs, longtable}",
        r"\begin{longtable}{llrll}",
        r"\caption{Chemical systems of the validation set. Each system is run "
        r"under every simulation protocol its chemistry admits, producing the "
        r"reference files the technical validation is computed on. Particle "
        r"counts are of the system as shipped; the hexABC entries are "
        r"solvent-stripped by the depositor and are simulated in vacuum, which "
        r"is why they carry no water model.}"
        r"\label{tab_validation_systems} \\",
        r"\toprule",
        header,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        header,
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    first = True
    for key, title in GROUPS:
        systems = by_group.get(key, [])
        if not systems:
            continue
        if not first:
            rows.append(r"\addlinespace")
        first = False
        rows.append(
            rf"\multicolumn{{5}}{{l}}{{\textit{{{title}}} "
            rf"($n={len(systems)}$)}} \\"
        )
        for system in sorted(systems, key=lambda s: n_particles(s, cache) or 0):
            n = n_particles(system, cache)
            rows.append(
                " & ".join(
                    [
                        rf"\texttt{{{tex_escape(system.id)}}}",
                        tex_escape(short_description(system)),
                        f"{n:,}" if n else "--",
                        tex_escape(
                            FF_LABEL.get(system.force_field, system.force_field)
                        ),
                        tex_escape(WATER_LABEL.get(system.water_model, "--")),
                    ]
                )
                + r" \\"
            )
    rows.append(r"\end{longtable}")
    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", help="write here instead of stdout")
    args = ap.parse_args()

    cache_path = Path(DATA_ROOT) / "natoms.json"
    cache = json.loads(cache_path.read_text()) if cache_path.is_file() else {}
    before = dict(cache)

    table = build(cache)

    if cache != before:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True))

    if args.out:
        Path(args.out).write_text(table + "\n")
        print(f"written to {args.out}")
    else:
        print(table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
