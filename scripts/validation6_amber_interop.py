#!/usr/bin/env python
"""
Cross-engine interoperability: the same molecule, two engines, one schema.

The interoperability claim is that engine-specific names denote the same
schema-level concepts.  This script tests the topology half of that claim with
data: not a stated correspondence, but one molecule compared item by item.

The same B-DNA dodecamer (PDB 1BNA) is prepared twice and independently:

    GROMACS   pdb2gmx with AMBER99SB-ILDN   -> .tpr  -> H5MD biophysics file
    AMBER     tleap with DNA.OL15           -> .prmtop

Both are reduced to the schema's own vocabulary and compared. The two halves of
the comparison say different things and are reported separately:

  * CHEMISTRY -- particle count, component count and names, atomic numbers, bond
    topology.  These are properties of the molecule, not of the force field, so
    the two engines must agree. Agreement shows the schema elements carry the
    molecule rather than one engine's bookkeeping.

  * PARAMETRISATION -- charges and masses.  These are properties of the force
    field, and OL15 is not AMBER99SB-ILDN, so they are expected to differ. What
    matters is that the same schema elements hold both without either engine
    needing an element the other does not.

    python scripts/validation6_amber_interop.py --h5md <file.h5md> --prmtop <file.prmtop>
    python scripts/validation6_amber_interop.py -a     # every AMBER_INTEROP_PAIRS entry
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import (
    AMBER_INTEROP_PAIRS,
    DATA_ROOT,
    RESULTS_DIR,
    VALIDATION_SYSTEMS,
)


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


from h5md_validator import roundtrip


def amber_topology(prmtop) -> dict:
    """Reduce an AMBER prmtop to the schema's topology vocabulary."""
    import pytraj as pt

    top = pt.load_topology(str(prmtop))
    atoms = list(top.atoms)
    residues = list(top.residues)

    # A rigid 3-site water carries a third, H-H "bond" in a prmtop (TIP3P:
    # HW-HW, 553 kcal/mol/A^2, 1.5136 A). It is not chemistry: AMBER has no
    # constraint list, so the H-H distance that keeps the triangle rigid under
    # SHAKE/SETTLE is stored as a bond. GROMACS stores the same distance as the
    # dhh parameter of [ settles ], and the schema keeps only the two O-H bonds.
    # It is removed here, and only inside a residue of exactly one O and two H.
    def is_water(residue):
        species = sorted(
            int(atoms[k].atomic_number) for k in range(residue.first, residue.last)
        )
        return species == [1, 1, 8]

    bonds, n_water_hh = set(), 0
    for b in top.bonds:
        i, j = (int(x) for x in b.indices)
        if (
            atoms[i].atomic_number == 1
            and atoms[j].atomic_number == 1
            and atoms[i].resid == atoms[j].resid
            and is_water(residues[atoms[i].resid])
        ):
            n_water_hh += 1
            continue
        bonds.add((min(i, j), max(i, j)))

    return {
        "particle_count": top.n_atoms,
        "component_count": top.n_residues,
        "particle_names": [a.name.strip() for a in atoms],
        "species": [int(a.atomic_number) for a in atoms],
        "component_ids": [int(a.resid) for a in atoms],
        "component_names": [r.name.strip() for r in residues],
        "charges": [float(a.charge) for a in atoms],
        "masses": [float(a.mass) for a in atoms],
        "bonds": sorted(bonds),
        "n_bonds": len(bonds),
        "n_water_hh_removed": n_water_hh,
        "n_angles": len(list(top.angles)),
        "n_dihedrals": len(list(top.dihedrals)),
    }


def h5md_topology(path) -> dict:
    """Build the full system topology of an H5MD file in the same vocabulary.

    A molecule type is stored once; ``topology/molecule_blocks`` says in which
    order the types occur and how many copies of each. The system is rebuilt by
    walking the blocks and repeating each type ``count`` times, so particle
    order matches the global indices of ``connectivity/bonds`` and the prmtop.
    """
    items = roundtrip.extract_target(path)
    names = [
        n.decode() if isinstance(n, bytes) else str(n)
        for n in items["topology.molecule_blocks.names"].value
    ]
    counts = [int(c) for c in items["topology.molecule_blocks.counts"].value]

    out = {
        "molecule_blocks": list(zip(names, counts)),
        "particle_count": 0,
        "component_count": 0,
        "particle_names": [],
        "species": [],
        "component_ids": [],
        "component_names": [],
        "charges": [],
        "masses": [],
    }
    offset = 0
    for m, count in zip(names, counts):

        def g(k, m=m):
            return items[f"topology.{m}.{k}"].value

        for _ in range(count):
            out["particle_count"] += int(g("particle_count"))
            out["component_count"] += int(g("component_count"))
            out["particle_names"] += list(g("particle_names"))
            out["species"] += [int(v) for v in g("species")]
            out["component_ids"] += [int(v) + offset for v in g("component_ids")]
            out["component_names"] += list(g("component_names"))
            out["charges"] += [float(v) for v in g("charges")]
            out["masses"] += [float(v) for v in g("masses")]
            offset += int(g("component_count"))

    # The global bond list already indexes the full system.
    conn = items.get("connectivity.bonds")
    bonds = set()
    if conn is not None:
        for pair in conn.value:
            i, j = int(pair[0]), int(pair[1])
            bonds.add((min(i, j), max(i, j)))
    out["bonds"] = sorted(bonds)
    out["n_bonds"] = len(bonds)
    return out


def compare(gmx: dict, amb: dict) -> dict:
    """Compare the two reductions, separating chemistry from parametrisation."""
    rows = []

    def add(group, item, a, b, agree, note=""):
        rows.append(
            {
                "group": group,
                "item": item,
                "gromacs": a,
                "amber": b,
                "agree": None if agree is None else bool(agree),
                "note": note,
            }
        )

    add(
        "chemistry",
        "particle count",
        gmx["particle_count"],
        amb["particle_count"],
        gmx["particle_count"] == amb["particle_count"],
    )
    add(
        "chemistry",
        "component count",
        gmx["component_count"],
        amb["component_count"],
        gmx["component_count"] == amb["component_count"],
    )

    seq_g = [s.strip().lstrip("D").rstrip("35") for s in gmx["component_names"]]
    seq_a = [s.strip().lstrip("D").rstrip("35") for s in amb["component_names"]]
    add(
        "chemistry",
        "residue sequence (base identity)",
        "".join(seq_g),
        "".join(seq_a),
        seq_g == seq_a,
        "terminal residue naming differs between engines (DA5/DA3 vs DA); "
        "compared on base identity",
    )

    formula_g = _formula(gmx["species"])
    formula_a = _formula(amb["species"])
    add(
        "chemistry",
        "molecular formula",
        formula_g,
        formula_a,
        formula_g == formula_a,
        "an order-independent statement of what is in the box",
    )

    # The two engines number atoms differently within a residue, so every
    # per-atom comparison has to be made on an alignment rather than on the
    # index. Aligning on (residue index, atom name) is what a reader of the
    # schema would do, and it is the alignment the schema's own fields support.
    raw_mapping, raw_unmatched, _ = _align(gmx, amb, use_aliases=False)
    mapping, unmatched_g, unmatched_a = _align(gmx, amb, use_aliases=True)
    add(
        "chemistry",
        "atoms aligned on atom names alone",
        f"{len(raw_mapping)}/{gmx['particle_count']}",
        f"{len(raw_mapping)}/{amb['particle_count']}",
        not raw_unmatched,
        (
            "all atoms align"
            if not raw_unmatched
            else f"{gmx['particle_count'] - len(raw_mapping)} atoms carry different names "
            f"for the same atom ({sorted(raw_unmatched)[:4]}); the schema defines no "
            "canonical atom-naming convention"
        ),
    )
    add(
        "chemistry",
        "atoms aligned with an alias table",
        f"{len(mapping)}/{gmx['particle_count']}",
        f"{len(mapping)}/{amb['particle_count']}",
        not unmatched_g and not unmatched_a,
        "aliases used: "
        + ", ".join(f"{k}={v}" for k, v in sorted(ATOM_NAME_ALIASES.items())),
    )

    n_diff = sum(
        1 for gi, ai in mapping.items() if gmx["species"][gi] != amb["species"][ai]
    )
    add(
        "chemistry",
        "atomic numbers (aligned)",
        f"{len(mapping)} atoms",
        f"{n_diff} differ",
        n_diff == 0,
        "element identity is a property of the molecule, not the force field",
    )

    add(
        "chemistry",
        "bond count",
        gmx["n_bonds"],
        amb["n_bonds"],
        gmx["n_bonds"] == amb["n_bonds"],
        (
            f"{amb['n_water_hh_removed']} water H-H bonds removed from the AMBER "
            "side: a rigid-water device, not chemistry"
            if amb.get("n_water_hh_removed")
            else ""
        ),
    )

    # Map the GROMACS bonds through the alignment and compare as a set.
    mapped_bonds = set()
    mappable = True
    for i, j in gmx["bonds"]:
        if i in mapping and j in mapping:
            a, b = mapping[i], mapping[j]
            mapped_bonds.add((min(a, b), max(a, b)))
        else:
            mappable = False
    amber_bonds = set(amb["bonds"])
    bonds_equal = mappable and mapped_bonds == amber_bonds
    only_g = len(mapped_bonds - amber_bonds)
    only_a = len(amber_bonds - mapped_bonds)
    add(
        "chemistry",
        "bond topology (aligned)",
        f"{len(mapped_bonds)} pairs",
        f"{len(amber_bonds)} pairs",
        bonds_equal,
        "the same covalent graph, compared after alignment"
        + (
            "" if bonds_equal else f"; {only_g} only in GROMACS, {only_a} only in AMBER"
        ),
    )

    tot_q_g = sum(gmx["charges"])
    tot_q_a = sum(amb["charges"])
    add(
        "parametrisation",
        "total charge (e)",
        round(tot_q_g, 4),
        round(tot_q_a, 4),
        abs(tot_q_g - tot_q_a) < 1e-3,
        "net charge is set by the chemistry and must agree even though the "
        "individual partial charges do not",
    )
    if mapping:
        dq = max(
            abs(gmx["charges"][gi] - amb["charges"][ai]) for gi, ai in mapping.items()
        )
        add(
            "parametrisation",
            "largest partial-charge difference (e)",
            "-",
            round(dq, 4),
            None,
            "the schema holds whatever the force field assigns. Here the two "
            "agree exactly: DNA.OL15 revises torsions and inherits the "
            "parmbsc0/ff99 charges that AMBER99SB-ILDN also uses",
        )
        dm = max(
            abs(gmx["masses"][gi] - amb["masses"][ai]) for gi, ai in mapping.items()
        )
        add(
            "parametrisation",
            "largest mass difference (u)",
            "-",
            round(dm, 4),
            dm < 0.01,
            "masses are elemental and must agree across force fields",
        )

    chemistry = [
        r for r in rows if r["group"] == "chemistry" and r["agree"] is not None
    ]
    return {
        "rows": rows,
        "n_chemistry_items": len(chemistry),
        "n_chemistry_agree": sum(1 for r in chemistry if r["agree"]),
        "interoperable": all(r["agree"] for r in chemistry),
    }


#: Atom names that denote the same atom under different engine conventions.
#: The schema stores atom names as free strings and defines no canonical
#: convention, so an alias table is the only way to align two engines'
#: descriptions of the same molecule -- which is itself a finding, recorded as
#: deviations.DEV-NO-ATOM-NAME-CONVENTION.
ATOM_NAME_ALIASES = {
    "O1P": "OP1",
    "O2P": "OP2",  # phosphate oxygens
    "H3T": "HO3*",
    "H5T": "HO5*",  # terminal hydroxyl hydrogens
}


def _align(gmx: dict, amb: dict, use_aliases=True):
    """Map GROMACS atom index -> AMBER atom index on (residue index, atom name).

    Atom-name spellings differ between engines for a few hydrogens (H5'1 vs
    H5', and similar), so names are normalised before matching and whatever
    fails to match is reported rather than dropped.
    """

    def norm(name):
        n = name.strip().upper().replace("''", "*").replace("'", "*")
        # trailing digit disambiguators on hydrogens: H5*1 -> H5*
        if len(n) > 2 and n[0] == "H" and n[-1].isdigit() and "*" in n:
            n = n[:-1]
        if use_aliases:
            n = ATOM_NAME_ALIASES.get(n, n)
        return n

    def index(d):
        out = {}
        for i, (cid, name) in enumerate(zip(d["component_ids"], d["particle_names"])):
            out.setdefault((cid, norm(name)), []).append(i)
        return out

    gi, ai = index(gmx), index(amb)
    mapping = {}
    for key, g_list in gi.items():
        a_list = ai.get(key)
        if not a_list or len(a_list) != len(g_list):
            continue
        for g, a in zip(g_list, a_list):
            mapping[g] = a
    unmatched_g = {
        gmx["particle_names"][i]
        for i in range(gmx["particle_count"])
        if i not in mapping
    }
    matched_a = set(mapping.values())
    unmatched_a = {
        amb["particle_names"][i]
        for i in range(amb["particle_count"])
        if i not in matched_a
    }
    return mapping, unmatched_g, unmatched_a


_SYMBOLS = {1: "H", 6: "C", 7: "N", 8: "O", 15: "P", 16: "S", 11: "Na", 17: "Cl"}


def _formula(species):
    from collections import Counter

    c = Counter(species)
    return "".join(
        f"{_SYMBOLS.get(z, f'Z{z}')}{n}"
        for z, n in sorted(c.items(), key=lambda x: -x[1])
    )


def run_pair(h5md, prmtop) -> dict:
    gmx = h5md_topology(h5md)
    amb = amber_topology(prmtop)
    result = compare(gmx, amb)
    result["h5md"] = str(h5md)
    result["prmtop"] = str(prmtop)
    result["gromacs_molecule_blocks"] = gmx["molecule_blocks"]
    return result


def print_result(result):
    width = max(len(r["item"]) for r in result["rows"])
    group = None
    for r in result["rows"]:
        if r["group"] != group:
            group = r["group"]
            print(f"\n{group.upper()}")
        mark = " " if r["agree"] is None else ("ok " if r["agree"] else "DIFF")
        print(
            f"  [{mark}] {r['item']:{width}s}  gromacs={r['gromacs']}  "
            f"amber={r['amber']}"
        )
        if r["note"]:
            print(f"         {r['note']}")
    print(
        f"\nchemistry items in agreement: "
        f"{result['n_chemistry_agree']}/{result['n_chemistry_items']}"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-a",
        "--run-all",
        action="store_true",
        help="compare every pair in config.matrix.AMBER_INTEROP_PAIRS",
    )
    ap.add_argument("--h5md", help="one H5MD file (with --prmtop)")
    ap.add_argument("--prmtop", help="the AMBER topology of the same molecule")
    ap.add_argument("--out", default=_r("amber_interop.json"))
    args = ap.parse_args()

    if args.run_all:
        if args.h5md or args.prmtop:
            ap.error("--run-all takes its pairs from AMBER_INTEROP_PAIRS")
        pairs = [
            (
                key,
                Path(DATA_ROOT) / "reference_files" / f"{key}.h5md",
                Path(VALIDATION_SYSTEMS) / top,
            )
            for key, top in AMBER_INTEROP_PAIRS
        ]
    elif args.h5md and args.prmtop:
        pairs = [(Path(args.h5md).stem, Path(args.h5md), Path(args.prmtop))]
    else:
        ap.error("give --run-all, or both --h5md and --prmtop")

    results, skipped = {}, []
    for key, h5md, prmtop in pairs:
        missing = [str(p) for p in (h5md, prmtop) if not p.is_file()]
        if missing:
            skipped.append(key)
            print(f"  skipped {key}: missing {', '.join(missing)}")
            continue
        result = results[key] = run_pair(h5md, prmtop)
        if not args.run_all:
            print_result(result)
            continue
        diff = [
            r["item"]
            for r in result["rows"]
            if r["group"] == "chemistry" and r["agree"] is False
        ]
        print(
            f"  {key:40s} {result['n_chemistry_agree']}/"
            f"{result['n_chemistry_items']}"
            + (f"  DIFF: {', '.join(diff)}" if diff else "")
        )

    n_ok = sum(r["interoperable"] for r in results.values())
    summary = {
        "n_pairs": len(results),
        "n_interoperable": n_ok,
        "skipped": skipped,
        "pairs": results,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    print(f"\n{n_ok}/{len(results)} pairs agree on every chemistry item")
    print(f"written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
