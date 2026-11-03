#!/usr/bin/env python
"""
Cross-engine interoperability: the same system, two engines, one schema.

The interoperability claim is that engine-specific names denote the same
schema-level concepts. This script tests the topology half of that claim with
data, item by item, on systems that were set up in AMBER.

Each pair in ``config.matrix.AMBER_INTEROP_PAIRS`` is a deposited AMBER topology
(the hexABC and MISATO groups of the validation dataset) and the H5MD file
GROMACS writes for the same system after ParmEd converted that topology:

    AMBER     deposited .prmtop                              -> read with pytraj
    GROMACS   .prmtop -> ParmEd -> topol.top -> .tpr -> H5MD biophysics file

The H5MD side is rebuilt to the full system from ``topology/molecule_blocks``.
On the AMBER side the H-H "bond" of rigid three-site water is dropped: it keeps
the water rigid under SHAKE and is not chemistry. Both sides are reduced to the
schema's own vocabulary and compared in two groups:

  * CHEMISTRY -- particle and component counts, residue sequence, molecular
    formula, atom names, atomic numbers and bond topology. These are properties
    of the molecule, so the two must agree.

  * PARAMETRISATION -- total charge, partial charges and masses, and every
    force-field term: bonds, angles, proper and improper dihedrals, 1-4 pairs
    and the Lennard-Jones type matrix. The same force field sits on both sides,
    so these agree too; what is tested is that the schema's elements carry the
    AMBER values without loss. For the terms the prmtop is read directly rather
    than through ParmEd, which did the conversion under test, and AMBER units
    and conventions are converted (kcal -> kJ, A -> nm, rad -> degree; bond and
    angle force constants x 2 for E = k/2 (x - x0)^2; 1-4 LJ / SCNB, Coulomb
    1 / SCEE). Compared per term at a relative tolerance of 1e-4:

      bonds       -> bond, or the distance of a constraint or SETTLE (water
                     O-H and H-H)
      angles      -> angle
      dihedrals   -> proper_dihedral / periodic_improper_dihedral, term by term;
                     a zero-force-constant AMBER term carries no energy, and
                     grompp removes it from the run input
      1-4 pairs   -> lennard_jones_14 and coulomb_14
      LJ A, B     -> the nonbonded T x T matrix; a 10-12 pair with zero
                     coefficients (water OW-HW in tleap) means no interaction

    python scripts/validation6_amber_interop.py -a     # every AMBER_INTEROP_PAIRS entry
    python scripts/validation6_amber_interop.py --h5md <file.h5md> --prmtop <file.prmtop>
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import pathlib
import re
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

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
        "residue_count": top.n_residues,
        "particle_names": [a.name.strip() for a in atoms],
        "species": [int(a.atomic_number) for a in atoms],
        "sequence_indices": [int(a.resid) for a in atoms],
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
        "residue_count": 0,
        "particle_names": [],
        "species": [],
        "sequence_indices": [],
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
            out["residue_count"] += int(g("residue_count"))
            out["particle_names"] += list(g("particle_names"))
            out["species"] += [int(v) for v in g("species")]
            out["sequence_indices"] += [int(v) + offset for v in g("sequence_indices")]
            out["component_names"] += list(g("component_names"))
            out["charges"] += [float(v) for v in g("charges")]
            out["masses"] += [float(v) for v in g("masses")]
            offset += int(g("residue_count"))

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
        "residue count",
        gmx["residue_count"],
        amb["residue_count"],
        gmx["residue_count"] == amb["residue_count"],
    )

    seq_g = [s.strip().lstrip("D").rstrip("35") for s in gmx["component_names"]]
    seq_a = [s.strip().lstrip("D").rstrip("35") for s in amb["component_names"]]
    add(
        "chemistry",
        "residue sequence (base identity)",
        "".join(seq_g),
        "".join(seq_a),
        seq_g == seq_a,
        "compared on residue identity, with terminal variants (e.g. DA5/DA3) "
        "folded into the base name",
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
        "net charge is set by the chemistry and must agree",
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
            "the same force field sits on both sides, so any difference is a "
            "loss in the conversion or in the writer",
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
#: descriptions of the same molecule.
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
        for i, (cid, name) in enumerate(
            zip(d["sequence_indices"], d["particle_names"])
        ):
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


# ==========================================================================
# Force-field terms
# ==========================================================================
KCAL = 4.184  # kJ per kcal
FF_RTOL = 1e-4


def read_prmtop(path) -> dict:
    """Every %FLAG section of a prmtop as a flat list of values."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    sections, flag, fmt, buf = {}, None, None, []
    for line in lines:
        if line.startswith("%FLAG"):
            if flag:
                sections[flag] = (fmt, buf)
            flag, fmt, buf = line.split()[1], None, []
        elif line.startswith("%FORMAT"):
            fmt = re.search(r"\((.*)\)", line).group(1)
        elif not line.startswith("%"):
            buf.append(line)
    if flag:
        sections[flag] = (fmt, buf)

    out = {}
    for flag, (fmt, body) in sections.items():
        m = re.match(r"(\d+)([aAiIeEfF])(\d+)", fmt or "")
        if not m:
            continue
        kind, width = m.group(2).upper(), int(m.group(3))
        conv = {"I": int, "E": float, "F": float}.get(kind, str.strip)
        out[flag] = [
            conv(line[i : i + width])
            for line in body
            for i in range(0, len(line), width)
            if line[i : i + width].strip()
        ]
    return out


def amber_terms(prmtop) -> dict:
    """Force-field terms of a prmtop, per atom index, in AMBER units."""
    p = read_prmtop(prmtop)
    natom, ntypes = p["POINTERS"][0], p["POINTERS"][1]
    nbidx = np.array(p["NONBONDED_PARM_INDEX"]).reshape(ntypes, ntypes)
    acoef, bcoef = p["LENNARD_JONES_ACOEF"], p["LENNARD_JONES_BCOEF"]
    hb_a, hb_b = p.get("HBOND_ACOEF", []), p.get("HBOND_BCOEF", [])

    def lj(a, b):
        """(A, B) of a type pair. A negative index points to a 10-12 term,
        which tleap uses with zero coefficients to switch a pair off (water
        OW-HW): no interaction. A non-zero 10-12 term -> None."""
        k = nbidx[a, b]
        if k > 0:
            return acoef[k - 1], bcoef[k - 1]
        if hb_a[-k - 1] == 0 and hb_b[-k - 1] == 0:
            return 0.0, 0.0
        return None

    ndih = len(p["DIHEDRAL_FORCE_CONSTANT"])
    scee = p.get("SCEE_SCALE_FACTOR") or [1.2] * ndih
    scnb = p.get("SCNB_SCALE_FACTOR") or [2.0] * ndih

    def entries(flags, width):
        for flag in flags:
            v = p.get(flag, [])
            for n in range(0, len(v), width):
                yield v[n : n + width]

    bonds = [
        (
            (i // 3, j // 3),
            p["BOND_FORCE_CONSTANT"][t - 1],
            p["BOND_EQUIL_VALUE"][t - 1],
        )
        for i, j, t in entries(("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"), 3)
    ]
    angles = [
        (
            (i // 3, j // 3, k // 3),
            p["ANGLE_FORCE_CONSTANT"][t - 1],
            p["ANGLE_EQUIL_VALUE"][t - 1],
        )
        for i, j, k, t in entries(("ANGLES_INC_HYDROGEN", "ANGLES_WITHOUT_HYDROGEN"), 4)
    ]
    # A negative third index: no 1-4 pair for this term; a negative fourth
    # index: an improper.
    dihedrals = [
        {
            "atoms": (i // 3, j // 3, abs(k) // 3, abs(m) // 3),
            "ignore_end": k < 0,
            "improper": m < 0,
            "pk": p["DIHEDRAL_FORCE_CONSTANT"][t - 1],
            "per": p["DIHEDRAL_PERIODICITY"][t - 1],
            "phase": p["DIHEDRAL_PHASE"][t - 1],
            "scee": scee[t - 1],
            "scnb": scnb[t - 1],
        }
        for i, j, k, m, t in entries(
            ("DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"), 5
        )
    ]
    return {
        "natom": natom,
        "names": p["ATOM_NAME"],
        "types": [t - 1 for t in p["ATOM_TYPE_INDEX"]],
        "lj": lj,
        "bonds": bonds,
        "angles": angles,
        "dihedrals": dihedrals,
    }


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def h5md_terms(path) -> dict:
    """Interaction terms of the whole system, with their parameter rows."""
    with h5py.File(path, "r") as f:
        tables = {}
        for cls in ("bonded", "nonbonded"):
            for t, g in f["force_field"].get(cls, {}).items():
                if "parameters" in g:
                    tables[t] = g["parameters"][()]
        topo = f["topology"]
        terms = defaultdict(list)  # type -> [(global atoms, parameter row)]
        names, type_id, offset = [], [], 0
        for row in topo["molecule_blocks"][()]:
            mol = topo[_s(row["name"])]
            pdata = mol["particle_data"][()]
            ntab = [_s(x) for x in mol["particle_name_table"][()]]
            inter = mol.get("interactions", {})
            local = {
                t: (g["particle_ids"][()], g["parameter_id"][()])
                for t, g in inter.items()
                if "particle_ids" in g
            }
            if "coulomb_14" in inter:  # its rows pair with lennard_jones_14's
                local["coulomb_14"] = (
                    local["lennard_jones_14"][0],
                    inter["coulomb_14"]["parameter_id"][()],
                )
            tid = (
                inter["lennard_jones_short_range"]["type_id"][()]
                if "lennard_jones_short_range" in inter
                else np.zeros(len(pdata), int)
            )
            for _ in range(int(row["count"])):
                for t, (ids, pid) in local.items():
                    terms[t] += [
                        (tuple(int(a) + offset for a in r), tables[t][q])
                        for r, q in zip(ids, pid)
                    ]
                names += [ntab[int(x)] for x in pdata["name"]]
                type_id += [int(x) for x in tid]
                offset += len(pdata)
        lj = tables.get("lennard_jones_short_range")
    T = math.isqrt(len(lj)) if lj is not None else 0
    return {
        "natom": offset,
        "names": names,
        "type_id": type_id,
        "terms": terms,
        "lj_matrix": lj.reshape(T, T) if lj is not None else None,
    }


class _Tally:
    """Term counts and the worst relative difference per compared field."""

    def __init__(self):
        self.n_amber = self.n_file = self.mismatched = self.missing = self.extra = 0
        self.worst = defaultdict(float)

    def value(self, field, a, b, atol=1e-6):
        # Both in double: with NumPy 2, a Python float minus a float32 is a float32.
        a, b = float(a), float(b)
        d = abs(a - b)
        scale = max(abs(a), abs(b))
        self.worst[field] = max(self.worst[field], d / scale if scale else 0.0)
        return d <= atol or d <= FF_RTOL * scale

    def record(self, ok):
        self.mismatched += not ok

    @property
    def agree(self):
        return self.mismatched == self.missing == self.extra == 0


def _canon(atoms):
    return min(tuple(atoms), tuple(reversed(atoms)))


def compare_forcefield(amb: dict, h5: dict) -> dict:
    """Every force-field term of the two sides, by category."""
    terms, out = h5["terms"], {}

    # Bonds: a bond, or the distance of a constraint or of a SETTLE.
    bond = {_canon(a): p for a, p in terms.get("bond", [])}
    constr = {_canon(a): p for a, p in terms.get("constraint", [])}
    settle = {}
    for (o, h1, h2), p in terms.get("settle", []):
        settle[_canon((o, h1))] = settle[_canon((o, h2))] = p["distance_OH"]
        settle[_canon((h1, h2))] = p["distance_HH"]
    t = _Tally()
    t.n_amber, t.n_file = len(amb["bonds"]), len(bond) + len(constr) + len(settle)
    seen = set()
    for atoms, k, r0 in amb["bonds"]:
        key = _canon(atoms)
        seen.add(key)
        if key in bond:
            p = bond[key]
            ok = t.value("k", 2 * k * KCAL * 100, p["force_constant"])
            ok &= t.value("r0", r0 * 0.1, p["equilibrium_distance"])
        elif key in constr:
            ok = t.value("r0", r0 * 0.1, constr[key]["distance"])
        elif key in settle:
            ok = t.value("r0", r0 * 0.1, settle[key])
        else:
            t.missing += 1
            continue
        t.record(ok)
    t.extra = len((set(bond) | set(constr) | set(settle)) - seen)
    out["bonds"] = (t, {"as constraint": len(constr), "as SETTLE": len(settle)})

    # Angles.
    ang = {_canon(a): p for a, p in terms.get("angle", [])}
    t = _Tally()
    t.n_amber, t.n_file = len(amb["angles"]), len(ang)
    seen = set()
    for atoms, k, th in amb["angles"]:
        key = _canon(atoms)
        seen.add(key)
        if key not in ang:
            t.missing += 1
            continue
        p = ang[key]
        ok = t.value("k", 2 * k * KCAL, p["force_constant"])
        ok &= t.value("theta0", math.degrees(th), p["equilibrium_angle"])
        t.record(ok)
    t.extra = len(set(ang) - seen)
    out["angles"] = (t, {})

    # Dihedrals, term by term.
    for label, improper, ftype in (
        ("proper dihedrals", False, "proper_dihedral"),
        ("improper dihedrals", True, "periodic_improper_dihedral"),
    ):
        fdih = defaultdict(list)
        for a, p in terms.get(ftype, []):
            fdih[_canon(a)].append(
                (int(p["multiplicity"]), float(p["force_constant"]), float(p["phase"]))
            )
        adih, n_zero = defaultdict(list), 0
        for d in amb["dihedrals"]:
            if d["improper"] != improper:
                continue
            if d["pk"] == 0:  # no energy; grompp removes it from the run input
                n_zero += 1
                continue
            adih[_canon(d["atoms"])].append(
                (round(d["per"]), d["pk"] * KCAL, math.degrees(d["phase"]))
            )
        # An improper written with its last two atoms swapped has the opposite
        # sign of angle; with phase 0 or 180 the energy is the same, so it is
        # the same term. Only then is the swapped order accepted.
        n_swapped = 0
        for key in [k for k in adih if improper and k not in fdih]:
            alt = _canon((key[0], key[1], key[3], key[2]))
            if (
                alt in fdih
                and alt not in adih
                and all(round(ph % 180, 3) in (0.0, 180.0) for _, _, ph in adih[key])
            ):
                adih[alt] = adih.pop(key)
                n_swapped += 1
        t = _Tally()
        t.n_amber = sum(map(len, adih.values()))
        t.n_file = sum(map(len, fdih.values()))
        for key, aterms in adih.items():
            fterms = sorted(fdih.get(key, []))
            for at in sorted(aterms):
                match = next((ft for ft in fterms if ft[0] == at[0]), None)
                if match is None:
                    t.missing += 1
                    continue
                fterms.remove(match)
                ok = t.value("k", at[1], match[1])
                ok &= t.value("phase", at[2] % 360, match[2] % 360, atol=1e-3)
                t.record(ok)
            t.extra += len(fterms)
        t.extra += sum(len(v) for k, v in fdih.items() if k not in adih)
        out[label] = (
            t,
            {"zero-force-constant terms": n_zero, "atom order swapped": n_swapped},
        )

    # 1-4 pairs, from the AMBER dihedrals that carry one.
    pairs = {}
    for d in amb["dihedrals"]:
        if not d["ignore_end"]:
            i, _, _, m = d["atoms"]
            pairs.setdefault(_canon((i, m)), (d["scee"], d["scnb"]))
    f14 = {_canon(a): p for a, p in terms.get("lennard_jones_14", [])}
    c14 = {_canon(a): p for a, p in terms.get("coulomb_14", [])}
    t = _Tally()
    t.n_amber, t.n_file = len(pairs), len(f14)
    for key, (scee, scnb) in pairs.items():
        if key not in f14:
            t.missing += 1
            continue
        ab = amb["lj"](amb["types"][key[0]], amb["types"][key[1]])
        ok = ab is not None
        if ok:
            ok &= t.value("c6", ab[1] / scnb * KCAL * 1e-6, f14[key]["c6"])
            ok &= t.value("c12", ab[0] / scnb * KCAL * 1e-12, f14[key]["c12"])
        if key in c14:
            ok &= t.value("scale", 1 / scee, c14[key]["scale_factor"])
        t.record(ok)
    t.extra = len(set(f14) - set(pairs))
    out["1-4 pairs"] = (t, {})

    # Lennard-Jones matrix, over the (AMBER type, file type) combinations present.
    M = h5["lj_matrix"]
    combos = sorted(set(zip(amb["types"], h5["type_id"])))
    t = _Tally()
    t.n_amber = t.n_file = len(combos) ** 2
    for a1, t1 in combos:
        for a2, t2 in combos:
            ab = amb["lj"](a1, a2)
            if ab is None:
                t.missing += 1
                continue
            ok = t.value("c6", ab[1] * KCAL * 1e-6, M[t1, t2]["c6"])
            ok &= t.value("c12", ab[0] * KCAL * 1e-12, M[t1, t2]["c12"])
            t.record(ok)
    out["Lennard-Jones matrix"] = (t, {})
    return out


def forcefield_rows(h5md, prmtop) -> list:
    """The force-field terms as parametrisation rows of the comparison."""
    amb, h5 = amber_terms(prmtop), h5md_terms(h5md)
    if amb["natom"] != h5["natom"] or amb["names"] != h5["names"]:
        return [
            {
                "group": "parametrisation",
                "item": "force-field terms",
                "gromacs": h5["natom"],
                "amber": amb["natom"],
                "agree": False,
                "note": "particle order differs; terms are compared per atom index",
            }
        ]
    rows = []
    for label, (t, extra) in compare_forcefield(amb, h5).items():
        worst = ", ".join(f"{k} {v:.1e}" for k, v in t.worst.items())
        notes = [f"worst relative difference: {worst}"] if worst else []
        notes += [f"{v} {k}" for k, v in extra.items() if v]
        notes += [
            f"{n} {what}"
            for n, what in (
                (t.mismatched, "differ"),
                (t.missing, "missing in file"),
                (t.extra, "only in file"),
            )
            if n
        ]
        rows.append(
            {
                "group": "parametrisation",
                "item": f"{label} (terms)",
                "gromacs": t.n_file,
                "amber": t.n_amber,
                "agree": t.agree,
                "note": "; ".join(notes),
                "worst_relative_difference": dict(t.worst),
            }
        )
    return rows


def run_pair(h5md, prmtop) -> dict:
    gmx = h5md_topology(h5md)
    amb = amber_topology(prmtop)
    result = compare(gmx, amb)
    result["rows"] += forcefield_rows(h5md, prmtop)
    param = [
        r
        for r in result["rows"]
        if r["group"] == "parametrisation" and r["agree"] is not None
    ]
    result["n_parametrisation_items"] = len(param)
    result["n_parametrisation_agree"] = sum(1 for r in param if r["agree"])
    result["interoperable"] = result["interoperable"] and all(r["agree"] for r in param)
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
        f"{result['n_chemistry_agree']}/{result['n_chemistry_items']}, "
        f"parametrisation items: "
        f"{result['n_parametrisation_agree']}/{result['n_parametrisation_items']}"
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
        diff = [r["item"] for r in result["rows"] if r["agree"] is False]
        print(
            f"  {key:40s} {result['n_chemistry_agree']}/"
            f"{result['n_chemistry_items']} "
            f"{result['n_parametrisation_agree']}/"
            f"{result['n_parametrisation_items']}"
            + (f"  DIFF: {', '.join(diff)}" if diff else ""),
            flush=True,
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
    print(
        f"\n{n_ok}/{len(results)} pairs agree on every chemistry and "
        "parametrisation item"
    )
    print(f"written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
