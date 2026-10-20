"""
Round-trip fidelity test: does the metadata come back unchanged?

    GROMACS run input (.tpr)
        -> gmx dump            [source side,  h5md_validator.gmxdump]
        -> serialize to H5MD   [GROMACS reference implementation]
        -> parse H5MD          [target side,  this module]
        -> compare

Both sides reduce to the same flat vocabulary of metadata items, so every item
is matched, mismatched, or absent on one side. That same reduction gives the
coverage matrix: an item in the source and not in the file is a coverage gap.

Numerical comparison uses a per-value tolerance rather than one global epsilon,
bounded by whichever is coarser -- ``gmx dump``'s field-dependent text
precision, or the file's float32 / float64 storage.

Used by:
    Validation 1: schema coverage
    Validation 2: round-trip test
    Validation 6: cross-engine mapping
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import h5py
import numpy as np

from . import vocabulary as vocab
from .gmxdump import MoleculeType, TprDump

#: float32 has ~7.2 decimal digits; leave a factor of a few for accumulated
#: conversion error in the writer.
_FLOAT32_EPS = 4 * float(np.finfo(np.float32).eps)


# --------------------------------------------------------------------------
@dataclass
class Value:
    kind: str  # categorical | integer | float | int_array | float_array | str_array | nested_int
    value: Any
    unit: str | None = None


@dataclass
class Comparison:
    key: str
    category: str
    kind: str
    status: str  # match | match_unordered | mismatch | missing_in_file | missing_in_source | absent_both
    source: Any = None
    target: Any = None
    detail: dict = field(default_factory=dict)

    def to_dict(self):
        d = asdict(self)
        for k in ("source", "target"):
            d[k] = _summarise(d[k])
        return d


def _summarise(v, limit=6):
    """Keep reports readable: long arrays are shown by head and length."""
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if isinstance(v, (list, tuple)):
        if len(v) > limit:
            return {"n": len(v), "head": [_summarise(x, 3) for x in v[:limit]]}
        return [_summarise(x, 3) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return float(v)
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    return v


def _decode(x):
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace")
    if isinstance(x, np.ndarray) and x.dtype.kind in ("S", "O"):
        return [_decode(v) for v in x.tolist()]
    if isinstance(x, np.generic):
        return x.item()
    return x


def _read(node):
    return _decode(node[()])


# ==========================================================================
# source side: the .tpr as GROMACS reports it
# ==========================================================================
def extract_source(d: TprDump) -> dict:
    ir = d.inputrec
    out: dict = {}

    def put(key, kind, value, unit=None):
        if value is not None:
            out[key] = Value(kind, value, unit)

    # ---------------- run control ----------------------------------------
    integrator = ir.get("integrator")
    simtype = vocab.SIMULATION_TYPE.get(integrator)
    put("parameters.run_control.simulation_type", "categorical", simtype)
    put("parameters.run_control.numsteps", "integer", ir.get("nsteps"))
    if simtype in ("molecular_dynamics", "stochastic_dynamics"):
        put(
            f"parameters.run_control.{simtype}.integrator",
            "categorical",
            vocab.INTEGRATOR.get(integrator),
        )
        put(f"parameters.run_control.{simtype}.timestep", "float", ir.get("dt"), "ps")
    elif simtype == "energy_minimization":
        put(
            "parameters.run_control.energy_minimization.algorithm",
            "categorical",
            vocab.INTEGRATOR.get(integrator),
        )
        put(
            "parameters.run_control.energy_minimization.step_size",
            "float",
            ir.get("emstep"),
            "nm",
        )
        put(
            "parameters.run_control.energy_minimization.tolerance",
            "float",
            ir.get("emtol"),
            "kJ mol-1 nm-1",
        )

    com_mode = vocab.COM_MODE.get(str(ir.get("comm-mode")))
    if com_mode is not None:
        put(
            "parameters.run_control.center_of_mass_motion_removal.mode",
            "categorical",
            com_mode,
        )
        put(
            "parameters.run_control.center_of_mass_motion_removal.interval",
            "integer",
            ir.get("nstcomm"),
        )
        put(
            "parameters.run_control.center_of_mass_motion_removal.groups",
            "str_array",
            d.groups.get("VCM", []),
        )
    factor = ir.get("mass-repartition-factor")
    if factor is not None and float(factor) != 1.0:
        put("parameters.run_control.mass_repartitioning.factor", "float", factor)

    # ---------------- output control -------------------------------------
    for key, mdp in (
        ("position", "nstxout"),
        ("velocity", "nstvout"),
        ("force", "nstfout"),
        ("energy", "nstenergy"),
        ("log", "nstlog"),
    ):
        put(f"parameters.output_control.{key}", "integer", ir.get(mdp))

    # ---------------- temperature coupling --------------------------------
    thermostat = vocab.THERMOSTAT.get(str(ir.get("tcoupl")))
    if thermostat:
        put("parameters.temperature_coupling.algorithm", "categorical", thermostat)
        put("parameters.temperature_coupling.interval", "integer", ir.get("nsttcouple"))
        put(
            "parameters.temperature_coupling.groups",
            "str_array",
            d.groups.get("T-Coupling", []),
        )
        put(
            "parameters.temperature_coupling.reference_temperature",
            "float_array",
            d.grpopts.get("ref-t"),
            "K",
        )
        tfield = vocab.THERMOSTAT_TIME_FIELD[thermostat]
        put(
            f"parameters.temperature_coupling.{thermostat}.{tfield}",
            "float_array",
            d.grpopts.get("tau-t"),
            "ps",
        )
        if thermostat == "nose_hoover":
            put(
                "parameters.temperature_coupling.nose_hoover.chain_length",
                "integer",
                ir.get("nh-chain-length"),
            )

    # ---------------- pressure coupling -----------------------------------
    barostat = vocab.BAROSTAT.get(str(ir.get("pcoupl")))
    if barostat:
        put("parameters.pressure_coupling.algorithm", "categorical", barostat)
        put("parameters.pressure_coupling.interval", "integer", ir.get("nstpcouple"))
        put(
            "parameters.pressure_coupling.reference_pressure",
            "float_array",
            _flat(ir.get("ref-p")),
            "10+5 Pa",
        )
        put(
            "parameters.pressure_coupling.compressibility",
            "float_array",
            _flat(ir.get("compressibility")),
            "10-5 Pa-1",
        )
        put(
            "parameters.pressure_coupling.reference_coordinate_scaling",
            "categorical",
            vocab.REFCOORD_SCALING.get(str(ir.get("refcoord-scaling"))),
        )
        pfield = vocab.BAROSTAT_TIME_FIELD[barostat]
        put(
            f"parameters.pressure_coupling.{barostat}.{pfield}",
            "float",
            ir.get("tau-p"),
            "ps",
        )

    # ---------------- coulomb ---------------------------------------------
    ctype = vocab.COULOMB_TYPE.get(str(ir.get("coulombtype")))
    put("parameters.coulomb_interactions.type", "categorical", ctype)
    put(
        "parameters.coulomb_interactions.cutoff_distance",
        "float",
        ir.get("rcoulomb"),
        "nm",
    )
    put(
        "parameters.coulomb_interactions.relative_permittivity",
        "float",
        ir.get("epsilon-r"),
    )
    if ctype == "cutoff":
        put(
            "parameters.coulomb_interactions.cutoff.potential_modifier",
            "categorical",
            vocab.POTENTIAL_MODIFIER.get(str(ir.get("coulomb-modifier"))),
        )
        if ir.get("rcoulomb-switch"):
            put(
                "parameters.coulomb_interactions.cutoff.switch_distance",
                "float",
                ir.get("rcoulomb-switch"),
                "nm",
            )
    if ctype == "reaction_field":
        eps_rf = ir.get("epsilon-rf")
        if isinstance(eps_rf, str):  # "inf" means infinite permittivity
            eps_rf = 0.0
        put(
            "parameters.coulomb_interactions.reaction_field.relative_permittivity",
            "float",
            eps_rf,
        )

    # ---------------- van der Waals ---------------------------------------
    vtype = vocab.VDW_TYPE.get(str(ir.get("vdw-type")))
    put("parameters.vdw_interactions.type", "categorical", vtype)
    put("parameters.vdw_interactions.cutoff_distance", "float", ir.get("rvdw"), "nm")
    put(
        "parameters.vdw_interactions.dispersion_correction",
        "categorical",
        vocab.DISPERSION_CORRECTION.get(str(ir.get("DispCorr"))),
    )
    if vtype == "cutoff":
        modifier = vocab.POTENTIAL_MODIFIER.get(str(ir.get("vdw-modifier")))
        put(
            "parameters.vdw_interactions.cutoff.potential_modifier",
            "categorical",
            modifier,
        )
        if modifier in ("potential_switch", "force_switch"):
            put(
                "parameters.vdw_interactions.cutoff.switch_distance",
                "float",
                ir.get("rvdw-switch"),
                "nm",
            )

    # ---------------- boundary conditions ---------------------------------
    pbc = str(ir.get("pbc", "")).lower()
    periodicity = {"xyz": [1, 1, 1], "xy": [1, 1, 0], "no": [0, 0, 0]}.get(pbc)
    put("parameters.boundary_conditions.periodicity", "int_array", periodicity)

    # ---------------- constraints -----------------------------------------
    algo = vocab.CONSTRAINT_ALGORITHM.get(str(ir.get("constraint-algorithm")))
    if algo and _has_constraints(d):
        put("parameters.constraints.algorithm", "categorical", algo)
        if algo == "lincs":
            put(
                "parameters.constraints.lincs.iterations",
                "integer",
                ir.get("lincs-iter"),
            )
            put("parameters.constraints.lincs.order", "integer", ir.get("lincs-order"))
        elif algo == "shake":
            put("parameters.constraints.shake.tolerance", "float", ir.get("shake-tol"))

    # ---------------- topology --------------------------------------------
    out.update(_source_topology(d))
    return out


def _flat(matrix):
    if matrix is None:
        return None
    if isinstance(matrix, list) and matrix and isinstance(matrix[0], list):
        return [float(v) for row in matrix for v in row]
    return [float(v) for v in np.atleast_1d(matrix)]


def _has_constraints(d: TprDump) -> bool:
    return any(k in ("CONSTR", "CONSTRNC") for m in d.moltypes for k in m.ilists)


def _source_topology(d: TprDump) -> dict:
    out: dict = {}
    out["topology.total_particle_count"] = Value("integer", d.n_atoms)

    total_residues = 0
    for block in d.molblocks:
        mol = d.moltypes[block["moltype"]]
        total_residues += block["nmol"] * mol.n_residues
    out["topology.total_residue_count"] = Value("integer", total_residues)
    out["topology.molecule_blocks.names"] = Value(
        "str_array", [b["name"] for b in d.molblocks]
    )
    out["topology.molecule_blocks.counts"] = Value(
        "int_array", [b["nmol"] for b in d.molblocks]
    )

    for mol in d.moltypes:
        p = f"topology.{mol.name}"
        out[f"{p}.particle_count"] = Value("integer", mol.n_atoms)
        out[f"{p}.residue_count"] = Value("integer", mol.n_residues)
        out[f"{p}.particle_names"] = Value("str_array", list(mol.atom_names))
        out[f"{p}.species"] = Value(
            "int_array", [a.get("atomnumber", -1) for a in mol.atoms]
        )
        out[f"{p}.charges"] = Value("float_array", [a.get("q") for a in mol.atoms], "e")
        out[f"{p}.masses"] = Value("float_array", [a.get("m") for a in mol.atoms], "u")
        out[f"{p}.sequence_indices"] = Value(
            "int_array", [a.get("resind", -1) for a in mol.atoms]
        )
        out[f"{p}.component_names"] = Value("str_array", [r[0] for r in mol.residues])
        out[f"{p}.residue_numbers"] = Value("int_array", [r[1] for r in mol.residues])
        # Each excluded pair once, (i, j) with i < j, as the file stores it. The
        # per-particle lists of gmx dump are symmetric and include the particle
        # itself, both implied by the pairs. Without a pair there is no
        # exclusions group in the file, so there is no entry here either.
        pairs = [[i, j] for i, e in enumerate(mol.excls) for j in sorted(e) if j > i]
        if pairs:
            out[f"{p}.exclusions"] = Value("nested_int", pairs)

        for schema_type, entries in _source_interactions(d, mol).items():
            ids, params = entries
            out[f"{p}.interactions.{schema_type}.particle_ids"] = Value(
                "int_array", ids
            )
            if params:
                out[f"{p}.interactions.{schema_type}.parameters"] = Value(
                    "float_array", params
                )
            # Else the dump parser cannot resolve this type's parameters -- CMAP
            # is the case: gmx dump prints the grids in their own section and
            # the tpr keeps every grid the force field defines, while the file
            # stores only the ones referenced, so the indices do not correspond.
            # Emitting an empty list would report a mismatch the comparison
            # cannot actually evaluate; leaving the key out reports the file's
            # values as missing_in_source, which is what they are.

    out["connectivity.bonds"] = Value("int_array", _source_bonds(d))
    return out


def _source_interactions(d: TprDump, mol: MoleculeType) -> dict:
    """Group a molecule type's interaction lists by schema interaction type.

    Returns ``{schema_type: (flattened particle ids, flattened resolved params)}``.
    Parameters are resolved through the GROMACS parameter index, so the test
    covers the index indirection end to end rather than only the tables.
    """
    grouped: dict = {}
    for functype, entries in sorted(mol.ilists.items()):
        schema_type = vocab.INTERACTION.get(functype)
        if schema_type is None or not entries:
            continue
        field_map = vocab.PARAMETER_FIELDS.get(functype, {})
        ids, params = grouped.setdefault(schema_type, ([], []))
        for param_index, atoms in entries:
            ids.extend(atoms)
            entry = d.functypes[param_index] if param_index < len(d.functypes) else None
            if entry is None:
                continue
            _, values = entry
            for gmx_key in field_map:
                if gmx_key in values:
                    params.append(values[gmx_key])
    return grouped


def _source_bonds(d: TprDump):
    """Expand molecule-local bond lists into the global 0-based bond list."""
    bonds = []
    offset = 0
    for block in d.molblocks:
        mol = d.moltypes[block["moltype"]]
        local = []
        for functype, entries in mol.ilists.items():
            if functype in vocab.BOND_LIKE:
                for _, atoms in entries:
                    if len(atoms) >= 2:
                        local.append((atoms[0], atoms[1]))
            elif functype == "SETTLE":
                # A rigid water: the two O-H bonds are chemical bonds.
                for _, atoms in entries:
                    if len(atoms) >= 3:
                        local.append((atoms[0], atoms[1]))
                        local.append((atoms[0], atoms[2]))
        for _ in range(block["nmol"]):
            for a, b in local:
                bonds.append((a + offset, b + offset))
            offset += mol.n_atoms
    return bonds


# ==========================================================================
# target side: the H5MD file
# ==========================================================================
def extract_target(path) -> dict:
    with h5py.File(path, "r") as f:
        out: dict = {}
        _target_parameters(f, out)
        _target_topology(f, out)
        return out


def _walk_parameters(group, prefix, out):
    for name in group:
        node = group[name]
        key = f"{prefix}.{name}"
        if isinstance(node, h5py.Group):
            _walk_parameters(node, key, out)
            continue
        value = _read(node)
        unit = _decode(node.attrs["unit"]) if "unit" in node.attrs else None
        if isinstance(value, str):
            out[key] = Value("categorical", value, unit)
        elif isinstance(value, list):
            if value and isinstance(value[0], str):
                out[key] = Value("str_array", value, unit)
            else:
                kind = "float_array" if node.dtype.kind == "f" else "int_array"
                out[key] = Value(kind, list(np.ravel(node[()])), unit)
        elif isinstance(node[()], np.ndarray) or node.shape not in ((), None):
            kind = "float_array" if node.dtype.kind == "f" else "int_array"
            out[key] = Value(kind, list(np.ravel(node[()])), unit)
        elif node.dtype.kind == "f":
            out[key] = Value("float", float(value), unit)
        else:
            out[key] = Value("integer", int(value), unit)


def _target_parameters(f, out):
    params = f.get("parameters")
    if params is not None:
        _walk_parameters(params, "parameters", out)


def _target_topology(f, out):
    topo = f.get("topology")
    if topo is None:
        return
    for attr, key in (
        ("total_particle_count", "topology.total_particle_count"),
        ("total_residue_count", "topology.total_residue_count"),
    ):
        if attr in topo.attrs:
            out[key] = Value("integer", int(np.atleast_1d(topo.attrs[attr])[0]))

    blocks = topo.get("molecule_blocks")
    if blocks is not None:
        data = blocks[()]
        out["topology.molecule_blocks.names"] = Value(
            "str_array", [_decode(v) for v in data["name"]]
        )
        out["topology.molecule_blocks.counts"] = Value(
            "int_array", [int(v) for v in data["count"]]
        )

    ff_params = _target_ff_parameters(f)

    for name in topo:
        mol = topo[name]
        if not isinstance(mol, h5py.Group) or "particle_data" not in mol:
            continue
        p = f"topology.{name}"
        for attr in ("particle_count", "residue_count"):
            if attr in mol.attrs:
                out[f"{p}.{attr}"] = Value(
                    "integer", int(np.atleast_1d(mol.attrs[attr])[0])
                )

        pdata = mol["particle_data"][()]
        pnames = (
            [_decode(v) for v in mol["particle_name_table"][()]]
            if "particle_name_table" in mol
            else []
        )
        name_field = pdata["name"]
        if name_field.dtype.kind in ("i", "u") and pnames:
            resolved = [pnames[i] for i in name_field]
        else:
            resolved = [_decode(v) for v in name_field]
        out[f"{p}.particle_names"] = Value("str_array", resolved)
        out[f"{p}.species"] = Value("int_array", [int(v) for v in pdata["species"]])
        units = [
            _decode(u)
            for u in np.atleast_1d(mol["particle_data"].attrs.get("units", []))
        ]
        out[f"{p}.charges"] = Value(
            "float_array",
            [float(v) for v in pdata["charge"]],
            units[2] if len(units) > 2 else None,
        )
        out[f"{p}.masses"] = Value(
            "float_array",
            [float(v) for v in pdata["mass"]],
            units[3] if len(units) > 3 else None,
        )
        out[f"{p}.sequence_indices"] = Value(
            "int_array", [int(v) for v in pdata["sequence_index"]]
        )

        if "sequence" in mol:
            seq = mol["sequence"][()]
            cnames = (
                [_decode(v) for v in mol["component_name_table"][()]]
                if "component_name_table" in mol
                else []
            )
            idx = seq["component_name"]
            if idx.dtype.kind in ("i", "u") and cnames:
                out[f"{p}.component_names"] = Value(
                    "str_array", [cnames[i] for i in idx]
                )
            else:
                out[f"{p}.component_names"] = Value(
                    "str_array", [_decode(v) for v in idx]
                )
            out[f"{p}.residue_numbers"] = Value(
                "int_array", [int(v) for v in seq["residue_number"]]
            )

        if "exclusions" in mol:
            out[f"{p}.exclusions"] = Value(
                "nested_int", _decode_exclusions(mol["exclusions"])
            )

        inter = mol.get("interactions")
        if isinstance(inter, h5py.Group):
            for tname in sorted(inter):
                grp = inter[tname]
                if "particle_ids" not in grp or "parameter_id" not in grp:
                    continue
                ids = np.asarray(grp["particle_ids"][()], dtype=np.int64)
                out[f"{p}.interactions.{tname}.particle_ids"] = Value(
                    "int_array", ids.ravel().tolist()
                )
                table = ff_params.get(tname)
                if table is not None:
                    pid = np.asarray(grp["parameter_id"][()], dtype=np.int64).ravel()
                    resolved = []
                    if table.dtype.names:
                        fields = list(table.dtype.names)
                        for i in pid:
                            if 0 <= i < table.shape[0]:
                                resolved.extend(float(table[i][fl]) for fl in fields)
                    else:
                        # A tabulated potential: the parameter table is a plain
                        # array of grids, one per type, not a compound of named
                        # fields. CMAP is the case in hand -- (n_types, 24, 24)
                        # of energies. Flatten the referenced grid so the values
                        # are compared like any other parameter set.
                        for i in pid:
                            if 0 <= i < table.shape[0]:
                                resolved.extend(
                                    float(x) for x in np.asarray(table[i]).ravel()
                                )
                    out[f"{p}.interactions.{tname}.parameters"] = Value(
                        "float_array", resolved
                    )

    conn = f.get("connectivity/bonds")
    if isinstance(conn, h5py.Dataset):
        data = np.asarray(conn[()], dtype=np.int64)
        out["connectivity.bonds"] = Value(
            "int_array", [tuple(int(x) for x in row) for row in data]
        )


def _decode_exclusions(grp):
    """The excluded pairs, one [i, j] per row of ``particle_ids``."""
    if "particle_ids" not in grp:
        return []
    pairs = np.atleast_2d(np.asarray(grp["particle_ids"][()], dtype=np.int64))
    return pairs.tolist()


def _target_ff_parameters(f):
    """interaction type -> the parameters compound dataset, wherever it lives."""
    tables = {}
    ff = f.get("force_field")
    if ff is None:
        return tables
    for layer in ("bonded", "nonbonded"):
        grp = ff.get(layer)
        if isinstance(grp, h5py.Group):
            for t in grp:
                if "parameters" in grp[t]:
                    tables[t] = grp[t]["parameters"][()]
    for t in ff:
        if t in ("bonded", "nonbonded"):
            continue
        node = ff[t]
        if isinstance(node, h5py.Group) and "parameters" in node:
            tables[t] = node["parameters"][()]
    return tables


# ==========================================================================
# observables: the .edr as the source, the file's observables group as target
# ==========================================================================
def extract_source_observables(edr) -> dict:
    """Energy terms from the .edr, keyed as the schema names them.

    The .edr is written by the same mdrun as the H5MD file, so this compares two
    writers over one set of values rather than two independent readers of a
    stored file -- weaker evidence than the .tpr round-trip, and enough to catch
    a wrong unit, a mis-mapped term or frames out of step.
    """
    out: dict = {}
    if not edr.steps:
        return out
    out["observables.step"] = Value("int_array", list(edr.steps))
    for term, element in vocab.OBSERVABLE.items():
        series = edr.terms.get(term)
        if series and None not in series:
            out[f"observables.{element}"] = Value("float_array", list(series))
    for element, pattern in vocab.OBSERVABLE_TENSOR.items():
        rows = []
        for i in vocab.TENSOR_AXES:
            for j in vocab.TENSOR_AXES:
                series = edr.terms.get(pattern.format(i, j))
                if not series or None in series:
                    rows = None
                    break
                rows.append(series)
            if rows is None:
                break
        if rows:
            # Frame-major, to match the file's (n_frames, 3, 3).
            out[f"observables.{element}"] = Value(
                "float_array",
                [rows[c][f] for f in range(edr.n_frames) for c in range(9)],
            )
    return out


def _target_observables(path, out):
    """Read the observables group, flattened to match the .edr side."""
    with h5py.File(path, "r") as f:
        obs = f.get("observables")
        if obs is None:
            return
        steps = None
        for name in obs:
            node = obs[name]
            if not isinstance(node, h5py.Group) or "value" not in node:
                continue  # per-group subgroups (observables/<group>/...)
            if (
                name not in vocab.OBSERVABLE.values()
                and name not in vocab.OBSERVABLE_TENSOR
            ):
                continue
            value = np.asarray(node["value"][()])
            out[f"observables.{name}"] = Value(
                "float_array", [float(x) for x in value.ravel()]
            )
            if steps is None and "step" in node:
                steps = [int(s) for s in np.atleast_1d(node["step"][()])]
        if steps is not None:
            out["observables.step"] = Value("int_array", steps)


# ==========================================================================
# comparison
# ==========================================================================
def _category_of(key: str) -> str:
    head = key.split(".")[0]
    if head == "parameters":
        return "md_parameters"
    if head == "connectivity":
        return "connectivity"
    if head == "observables":
        return "observables"
    if head == "topology":
        if ".outside_predefined." in key:
            return "outside_predefined"
        if ".interactions." in key:
            return "force_field"
        return "topology"
    return head


def _tolerance(source_value) -> float:
    """Relative tolerance for one float, from printed and stored precision."""
    sig = getattr(source_value, "sig_digits", None)
    eps_print = 10.0 ** (-(sig - 1)) if sig else 0.0
    return max(eps_print, _FLOAT32_EPS)


def _floats_equal(a, b, tol) -> bool:
    a, b = float(a), float(b)
    if a == b:
        return True
    return abs(a - b) <= tol * max(abs(a), abs(b))


def compare(source: dict, target: dict, *, ignore=()) -> list:
    """Compare the two extractions item by item."""
    results = []
    for key in sorted(set(source) | set(target)):
        if key in ignore:
            continue
        src, tgt = source.get(key), target.get(key)
        category = _category_of(key)
        if src is None:
            results.append(
                Comparison(
                    key, category, tgt.kind, "missing_in_source", None, tgt.value
                )
            )
            continue
        if tgt is None:
            results.append(
                Comparison(key, category, src.kind, "missing_in_file", src.value, None)
            )
            continue
        results.append(_compare_one(key, category, src, tgt))
    return results


def _compare_one(key, category, src: Value, tgt: Value) -> Comparison:
    kind = src.kind
    detail = {}
    if src.unit and tgt.unit and src.unit != tgt.unit:
        detail["unit_source"] = src.unit
        detail["unit_file"] = tgt.unit

    if kind == "categorical":
        ok = str(src.value) == str(tgt.value)
        return Comparison(
            key,
            category,
            kind,
            "match" if ok else "mismatch",
            src.value,
            tgt.value,
            detail,
        )

    if kind == "integer":
        ok = int(src.value) == int(tgt.value)
        return Comparison(
            key,
            category,
            kind,
            "match" if ok else "mismatch",
            src.value,
            tgt.value,
            detail,
        )

    if kind == "float":
        tol = _tolerance(src.value)
        ok = _floats_equal(src.value, tgt.value, tol)
        detail["tolerance"] = tol
        if not ok:
            detail["abs_diff"] = abs(float(src.value) - float(tgt.value))
        return Comparison(
            key,
            category,
            kind,
            "match" if ok else "mismatch",
            src.value,
            tgt.value,
            detail,
        )

    if kind == "str_array":
        a = [str(x) for x in src.value]
        b = [str(x) for x in tgt.value]
        if a == b:
            return Comparison(key, category, kind, "match", a, b, detail)
        if sorted(a) == sorted(b):
            return Comparison(key, category, kind, "match_unordered", a, b, detail)
        detail.update(_first_difference(a, b))
        return Comparison(key, category, kind, "mismatch", a, b, detail)

    if kind == "int_array":
        a = _as_int_list(src.value)
        b = _as_int_list(tgt.value)
        if a == b:
            return Comparison(key, category, kind, "match", a, b, detail)
        if sorted(a) == sorted(b):
            return Comparison(key, category, kind, "match_unordered", a, b, detail)
        detail.update(_first_difference(a, b))
        return Comparison(key, category, kind, "mismatch", a, b, detail)

    if kind == "float_array":
        a = list(src.value)
        b = [float(x) for x in tgt.value]
        if len(a) != len(b):
            detail["n_source"], detail["n_file"] = len(a), len(b)
            return Comparison(key, category, kind, "mismatch", a, b, detail)
        tol = max((_tolerance(x) for x in a), default=_FLOAT32_EPS)
        detail["tolerance"] = tol
        bad = [i for i, (x, y) in enumerate(zip(a, b)) if not _floats_equal(x, y, tol)]
        if not bad:
            return Comparison(key, category, kind, "match", a, b, detail)
        detail["n_outside_tolerance"] = len(bad)
        detail["first_offsets"] = bad[:8]
        detail["examples"] = [
            {"i": i, "source": float(a[i]), "file": float(b[i])} for i in bad[:4]
        ]
        return Comparison(key, category, kind, "mismatch", a, b, detail)

    if kind == "nested_int":
        a = [list(map(int, x)) for x in src.value]
        b = [list(map(int, x)) for x in tgt.value]
        if a == b:
            return Comparison(key, category, kind, "match", len(a), len(b), detail)
        detail["n_source"], detail["n_file"] = len(a), len(b)
        bad = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        detail["n_differing_rows"] = len(bad)
        if bad:
            i = bad[0]
            detail["first_difference"] = {
                "row": i,
                "source": a[i][:12],
                "file": b[i][:12],
            }
        return Comparison(key, category, kind, "mismatch", len(a), len(b), detail)

    return Comparison(key, category, kind, "mismatch", src.value, tgt.value, detail)


def _as_int_list(values):
    out = []
    for v in values:
        if isinstance(v, (tuple, list)):
            out.append(tuple(int(x) for x in v))
        else:
            out.append(int(v))
    return out


def _first_difference(a, b):
    if len(a) != len(b):
        return {"n_source": len(a), "n_file": len(b)}
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return {"first_difference": {"index": i, "source": x, "file": y}}
    return {}


# ==========================================================================
@dataclass
class RoundTripReport:
    file: str
    tpr: str
    comparisons: list = field(default_factory=list)

    def _by(self, status):
        return [c for c in self.comparisons if c.status == status]

    @property
    def matches(self):
        return [c for c in self.comparisons if c.status in ("match", "match_unordered")]

    @property
    def mismatches(self):
        return self._by("mismatch")

    @property
    def missing_in_file(self):
        """Items of a predefined kind that the source had and the file does not."""
        return [
            c for c in self._by("missing_in_file") if c.category != "outside_predefined"
        ]

    @property
    def outside_predefined(self):
        """Source terms whose functional form is outside Table SI.

        Their presence here is not itself a finding -- the predefined set is
        deliberately not exhaustive. What the comparison shows is whether each
        one reached the file at all.
        """
        return [c for c in self.comparisons if c.category == "outside_predefined"]

    @property
    def dropped_terms(self):
        """Terms outside the predefined set that reached the file nowhere."""
        return [c for c in self.outside_predefined if c.status == "missing_in_file"]

    @property
    def extra_in_file(self):
        return self._by("missing_in_source")

    @property
    def faithful(self) -> bool:
        """Fidelity: everything the file does carry came back unchanged.

        Deliberately independent of what happens to terms outside the
        predefined set, which is reported by ``all_terms_written`` and by
        scripts/validation1_interaction_coverage.py.
        """
        return not self.mismatches and not self.missing_in_file

    @property
    def all_terms_written(self) -> bool:
        """Every source term reached the file, predefined or not."""
        return not self.dropped_terms

    def counts_by_category(self):
        out = {}
        for c in self.comparisons:
            e = out.setdefault(
                c.category,
                {
                    "compared": 0,
                    "match": 0,
                    "mismatch": 0,
                    "missing_in_file": 0,
                    "extra_in_file": 0,
                },
            )
            if c.status in ("match", "match_unordered"):
                e["compared"] += 1
                e["match"] += 1
            elif c.status == "mismatch":
                e["compared"] += 1
                e["mismatch"] += 1
            elif c.status == "missing_in_file":
                e["missing_in_file"] += 1
            elif c.status == "missing_in_source":
                e["extra_in_file"] += 1
        return out

    def to_dict(self):
        return {
            "file": self.file,
            "tpr": self.tpr,
            "faithful": self.faithful,
            "all_terms_written": self.all_terms_written,
            "n_items": len(self.comparisons),
            "n_match": len(self.matches),
            "n_mismatch": len(self.mismatches),
            "n_missing_in_file": len(self.missing_in_file),
            "n_outside_predefined": len(self.outside_predefined),
            "n_dropped_terms": len(self.dropped_terms),
            "dropped_terms": [
                {"term": c.key.split(".")[-1], "entries": c.source}
                for c in self.dropped_terms
            ],
            "n_extra_in_file": len(self.extra_in_file),
            "by_category": self.counts_by_category(),
            "comparisons": [c.to_dict() for c in self.comparisons],
        }

    def summary(self):
        verdict = "FAITHFUL" if self.faithful else "DIVERGENT"
        gap = "" if self.all_terms_written else f" dropped={len(self.dropped_terms)}"
        return (
            f"{verdict}  {self.file}  items={len(self.comparisons)} "
            f"match={len(self.matches)} mismatch={len(self.mismatches)} "
            f"missing={len(self.missing_in_file)} "
            f"extra={len(self.extra_in_file)}{gap}"
        )


def roundtrip(
    dump: TprDump, h5md_path, tpr_path="", ignore=(), edr=None
) -> RoundTripReport:
    source = extract_source(dump)
    target = extract_target(h5md_path)
    if edr is not None:
        source.update(extract_source_observables(edr))
        _target_observables(h5md_path, target)
    return RoundTripReport(
        str(h5md_path), str(tpr_path), compare(source, target, ignore=ignore)
    )
