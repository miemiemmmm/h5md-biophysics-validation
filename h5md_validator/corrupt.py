"""
Deliberately corrupted files with known expected violations.

Proving that valid files pass is only half of a schema claim; the other half is
that invalid metadata is caught. Each :class:`Corruption` injects exactly one
defect into a copy of a valid file and declares the finding code the validator
must raise. An expected code that is not raised is a detection failure.

A corruption that does not apply to a base file (no barostat group, no
selections, ...) returns False and is skipped, not counted as a failure.

Used by:
    Validation 4: negative tests
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import h5py
import numpy as np


@dataclass(frozen=True)
class Corruption:
    code: str
    category: str
    description: str
    #: Finding code the validator must raise for this file.
    expect: str
    apply: Callable

    def __call__(self, f) -> bool:
        return self.apply(f)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _replace(group, name, data, attrs=None):
    """Delete and rewrite a dataset, preserving its attributes by default."""
    keep = dict(group[name].attrs) if attrs is None and name in group else (attrs or {})
    if name in group:
        del group[name]
    ds = group.create_dataset(name, data=data)
    for k, v in keep.items():
        ds.attrs[k] = v
    return ds


def _mol_types(f):
    topo = f.get("topology")
    if topo is None:
        return []
    return [n for n in topo if isinstance(topo[n], h5py.Group)]


def _largest_mol(f):
    """The molecule type with the most particles -- the most informative target."""
    best, best_n = None, -1
    topo = f.get("topology")
    if topo is None:
        return None
    for n in _mol_types(f):
        count = int(np.atleast_1d(topo[n].attrs.get("particle_count", 0))[0])
        if count > best_n:
            best, best_n = n, count
    return None if best is None else f["topology"][best]


def _first_interaction(f):
    mol = _largest_mol(f)
    if mol is None or "interactions" not in mol:
        return None, None
    inter = mol["interactions"]
    for t in inter:
        if "parameter_id" in inter[t] and "particle_ids" in inter[t]:
            return inter[t], t
    return None, None


def _first_observable(f):
    obs = f.get("observables")
    if obs is None:
        return None
    for name in obs:
        node = obs[name]
        if isinstance(node, h5py.Group) and "step" in node and "value" in node:
            return node
    return None


def _ff_parameters(f):
    """Yield (group, interaction_name) for every force-field parameters dataset."""
    ff = f.get("force_field")
    if ff is None:
        return
    for layer in ("bonded", "nonbonded"):
        if layer in ff:
            for t in ff[layer]:
                if "parameters" in ff[layer][t]:
                    yield ff[layer][t], t
    for t in ff:
        if t in ("bonded", "nonbonded"):
            continue
        if isinstance(ff[t], h5py.Group) and "parameters" in ff[t]:
            yield ff[t], t


# --------------------------------------------------------------------------
# mutators
# --------------------------------------------------------------------------
def _m_module_version_unsupported(f):
    f["h5md/modules/biophysics"].attrs["version"] = np.array([9, 9], dtype=np.int32)
    return True


def _m_module_version_missing(f):
    del f["h5md/modules/biophysics"].attrs["version"]
    return True


def _m_module_missing(f):
    del f["h5md/modules/biophysics"]
    return True


def _m_total_particle_count(f):
    t = f["topology"]
    t.attrs["total_particle_count"] = (
        int(np.atleast_1d(t.attrs["total_particle_count"])[0]) + 1
    )
    return True


def _m_total_residue_count(f):
    t = f["topology"]
    t.attrs["total_residue_count"] = (
        int(np.atleast_1d(t.attrs["total_residue_count"])[0]) + 1
    )
    return True


def _m_block_count(f):
    blocks = f["topology/molecule_blocks"]
    data = blocks[()]
    data["count"][0] += 1
    _replace(f["topology"], "molecule_blocks", data)
    return True


def _m_block_dangling(f):
    blocks = f["topology/molecule_blocks"]
    data = blocks[()]
    data["name"][0] = (
        b"NO_SUCH_MOLECULE_TYPE"
        if data.dtype["name"].kind == "S"
        else "NO_SUCH_MOLECULE_TYPE"
    )
    _replace(f["topology"], "molecule_blocks", data)
    return True


def _m_blocks_missing(f):
    del f["topology/molecule_blocks"]
    return True


def _m_particle_data_missing(f):
    mol = _largest_mol(f)
    del mol["particle_data"]
    return True


def _m_particle_data_truncated(f):
    mol = _largest_mol(f)
    data = mol["particle_data"][()]
    if data.shape[0] < 2:
        return False
    _replace(mol, "particle_data", data[:-1])
    return True


def _m_particle_count_attr(f):
    mol = _largest_mol(f)
    mol.attrs["particle_count"] = int(np.atleast_1d(mol.attrs["particle_count"])[0]) + 1
    return True


def _m_sequence_truncated(f):
    mol = _largest_mol(f)
    data = mol["sequence"][()]
    if data.shape[0] < 2:
        return False
    _replace(mol, "sequence", data[:-1])
    return True


def _m_sequence_index_out_of_range(f):
    mol = _largest_mol(f)
    data = mol["particle_data"][()]
    n = int(np.atleast_1d(mol.attrs["residue_count"])[0])
    data["sequence_index"][0] = n + 5
    _replace(mol, "particle_data", data)
    return True


def _m_name_index_out_of_range(f):
    mol = _largest_mol(f)
    if mol["particle_data"].dtype["name"].kind not in ("i", "u"):
        return False
    data = mol["particle_data"][()]
    data["name"][0] = mol["particle_name_table"].shape[0] + 3
    _replace(mol, "particle_data", data)
    return True


def _m_component_name_out_of_range(f):
    mol = _largest_mol(f)
    data = mol["sequence"][()]
    data["component_name"][0] = mol["component_name_table"].shape[0] + 3
    _replace(mol, "sequence", data)
    return True


def _m_negative_mass(f):
    mol = _largest_mol(f)
    data = mol["particle_data"][()]
    data["mass"][0] = -1.0
    _replace(mol, "particle_data", data)
    return True


def _m_name_table_missing(f):
    mol = _largest_mol(f)
    del mol["particle_name_table"]
    return True


def _m_particle_data_field_renamed(f):
    mol = _largest_mol(f)
    data = mol["particle_data"][()]
    names = list(data.dtype.names)
    names[names.index("charge")] = "partial_charge"
    data.dtype = np.dtype([(n, data.dtype[i]) for i, n in enumerate(names)])
    _replace(mol, "particle_data", data)
    return True


def _m_parameter_id_dangling(f):
    grp, _ = _first_interaction(f)
    if grp is None:
        return False
    data = grp["parameter_id"][()]
    data[0] = 10**6
    _replace(grp, "parameter_id", data)
    return True


def _m_parameter_id_negative(f):
    grp, _ = _first_interaction(f)
    if grp is None:
        return False
    data = grp["parameter_id"][()]
    data[0] = -1
    _replace(grp, "parameter_id", data)
    return True


def _m_parameter_table_missing(f):
    grp, name = _first_interaction(f)
    if grp is None:
        return False
    for pgrp, tname in _ff_parameters(f):
        if tname == name:
            del pgrp["parameters"]
            return True
    return False


def _m_ff_field_renamed(f):
    for pgrp, _ in _ff_parameters(f):
        data = pgrp["parameters"][()]
        names = list(data.dtype.names)
        names[0] = "renamed_parameter"
        attrs = dict(pgrp["parameters"].attrs)
        data.dtype = np.dtype([(n, data.dtype[i]) for i, n in enumerate(names)])
        _replace(pgrp, "parameters", data, attrs=attrs)
        return True
    return False


def _m_units_missing(f):
    for pgrp, _ in _ff_parameters(f):
        if "units" in pgrp["parameters"].attrs:
            del pgrp["parameters"].attrs["units"]
            return True
    return False


def _m_units_wrong_length(f):
    # Only a table with more than one compound field can be shortened to one
    # entry. ``constraint`` has a single field and is the first table several
    # force fields emit, so taking the first table unconditionally rewrote
    # ['nm'] over ['nm']: the corruption reported success, the file stayed
    # valid, and the validator was marked as having missed a defect that was
    # never introduced.
    for pgrp, _ in _ff_parameters(f):
        ds = pgrp["parameters"]
        n_fields = len(ds.dtype.names) if ds.dtype.names else 1
        if "units" in ds.attrs and n_fields > 1:
            ds.attrs["units"] = np.array([b"nm"], dtype="S8")
            return True
    return False


def _m_unit_malformed(f):
    for pgrp, _ in _ff_parameters(f):
        ds = pgrp["parameters"]
        if "units" in ds.attrs:
            vals = list(np.atleast_1d(ds.attrs["units"]))
            vals[0] = b"nm^3"
            ds.attrs["units"] = np.array(vals, dtype="S16")
            return True
    return False


def _m_particle_ids_out_of_range(f):
    grp, _ = _first_interaction(f)
    if grp is None:
        return False
    data = grp["particle_ids"][()]
    data[0, 0] = 10**6
    _replace(grp, "particle_ids", data)
    return True


def _m_interaction_arity(f):
    grp, _ = _first_interaction(f)
    if grp is None:
        return False
    data = grp["particle_ids"][()]
    if data.ndim != 2 or data.shape[1] < 2:
        return False
    _replace(grp, "particle_ids", data[:, :-1])
    return True


def _m_interaction_length_mismatch(f):
    grp, _ = _first_interaction(f)
    if grp is None or grp["parameter_id"].shape[0] < 2:
        return False
    data = grp["parameter_id"][()]
    _replace(grp, "parameter_id", data[:-1])
    return True


def _m_exclusion_order(f):
    mol = _largest_mol(f)
    if mol is None or "exclusions" not in mol:
        return False
    pairs = mol["exclusions/particle_ids"][()]
    pairs[0] = pairs[0, ::-1]  # store the pair higher index first
    _replace(mol["exclusions"], "particle_ids", pairs)
    return True


def _m_exclusion_index(f):
    mol = _largest_mol(f)
    if mol is None or "exclusions" not in mol:
        return False
    pairs = mol["exclusions/particle_ids"][()]
    pairs[0, 1] = 10**6  # the second index, so the pair still has i < j
    _replace(mol["exclusions"], "particle_ids", pairs)
    return True


def _m_exclusion_duplicate(f):
    mol = _largest_mol(f)
    if mol is None or "exclusions" not in mol:
        return False
    pairs = mol["exclusions/particle_ids"][()]
    _replace(mol["exclusions"], "particle_ids", np.vstack([pairs, pairs[:1]]))
    return True


def _m_bond_out_of_range(f):
    if "connectivity/bonds" not in f:
        return False
    data = f["connectivity/bonds"][()]
    data[0, 0] = 10**9
    _replace(f["connectivity"], "bonds", data)
    return True


def _m_bond_self(f):
    if "connectivity/bonds" not in f:
        return False
    data = f["connectivity/bonds"][()]
    data[0, 1] = data[0, 0]
    _replace(f["connectivity"], "bonds", data)
    return True


def _m_bond_shape(f):
    if "connectivity/bonds" not in f:
        return False
    data = f["connectivity/bonds"][()]
    _replace(f["connectivity"], "bonds", np.column_stack([data, data[:, :1]]))
    return True


def _m_step_value_mismatch(f):
    obs = _first_observable(f)
    if obs is None or obs["step"].shape[0] < 2:
        return False
    _replace(obs, "step", obs["step"][()][:-1])
    return True


def _m_time_not_monotonic(f):
    obs = _first_observable(f)
    if obs is None or "time" not in obs or obs["time"].shape[0] < 3:
        return False
    data = obs["time"][()]
    data[1], data[2] = data[2], data[1]
    _replace(obs, "time", data)
    return True


def _m_time_unit_missing(f):
    obs = _first_observable(f)
    if obs is None or "time" not in obs or "unit" not in obs["time"].attrs:
        return False
    del obs["time"].attrs["unit"]
    return True


def _m_observable_shape(f):
    obs = f.get("observables")
    if obs is None or "pressure_tensor" not in obs:
        return False
    node = obs["pressure_tensor"]
    data = node["value"][()]
    if data.ndim != 3:
        return False
    _replace(node, "value", data[:, :, :2])
    return True


def _m_observable_unit_missing(f):
    obs = _first_observable(f)
    if obs is None or "unit" not in obs["value"].attrs:
        return False
    del obs["value"].attrs["unit"]
    return True


def _m_observable_unknown(f):
    obs = f.get("observables")
    if obs is None:
        return False
    obs.create_dataset("magnetisation", data=np.zeros(3))
    obs["magnetisation"].attrs["unit"] = "A m-1"
    return True


def _m_enum_invalid(f):
    rc = f.get("parameters/run_control")
    if rc is None or "simulation_type" not in rc:
        return False
    _replace(rc, "simulation_type", np.bytes_("quantum_dynamics"))
    return True


def _m_conditional_subgroup(f):
    tc = f.get("parameters/temperature_coupling")
    if tc is None or "algorithm" not in tc:
        return False
    current = tc["algorithm"][()]
    current = current.decode() if isinstance(current, bytes) else str(current)
    other = "nose_hoover" if current != "nose_hoover" else "berendsen"
    _replace(tc, "algorithm", np.bytes_(other))
    return True


def _m_per_group_length(f):
    tc = f.get("parameters/temperature_coupling")
    if tc is None or "reference_temperature" not in tc:
        return False
    data = np.atleast_1d(tc["reference_temperature"][()])
    _replace(tc, "reference_temperature", np.concatenate([data, data[:1]]))
    return True


def _m_periodicity_shape(f):
    bc = f.get("parameters/boundary_conditions")
    if bc is None or "periodicity" not in bc:
        return False
    _replace(bc, "periodicity", np.array([1, 1], dtype=np.int64))
    return True


def _m_timestep_dtype(f):
    md = f.get("parameters/run_control/molecular_dynamics")
    if md is None or "timestep" not in md:
        return False
    attrs = dict(md["timestep"].attrs)
    _replace(md, "timestep", np.int64(2), attrs=attrs)
    return True


def _m_timestep_unit_missing(f):
    md = f.get("parameters/run_control/molecular_dynamics")
    if md is None or "timestep" not in md or "unit" not in md["timestep"].attrs:
        return False
    del md["timestep"].attrs["unit"]
    return True


def _m_simulation_type_missing(f):
    rc = f.get("parameters/run_control")
    if rc is None or "simulation_type" not in rc:
        return False
    del rc["simulation_type"]
    return True


def _m_selection_index(f):
    sel = f.get("selections")
    if sel is None:
        return False
    for name in sel:
        grp = sel[name]
        if "indices" in grp and grp["indices"].size:
            data = grp["indices"][()]
            data[0] = 10**9
            _replace(grp, "indices", data)
            return True
    return False


def _m_selection_range(f):
    sel = f.get("selections")
    if sel is None:
        return False
    for name in sel:
        grp = sel[name]
        if "ranges" in grp and grp["ranges"].size:
            data = grp["ranges"][()]
            data[0, 1] = 10**9
            _replace(grp, "ranges", data)
            return True
    return False


def _m_selection_empty(f):
    sel = f.get("selections")
    if sel is None:
        return False
    sel.create_group("empty_selection")
    return True


def _m_trajectory_mismatch(f):
    particles = f.get("particles")
    if particles is None:
        return False
    for sysname in particles:
        pos = particles[sysname].get("position")
        if pos is None or "value" not in pos:
            continue
        data = pos["value"][()]
        if data.ndim != 3 or data.shape[1] < 2:
            continue
        attrs = dict(pos["value"].attrs)
        _replace(pos, "value", data[:, :-1, :], attrs=attrs)
        return True
    return False


# --------------------------------------------------------------------------
# catalogue
# --------------------------------------------------------------------------
CORRUPTIONS = [
    Corruption(
        "C-MOD-VERSION-UNSUPPORTED",
        "module",
        "declare module version 9.9",
        "MOD-VERSION-UNSUPPORTED",
        _m_module_version_unsupported,
    ),
    Corruption(
        "C-MOD-VERSION-MISSING",
        "module",
        "remove the module version attribute",
        "MOD-VERSION-MISSING",
        _m_module_version_missing,
    ),
    Corruption(
        "C-MOD-MISSING",
        "module",
        "remove the biophysics module declaration",
        "MOD-MISSING",
        _m_module_missing,
    ),
    # Topology related corruptions
    Corruption(
        "C-TOPO-TOTAL-PARTICLES",
        "topology invariant",
        "inflate total_particle_count by one",
        "TOPO-INVARIANT-PARTICLES",
        _m_total_particle_count,
    ),
    Corruption(
        "C-TOPO-TOTAL-RESIDUES",
        "topology invariant",
        "inflate total_residue_count by one",
        "TOPO-INVARIANT-RESIDUES",
        _m_total_residue_count,
    ),
    Corruption(
        "C-TOPO-BLOCK-COUNT",
        "topology invariant",
        "inflate the molecule count of the first block",
        "TOPO-INVARIANT-PARTICLES",
        _m_block_count,
    ),
    Corruption(
        "C-TOPO-BLOCK-DANGLING",
        "referential integrity",
        "point a molecule block at a molecule type that does not exist",
        "REF-BLOCK-DANGLING",
        _m_block_dangling,
    ),
    Corruption(
        "C-TOPO-BLOCKS-MISSING",
        "required element",
        "remove the molecule_blocks dataset",
        "STRUCT-MISSING-DATASET",
        _m_blocks_missing,
    ),
    Corruption(
        "C-TOPO-PARTICLE-DATA-MISSING",
        "required element",
        "remove particle_data from a molecule type",
        "STRUCT-MISSING-DATASET",
        _m_particle_data_missing,
    ),
    Corruption(
        "C-TOPO-PARTICLE-DATA-SHORT",
        "topology invariant",
        "drop the last row of particle_data",
        "TOPO-PARTICLE-DATA-LENGTH",
        _m_particle_data_truncated,
    ),
    Corruption(
        "C-TOPO-PARTICLE-COUNT-ATTR",
        "topology invariant",
        "inflate particle_count of a molecule type",
        "TOPO-PARTICLE-DATA-LENGTH",
        _m_particle_count_attr,
    ),
    Corruption(
        "C-TOPO-SEQUENCE-SHORT",
        "topology invariant",
        "drop the last row of sequence",
        "TOPO-SEQUENCE-LENGTH",
        _m_sequence_truncated,
    ),
    Corruption(
        "C-TOPO-SEQUENCE-INDEX",
        "referential integrity",
        "set a sequence_index beyond residue_count",
        "REF-SEQUENCE-INDEX",
        _m_sequence_index_out_of_range,
    ),
    Corruption(
        "C-TOPO-NAME-INDEX",
        "referential integrity",
        "set a particle name index beyond particle_name_table",
        "REF-PARTICLE-NAME-INDEX",
        _m_name_index_out_of_range,
    ),
    Corruption(
        "C-TOPO-COMPONENT-NAME-INDEX",
        "referential integrity",
        "set a component name index beyond component_name_table",
        "REF-COMPONENT-NAME-INDEX",
        _m_component_name_out_of_range,
    ),
    Corruption(
        "C-TOPO-NEGATIVE-MASS",
        "physical validity",
        "give a particle a negative mass",
        "TOPO-NEGATIVE-MASS",
        _m_negative_mass,
    ),
    Corruption(
        "C-TOPO-NAME-TABLE-MISSING",
        "required element",
        "remove particle_name_table",
        "STRUCT-MISSING-DATASET",
        _m_name_table_missing,
    ),
    Corruption(
        "C-TOPO-FIELD-RENAMED",
        "structure",
        "rename the charge field of particle_data",
        "STRUCT-COMPOUND-FIELDS",
        _m_particle_data_field_renamed,
    ),
    # Force-field related corruptions
    Corruption(
        "C-FF-PARAM-ID-DANGLING",
        "referential integrity",
        "reference a force-field parameter index that does not exist",
        "REF-PARAMETER-ID",
        _m_parameter_id_dangling,
    ),
    Corruption(
        "C-FF-PARAM-ID-NEGATIVE",
        "referential integrity",
        "reference a negative force-field parameter index",
        "REF-PARAMETER-ID",
        _m_parameter_id_negative,
    ),
    Corruption(
        "C-FF-PARAM-TABLE-MISSING",
        "referential integrity",
        "remove a parameters dataset that an interaction list references",
        "REF-PARAMETER-TABLE-MISSING",
        _m_parameter_table_missing,
    ),
    Corruption(
        "C-FF-FIELD-RENAMED",
        "structure",
        "rename a force-field compound parameter field",
        "STRUCT-COMPOUND-FIELDS",
        _m_ff_field_renamed,
    ),
    Corruption(
        "C-FF-UNITS-MISSING",
        "units",
        "remove the units attribute of a compound parameters dataset",
        "UNITS-MISSING",
        _m_units_missing,
    ),
    Corruption(
        "C-FF-UNITS-LENGTH",
        "units",
        "shorten the units attribute to one entry",
        "UNITS-LENGTH",
        _m_units_wrong_length,
    ),
    Corruption(
        "C-FF-UNIT-MALFORMED",
        "units",
        "write a unit as 'nm^3' instead of 'nm+3'",
        "UNIT-SYNTAX",
        _m_unit_malformed,
    ),
    Corruption(
        "C-INT-PARTICLE-ID",
        "referential integrity",
        "set an interaction particle_id beyond particle_count",
        "REF-PARTICLE-ID",
        _m_particle_ids_out_of_range,
    ),
    Corruption(
        "C-INT-ARITY",
        "structure",
        "drop a column from particle_ids",
        "STRUCT-INTERACTION-ARITY",
        _m_interaction_arity,
    ),
    Corruption(
        "C-INT-LENGTH",
        "structure",
        "make parameter_id and particle_ids differ in length",
        "STRUCT-INTERACTION-LENGTH",
        _m_interaction_length_mismatch,
    ),
    Corruption(
        "C-EXCL-ORDER",
        "referential integrity",
        "store an exclusion pair higher index first",
        "REF-EXCL-ORDER",
        _m_exclusion_order,
    ),
    Corruption(
        "C-EXCL-INDEX",
        "referential integrity",
        "set an exclusion index beyond particle_count",
        "REF-EXCL-INDEX",
        _m_exclusion_index,
    ),
    Corruption(
        "C-EXCL-DUPLICATE",
        "structure",
        "list an exclusion pair twice",
        "REF-EXCL-DUPLICATE",
        _m_exclusion_duplicate,
    ),
    # Connectivity related corruptions
    Corruption(
        "C-CONN-INDEX",
        "referential integrity",
        "set a bond particle index beyond the system size",
        "CONN-PARTICLE-ID",
        _m_bond_out_of_range,
    ),
    Corruption(
        "C-CONN-SELF",
        "physical validity",
        "bond a particle to itself",
        "CONN-SELF-BOND",
        _m_bond_self,
    ),
    Corruption(
        "C-CONN-SHAPE",
        "structure",
        "store bonds as Integer[N][3]",
        "STRUCT-SHAPE",
        _m_bond_shape,
    ),
    # Time-dependent dataset related corruptions
    Corruption(
        "C-TIME-LENGTH",
        "time alignment",
        "truncate the step dataset of a time-dependent observable",
        "TIME-LENGTH-MISMATCH",
        _m_step_value_mismatch,
    ),
    Corruption(
        "C-TIME-MONOTONIC",
        "time alignment",
        "swap two time values so the series is not monotonic",
        "TIME-NOT-MONOTONIC",
        _m_time_not_monotonic,
    ),
    Corruption(
        "C-TIME-UNIT-MISSING",
        "units",
        "remove the unit attribute of a time dataset",
        "UNIT-MISSING",
        _m_time_unit_missing,
    ),
    # Observable related corruptions
    Corruption(
        "C-OBS-SHAPE",
        "structure",
        "store the pressure tensor as 3x2",
        "STRUCT-SHAPE",
        _m_observable_shape,
    ),
    Corruption(
        "C-OBS-UNIT-MISSING",
        "units",
        "remove the unit attribute of an observable value",
        "UNIT-MISSING",
        _m_observable_unit_missing,
    ),
    Corruption(
        "C-OBS-UNKNOWN",
        "structure",
        "add an observable that the specification does not define",
        "STRUCT-UNDOCUMENTED",
        _m_observable_unknown,
    ),
    # MD Parameter related corruptions
    Corruption(
        "C-PARAM-ENUM",
        "enumeration",
        "set simulation_type to a value outside the defined set",
        "ENUM-UNDEFINED-VALUE",
        _m_enum_invalid,
    ),
    Corruption(
        "C-PARAM-CONDITIONAL",
        "conditional structure",
        "name a thermostat whose parameter subgroup is absent",
        "COND-SUBGROUP-MISSING",
        _m_conditional_subgroup,
    ),
    Corruption(
        "C-PARAM-PER-GROUP",
        "structure",
        "give reference_temperature one more entry than there are groups",
        "STRUCT-PER-GROUP-LENGTH",
        _m_per_group_length,
    ),
    Corruption(
        "C-PARAM-PERIODICITY-SHAPE",
        "structure",
        "store periodicity as a 2-vector",
        "STRUCT-SHAPE",
        _m_periodicity_shape,
    ),
    Corruption(
        "C-PARAM-DTYPE",
        "structure",
        "store the timestep as an Integer",
        "STRUCT-DTYPE",
        _m_timestep_dtype,
    ),
    Corruption(
        "C-PARAM-UNIT-MISSING",
        "units",
        "remove the unit attribute of the timestep",
        "UNIT-MISSING",
        _m_timestep_unit_missing,
    ),
    Corruption(
        "C-PARAM-REQUIRED-MISSING",
        "required element",
        "remove run_control/simulation_type",
        "STRUCT-MISSING-DATASET",
        _m_simulation_type_missing,
    ),
    # Selection related corruptions
    Corruption(
        "C-SEL-INDEX",
        "referential integrity",
        "set a selection index beyond the system size",
        "REF-SEL-INDEX",
        _m_selection_index,
    ),
    Corruption(
        "C-SEL-RANGE",
        "referential integrity",
        "set a selection range end beyond the system size",
        "REF-SEL-RANGE",
        _m_selection_range,
    ),
    Corruption(
        "C-SEL-EMPTY",
        "structure",
        "add a selection group that defines nothing",
        "SEL-EMPTY",
        _m_selection_empty,
    ),
    # Trajectory related corruptions
    Corruption(
        "C-TRAJ-MISMATCH",
        "cross-group consistency",
        "drop a particle from the trajectory but not from the topology",
        "TOPO-TRAJECTORY-MISMATCH",
        _m_trajectory_mismatch,
    ),
]

CORRUPTIONS_BY_CODE = {c.code: c for c in CORRUPTIONS}


def generate(base_file, out_dir, corruptions=None):
    """Write one corrupted copy of ``base_file`` per applicable corruption.

    Returns a list of ``(Corruption, path)`` for the files actually produced.
    """
    base_file = Path(base_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    produced = []
    for corruption in corruptions or CORRUPTIONS:
        target = out_dir / f"{base_file.stem}__{corruption.code}.h5md"
        shutil.copyfile(base_file, target)
        try:
            with h5py.File(target, "r+") as f:
                applied = corruption(f)
        except Exception as exc:  # a mutator that cannot run is not a test failure
            target.unlink(missing_ok=True)
            produced.append((corruption, None, f"{type(exc).__name__}: {exc}"))
            continue
        if not applied:
            target.unlink(missing_ok=True)
            produced.append((corruption, None, "not applicable to this base file"))
        else:
            produced.append((corruption, target, None))
    return produced
