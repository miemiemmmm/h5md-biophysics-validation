"""
Conformance and internal-consistency validator for the H5MD ``biophysics`` module.

Implements the checks of Validation 3 (internal consistency) of the technical
validation plan, on top of a structural conformance pass against the normative
schema model in :mod:`h5md_validator.spec`.

Check families, by code prefix:

    MOD     module declaration and version
    STRUCT  presence, type and shape of specified elements
    ENUM    values of enumerated elements
    UNIT    presence and syntax of unit / units attributes
    COND    conditional subgroup instantiation (S7)
    TIME    step / time / value first-dimension alignment
    TOPO    topology invariants (particle and component counts)
    REF     referential integrity of every index in the file
    SEL     selection bounds
    CONN    connectivity bounds

Used by:
    Validation 3: conformance and internal consistency
    Validation 4: negative / corruption tests -- the same checks are what the
                  deliberately broken files of :mod:`h5md_validator.corrupt`
                  have to fail
"""

from __future__ import annotations

import numpy as np
import h5py

from . import spec, units as unitmod
from .deviations import DEVIATIONS
from .findings import Report, Severity


# --------------------------------------------------------------------------
# small HDF5 helpers
# --------------------------------------------------------------------------
def _decode(value):
    """Normalise an HDF5 scalar/array of strings into str / list[str]."""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.ndarray):
        if value.dtype.kind in ("S", "O"):
            return [_decode(v) for v in value.tolist()]
        return value
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read(ds):
    return _decode(ds[()])


def _dtype_class(dtype) -> str:
    """Map an HDF5 dtype onto the specification's type vocabulary."""
    if dtype.names is not None:
        return spec.COMPOUND
    kind = dtype.kind
    if kind in ("i", "u"):
        return spec.INTEGER
    if kind == "f":
        return spec.FLOAT
    if kind in ("S", "U"):
        return spec.STRING
    if kind == "O":
        # h5py stores variable-length strings as object arrays.
        return spec.STRING
    if kind == "b":
        return spec.BOOLEAN
    return f"unknown({kind})"


def _is_time_dependent(obj) -> bool:
    """H5MD time-dependent element: a group holding 'step' and 'value'."""
    return isinstance(obj, h5py.Group) and "step" in obj and "value" in obj


class Validator:
    """Validates one H5MD file against the ``biophysics`` v1.0 specification."""

    def __init__(self, *, strict_units: bool = True, report_deviations: bool = True):
        self.strict_units = strict_units
        self.report_deviations = report_deviations

    # ----------------------------------------------------------------- entry
    def validate(self, path) -> Report:
        report = Report(file=str(path))
        try:
            handle = h5py.File(path, "r")
        except OSError as exc:
            # A file the HDF5 library cannot open is a validation verdict, not a
            # crash: a writer that aborts mid-write can leave exactly this.
            report.add("FILE-UNREADABLE", Severity.ERROR, "",
                       f"file cannot be opened as HDF5: {exc}",
                       error=str(exc))
            return report
        with handle as f:
            self._check_module(f, report)
            self._check_parameters(f, report)
            self._check_observables(f, report)
            self._check_force_field(f, report)
            topo = self._check_topology(f, report)
            self._check_connectivity(f, report, topo)
            self._check_selections(f, report, topo)
            self._check_cross_group(f, report)
        return report

    def _check_cross_group(self, f, rep):
        """Checks that need two parts of the file at once."""
        constraints = f.get("parameters/constraints")
        if constraints is None:
            return
        if "mode" not in constraints and self.report_deviations:
            rep.add("DEV-CONSTRAINTS-MODE-UNWRITTEN", Severity.DEVIATION,
                    "parameters/constraints",
                    DEVIATIONS["DEV-CONSTRAINTS-MODE-UNWRITTEN"].title,
                    deviation="DEV-CONSTRAINTS-MODE-UNWRITTEN")
        if "algorithm" not in constraints:
            return
        # An algorithm is declared; check that something is actually constrained.
        n_constraints = 0
        topo = f.get("topology")
        if topo is not None:
            for name in topo:
                mol = topo[name]
                if not isinstance(mol, h5py.Group):
                    continue
                inter = mol.get("interactions")
                if not isinstance(inter, h5py.Group):
                    continue
                # Only 'constraint' entries are solved by LINCS or SHAKE.
                # SETTLE is a separate analytic algorithm for rigid water and
                # says nothing about the configured constraint solver (S7.8).
                node = inter.get("constraint")
                if node is not None and "parameter_id" in node:
                    n_constraints += node["parameter_id"].shape[0]
        if n_constraints == 0 and self.report_deviations:
            algorithm = _read(constraints["algorithm"])
            rep.add("DEV-CONSTRAINTS-ALGORITHM-UNCONDITIONAL", Severity.DEVIATION,
                    "parameters/constraints/algorithm",
                    DEVIATIONS["DEV-CONSTRAINTS-ALGORITHM-UNCONDITIONAL"].title,
                    deviation="DEV-CONSTRAINTS-ALGORITHM-UNCONDITIONAL",
                    algorithm=algorithm)

    # ------------------------------------------------------------- S4.1 root
    def _check_module(self, f, rep):
        modpath = "h5md/modules/biophysics"
        mod = f.get(modpath)
        if mod is None:
            rep.add(
                "MOD-MISSING", Severity.ERROR, modpath,
                "file does not declare the biophysics module",
            )
            return
        if "version" not in mod.attrs:
            rep.add(
                "MOD-VERSION-MISSING", Severity.ERROR, modpath,
                "biophysics module has no version attribute",
            )
            return
        version = tuple(int(v) for v in np.atleast_1d(mod.attrs["version"]))
        if len(version) != 2:
            rep.add(
                "MOD-VERSION-SHAPE", Severity.ERROR, modpath,
                f"version must be Integer[2], got {len(version)} element(s)",
                actual=list(version),
            )
        elif version != spec.MODULE_VERSION:
            rep.add(
                "MOD-VERSION-UNSUPPORTED", Severity.ERROR, modpath,
                f"unsupported module version {version[0]}.{version[1]}, "
                f"validator implements {spec.MODULE_VERSION[0]}.{spec.MODULE_VERSION[1]}",
                actual=list(version), expected=list(spec.MODULE_VERSION),
            )

    # ------------------------------------------------------------ S4.2 units
    def _check_unit_attr(self, ds, path, rep, required: bool):
        if "unit" not in ds.attrs:
            if required:
                rep.add(
                    "UNIT-MISSING", Severity.ERROR, path,
                    "dataset holds a physical quantity but has no unit attribute",
                )
            return
        unit = _decode(ds.attrs["unit"])
        if not isinstance(unit, str):
            rep.add("UNIT-TYPE", Severity.ERROR, path,
                    "unit attribute must be a String", actual=str(type(unit)))
            return
        if not unitmod.is_valid(unit):
            rep.add("UNIT-SYNTAX", Severity.ERROR, path,
                    f"unit {unit!r} is not valid H5MD unit notation", unit=unit)
        elif self.strict_units:
            unknown = unitmod.unknown_symbols(unit)
            if unknown:
                rep.add("UNIT-SYMBOL-UNKNOWN", Severity.WARNING, path,
                        f"unit {unit!r} uses unrecognised symbol(s) {unknown}",
                        unit=unit, symbols=unknown)

    def _check_units_attr(self, ds, path, rep, n_fields: int):
        """S4.2: a compound dataset carries a ``units`` list, one per field."""
        if "units" not in ds.attrs:
            rep.add("UNITS-MISSING", Severity.ERROR, path,
                    "compound dataset has no units attribute")
            return
        vals = _decode(ds.attrs["units"])
        vals = [vals] if isinstance(vals, str) else list(vals)
        if len(vals) != n_fields:
            rep.add("UNITS-LENGTH", Severity.ERROR, path,
                    f"units has {len(vals)} entries for {n_fields} compound fields",
                    expected=n_fields, actual=len(vals))
        for u in vals:
            if u and not unitmod.is_valid(u):
                rep.add("UNIT-SYNTAX", Severity.ERROR, path,
                        f"unit {u!r} is not valid H5MD unit notation", unit=u)

    # ------------------------------------------------- generic tree conformance
    def _check_group(self, parent, node: spec.Group, path: str, rep):
        """Recursively check a specified group against the file."""
        obj = parent.get(node.name)
        if obj is None:
            if node.required:
                rep.add("STRUCT-MISSING-GROUP", Severity.ERROR, path,
                        f"required group {node.name!r} is absent")
            return
        if not isinstance(obj, h5py.Group):
            rep.add("STRUCT-NOT-GROUP", Severity.ERROR, path,
                    f"{node.name!r} is specified as a group but is a dataset")
            return

        for attr in node.attrs:
            if attr.name not in obj.attrs and attr.required:
                rep.add("STRUCT-MISSING-ATTR", Severity.ERROR, path,
                        f"required attribute {attr.name!r} is absent")

        for dsnode in node.datasets:
            self._check_dataset(obj, dsnode, f"{path}/{dsnode.name}", rep)

        for sub in node.groups:
            self._check_group(obj, sub, f"{path}/{sub.name}", rep)

        # Elements present in the file that the specification does not define.
        specified = {d.name for d in node.datasets} | {g.name for g in node.groups}
        for key in obj:
            if key not in specified:
                self._report_undocumented(f"{path}/{key}", rep, node.name)

    def _report_undocumented(self, path, rep, parent_name):
        """An element present in the file that the specification does not define."""
        rep.add("STRUCT-UNDOCUMENTED", Severity.WARNING, path,
                f"element is not defined by the specification for {parent_name!r}")

    def _check_dataset(self, parent, node: spec.Dataset, path: str, rep):
        obj = parent.get(node.name)
        if obj is None:
            if node.required:
                rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, path,
                        f"required dataset {node.name!r} is absent")
            return
        if isinstance(obj, h5py.Group):
            # Permitted when the element is stored time-dependently.
            if _is_time_dependent(obj):
                self._check_time_alignment(obj, path, rep)
                return
            rep.add("STRUCT-NOT-DATASET", Severity.ERROR, path,
                    f"{node.name!r} is specified as a dataset but is a group")
            return

        actual = _dtype_class(obj.dtype)
        expected = spec.INTEGER if node.dtype == spec.BOOLEAN else node.dtype
        if actual != expected:
            rep.add("STRUCT-DTYPE", Severity.ERROR, path,
                    f"expected {node.dtype}, found {actual}",
                    expected=node.dtype, actual=actual)
        else:
            self._check_shape(obj, node, path, rep)
            if node.dtype == spec.BOOLEAN:
                bad = sorted({int(v) for v in np.atleast_1d(obj[()]).ravel()} - {0, 1})
                if bad:
                    # S3.3 admits any non-zero as true, so this is informational.
                    rep.add("ENUM-BOOLEAN-NONCANONICAL", Severity.INFO, path,
                            f"Boolean stored with non-canonical value(s) {bad}",
                            values=bad)
            if node.enum is not None:
                self._check_enum(obj, node, path, rep)

        self._check_unit_attr(obj, path, rep, required=node.unit)

    def _check_shape(self, obj, node: spec.Dataset, path, rep):
        want = node.shape
        got = tuple(obj.shape)
        if want is None:
            if got != ():
                rep.add("STRUCT-SHAPE", Severity.ERROR, path,
                        f"expected a scalar, found shape {list(got)}",
                        expected="scalar", actual=list(got))
        elif want == ("scalar_or_tensor",):
            if got not in ((), (3, 3)):
                rep.add("STRUCT-SHAPE", Severity.ERROR, path,
                        f"expected a scalar or a 3x3 tensor, found shape {list(got)}",
                        expected="scalar or [3,3]", actual=list(got))
        elif all(isinstance(d, int) for d in want):
            if got != want:
                rep.add("STRUCT-SHAPE", Severity.ERROR, path,
                        f"expected shape {list(want)}, found {list(got)}",
                        expected=list(want), actual=list(got))
        else:
            if len(got) != len(want):
                rep.add("STRUCT-SHAPE", Severity.ERROR, path,
                        f"expected a {len(want)}-D dataset, found {len(got)}-D",
                        expected=len(want), actual=len(got))

    def _check_enum(self, obj, node: spec.Dataset, path, rep):
        value = _read(obj)
        vals = value if isinstance(value, list) else [value]
        for v in vals:
            if not isinstance(v, str):
                continue
            if v not in node.enum:
                # S7 permits algorithms outside the predefined set, provided a
                # doi/url points at a definition.  That is a WARNING, not an error.
                rep.add("ENUM-UNDEFINED-VALUE", Severity.WARNING, path,
                        f"value {v!r} is outside the predefined set "
                        f"{sorted(node.enum)}; a doi or url attribute is recommended",
                        value=v, predefined=sorted(node.enum))

    # ------------------------------------------------- S3.2 time-dependent data
    def _check_time_alignment(self, grp, path, rep):
        """N_step == N_time == N_value on the first dimension (S3.2)."""
        lengths = {}
        for key in ("step", "time", "value"):
            if key in grp and isinstance(grp[key], h5py.Dataset):
                shape = grp[key].shape
                lengths[key] = shape[0] if shape else None
        if "step" not in lengths:
            rep.add("TIME-STEP-MISSING", Severity.ERROR, path,
                    "time-dependent element has no step dataset")
        if "value" not in lengths:
            rep.add("TIME-VALUE-MISSING", Severity.ERROR, path,
                    "time-dependent element has no value dataset")
        present = {k: v for k, v in lengths.items() if v is not None}
        if len(set(present.values())) > 1:
            rep.add("TIME-LENGTH-MISMATCH", Severity.ERROR, path,
                    "first dimension of step/time/value datasets disagree: "
                    + ", ".join(f"{k}={v}" for k, v in sorted(present.items())),
                    lengths=present)
        for key in ("step", "time"):
            if key in grp and isinstance(grp[key], h5py.Dataset):
                data = np.asarray(grp[key][()]).ravel()
                if data.size > 1 and np.any(np.diff(data) < 0):
                    rep.add("TIME-NOT-MONOTONIC", Severity.ERROR, f"{path}/{key}",
                            f"{key} values are not monotonically increasing")
        if "time" in grp and isinstance(grp["time"], h5py.Dataset):
            self._check_unit_attr(grp["time"], f"{path}/time", rep, required=True)

    # ---------------------------------------------------------- S7 parameters
    def _check_parameters(self, f, rep):
        params = f.get("parameters")
        if params is None:
            return
        self._check_group(f, spec.PARAMETERS, "parameters", rep)

        # S7: the subgroup named by the selector dataset must be instantiated.
        for gpath, (selector, owning) in spec.CONDITIONAL_SELECTORS.items():
            grp = f.get(gpath)
            if grp is None or selector not in grp:
                continue
            value = _read(grp[selector])
            if not isinstance(value, str):
                continue
            if value in owning and value not in grp:
                rep.add("COND-SUBGROUP-MISSING", Severity.ERROR, gpath,
                        f"{selector} is {value!r} but no {value!r} subgroup is present",
                        selector=selector, value=value)
            # A subgroup for an algorithm that was not selected is contradictory.
            for name in owning:
                if name != value and name in grp and isinstance(grp[name], h5py.Group):
                    rep.add("COND-SUBGROUP-ORPHAN", Severity.ERROR, f"{gpath}/{name}",
                            f"subgroup {name!r} is present but {selector} is {value!r}",
                            selector=selector, value=value, orphan=name)

        # S7.3: per-group arrays must match the length of 'groups'.
        for gpath, names in spec.PER_GROUP_DATASETS.items():
            grp = f.get(gpath)
            if grp is None or "groups" not in grp:
                continue
            n = len(np.atleast_1d(grp["groups"][()]))
            for rel in names:
                ds = grp.get(rel)
                if ds is None or not isinstance(ds, h5py.Dataset):
                    continue
                got = ds.shape[0] if ds.shape else 1
                if got != n:
                    rep.add("STRUCT-PER-GROUP-LENGTH", Severity.ERROR, f"{gpath}/{rel}",
                            f"has {got} value(s) for {n} coupling group(s)",
                            expected=n, actual=got)

        # S7: an algorithm outside the predefined set should carry a doi or url.
        for gpath, (selector, predefined) in spec.SELECTOR_ENUMS.items():
            grp = f.get(gpath)
            if grp is None or selector not in grp:
                continue
            value = _read(grp[selector])
            if isinstance(value, str) and value not in predefined:
                sub = grp.get(value)
                holder = sub if sub is not None else grp.get(selector)
                if holder is not None and not ({"doi", "url"} & set(holder.attrs)):
                    rep.add("REF-PROVENANCE-MISSING", Severity.WARNING,
                            f"{gpath}/{value}" if sub is not None else f"{gpath}/{selector}",
                            f"{value!r} is not a predefined value; a doi or url "
                            "attribute is recommended", value=value)

    # --------------------------------------------------------- S6 observables
    def _check_observables(self, f, rep):
        obs = f.get("observables")
        if obs is None:
            return
        for name in obs:
            node = obs[name]
            path = f"observables/{name}"
            if _is_time_dependent(node):
                self._check_time_alignment(node, path, rep)
                self._check_observable_shape(name, node["value"], path, rep)
                self._check_unit_attr(node["value"], f"{path}/value", rep, required=True)
            elif isinstance(node, h5py.Dataset):
                self._check_observable_shape(name, node, path, rep)
                self._check_unit_attr(node, path, rep, required=True)
            elif isinstance(node, h5py.Group):
                # Per-selection-group nesting: observables/<group>/<observable>.
                if self.report_deviations:
                    rep.add("DEV-OBS-GROUP-NESTING", Severity.DEVIATION, path,
                            DEVIATIONS["DEV-OBS-GROUP-NESTING"].title,
                            deviation="DEV-OBS-GROUP-NESTING")
                for sub in node:
                    subnode = node[sub]
                    subpath = f"{path}/{sub}"
                    if _is_time_dependent(subnode):
                        self._check_time_alignment(subnode, subpath, rep)
                        self._check_observable_shape(sub, subnode["value"], subpath, rep)
                        self._check_unit_attr(subnode["value"], f"{subpath}/value",
                                              rep, required=True)

    def _check_observable_shape(self, name, ds, path, rep):
        want = spec.OBSERVABLES.get(name)
        if want is None:
            rep.add("STRUCT-UNDOCUMENTED", Severity.WARNING, path,
                    f"{name!r} is not one of the predefined observables of S6",
                    name=name)
            return
        got = tuple(ds.shape)
        # Drop the leading frame axis for time-dependent storage.
        record = got[1:] if len(got) > len(want) else got
        if tuple(record) != tuple(want):
            rep.add("STRUCT-SHAPE", Severity.ERROR, path,
                    f"observable {name!r} expects record shape {list(want)}, "
                    f"found {list(record)}", expected=list(want), actual=list(record))

    # --------------------------------------------------------- S8 force field
    def _check_force_field(self, f, rep):
        ff = f.get("force_field")
        if ff is None:
            return {}

        # DEV-FF-BONDED-LAYER: the spec puts <type> directly under force_field.
        layered = {k for k in ff if k in ("bonded", "nonbonded")}
        if layered and self.report_deviations:
            rep.add("DEV-FF-BONDED-LAYER", Severity.DEVIATION, "force_field",
                    DEVIATIONS["DEV-FF-BONDED-LAYER"].title,
                    deviation="DEV-FF-BONDED-LAYER", subgroups=sorted(layered))

        for container, base in self._force_field_containers(ff):
            for tname in container:
                tgrp = container[tname]
                path = f"{base}/{tname}"
                itype = spec.INTERACTION_TYPES.get(tname)
                if itype is None:
                    # Table SI is not exhaustive by design. Treat an unlisted
                    # type the way S7 treats an algorithm outside its predefined
                    # set: permitted, and expected to point at a definition.
                    # See deviations.DEV-FF-NO-EXTENSION-CONVENTION -- S8 does
                    # not yet state this, so the check is advisory.
                    if not ({"doi", "url"} & set(tgrp.attrs)):
                        rep.add("REF-PROVENANCE-MISSING", Severity.WARNING, path,
                                f"{tname!r} is outside the predefined set of "
                                "Table SI; an extension is permitted, but a doi or "
                                "url attribute naming the functional form is "
                                "recommended", name=tname)
                    else:
                        rep.add("STRUCT-EXTENSION", Severity.INFO, path,
                                f"{tname!r} is an extension outside Table SI, with "
                                "a reference attached", name=tname)
                    continue
                if itype.bonded and base.endswith("nonbonded"):
                    rep.add("STRUCT-FF-CLASSIFICATION", Severity.ERROR, path,
                            f"{tname!r} is a bonded interaction type but is stored "
                            "under nonbonded")
                if not itype.bonded and base.endswith("/bonded"):
                    rep.add("STRUCT-FF-CLASSIFICATION", Severity.ERROR, path,
                            f"{tname!r} is a nonbonded interaction type but is stored "
                            "under bonded")
                for pname in ("parameters", "parameters_B"):
                    if pname not in tgrp:
                        if pname == "parameters":
                            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR,
                                    f"{path}/parameters",
                                    "interaction type has no parameters dataset")
                        continue
                    self._check_parameter_compound(tgrp[pname], f"{path}/{pname}",
                                                   itype, rep)

    def _force_field_containers(self, ff):
        """Yield (group, path) pairs that hold interaction-type subgroups."""
        found = False
        for layer in ("bonded", "nonbonded"):
            if layer in ff and isinstance(ff[layer], h5py.Group):
                found = True
                yield ff[layer], f"force_field/{layer}"
        if not found:
            yield ff, "force_field"

    def _check_parameter_compound(self, ds, path, itype, rep):
        if ds.dtype.names is None:
            rep.add("STRUCT-DTYPE", Severity.ERROR, path,
                    "parameters must be a compound dataset",
                    expected=spec.COMPOUND, actual=_dtype_class(ds.dtype))
            return
        got = tuple(ds.dtype.names)
        accepted = [tuple(fl.name for fl in itype.params)]
        if itype.alt_params:
            accepted.append(tuple(fl.name for fl in itype.alt_params))
        if got not in accepted:
            rep.add("STRUCT-COMPOUND-FIELDS", Severity.ERROR, path,
                    f"compound fields {list(got)} do not match any layout "
                    f"specified for {itype.name!r}: {[list(a) for a in accepted]}",
                    expected=[list(a) for a in accepted], actual=list(got))
        else:
            layout = itype.params if got == accepted[0] else itype.alt_params
            for fl in layout:
                actual = _dtype_class(ds.dtype[fl.name])
                if actual != fl.dtype:
                    rep.add("STRUCT-COMPOUND-DTYPE", Severity.ERROR, f"{path}.{fl.name}",
                            f"field {fl.name!r} expected {fl.dtype}, found {actual}",
                            expected=fl.dtype, actual=actual)
        self._check_units_attr(ds, path, rep, len(got))

    # ------------------------------------------------------------ S10 topology
    def _check_topology(self, f, rep):
        """Validate the topology group; return the resolved topology summary."""
        topo = f.get("topology")
        if topo is None:
            return None

        info = {"molecule_types": {}, "total_particle_count": None,
                "total_component_count": None, "n_particles": None}

        for attr in ("total_particle_count", "total_component_count"):
            if attr not in topo.attrs:
                rep.add("STRUCT-MISSING-ATTR", Severity.ERROR, "topology",
                        f"required attribute {attr!r} is absent")
            else:
                info[attr] = int(np.atleast_1d(topo.attrs[attr])[0])
        if "system_name" in topo.attrs and self.report_deviations:
            rep.add("DEV-TOPO-SYSTEM-NAME", Severity.DEVIATION, "topology",
                    DEVIATIONS["DEV-TOPO-SYSTEM-NAME"].title,
                    deviation="DEV-TOPO-SYSTEM-NAME")

        ff_params = self._collect_parameter_counts(f)

        # --- molecule types -------------------------------------------------
        for name in topo:
            if name == "molecule_blocks":
                continue
            mol = topo[name]
            if not isinstance(mol, h5py.Group):
                continue
            info["molecule_types"][name] = self._check_molecule_type(
                mol, f"topology/{name}", rep, ff_params)

        # --- molecule_blocks and the S10 count invariants ---------------------
        blocks = topo.get("molecule_blocks")
        if blocks is None:
            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, "topology/molecule_blocks",
                    "required dataset molecule_blocks is absent")
            return info
        if blocks.dtype.names is None:
            rep.add("STRUCT-DTYPE", Severity.ERROR, "topology/molecule_blocks",
                    "molecule_blocks must be a compound dataset")
            return info
        expected_fields = tuple(fl.name for fl in spec.MOLECULE_BLOCK_FIELDS)
        if tuple(blocks.dtype.names) != expected_fields:
            rep.add("STRUCT-COMPOUND-FIELDS", Severity.ERROR, "topology/molecule_blocks",
                    f"compound fields {list(blocks.dtype.names)} do not match "
                    f"{list(expected_fields)}",
                    expected=list(expected_fields), actual=list(blocks.dtype.names))
            return info

        data = blocks[()]
        total_p = total_c = 0
        for i, row in enumerate(data):
            bname = _decode(row["name"])
            count = int(row["count"])
            if count < 0:
                rep.add("TOPO-BLOCK-COUNT-NEGATIVE", Severity.ERROR,
                        f"topology/molecule_blocks[{i}]",
                        f"molecule block {bname!r} has negative count {count}",
                        index=i, name=bname, count=count)
            mol = info["molecule_types"].get(bname)
            if mol is None:
                rep.add("REF-BLOCK-DANGLING", Severity.ERROR,
                        f"topology/molecule_blocks[{i}]",
                        f"molecule block references molecule type {bname!r}, "
                        "which has no group under topology/",
                        index=i, name=bname)
                continue
            total_p += count * mol["particle_count"]
            total_c += count * mol["component_count"]

        info["n_particles"] = total_p
        # sum_i n_i N_i == total_particle_count
        if info["total_particle_count"] is not None and total_p != info["total_particle_count"]:
            rep.add("TOPO-INVARIANT-PARTICLES", Severity.ERROR, "topology",
                    f"sum over molecule blocks of count x particle_count is {total_p}, "
                    f"but total_particle_count is {info['total_particle_count']}",
                    computed=total_p, declared=info["total_particle_count"])
        if info["total_component_count"] is not None and total_c != info["total_component_count"]:
            rep.add("TOPO-INVARIANT-COMPONENTS", Severity.ERROR, "topology",
                    f"sum over molecule blocks of count x component_count is {total_c}, "
                    f"but total_component_count is {info['total_component_count']}",
                    computed=total_c, declared=info["total_component_count"])

        # The trajectory must describe the same number of particles.
        self._check_particle_count_vs_trajectory(f, info, rep)
        return info

    def _check_particle_count_vs_trajectory(self, f, info, rep):
        particles = f.get("particles")
        if particles is None or info["n_particles"] is None:
            return
        for sysname in particles:
            pos = particles[sysname].get("position")
            if pos is None:
                continue
            value = pos["value"] if _is_time_dependent(pos) else pos
            if not isinstance(value, h5py.Dataset) or value.ndim < 2:
                continue
            n = value.shape[-2]
            if n != info["n_particles"]:
                rep.add("TOPO-TRAJECTORY-MISMATCH", Severity.ERROR,
                        f"particles/{sysname}/position",
                        f"trajectory stores {n} particles but the topology describes "
                        f"{info['n_particles']}",
                        trajectory=n, topology=info["n_particles"])

    def _collect_parameter_counts(self, f):
        """Map interaction type name -> number of force-field parameter entries."""
        counts = {}
        ff = f.get("force_field")
        if ff is None:
            return counts
        for layer in ("bonded", "nonbonded"):
            grp = ff.get(layer)
            if isinstance(grp, h5py.Group):
                for tname in grp:
                    ds = grp[tname].get("parameters")
                    if isinstance(ds, h5py.Dataset):
                        counts[tname] = ds.shape[0]
        for tname in ff:
            if tname in ("bonded", "nonbonded"):
                continue
            node = ff[tname]
            if isinstance(node, h5py.Group):
                ds = node.get("parameters")
                if isinstance(ds, h5py.Dataset):
                    counts[tname] = ds.shape[0]
        return counts

    def _check_molecule_type(self, mol, path, rep, ff_params):
        summary = {"particle_count": 0, "component_count": 0}

        for attr in ("particle_count", "component_count"):
            if attr not in mol.attrs:
                rep.add("STRUCT-MISSING-ATTR", Severity.ERROR, path,
                        f"required attribute {attr!r} is absent")
            else:
                summary[attr] = int(np.atleast_1d(mol.attrs[attr])[0])

        n_particles = summary["particle_count"]
        n_components = summary["component_count"]

        pnames = mol.get("particle_name_table")
        cnames = mol.get("component_name_table")
        n_pnames = pnames.shape[0] if isinstance(pnames, h5py.Dataset) else 0
        n_cnames = cnames.shape[0] if isinstance(cnames, h5py.Dataset) else 0
        for nm, ds in (("particle_name_table", pnames), ("component_name_table", cnames)):
            if ds is None:
                rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/{nm}",
                        f"required dataset {nm!r} is absent")
            elif _dtype_class(ds.dtype) != spec.STRING:
                rep.add("STRUCT-DTYPE", Severity.ERROR, f"{path}/{nm}",
                        f"expected String, found {_dtype_class(ds.dtype)}")

        # --- particle_data ---------------------------------------------------
        pdata = mol.get("particle_data")
        if pdata is None:
            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/particle_data",
                    "required dataset particle_data is absent")
        else:
            self._check_particle_data(pdata, f"{path}/particle_data", rep,
                                      n_particles, n_components, n_pnames)

        # --- sequence --------------------------------------------------------
        seq = mol.get("sequence")
        if seq is None:
            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/sequence",
                    "required dataset sequence is absent")
        else:
            self._check_sequence(seq, f"{path}/sequence", rep, n_components, n_cnames)

        # --- interactions ----------------------------------------------------
        inter = mol.get("interactions")
        if inter is None and self.report_deviations:
            # S10 marks the group mandatory; monatomic species legitimately have
            # no bonded interactions and the writer omits it.
            rep.add("DEV-TOPO-INTERACTIONS-REQUIRED", Severity.DEVIATION, path,
                    DEVIATIONS["DEV-TOPO-INTERACTIONS-REQUIRED"].title,
                    deviation="DEV-TOPO-INTERACTIONS-REQUIRED",
                    particle_count=n_particles)
        if isinstance(inter, h5py.Group):
            for tname in inter:
                self._check_interaction(inter[tname], f"{path}/interactions/{tname}",
                                        tname, rep, n_particles, ff_params)

        # --- exclusions ------------------------------------------------------
        exc = mol.get("exclusions")
        if isinstance(exc, h5py.Group):
            self._check_exclusions(exc, f"{path}/exclusions", rep, n_particles)
        return summary

    def _check_exclusions(self, grp, path, rep, n_particles):
        """Validate the per-particle exclusion list.

        S10 shows the two datasets but defines neither (DEV-EXCL-UNDOCUMENTED).
        The encoding in use is a compressed per-particle list: ``ranges`` has one
        row per particle, and row i is a 1-based, right-closed span of offsets
        into ``indices``; ``indices`` holds 0-based particle indices.
        """
        idx = grp.get("indices")
        rng = grp.get("ranges")
        if idx is None or rng is None:
            for nm, ds in (("indices", idx), ("ranges", rng)):
                if ds is None:
                    rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/{nm}",
                            f"required dataset {nm!r} is absent")
            return
        if self.report_deviations:
            rep.add("DEV-EXCL-UNDOCUMENTED", Severity.DEVIATION, path,
                    DEVIATIONS["DEV-EXCL-UNDOCUMENTED"].title,
                    deviation="DEV-EXCL-UNDOCUMENTED")

        indices = np.asarray(idx[()], dtype=np.int64).ravel()
        # Exclusions name particles of this molecule type.
        self._check_index_range(indices, n_particles, f"{path}/indices", rep,
                                "REF-EXCL-INDEX", "excluded particle index")

        ranges = np.asarray(rng[()], dtype=np.int64)
        if ranges.ndim != 2 or ranges.shape[1] != 3:
            rep.add("STRUCT-SHAPE", Severity.ERROR, f"{path}/ranges",
                    f"ranges must be Integer[N][3], found shape {list(ranges.shape)}",
                    expected="[N,3]", actual=list(ranges.shape))
            return
        if ranges.shape[0] != n_particles:
            rep.add("TOPO-EXCLUSION-ROWS", Severity.ERROR, f"{path}/ranges",
                    f"ranges has {ranges.shape[0]} rows for {n_particles} particle(s); "
                    "one row per particle is expected",
                    expected=n_particles, actual=int(ranges.shape[0]))
            return

        first, last, stride = ranges[:, 0], ranges[:, 1], ranges[:, 2]
        n_idx = indices.size
        # 1-based, right-closed offsets into ``indices``.
        bad = np.flatnonzero((first < 1) | (last > n_idx) | (first > last + 1))
        if bad.size:
            rep.add("REF-EXCL-RANGE", Severity.ERROR, f"{path}/ranges",
                    f"{bad.size} row(s) span outside the 1-based inclusive bounds "
                    f"[1, {n_idx}] of the indices dataset",
                    count=int(bad.size), bound=int(n_idx),
                    first_rows=bad[:10].tolist())
        bad = np.flatnonzero(stride != 1)
        if bad.size:
            rep.add("REF-EXCL-STRIDE", Severity.ERROR, f"{path}/ranges",
                    f"{bad.size} row(s) have a stride other than 1; the exclusion "
                    "list is a contiguous span per particle",
                    count=int(bad.size), first_rows=bad[:10].tolist())
        # The per-particle spans must tile the indices array exactly once.
        if n_particles:
            if int(first[0]) != 1:
                rep.add("REF-EXCL-COVERAGE", Severity.ERROR, f"{path}/ranges",
                        f"first span starts at {int(first[0])}, expected 1",
                        actual=int(first[0]))
            gaps = np.flatnonzero(first[1:] != last[:-1] + 1)
            if gaps.size:
                rep.add("REF-EXCL-COVERAGE", Severity.ERROR, f"{path}/ranges",
                        f"{gaps.size} discontinuity/ies between consecutive spans; "
                        "the spans must tile the indices dataset",
                        count=int(gaps.size), first_rows=gaps[:10].tolist())
            if int(last[-1]) != n_idx:
                rep.add("REF-EXCL-COVERAGE", Severity.ERROR, f"{path}/ranges",
                        f"last span ends at {int(last[-1])} but indices has {n_idx} "
                        "entries", actual=int(last[-1]), expected=int(n_idx))

    def _check_particle_data(self, ds, path, rep, n_particles, n_components, n_pnames):
        if ds.dtype.names is None:
            rep.add("STRUCT-DTYPE", Severity.ERROR, path,
                    "particle_data must be a compound dataset")
            return
        got = tuple(ds.dtype.names)
        expected = tuple(fl.name for fl in spec.PARTICLE_DATA_FIELDS)
        if got != expected:
            rep.add("STRUCT-COMPOUND-FIELDS", Severity.ERROR, path,
                    f"compound fields {list(got)} do not match the specified "
                    f"layout {list(expected)}",
                    expected=list(expected), actual=list(got))
            return

        for fl in spec.PARTICLE_DATA_FIELDS:
            actual = _dtype_class(ds.dtype[fl.name])
            if actual == fl.dtype:
                continue
            if fl.name == "name" and actual == spec.INTEGER:
                # Registered mismatch: the writer stores an index into
                # particle_name_table where the spec declares a String.
                if self.report_deviations:
                    rep.add("DEV-TOPO-NAME-TYPE", Severity.DEVIATION, f"{path}.name",
                            DEVIATIONS["DEV-TOPO-NAME-TYPE"].title,
                            deviation="DEV-TOPO-NAME-TYPE")
                continue
            rep.add("STRUCT-COMPOUND-DTYPE", Severity.ERROR, f"{path}.{fl.name}",
                    f"field {fl.name!r} expected {fl.dtype}, found {actual}",
                    expected=fl.dtype, actual=actual)

        if ds.shape[0] != n_particles:
            rep.add("TOPO-PARTICLE-DATA-LENGTH", Severity.ERROR, path,
                    f"particle_data has {ds.shape[0]} rows but particle_count is "
                    f"{n_particles}", expected=n_particles, actual=ds.shape[0])

        self._check_units_attr(ds, path, rep, len(got))

        data = ds[()]
        # component_id must index a component of this molecule type.
        cid = np.asarray(data["component_id"], dtype=np.int64)
        self._check_index_range(cid, n_components, f"{path}.component_id", rep,
                                "REF-COMPONENT-ID", "component_id")
        # When name is an index it must address particle_name_table.
        if _dtype_class(ds.dtype["name"]) == spec.INTEGER and n_pnames:
            nm = np.asarray(data["name"], dtype=np.int64)
            self._check_index_range(nm, n_pnames, f"{path}.name", rep,
                                    "REF-PARTICLE-NAME-INDEX", "particle name index")
        # Physical sanity: masses are non-negative.
        mass = np.asarray(data["mass"], dtype=np.float64)
        bad = np.flatnonzero(mass < 0)
        if bad.size:
            rep.add("TOPO-NEGATIVE-MASS", Severity.ERROR, f"{path}.mass",
                    f"{bad.size} particle(s) have negative mass",
                    count=int(bad.size), first_indices=bad[:10].tolist())

    def _check_sequence(self, ds, path, rep, n_components, n_cnames):
        if ds.dtype.names is None:
            rep.add("STRUCT-DTYPE", Severity.ERROR, path,
                    "sequence must be a compound dataset")
            return
        got = tuple(ds.dtype.names)
        expected = tuple(fl.name for fl in spec.SEQUENCE_FIELDS)
        if got != expected:
            rep.add("STRUCT-COMPOUND-FIELDS", Severity.ERROR, path,
                    f"compound fields {list(got)} do not match {list(expected)}",
                    expected=list(expected), actual=list(got))
            return
        if ds.shape[0] != n_components:
            rep.add("TOPO-SEQUENCE-LENGTH", Severity.ERROR, path,
                    f"sequence has {ds.shape[0]} rows but component_count is "
                    f"{n_components}", expected=n_components, actual=ds.shape[0])
        data = ds[()]
        if n_cnames:
            idx = np.asarray(data["component_name"], dtype=np.int64)
            self._check_index_range(idx, n_cnames, f"{path}.component_name", rep,
                                    "REF-COMPONENT-NAME-INDEX", "component name index")

    def _check_interaction(self, grp, path, tname, rep, n_particles, ff_params):
        itype = spec.INTERACTION_TYPES.get(tname)
        if itype is None:
            # An interaction list for a type outside Table SI is permitted as an
            # extension; its arity and parameter layout cannot be checked here,
            # only that it references a parameter table that exists.
            rep.add("STRUCT-EXTENSION", Severity.INFO, path,
                    f"{tname!r} is outside the predefined set of Table SI; its "
                    "layout is not checked", name=tname)
            if "parameter_id" in grp and tname not in ff_params:
                rep.add("REF-PARAMETER-TABLE-MISSING", Severity.ERROR,
                        f"{path}/parameter_id",
                        f"extension type {tname!r} has an interaction list but no "
                        "force_field parameters dataset", name=tname)
            return
        if not itype.bonded:
            rep.add("STRUCT-FF-CLASSIFICATION", Severity.ERROR, path,
                    f"{tname!r} is a nonbonded type and has no per-molecule "
                    "interaction list")
        pid = grp.get("parameter_id")
        ids = grp.get("particle_ids")
        if pid is None:
            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/parameter_id",
                    "required dataset parameter_id is absent")
        if ids is None:
            rep.add("STRUCT-MISSING-DATASET", Severity.ERROR, f"{path}/particle_ids",
                    "required dataset particle_ids is absent")
        if pid is None or ids is None:
            return

        n_entries = pid.shape[0]
        if ids.shape[0] != n_entries:
            rep.add("STRUCT-INTERACTION-LENGTH", Severity.ERROR, path,
                    f"parameter_id has {n_entries} entries but particle_ids has "
                    f"{ids.shape[0]}", expected=n_entries, actual=ids.shape[0])
        if itype.n_particles is not None and ids.ndim == 2:
            if ids.shape[1] != itype.n_particles:
                rep.add("STRUCT-INTERACTION-ARITY", Severity.ERROR,
                        f"{path}/particle_ids",
                        f"{tname!r} relates {itype.n_particles} particles but "
                        f"particle_ids has {ids.shape[1]} columns "
                        "(arity assumed; see DEV-SPEC-GAP-ARITY)",
                        expected=itype.n_particles, actual=ids.shape[1])

        # parameter_id must address an existing force-field parameter entry.
        n_params = ff_params.get(tname)
        pid_data = np.asarray(pid[()], dtype=np.int64).ravel()
        if n_params is None:
            if pid_data.size:
                rep.add("REF-PARAMETER-TABLE-MISSING", Severity.ERROR, f"{path}/parameter_id",
                        f"interaction list references {tname!r} parameters but no "
                        f"force_field parameters dataset exists for it", name=tname)
        else:
            self._check_index_range(pid_data, n_params, f"{path}/parameter_id", rep,
                                    "REF-PARAMETER-ID", "parameter_id")

        # particle_ids are molecule-local (see DEV-TOPO-PARTICLE-ID-SCOPE).
        id_data = np.asarray(ids[()], dtype=np.int64).ravel()
        self._check_index_range(id_data, n_particles, f"{path}/particle_ids", rep,
                                "REF-PARTICLE-ID", "particle_id")

    def _check_index_container(self, grp, path, rep, n_particles, prefix):
        """Validate a (ranges, indices) pair against an upper bound."""
        idx = grp.get("indices")
        if isinstance(idx, h5py.Dataset) and idx.size:
            data = np.asarray(idx[()], dtype=np.int64).ravel()
            self._check_index_range(data, n_particles, f"{path}/indices", rep,
                                    f"REF-{prefix}-INDEX", "index")
        rng = grp.get("ranges")
        if isinstance(rng, h5py.Dataset) and rng.size:
            data = np.asarray(rng[()], dtype=np.int64)
            if data.ndim != 2 or data.shape[1] != 3:
                rep.add("STRUCT-SHAPE", Severity.ERROR, f"{path}/ranges",
                        f"ranges must be Integer[N][3], found shape {list(data.shape)}",
                        expected="[N,3]", actual=list(data.shape))
                return
            first, last, stride = data[:, 0], data[:, 1], data[:, 2]
            bad = np.flatnonzero((first < 0) | (first >= n_particles))
            if bad.size:
                rep.add(f"REF-{prefix}-RANGE", Severity.ERROR, f"{path}/ranges",
                        f"{bad.size} range(s) start outside [0, {n_particles})",
                        count=int(bad.size), first_rows=bad[:10].tolist())
            # S9 says 'last' is exclusive; Algorithm S1 says inclusive.  Accept
            # both bounds here and report the ambiguity once, per the registry.
            bad = np.flatnonzero((last < 0) | (last > n_particles))
            if bad.size:
                rep.add(f"REF-{prefix}-RANGE", Severity.ERROR, f"{path}/ranges",
                        f"{bad.size} range(s) end outside [0, {n_particles}]",
                        count=int(bad.size), first_rows=bad[:10].tolist())
            bad = np.flatnonzero(stride == 0)
            if bad.size:
                rep.add(f"REF-{prefix}-STRIDE", Severity.ERROR, f"{path}/ranges",
                        f"{bad.size} range(s) have a zero stride",
                        count=int(bad.size), first_rows=bad[:10].tolist())

    def _check_index_range(self, data, upper, path, rep, code, label):
        if data.size == 0 or upper is None:
            return
        lo, hi = int(data.min()), int(data.max())
        if lo < 0 or hi >= upper:
            bad = np.flatnonzero((data < 0) | (data >= upper))
            rep.add(code, Severity.ERROR, path,
                    f"{bad.size} {label} value(s) outside [0, {upper}); "
                    f"observed range [{lo}, {hi}]",
                    count=int(bad.size), bound=int(upper),
                    observed_min=lo, observed_max=hi,
                    first_offsets=bad[:10].tolist())

    # -------------------------------------------------------- S5 connectivity
    def _check_connectivity(self, f, rep, topo):
        conn = f.get("connectivity")
        if conn is None:
            return
        bonds = conn.get("bonds")
        if bonds is None:
            return
        if _is_time_dependent(bonds):
            self._check_time_alignment(bonds, "connectivity/bonds", rep)
            bonds = bonds["value"]
        if not isinstance(bonds, h5py.Dataset):
            return
        if _dtype_class(bonds.dtype) != spec.INTEGER:
            rep.add("STRUCT-DTYPE", Severity.ERROR, "connectivity/bonds",
                    f"expected Integer, found {_dtype_class(bonds.dtype)}")
            return
        if bonds.ndim != 2 or bonds.shape[-1] != 2:
            rep.add("STRUCT-SHAPE", Severity.ERROR, "connectivity/bonds",
                    f"bonds must be Integer[N][2], found shape {list(bonds.shape)}",
                    expected="[N,2]", actual=list(bonds.shape))
            return
        if topo and topo.get("n_particles"):
            data = np.asarray(bonds[()], dtype=np.int64).ravel()
            self._check_index_range(data, topo["n_particles"], "connectivity/bonds",
                                    rep, "CONN-PARTICLE-ID", "global particle index")
        data2 = np.asarray(bonds[()], dtype=np.int64)
        self_bonded = np.flatnonzero(data2[:, 0] == data2[:, 1])
        if self_bonded.size:
            rep.add("CONN-SELF-BOND", Severity.ERROR, "connectivity/bonds",
                    f"{self_bonded.size} bond(s) connect a particle to itself",
                    count=int(self_bonded.size),
                    first_rows=self_bonded[:10].tolist())

    # ------------------------------------------------------------ S9 selections
    def _check_selections(self, f, rep, topo):
        sel = f.get("selections")
        if sel is None:
            return
        n = topo.get("n_particles") if topo else None
        reported_convention = False
        for name in sel:
            grp = sel[name]
            path = f"selections/{name}"
            if not isinstance(grp, h5py.Group):
                rep.add("STRUCT-NOT-GROUP", Severity.ERROR, path,
                        "a selection must be a group")
                continue
            keys = set(grp)
            index_based = bool(keys & {"ranges", "indices"})
            string_based = bool(keys & {"selection", "software"})
            if not index_based and not string_based:
                rep.add("SEL-EMPTY", Severity.ERROR, path,
                        "selection defines neither ranges/indices nor a selection string")
            if index_based and string_based:
                rep.add("SEL-MIXED", Severity.WARNING, path,
                        "selection mixes index-based and string-based definitions")
            if index_based:
                if n is not None:
                    self._check_index_container(grp, path, rep, n, "SEL")
                if "ranges" in grp and grp["ranges"].size and self.report_deviations \
                        and not reported_convention:
                    rep.add("DEV-SEL-RANGE-CONVENTION", Severity.DEVIATION, path,
                            DEVIATIONS["DEV-SEL-RANGE-CONVENTION"].title,
                            deviation="DEV-SEL-RANGE-CONVENTION")
                    reported_convention = True
            if string_based and "selection" in keys and "software" not in keys:
                rep.add("SEL-NO-SOFTWARE", Severity.WARNING, path,
                        "string-based selection does not name the software that "
                        "can parse it")


def validate(path, **kwargs) -> Report:
    """Convenience wrapper: validate a single file."""
    return Validator(**kwargs).validate(path)
