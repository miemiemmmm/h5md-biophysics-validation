"""
Conformance and internal-consistency validator for the ``biophysics`` module.

Checks run on top of a structural conformance pass against the normative schema
model in :mod:`h5md_validator.spec`.

Check families, by code prefix:

    MOD     module declaration and version
    STRUCT  presence, type and shape of specified elements
    ENUM    values of enumerated elements
    UNIT    presence and syntax of unit / units attributes
    COND    conditional subgroup instantiation (S7)
    TIME    step / time / value first-dimension alignment
    TOPO    topology invariants (particle and residue counts, 1-4 pairing)
    REF     referential integrity of every index and selection name
    SEL     selection bounds and coverage
    CONN    connectivity bounds
    EXT     extensions outside a predefined set (S7, S8)

Used by:
    Validation 3: internal consistency
    Validation 4: negative tests
"""

from __future__ import annotations

import math

import h5py
import numpy as np

from . import spec
from . import units as unitmod
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
        # h5py writes NumPy bools as an HDF5 enum {FALSE, TRUE}; the
        # specification's Boolean is an Integer (S3.3), so this is not one.
        return "Enum (h5py bool)"
    return f"unknown({kind})"


def _is_time_dependent(obj) -> bool:
    """H5MD time-dependent element: a group of step, (time) and value (S3.2).

    A group holding any of the three counts, so that an element missing its step
    or value is checked as a broken time-dependent element, not mistaken for
    something else.
    """
    return isinstance(obj, h5py.Group) and bool({"step", "time", "value"} & set(obj))


class Validator:
    """Validates one H5MD file against the ``biophysics`` v1.0 specification."""

    def __init__(self, *, strict_units: bool = True):
        self.strict_units = strict_units

    # ----------------------------------------------------------------- entry
    def validate(self, path) -> Report:
        report = Report(file=str(path))
        try:
            handle = h5py.File(path, "r")
        except OSError as exc:
            # A file the HDF5 library cannot open is a validation verdict, not a
            # crash: a writer that aborts mid-write can leave exactly this.
            report.add(
                "FILE-UNREADABLE",
                Severity.ERROR,
                "",
                f"file cannot be opened as HDF5: {exc}",
                error=str(exc),
            )
            return report
        with handle as f:
            self._check_module(f, report)
            self._check_parameters(f, report)
            self._check_observables(f, report)
            self._check_force_field(f, report)
            topo = self._check_topology(f, report)
            self._check_connectivity(f, report, topo)
            self._check_selections(f, report, topo)
            self._check_selection_references(f, report, topo)
        return report

    # ------------------------------------------------------------- S4.1 root
    def _check_module(self, f, rep):
        modpath = f"h5md/modules/{spec.MODULE_NAME}"
        mod = f.get(modpath)
        if mod is None:
            rep.add(
                "MOD-MISSING",
                Severity.ERROR,
                modpath,
                f"file does not declare the {spec.MODULE_NAME} module",
            )
            return
        if "version" not in mod.attrs:
            rep.add(
                "MOD-VERSION-MISSING",
                Severity.ERROR,
                modpath,
                f"{spec.MODULE_NAME} module has no version attribute",
            )
            return
        version = tuple(int(v) for v in np.atleast_1d(mod.attrs["version"]))
        if len(version) != 2:
            rep.add(
                "MOD-VERSION-SHAPE",
                Severity.ERROR,
                modpath,
                f"version must be Integer[2], got {len(version)} element(s)",
                actual=list(version),
            )
        elif version != spec.MODULE_VERSION:
            rep.add(
                "MOD-VERSION-UNSUPPORTED",
                Severity.ERROR,
                modpath,
                f"unsupported module version {version[0]}.{version[1]}, "
                f"validator implements {spec.MODULE_VERSION[0]}.{spec.MODULE_VERSION[1]}",
                actual=list(version),
                expected=list(spec.MODULE_VERSION),
            )

    # ------------------------------------------------------------ S4.2 units
    def _check_unit_attr(self, ds, path, rep, required: bool):
        if "unit" not in ds.attrs:
            if required:
                rep.add(
                    "UNIT-MISSING",
                    Severity.ERROR,
                    path,
                    "dataset holds a physical quantity but has no unit attribute",
                )
            return
        unit = _decode(ds.attrs["unit"])
        if not isinstance(unit, str):
            rep.add(
                "UNIT-TYPE",
                Severity.ERROR,
                path,
                "unit attribute must be a String",
                actual=str(type(unit)),
            )
            return
        if not unitmod.is_valid(unit):
            rep.add(
                "UNIT-SYNTAX",
                Severity.ERROR,
                path,
                f"unit {unit!r} is not valid H5MD unit notation",
                unit=unit,
            )
        elif self.strict_units:
            unknown = unitmod.unknown_symbols(unit)
            if unknown:
                rep.add(
                    "UNIT-SYMBOL-UNKNOWN",
                    Severity.WARNING,
                    path,
                    f"unit {unit!r} uses unrecognised symbol(s) {unknown}",
                    unit=unit,
                    symbols=unknown,
                )

    def _check_units_attr(self, ds, path, rep, n_fields: int):
        """S4.2: a compound dataset carries a ``units`` list, one per field."""
        if "units" not in ds.attrs:
            rep.add(
                "UNITS-MISSING",
                Severity.ERROR,
                path,
                "compound dataset has no units attribute",
            )
            return
        vals = _decode(ds.attrs["units"])
        vals = [vals] if isinstance(vals, str) else list(vals)
        if len(vals) != n_fields:
            rep.add(
                "UNITS-LENGTH",
                Severity.ERROR,
                path,
                f"units has {len(vals)} entries for {n_fields} compound fields",
                expected=n_fields,
                actual=len(vals),
            )
        for u in vals:
            if u and not unitmod.is_valid(u):
                rep.add(
                    "UNIT-SYNTAX",
                    Severity.ERROR,
                    path,
                    f"unit {u!r} is not valid H5MD unit notation",
                    unit=u,
                )

    # ------------------------------------------------- generic tree conformance
    def _check_against_model(self, parent, node: spec.Group, path: str, rep):
        """Recursively check a file group against its node in the spec.py model."""
        obj = parent.get(node.name)
        if obj is None:
            if node.required:
                rep.add(
                    "STRUCT-MISSING-GROUP",
                    Severity.ERROR,
                    path,
                    f"required group {node.name!r} is absent",
                )
            return
        if not isinstance(obj, h5py.Group):
            rep.add(
                "STRUCT-NOT-GROUP",
                Severity.ERROR,
                path,
                f"{node.name!r} is specified as a group but is a dataset",
            )
            return

        for attr in node.attrs:
            if attr.name not in obj.attrs and attr.required:
                rep.add(
                    "STRUCT-MISSING-ATTR",
                    Severity.ERROR,
                    path,
                    f"required attribute {attr.name!r} is absent",
                )

        for dsnode in node.datasets:
            self._check_dataset(obj, dsnode, f"{path}/{dsnode.name}", rep)

        for sub in node.groups:
            self._check_against_model(obj, sub, f"{path}/{sub.name}", rep)

        # Elements present in the file that the specification does not define.
        specified = {d.name for d in node.datasets} | {g.name for g in node.groups}
        for key in obj:
            if key not in specified:
                self._report_undocumented(f"{path}/{key}", rep, node.name)

    def _report_undocumented(self, path, rep, parent_name):
        """An element present in the file that the specification does not define."""
        rep.add(
            "STRUCT-UNDOCUMENTED",
            Severity.WARNING,
            path,
            f"element is not defined by the specification for {parent_name!r}",
        )

    def _check_dataset(self, parent, node: spec.Dataset, path: str, rep):
        obj = parent.get(node.name)
        if obj is None:
            if node.required:
                rep.add(
                    "STRUCT-MISSING-DATASET",
                    Severity.ERROR,
                    path,
                    f"required dataset {node.name!r} is absent",
                )
            return
        if isinstance(obj, h5py.Group):
            # Permitted when the element is stored time-dependently.
            if _is_time_dependent(obj):
                self._check_time_alignment(obj, path, rep)
                return
            rep.add(
                "STRUCT-NOT-DATASET",
                Severity.ERROR,
                path,
                f"{node.name!r} is specified as a dataset but is a group",
            )
            return

        actual = _dtype_class(obj.dtype)
        expected = spec.INTEGER if node.dtype == spec.BOOLEAN else node.dtype
        want = "Boolean (Integer 0/1)" if node.dtype == spec.BOOLEAN else node.dtype
        if actual != expected:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                path,
                f"expected {want}, found {actual}",
                expected=node.dtype,
                actual=actual,
            )
        else:
            self._check_shape(obj, node, path, rep)
            if node.dtype == spec.BOOLEAN:
                bad = sorted({int(v) for v in np.atleast_1d(obj[()]).ravel()} - {0, 1})
                if bad:
                    # S3.3 admits any non-zero as true, so this is informational.
                    rep.add(
                        "ENUM-BOOLEAN-NONCANONICAL",
                        Severity.INFO,
                        path,
                        f"Boolean stored with non-canonical value(s) {bad}",
                        values=bad,
                    )
            if node.enum is not None:
                self._check_enum(obj, node, path, rep)

        self._check_unit_attr(obj, path, rep, required=node.unit)

    def _check_shape(self, obj, node: spec.Dataset, path, rep):
        want = node.shape
        got = tuple(obj.shape)
        if want is None:
            if got != ():
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    path,
                    f"expected a scalar, found shape {list(got)}",
                    expected="scalar",
                    actual=list(got),
                )
        elif want == ("scalar_or_tensor",):
            if got not in ((), (3, 3)):
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    path,
                    f"expected a scalar or a 3x3 tensor, found shape {list(got)}",
                    expected="scalar or [3,3]",
                    actual=list(got),
                )
        elif all(isinstance(d, int) for d in want):
            if got != want:
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    path,
                    f"expected shape {list(want)}, found {list(got)}",
                    expected=list(want),
                    actual=list(got),
                )
        else:
            if len(got) != len(want):
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    path,
                    f"expected a {len(want)}-D dataset, found {len(got)}-D",
                    expected=len(want),
                    actual=len(got),
                )

    def _check_enum(self, obj, node: spec.Dataset, path, rep):
        value = _read(obj)
        vals = value if isinstance(value, list) else [value]
        for v in vals:
            if not isinstance(v, str):
                continue
            if v not in node.enum:
                # S7 permits algorithms outside the predefined set, provided a
                # doi/url points at a definition.  That is a WARNING, not an error.
                rep.add(
                    "ENUM-UNDEFINED-VALUE",
                    Severity.WARNING,
                    path,
                    f"value {v!r} is outside the predefined set "
                    f"{sorted(node.enum)}; a doi or url attribute is recommended",
                    value=v,
                    predefined=sorted(node.enum),
                )

    # ------------------------------------------------- S3.2 time-dependent data
    def _check_time_alignment(self, grp, path, rep):
        """N_step == N_time == N_value on the first dimension (S3.2)."""
        lengths = {}
        for key in ("step", "time", "value"):
            if key in grp and isinstance(grp[key], h5py.Dataset):
                shape = grp[key].shape
                lengths[key] = shape[0] if shape else None
        if "step" not in lengths:
            rep.add(
                "TIME-STEP-MISSING",
                Severity.ERROR,
                path,
                "time-dependent element has no step dataset",
            )
        if "value" not in lengths:
            rep.add(
                "TIME-VALUE-MISSING",
                Severity.ERROR,
                path,
                "time-dependent element has no value dataset",
            )
        present = {k: v for k, v in lengths.items() if v is not None}
        if len(set(present.values())) > 1:
            rep.add(
                "TIME-LENGTH-MISMATCH",
                Severity.ERROR,
                path,
                "first dimension of step/time/value datasets disagree: "
                + ", ".join(f"{k}={v}" for k, v in sorted(present.items())),
                lengths=present,
            )
        # S3.2: step is Integer[variable] and time is Float[variable].
        for key, want in (("step", spec.INTEGER), ("time", spec.FLOAT)):
            ds = grp.get(key)
            if not isinstance(ds, h5py.Dataset):
                continue
            got = _dtype_class(ds.dtype)
            if got != want:
                rep.add(
                    "STRUCT-DTYPE",
                    Severity.ERROR,
                    f"{path}/{key}",
                    f"{key} must be {want}[variable], found {got}",
                    expected=want,
                    actual=got,
                )
            if ds.ndim != 1:
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    f"{path}/{key}",
                    f"{key} must be one-dimensional, found shape {list(ds.shape)}",
                    expected="[variable]",
                    actual=list(ds.shape),
                )
        for key in ("step", "time"):
            if key in grp and isinstance(grp[key], h5py.Dataset):
                data = np.asarray(grp[key][()]).ravel()
                if data.size > 1 and np.any(np.diff(data) < 0):
                    rep.add(
                        "TIME-NOT-MONOTONIC",
                        Severity.ERROR,
                        f"{path}/{key}",
                        f"{key} values are not monotonically increasing",
                    )
        if "time" in grp and isinstance(grp["time"], h5py.Dataset):
            self._check_unit_attr(grp["time"], f"{path}/time", rep, required=True)

    # ---------------------------------------------------------- S7 parameters
    def _check_parameters(self, f, rep):
        params = f.get("parameters")
        if params is None:
            return
        self._check_against_model(f, spec.PARAMETERS, "parameters", rep)

        # S7: the subgroup named by the selector dataset must be instantiated.
        for gpath, (selector, owning) in spec.CONDITIONAL_SELECTORS.items():
            grp = f.get(gpath)
            if grp is None or selector not in grp:
                continue
            value = _read(grp[selector])
            if not isinstance(value, str):
                continue
            if value in owning and value not in grp:
                rep.add(
                    "COND-SUBGROUP-MISSING",
                    Severity.ERROR,
                    gpath,
                    f"{selector} is {value!r} but no {value!r} subgroup is present",
                    selector=selector,
                    value=value,
                )
            # A subgroup for an algorithm that was not selected is contradictory.
            for name in owning:
                if name != value and name in grp and isinstance(grp[name], h5py.Group):
                    rep.add(
                        "COND-SUBGROUP-ORPHAN",
                        Severity.ERROR,
                        f"{gpath}/{name}",
                        f"subgroup {name!r} is present but {selector} is {value!r}",
                        selector=selector,
                        value=value,
                        orphan=name,
                    )

        # S7.1, S7.3: per-group arrays must match the length of 'groups'; where
        # 'groups' is optional and absent, one group of every particle is implied.
        for gpath, names in spec.PER_GROUP_DATASETS.items():
            grp = f.get(gpath)
            if grp is None:
                continue
            if "groups" in grp:
                n = len(np.atleast_1d(grp["groups"][()]))
            elif gpath in spec.GROUPS_OPTIONAL:
                n = 1
            else:
                continue
            for rel in names:
                ds = grp.get(rel)
                if ds is None or not isinstance(ds, h5py.Dataset):
                    continue
                got = ds.shape[0] if ds.shape else 1
                if got != n:
                    rep.add(
                        "STRUCT-PER-GROUP-LENGTH",
                        Severity.ERROR,
                        f"{gpath}/{rel}",
                        f"has {got} value(s) for {n} coupling group(s)",
                        expected=n,
                        actual=got,
                    )

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
                    rep.add(
                        "EXT-NO-REFERENCE",
                        Severity.WARNING,
                        f"{gpath}/{value}"
                        if sub is not None
                        else f"{gpath}/{selector}",
                        f"{value!r} is not a predefined value; a doi or url "
                        "attribute is recommended",
                        value=value,
                    )

        # S7.1: stochastic dynamics controls the temperature itself.
        rc = f.get("parameters/run_control")
        if (
            rc is not None
            and "simulation_type" in rc
            and _read(rc["simulation_type"]) == "stochastic_dynamics"
            and "temperature_coupling" in params
        ):
            rep.add(
                "COND-SD-TEMPERATURE-COUPLING",
                Severity.ERROR,
                "parameters/temperature_coupling",
                "temperature_coupling is present, but stochastic dynamics controls "
                "the temperature itself",
            )
        # S7.1: Brownian friction is given in exactly one of two ways.
        bd = f.get("parameters/run_control/stochastic_dynamics/brownian")
        if isinstance(bd, h5py.Group):
            given = [k for k in ("coupling_time", "friction_coefficient") if k in bd]
            if len(given) != 1:
                rep.add(
                    "COND-BROWNIAN-FRICTION",
                    Severity.ERROR,
                    "parameters/run_control/stochastic_dynamics/brownian",
                    "exactly one of coupling_time and friction_coefficient must be "
                    f"present, found {given or 'neither'}",
                    present=given,
                )

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
                if isinstance(node.get("value"), h5py.Dataset):
                    self._check_observable_shape(name, node["value"], path, rep)
                    self._check_unit_attr(
                        node["value"], f"{path}/value", rep, required=True
                    )
            elif isinstance(node, h5py.Dataset):
                self._check_observable_shape(name, node, path, rep)
                self._check_unit_attr(node, path, rep, required=True)
            elif isinstance(node, h5py.Group):
                # S6: observables of a subset of particles sit in a subgroup named
                # after its selection group, which must not shadow an observable.
                # That the selection exists is checked with the other references.
                if name in spec.OBSERVABLES:
                    rep.add(
                        "STRUCT-OBS-SUBGROUP-NAME",
                        Severity.ERROR,
                        path,
                        f"{name!r} is a predefined observable, so it cannot name a "
                        "selection subgroup",
                        name=name,
                    )
                for sub in node:
                    subnode = node[sub]
                    subpath = f"{path}/{sub}"
                    if _is_time_dependent(subnode):
                        self._check_time_alignment(subnode, subpath, rep)
                        if isinstance(subnode.get("value"), h5py.Dataset):
                            self._check_observable_shape(
                                sub, subnode["value"], subpath, rep
                            )
                            self._check_unit_attr(
                                subnode["value"],
                                f"{subpath}/value",
                                rep,
                                required=True,
                            )
                    elif isinstance(subnode, h5py.Dataset):
                        self._check_observable_shape(sub, subnode, subpath, rep)
                        self._check_unit_attr(subnode, subpath, rep, required=True)

    def _check_observable_shape(self, name, ds, path, rep):
        want = spec.OBSERVABLES.get(name)
        if want is None:
            rep.add(
                "STRUCT-UNDOCUMENTED",
                Severity.WARNING,
                path,
                f"{name!r} is not one of the predefined observables of S6",
                name=name,
            )
            return
        got = tuple(ds.shape)
        # Drop the leading frame axis for time-dependent storage.
        record = got[1:] if len(got) > len(want) else got
        if tuple(record) != tuple(want):
            rep.add(
                "STRUCT-SHAPE",
                Severity.ERROR,
                path,
                f"observable {name!r} expects record shape {list(want)}, "
                f"found {list(record)}",
                expected=list(want),
                actual=list(record),
            )

    # --------------------------------------------------------- S8 force field
    def _check_force_field(self, f, rep):
        ff = f.get("force_field")
        if ff is None:
            return {}

        # S8: every interaction type sits under bonded or nonbonded.
        for key in ff:
            if key in spec.FORCE_FIELD_LAYERS:
                continue
            if isinstance(ff[key], h5py.Group):
                rep.add(
                    "STRUCT-FF-CLASSIFICATION",
                    Severity.ERROR,
                    f"force_field/{key}",
                    f"{key!r} is stored directly under force_field, outside "
                    "bonded and nonbonded",
                    name=key,
                )
            else:
                self._report_undocumented(f"force_field/{key}", rep, "force_field")

        for container, base in self._force_field_containers(ff):
            for tname in container:
                tgrp = container[tname]
                path = f"{base}/{tname}"
                itype = spec.INTERACTION_TYPES.get(tname)
                if itype is None:
                    # S8: a type outside Table SI is permitted in the same layout,
                    # and should carry a doi or url defining its functional form.
                    if not ({"doi", "url"} & set(tgrp.attrs)):
                        rep.add(
                            "EXT-NO-REFERENCE",
                            Severity.WARNING,
                            path,
                            f"{tname!r} is outside the predefined set of "
                            "Table SI; an extension is permitted, but a doi or "
                            "url attribute naming the functional form is "
                            "recommended",
                            name=tname,
                        )
                    else:
                        rep.add(
                            "EXT-REFERENCED",
                            Severity.INFO,
                            path,
                            f"{tname!r} is an extension outside Table SI, with "
                            "a reference attached",
                            name=tname,
                        )
                    continue
                if itype.bonded and base.endswith("nonbonded"):
                    rep.add(
                        "STRUCT-FF-CLASSIFICATION",
                        Severity.ERROR,
                        path,
                        f"{tname!r} is a bonded interaction type but is stored "
                        "under nonbonded",
                    )
                if not itype.bonded and base.endswith("/bonded"):
                    rep.add(
                        "STRUCT-FF-CLASSIFICATION",
                        Severity.ERROR,
                        path,
                        f"{tname!r} is a nonbonded interaction type but is stored "
                        "under bonded",
                    )
                for pname in ("parameters", "parameters_end_state"):
                    if pname not in tgrp:
                        if pname == "parameters":
                            rep.add(
                                "STRUCT-MISSING-DATASET",
                                Severity.ERROR,
                                f"{path}/parameters",
                                "interaction type has no parameters dataset",
                            )
                        continue
                    if itype.grid:
                        self._check_parameter_grid(tgrp[pname], f"{path}/{pname}", rep)
                    else:
                        self._check_parameter_compound(
                            tgrp[pname], f"{path}/{pname}", itype, rep
                        )
                    if not itype.bonded:
                        self._check_type_matrix(tgrp[pname], f"{path}/{pname}", rep)

        # S8: coulomb_14 is present exactly when lennard_jones_14 is.
        bonded = ff.get("bonded")
        if isinstance(bonded, h5py.Group) and (
            ("lennard_jones_14" in bonded) != ("coulomb_14" in bonded)
        ):
            rep.add(
                "TOPO-COULOMB14-PAIRING",
                Severity.ERROR,
                "force_field/bonded",
                "coulomb_14 must be present exactly when lennard_jones_14 is",
                lennard_jones_14="lennard_jones_14" in bonded,
                coulomb_14="coulomb_14" in bonded,
            )

    def _force_field_containers(self, ff):
        """Yield (group, path) pairs that hold interaction-type subgroups."""
        for layer in spec.FORCE_FIELD_LAYERS:
            if layer in ff and isinstance(ff[layer], h5py.Group):
                yield ff[layer], f"force_field/{layer}"

    def _check_parameter_grid(self, ds, path, rep):
        """S8: a tabulated potential (cmap) is Float[K][G][G] with one unit."""
        if _dtype_class(ds.dtype) != spec.FLOAT:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                path,
                f"expected Float, found {_dtype_class(ds.dtype)}",
                expected=spec.FLOAT,
                actual=_dtype_class(ds.dtype),
            )
        if ds.ndim != 3 or ds.shape[1] != ds.shape[2]:
            rep.add(
                "STRUCT-SHAPE",
                Severity.ERROR,
                path,
                f"expected Float[K][G][G], found shape {list(ds.shape)}",
                expected="[K,G,G]",
                actual=list(ds.shape),
            )
        self._check_unit_attr(ds, path, rep, required=True)

    def _check_type_matrix(self, ds, path, rep):
        """S8: nonbonded parameters form a T x T matrix, K = T^2 records."""
        k = ds.shape[0] if ds.shape else 0
        t = math.isqrt(k)
        if t * t != k:
            rep.add(
                "STRUCT-NONBONDED-MATRIX",
                Severity.ERROR,
                path,
                f"{k} records do not form a square matrix of particle types",
                records=k,
            )

    def _check_parameter_compound(self, ds, path, itype, rep):
        if ds.dtype.names is None:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                path,
                "parameters must be a compound dataset",
                expected=spec.COMPOUND,
                actual=_dtype_class(ds.dtype),
            )
            return
        got = tuple(ds.dtype.names)
        accepted = [tuple(fl.name for fl in itype.params)]
        if itype.alt_params:
            accepted.append(tuple(fl.name for fl in itype.alt_params))
        if got not in accepted:
            rep.add(
                "STRUCT-COMPOUND-FIELDS",
                Severity.ERROR,
                path,
                f"compound fields {list(got)} do not match any layout "
                f"specified for {itype.name!r}: {[list(a) for a in accepted]}",
                expected=[list(a) for a in accepted],
                actual=list(got),
            )
        else:
            layout = itype.params if got == accepted[0] else itype.alt_params
            for fl in layout:
                actual = _dtype_class(ds.dtype[fl.name])
                if actual != fl.dtype:
                    rep.add(
                        "STRUCT-COMPOUND-DTYPE",
                        Severity.ERROR,
                        f"{path}.{fl.name}",
                        f"field {fl.name!r} expected {fl.dtype}, found {actual}",
                        expected=fl.dtype,
                        actual=actual,
                    )
        self._check_units_attr(ds, path, rep, len(got))

    # ------------------------------------------------------------ S10 topology
    def _check_topology(self, f, rep):
        """Validate the topology group; return the resolved topology summary."""
        topo = f.get("topology")
        if topo is None:
            return None

        info = {
            "molecule_types": {},
            "total_particle_count": None,
            "total_residue_count": None,
            "n_particles": None,
        }

        for attr in ("total_particle_count", "total_residue_count"):
            if attr not in topo.attrs:
                rep.add(
                    "STRUCT-MISSING-ATTR",
                    Severity.ERROR,
                    "topology",
                    f"required attribute {attr!r} is absent",
                )
            else:
                info[attr] = int(np.atleast_1d(topo.attrs[attr])[0])
        # system_name is an optional attribute (S10); nothing to check when absent.
        if "system_name" in topo.attrs and not isinstance(
            _decode(topo.attrs["system_name"]), str
        ):
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                "topology/@system_name",
                "system_name must be a String",
                expected=spec.STRING,
            )

        ff_params = self._collect_parameter_counts(f)
        nb = f.get("force_field/nonbonded")
        nonbonded = sorted(nb) if isinstance(nb, h5py.Group) else []

        # --- molecule types -------------------------------------------------
        for name in topo:
            if name == "molecule_blocks":
                continue
            mol = topo[name]
            if not isinstance(mol, h5py.Group):
                continue
            info["molecule_types"][name] = self._check_molecule_type(
                mol, f"topology/{name}", rep, ff_params, nonbonded
            )

        # --- molecule_blocks and the S10 count invariants ---------------------
        blocks = topo.get("molecule_blocks")
        if blocks is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                "topology/molecule_blocks",
                "required dataset molecule_blocks is absent",
            )
            return info
        if blocks.dtype.names is None:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                "topology/molecule_blocks",
                "molecule_blocks must be a compound dataset",
            )
            return info
        expected_fields = tuple(fl.name for fl in spec.MOLECULE_BLOCK_FIELDS)
        if tuple(blocks.dtype.names) != expected_fields:
            rep.add(
                "STRUCT-COMPOUND-FIELDS",
                Severity.ERROR,
                "topology/molecule_blocks",
                f"compound fields {list(blocks.dtype.names)} do not match "
                f"{list(expected_fields)}",
                expected=list(expected_fields),
                actual=list(blocks.dtype.names),
            )
            return info

        data = blocks[()]
        total_p = total_r = 0
        for i, row in enumerate(data):
            bname = _decode(row["name"])
            count = int(row["count"])
            if count < 0:
                rep.add(
                    "TOPO-BLOCK-COUNT-NEGATIVE",
                    Severity.ERROR,
                    f"topology/molecule_blocks[{i}]",
                    f"molecule block {bname!r} has negative count {count}",
                    index=i,
                    name=bname,
                    count=count,
                )
            mol = info["molecule_types"].get(bname)
            if mol is None:
                rep.add(
                    "REF-BLOCK-DANGLING",
                    Severity.ERROR,
                    f"topology/molecule_blocks[{i}]",
                    f"molecule block references molecule type {bname!r}, "
                    "which has no group under topology/",
                    index=i,
                    name=bname,
                )
                continue
            total_p += count * mol["particle_count"]
            total_r += count * mol["residue_count"]

        info["n_particles"] = total_p
        # sum_i n_i N_i == total_particle_count
        if (
            info["total_particle_count"] is not None
            and total_p != info["total_particle_count"]
        ):
            rep.add(
                "TOPO-INVARIANT-PARTICLES",
                Severity.ERROR,
                "topology",
                f"sum over molecule blocks of count x particle_count is {total_p}, "
                f"but total_particle_count is {info['total_particle_count']}",
                computed=total_p,
                declared=info["total_particle_count"],
            )
        if (
            info["total_residue_count"] is not None
            and total_r != info["total_residue_count"]
        ):
            rep.add(
                "TOPO-INVARIANT-RESIDUES",
                Severity.ERROR,
                "topology",
                f"sum over molecule blocks of count x residue_count is {total_r}, "
                f"but total_residue_count is {info['total_residue_count']}",
                computed=total_r,
                declared=info["total_residue_count"],
            )

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
            value = pos.get("value") if _is_time_dependent(pos) else pos
            if not isinstance(value, h5py.Dataset) or value.ndim < 2:
                continue
            n = value.shape[-2]
            if n != info["n_particles"]:
                rep.add(
                    "TOPO-TRAJECTORY-MISMATCH",
                    Severity.ERROR,
                    f"particles/{sysname}/position",
                    f"trajectory stores {n} particles but the topology describes "
                    f"{info['n_particles']}",
                    trajectory=n,
                    topology=info["n_particles"],
                )

    def _collect_parameter_counts(self, f):
        """Map interaction type name -> K, its number of parameter records."""
        counts = {}
        ff = f.get("force_field")
        if ff is None:
            return counts
        for layer in spec.FORCE_FIELD_LAYERS:
            grp = ff.get(layer)
            if isinstance(grp, h5py.Group):
                for tname in grp:
                    ds = grp[tname].get("parameters")
                    if isinstance(ds, h5py.Dataset) and ds.shape:
                        counts[tname] = ds.shape[0]
        return counts

    def _check_molecule_type(self, mol, path, rep, ff_params, nonbonded=()):
        summary = {"particle_count": 0, "residue_count": 0}

        for attr in ("particle_count", "residue_count"):
            if attr not in mol.attrs:
                rep.add(
                    "STRUCT-MISSING-ATTR",
                    Severity.ERROR,
                    path,
                    f"required attribute {attr!r} is absent",
                )
            else:
                summary[attr] = int(np.atleast_1d(mol.attrs[attr])[0])

        n_particles = summary["particle_count"]
        n_residues = summary["residue_count"]

        pnames = mol.get("particle_name_table")
        cnames = mol.get("component_name_table")
        n_pnames = pnames.shape[0] if isinstance(pnames, h5py.Dataset) else 0
        n_cnames = cnames.shape[0] if isinstance(cnames, h5py.Dataset) else 0
        for nm, ds in (
            ("particle_name_table", pnames),
            ("component_name_table", cnames),
        ):
            if ds is None:
                rep.add(
                    "STRUCT-MISSING-DATASET",
                    Severity.ERROR,
                    f"{path}/{nm}",
                    f"required dataset {nm!r} is absent",
                )
            elif _dtype_class(ds.dtype) != spec.STRING:
                rep.add(
                    "STRUCT-DTYPE",
                    Severity.ERROR,
                    f"{path}/{nm}",
                    f"expected String, found {_dtype_class(ds.dtype)}",
                )

        # --- particle_data ---------------------------------------------------
        pdata = mol.get("particle_data")
        if pdata is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/particle_data",
                "required dataset particle_data is absent",
            )
        else:
            self._check_particle_data(
                pdata, f"{path}/particle_data", rep, n_particles, n_residues, n_pnames
            )

        # --- sequence --------------------------------------------------------
        seq = mol.get("sequence")
        if seq is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/sequence",
                "required dataset sequence is absent",
            )
        else:
            self._check_sequence(seq, f"{path}/sequence", rep, n_residues, n_cnames)

        # --- interactions ----------------------------------------------------
        # Optional (S10): absent when the molecule type has no interaction.
        inter = mol.get("interactions")
        if isinstance(inter, h5py.Group):
            for tname in inter:
                self._check_interaction(
                    inter[tname],
                    f"{path}/interactions/{tname}",
                    tname,
                    rep,
                    n_particles,
                    ff_params,
                )
        # S8: every particle carries a type for each nonbonded interaction.
        for tname in nonbonded:
            if not (isinstance(inter, h5py.Group) and tname in inter):
                rep.add(
                    "STRUCT-MISSING-DATASET",
                    Severity.ERROR,
                    f"{path}/interactions/{tname}/type_id",
                    f"no type_id for the nonbonded interaction {tname!r}",
                    name=tname,
                )
        # S8: coulomb_14 row k applies to lennard_jones_14 row k.
        lj14 = inter.get("lennard_jones_14") if isinstance(inter, h5py.Group) else None
        c14 = inter.get("coulomb_14") if isinstance(inter, h5py.Group) else None
        if (lj14 is None) != (c14 is None):
            rep.add(
                "TOPO-COULOMB14-PAIRING",
                Severity.ERROR,
                f"{path}/interactions",
                "coulomb_14 must be present exactly when lennard_jones_14 is",
                lennard_jones_14=lj14 is not None,
                coulomb_14=c14 is not None,
            )
        elif lj14 is not None and "parameter_id" in lj14 and "parameter_id" in c14:
            n_lj, n_c = lj14["parameter_id"].shape[0], c14["parameter_id"].shape[0]
            if n_lj != n_c:
                rep.add(
                    "TOPO-COULOMB14-PAIRING",
                    Severity.ERROR,
                    f"{path}/interactions/coulomb_14/parameter_id",
                    f"coulomb_14 has {n_c} entries for {n_lj} lennard_jones_14 pairs",
                    expected=n_lj,
                    actual=n_c,
                )

        # --- exclusions ------------------------------------------------------
        exc = mol.get("exclusions")
        if isinstance(exc, h5py.Group):
            self._check_exclusions(exc, f"{path}/exclusions", rep, n_particles)
        return summary

    def _check_exclusions(self, grp, path, rep, n_particles):
        """Validate the excluded particle pairs.

        S10: one dataset, ``particle_ids`` of shape [N_e][2], holding every
        excluded pair once as 0-based particle indices of this molecule type
        with i < j; the exclusion of a particle from itself is implied. The
        group is absent when the molecule type has no excluded pair.
        """
        ds = grp.get("particle_ids")
        if ds is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/particle_ids",
                "required dataset 'particle_ids' is absent",
            )
            return

        pairs = np.asarray(ds[()], dtype=np.int64)
        if pairs.ndim != 2 or pairs.shape[1] != 2:
            rep.add(
                "STRUCT-SHAPE",
                Severity.ERROR,
                f"{path}/particle_ids",
                f"particle_ids must be Integer[P][2], found shape {list(pairs.shape)}",
                expected="[P,2]",
                actual=list(pairs.shape),
            )
            return
        if pairs.shape[0] == 0:
            rep.add(
                "TOPO-EXCLUSION-EMPTY",
                Severity.WARNING,
                f"{path}/particle_ids",
                "the exclusions group holds no pair; a molecule type without "
                "excluded pairs is expected to have no exclusions group",
            )
            return

        # Exclusions name particles of this molecule type.
        self._check_index_range(
            pairs.ravel(),
            n_particles,
            f"{path}/particle_ids",
            rep,
            "REF-EXCL-INDEX",
            "excluded particle index",
        )
        # Each pair is stored once, lower index first, which also rules out a
        # particle excluding itself.
        bad = np.flatnonzero(pairs[:, 0] >= pairs[:, 1])
        if bad.size:
            rep.add(
                "REF-EXCL-ORDER",
                Severity.ERROR,
                f"{path}/particle_ids",
                f"{bad.size} pair(s) do not satisfy i < j",
                count=int(bad.size),
                first_rows=bad[:10].tolist(),
            )
        # Compare unordered, so (i, j) and (j, i) count as the same pair.
        n_unique = np.unique(np.sort(pairs, axis=1), axis=0).shape[0]
        if n_unique != pairs.shape[0]:
            rep.add(
                "REF-EXCL-DUPLICATE",
                Severity.ERROR,
                f"{path}/particle_ids",
                f"{pairs.shape[0] - n_unique} pair(s) are listed more than once",
                count=int(pairs.shape[0] - n_unique),
            )

    def _check_particle_data(self, ds, path, rep, n_particles, n_residues, n_pnames):
        if ds.dtype.names is None:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                path,
                "particle_data must be a compound dataset",
            )
            return
        got = tuple(ds.dtype.names)
        expected = tuple(fl.name for fl in spec.PARTICLE_DATA_FIELDS)
        if got != expected:
            rep.add(
                "STRUCT-COMPOUND-FIELDS",
                Severity.ERROR,
                path,
                f"compound fields {list(got)} do not match the specified "
                f"layout {list(expected)}",
                expected=list(expected),
                actual=list(got),
            )
            return

        for fl in spec.PARTICLE_DATA_FIELDS:
            actual = _dtype_class(ds.dtype[fl.name])
            if actual == fl.dtype:
                continue
            rep.add(
                "STRUCT-COMPOUND-DTYPE",
                Severity.ERROR,
                f"{path}.{fl.name}",
                f"field {fl.name!r} expected {fl.dtype}, found {actual}",
                expected=fl.dtype,
                actual=actual,
            )

        if ds.shape[0] != n_particles:
            rep.add(
                "TOPO-PARTICLE-DATA-LENGTH",
                Severity.ERROR,
                path,
                f"particle_data has {ds.shape[0]} rows but particle_count is "
                f"{n_particles}",
                expected=n_particles,
                actual=ds.shape[0],
            )

        self._check_units_attr(ds, path, rep, len(got))

        data = ds[()]
        # sequence_index must address a row of this molecule type's sequence.
        cid = np.asarray(data["sequence_index"], dtype=np.int64)
        self._check_index_range(
            cid,
            n_residues,
            f"{path}.sequence_index",
            rep,
            "REF-SEQUENCE-INDEX",
            "sequence_index",
        )
        # name is an index into particle_name_table (S10).
        if _dtype_class(ds.dtype["name"]) == spec.INTEGER and n_pnames:
            nm = np.asarray(data["name"], dtype=np.int64)
            self._check_index_range(
                nm,
                n_pnames,
                f"{path}.name",
                rep,
                "REF-PARTICLE-NAME-INDEX",
                "particle name index",
            )
        # Physical sanity: masses are non-negative.
        mass = np.asarray(data["mass"], dtype=np.float64)
        bad = np.flatnonzero(mass < 0)
        if bad.size:
            rep.add(
                "TOPO-NEGATIVE-MASS",
                Severity.ERROR,
                f"{path}.mass",
                f"{bad.size} particle(s) have negative mass",
                count=int(bad.size),
                first_indices=bad[:10].tolist(),
            )

    def _check_sequence(self, ds, path, rep, n_residues, n_cnames):
        if ds.dtype.names is None:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                path,
                "sequence must be a compound dataset",
            )
            return
        got = tuple(ds.dtype.names)
        expected = tuple(fl.name for fl in spec.SEQUENCE_FIELDS)
        if got != expected:
            rep.add(
                "STRUCT-COMPOUND-FIELDS",
                Severity.ERROR,
                path,
                f"compound fields {list(got)} do not match {list(expected)}",
                expected=list(expected),
                actual=list(got),
            )
            return
        if ds.shape[0] != n_residues:
            rep.add(
                "TOPO-SEQUENCE-LENGTH",
                Severity.ERROR,
                path,
                f"sequence has {ds.shape[0]} rows but residue_count is {n_residues}",
                expected=n_residues,
                actual=ds.shape[0],
            )
        data = ds[()]
        if n_cnames:
            idx = np.asarray(data["component_name"], dtype=np.int64)
            self._check_index_range(
                idx,
                n_cnames,
                f"{path}.component_name",
                rep,
                "REF-COMPONENT-NAME-INDEX",
                "component name index",
            )

    def _check_interaction(self, grp, path, tname, rep, n_particles, ff_params):
        itype = spec.INTERACTION_TYPES.get(tname)
        if itype is None:
            # An interaction list for a type outside Table SI is permitted as an
            # extension; its arity and parameter layout cannot be checked here,
            # only that it references a parameter table that exists.
            rep.add(
                "EXT-LAYOUT-UNCHECKED",
                Severity.INFO,
                path,
                f"{tname!r} is outside the predefined set of Table SI; its "
                "layout is not checked",
                name=tname,
            )
            if "parameter_id" in grp and tname not in ff_params:
                rep.add(
                    "REF-PARAMETER-TABLE-MISSING",
                    Severity.ERROR,
                    f"{path}/parameter_id",
                    f"extension type {tname!r} has an interaction list but no "
                    "force_field parameters dataset",
                    name=tname,
                )
            return
        if not itype.bonded:
            self._check_type_ids(grp, path, tname, rep, n_particles, ff_params)
            return
        pid = grp.get("parameter_id")
        ids = grp.get("particle_ids")
        if pid is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/parameter_id",
                "required dataset parameter_id is absent",
            )
        if not itype.has_particle_ids:
            # coulomb_14: only parameter_id, aligned with lennard_jones_14.
            if pid is not None:
                self._check_parameter_ids(pid, path, tname, rep, ff_params)
            return
        if ids is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/particle_ids",
                "required dataset particle_ids is absent",
            )
        if pid is None or ids is None:
            return

        n_entries = pid.shape[0]
        if ids.shape[0] != n_entries:
            rep.add(
                "STRUCT-INTERACTION-LENGTH",
                Severity.ERROR,
                path,
                f"parameter_id has {n_entries} entries but particle_ids has "
                f"{ids.shape[0]}",
                expected=n_entries,
                actual=ids.shape[0],
            )
        if (
            itype.n_particles is not None
            and ids.ndim == 2
            and ids.shape[1] != itype.n_particles
        ):
            rep.add(
                "STRUCT-INTERACTION-ARITY",
                Severity.ERROR,
                f"{path}/particle_ids",
                f"{tname!r} relates {itype.n_particles} particles but "
                f"particle_ids has {ids.shape[1]} columns",
                expected=itype.n_particles,
                actual=ids.shape[1],
            )

        self._check_parameter_ids(pid, path, tname, rep, ff_params)

        # particle_ids are 0-based and local to the molecule type (S10).
        id_data = np.asarray(ids[()], dtype=np.int64).ravel()
        self._check_index_range(
            id_data,
            n_particles,
            f"{path}/particle_ids",
            rep,
            "REF-PARTICLE-ID",
            "particle_id",
        )

    def _check_parameter_ids(self, pid, path, tname, rep, ff_params):
        """parameter_id must address a record of the type's parameters (S8)."""
        n_params = ff_params.get(tname)
        pid_data = np.asarray(pid[()], dtype=np.int64).ravel()
        if n_params is None:
            if pid_data.size:
                rep.add(
                    "REF-PARAMETER-TABLE-MISSING",
                    Severity.ERROR,
                    f"{path}/parameter_id",
                    f"interaction list references {tname!r} parameters but no "
                    f"force_field parameters dataset exists for it",
                    name=tname,
                )
            return
        self._check_index_range(
            pid_data,
            n_params,
            f"{path}/parameter_id",
            rep,
            "REF-PARAMETER-ID",
            "parameter_id",
        )

    def _check_type_ids(self, grp, path, tname, rep, n_particles, ff_params):
        """S8: a nonbonded type stores type_id, one per particle, in [0, T)."""
        ds = grp.get("type_id")
        if ds is None:
            rep.add(
                "STRUCT-MISSING-DATASET",
                Severity.ERROR,
                f"{path}/type_id",
                "required dataset type_id is absent",
            )
            return
        if _dtype_class(ds.dtype) != spec.INTEGER:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                f"{path}/type_id",
                f"expected Integer, found {_dtype_class(ds.dtype)}",
            )
            return
        if ds.ndim != 1:
            rep.add(
                "STRUCT-SHAPE",
                Severity.ERROR,
                f"{path}/type_id",
                f"type_id must be one-dimensional, found shape {list(ds.shape)}",
                expected="[N_p]",
                actual=list(ds.shape),
            )
        elif ds.shape[0] != n_particles:
            rep.add(
                "TOPO-TYPE-ID-LENGTH",
                Severity.ERROR,
                f"{path}/type_id",
                f"type_id has {ds.shape[0]} entries but particle_count is "
                f"{n_particles}",
                expected=n_particles,
                actual=ds.shape[0],
            )
        k = ff_params.get(tname)
        if k is None:
            rep.add(
                "REF-PARAMETER-TABLE-MISSING",
                Severity.ERROR,
                f"{path}/type_id",
                f"type_id references {tname!r} parameters but no force_field "
                "parameters dataset exists for it",
                name=tname,
            )
            return
        t = math.isqrt(k)
        if t * t == k:
            self._check_index_range(
                np.asarray(ds[()], dtype=np.int64).ravel(),
                t,
                f"{path}/type_id",
                rep,
                "REF-TYPE-ID",
                "type_id",
            )

    def _check_index_container(self, grp, path, rep, n_particles, prefix):
        """Validate a (ranges, indices) pair against an upper bound."""
        idx = grp.get("indices")
        if isinstance(idx, h5py.Dataset) and idx.size:
            data = np.asarray(idx[()], dtype=np.int64).ravel()
            self._check_index_range(
                data,
                n_particles,
                f"{path}/indices",
                rep,
                f"REF-{prefix}-INDEX",
                "index",
            )
        rng = grp.get("ranges")
        if isinstance(rng, h5py.Dataset) and rng.size:
            data = np.asarray(rng[()], dtype=np.int64)
            if data.ndim != 2 or data.shape[1] != 3:
                rep.add(
                    "STRUCT-SHAPE",
                    Severity.ERROR,
                    f"{path}/ranges",
                    f"ranges must be Integer[N][3], found shape {list(data.shape)}",
                    expected="[N,3]",
                    actual=list(data.shape),
                )
                return
            first, last, stride = data[:, 0], data[:, 1], data[:, 2]
            bad = np.flatnonzero((first < 0) | (first >= n_particles))
            if bad.size:
                rep.add(
                    f"REF-{prefix}-RANGE",
                    Severity.ERROR,
                    f"{path}/ranges",
                    f"{bad.size} range(s) start outside [0, {n_particles})",
                    count=int(bad.size),
                    first_rows=bad[:10].tolist(),
                )
            # S9: stop is an exclusive bound, so it may equal the system size.
            bad = np.flatnonzero((last < 0) | (last > n_particles))
            if bad.size:
                rep.add(
                    f"REF-{prefix}-RANGE",
                    Severity.ERROR,
                    f"{path}/ranges",
                    f"{bad.size} range(s) end outside [0, {n_particles}]",
                    count=int(bad.size),
                    first_rows=bad[:10].tolist(),
                )
            bad = np.flatnonzero(stride <= 0)
            if bad.size:
                rep.add(
                    f"REF-{prefix}-STRIDE",
                    Severity.ERROR,
                    f"{path}/ranges",
                    f"{bad.size} range(s) have a stride that is not positive",
                    count=int(bad.size),
                    first_rows=bad[:10].tolist(),
                )

    def _check_index_range(self, data, upper, path, rep, code, label):
        if data.size == 0 or upper is None:
            return
        lo, hi = int(data.min()), int(data.max())
        if lo < 0 or hi >= upper:
            bad = np.flatnonzero((data < 0) | (data >= upper))
            rep.add(
                code,
                Severity.ERROR,
                path,
                f"{bad.size} {label} value(s) outside [0, {upper}); "
                f"observed range [{lo}, {hi}]",
                count=int(bad.size),
                bound=int(upper),
                observed_min=lo,
                observed_max=hi,
                first_offsets=bad[:10].tolist(),
            )

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
            bonds = bonds.get("value")
        if not isinstance(bonds, h5py.Dataset):
            return
        if _dtype_class(bonds.dtype) != spec.INTEGER:
            rep.add(
                "STRUCT-DTYPE",
                Severity.ERROR,
                "connectivity/bonds",
                f"expected Integer, found {_dtype_class(bonds.dtype)}",
            )
            return
        if bonds.ndim != 2 or bonds.shape[-1] != 2:
            rep.add(
                "STRUCT-SHAPE",
                Severity.ERROR,
                "connectivity/bonds",
                f"bonds must be Integer[N][2], found shape {list(bonds.shape)}",
                expected="[N,2]",
                actual=list(bonds.shape),
            )
            return
        if topo and topo.get("n_particles"):
            data = np.asarray(bonds[()], dtype=np.int64).ravel()
            self._check_index_range(
                data,
                topo["n_particles"],
                "connectivity/bonds",
                rep,
                "CONN-PARTICLE-ID",
                "global particle index",
            )
        data2 = np.asarray(bonds[()], dtype=np.int64)
        self_bonded = np.flatnonzero(data2[:, 0] == data2[:, 1])
        if self_bonded.size:
            rep.add(
                "CONN-SELF-BOND",
                Severity.ERROR,
                "connectivity/bonds",
                f"{self_bonded.size} bond(s) connect a particle to itself",
                count=int(self_bonded.size),
                first_rows=self_bonded[:10].tolist(),
            )

    # ------------------------------------------------------------ S9 selections
    def _check_selections(self, f, rep, topo):
        sel = f.get("selections")
        if sel is None:
            return
        n = topo.get("n_particles") if topo else None
        for name in sel:
            grp = sel[name]
            path = f"selections/{name}"
            if not isinstance(grp, h5py.Group):
                rep.add(
                    "STRUCT-NOT-GROUP",
                    Severity.ERROR,
                    path,
                    "a selection must be a group",
                )
                continue
            keys = set(grp)
            index_based = bool(keys & {"ranges", "indices"})
            string_based = bool(keys & {"selection", "software"})
            if not index_based and not string_based:
                rep.add(
                    "SEL-EMPTY",
                    Severity.ERROR,
                    path,
                    "selection defines neither ranges/indices nor a selection string",
                )
            if index_based and string_based:
                rep.add(
                    "SEL-MIXED",
                    Severity.WARNING,
                    path,
                    "selection mixes index-based and string-based definitions",
                )
            if index_based and n is not None:
                self._check_index_container(grp, path, rep, n, "SEL")
            # A selection string is unreadable without the software that
            # parses it, so its absence is an error, not advice.
            if string_based and "selection" in keys and "software" not in keys:
                rep.add(
                    "SEL-NO-SOFTWARE",
                    Severity.ERROR,
                    path,
                    "string-based selection does not name the software that "
                    "can parse it",
                )

    def _check_selection_references(self, f, rep, topo):
        """S9: every selection name used elsewhere must exist in selections."""
        sel = f.get("selections")
        names = set(sel) if isinstance(sel, h5py.Group) else set()
        refs = []
        for dpath in spec.SELECTION_REFERENCES:
            ds = f.get(dpath)
            if isinstance(ds, h5py.Dataset):
                vals = _read(ds)
                refs += [
                    (dpath, g) for g in (vals if isinstance(vals, list) else [vals])
                ]
        obs = f.get("observables")
        if isinstance(obs, h5py.Group):
            for name in obs:
                node = obs[name]
                if isinstance(node, h5py.Group) and not _is_time_dependent(node):
                    refs.append((f"observables/{name}", name))
        for path, group in refs:
            if group not in names:
                rep.add(
                    "REF-SELECTION-MISSING",
                    Severity.ERROR,
                    path,
                    f"selection {group!r} is referenced but not defined in selections",
                    name=group,
                )

        # S7.3: the coupling groups should together cover the whole system.
        n = topo.get("n_particles") if topo else None
        tc = f.get("parameters/temperature_coupling/groups")
        if not (n and isinstance(tc, h5py.Dataset) and isinstance(sel, h5py.Group)):
            return
        vals = _read(tc)
        covered = np.zeros(n, dtype=bool)
        for g in vals if isinstance(vals, list) else [vals]:
            idx = self._selection_indices(sel.get(g), n)
            if idx is None:
                return  # missing or string-based: coverage cannot be decided
            covered[idx] = True
        if not covered.all():
            rep.add(
                "SEL-COUPLING-COVERAGE",
                Severity.WARNING,
                "parameters/temperature_coupling/groups",
                f"the temperature-coupling groups leave {int((~covered).sum())} of "
                f"{n} particles uncoupled",
                uncovered=int((~covered).sum()),
            )

    @staticmethod
    def _selection_indices(grp, n):
        """The particle indices of an index-based selection, or None."""
        if not isinstance(grp, h5py.Group) or not ({"ranges", "indices"} & set(grp)):
            return None
        parts = []
        if "indices" in grp:
            parts.append(np.asarray(grp["indices"][()], dtype=np.int64).ravel())
        if "ranges" in grp and grp["ranges"].size:
            for start, stop, stride in np.asarray(grp["ranges"][()], dtype=np.int64):
                if stride > 0:
                    parts.append(np.arange(start, stop, stride))
        idx = np.concatenate(parts) if parts else np.zeros(0, dtype=np.int64)
        return idx[(idx >= 0) & (idx < n)]


def validate(path, **kwargs) -> Report:
    """Convenience wrapper: validate a single file."""
    return Validator(**kwargs).validate(path)
