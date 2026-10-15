"""
Declarative model of the H5MD ``biophysics`` module, version 1.0.

Normative source: the module specification, version 1.0 (``SI-Specification.pdf``
in the paper's supplementary material); section numbers in the comments below
(S4, S5, ...) refer to it. Nothing here is derived from the GROMACS
reference implementation -- disagreements are recorded in
:mod:`h5md_validator.deviations` rather than absorbed into the model.

Used by:
    Validation 1: schema coverage
    Validation 3: internal consistency
"""

from __future__ import annotations

from dataclasses import dataclass

# --------------------------------------------------------------------------
# Type vocabulary (S3)
# --------------------------------------------------------------------------
# The specification names four scalar types plus compounds.  ``Boolean`` is not
# a native HDF5 type: S3.3 defines it as an Integer where 0 is false and any
# non-zero value is true, so it validates as an Integer with a value check.
INTEGER = "Integer"
FLOAT = "Float"
STRING = "String"
BOOLEAN = "Boolean"
COMPOUND = "Compound"

#: Shape sentinel for a dimension whose extent the specification leaves free.
VAR = "var"


@dataclass(frozen=True)
class Attr:
    """An HDF5 attribute (``\\++ name`` in the specification's notation)."""

    name: str
    dtype: str
    shape: tuple | None = None
    required: bool = True
    enum: frozenset | None = None
    doc: str = ""


@dataclass(frozen=True)
class Field:
    """One field of a compound dataset."""

    name: str
    dtype: str


@dataclass(frozen=True)
class Dataset:
    """An HDF5 dataset (``\\-- name: <type>[N]`` in the specification)."""

    name: str
    dtype: str
    shape: tuple | None = None
    required: bool = True
    enum: frozenset | None = None
    #: True when S4.2 requires a ``unit`` attribute (physical quantity).
    unit: bool = False
    #: Compound layout, for ``dtype == COMPOUND``.
    fields: tuple | None = None
    #: True when the compound requires a ``units`` attribute (S4.2).
    units: bool = False
    attrs: tuple = ()
    doc: str = ""


@dataclass(frozen=True)
class Group:
    """An HDF5 group (``\\-- name`` in the specification)."""

    name: str
    required: bool = True
    attrs: tuple = ()
    datasets: tuple = ()
    groups: tuple = ()
    #: True for ``<item>``: a group whose *name* is not fixed by the spec.
    variadic: bool = False
    doc: str = ""


# ==========================================================================
# S7.1  Run control
# ==========================================================================
SIMULATION_TYPES = frozenset(
    {"molecular_dynamics", "stochastic_dynamics", "energy_minimization"}
)
MD_INTEGRATORS = frozenset({"leap_frog", "velocity_verlet", "euler"})
SD_INTEGRATORS = frozenset({"leap_frog", "langevin", "brownian"})
EM_ALGORITHMS = frozenset({"steepest_descent", "conjugate_gradient", "l_bfgs"})
COM_MODES = frozenset({"translation", "rotation", "translation_and_rotation"})

RUN_CONTROL = Group(
    "run_control",
    required=False,
    datasets=(
        Dataset("simulation_type", STRING, enum=SIMULATION_TYPES),
        Dataset("numsteps", INTEGER),
    ),
    groups=(
        Group(
            "molecular_dynamics",
            required=False,
            datasets=(
                Dataset("integrator", STRING, enum=MD_INTEGRATORS),
                Dataset("timestep", FLOAT, unit=True),
            ),
        ),
        Group(
            "stochastic_dynamics",
            required=False,
            datasets=(
                Dataset("integrator", STRING, enum=SD_INTEGRATORS),
                Dataset("timestep", FLOAT, unit=True),
            ),
        ),
        Group(
            "center_of_mass_motion_removal",
            required=False,
            datasets=(
                Dataset("groups", STRING, shape=("N",)),
                Dataset("mode", STRING, enum=COM_MODES),
                Dataset("interval", INTEGER),
            ),
        ),
        Group(
            "energy_minimization",
            required=False,
            datasets=(
                Dataset("algorithm", STRING, enum=EM_ALGORITHMS),
                Dataset("step_size", FLOAT, unit=True),
                Dataset("tolerance", FLOAT, unit=True),
            ),
        ),
        Group(
            "mass_repartitioning",
            required=False,
            datasets=(Dataset("factor", FLOAT),),
        ),
    ),
)

# ==========================================================================
# S7.2  Output control
# ==========================================================================
OUTPUT_CONTROL = Group(
    "output_control",
    required=False,
    datasets=(
        Dataset("position", INTEGER),
        Dataset("velocity", INTEGER),
        Dataset("force", INTEGER),
        Dataset("energy", INTEGER),
        Dataset("log", INTEGER),
    ),
)

# ==========================================================================
# S7.3  Temperature coupling
# ==========================================================================
THERMOSTATS = frozenset(
    {"stochastic_velocity_rescaling", "andersen", "nose_hoover", "berendsen"}
)

TEMPERATURE_COUPLING = Group(
    "temperature_coupling",
    required=False,
    datasets=(
        Dataset("algorithm", STRING, enum=THERMOSTATS),
        Dataset("interval", INTEGER),
        Dataset("groups", STRING, shape=("N",)),
        Dataset("reference_temperature", FLOAT, shape=("N",), unit=True),
    ),
    groups=(
        Group(
            "stochastic_velocity_rescaling",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("coupling_time", FLOAT, shape=("N",), unit=True),),
        ),
        Group(
            "andersen",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("coupling_time", FLOAT, shape=("N",), unit=True),),
        ),
        Group(
            "nose_hoover",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(
                Dataset("period_time", FLOAT, shape=("N",), unit=True),
                Dataset("chain_length", INTEGER),
            ),
        ),
        Group(
            "berendsen",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("coupling_time", FLOAT, shape=("N",), unit=True),),
        ),
    ),
)

# ==========================================================================
# S7.4  Pressure coupling
# ==========================================================================
BAROSTATS = frozenset({"stochastic_cell_rescaling", "parrinello_rahman", "berendsen"})

PRESSURE_COUPLING = Group(
    "pressure_coupling",
    required=False,
    datasets=(
        Dataset("algorithm", STRING, enum=BAROSTATS),
        # S7.4: scalar if isotropic, 3x3 tensor if anisotropic.
        Dataset("reference_pressure", FLOAT, shape=("scalar_or_tensor",), unit=True),
        Dataset("compressibility", FLOAT, shape=("scalar_or_tensor",), unit=True),
        Dataset("interval", INTEGER),
        Dataset("reference_coordinate_scaling", STRING),
    ),
    groups=(
        Group(
            "stochastic_cell_rescaling",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("coupling_time", FLOAT, unit=True),),
        ),
        Group(
            "parrinello_rahman",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("period_time", FLOAT, unit=True),),
        ),
        Group(
            "berendsen",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("coupling_time", FLOAT, unit=True),),
        ),
    ),
)

# ==========================================================================
# S7.5  Coulomb interactions
# ==========================================================================
COULOMB_TYPES = frozenset({"pme", "fmm", "pppm", "ewald", "reaction_field", "cutoff"})
COULOMB_MODIFIERS = frozenset({"none", "potential_shift"})

COULOMB_INTERACTIONS = Group(
    "coulomb_interactions",
    required=False,
    datasets=(
        Dataset("type", STRING, enum=COULOMB_TYPES),
        Dataset("cutoff_distance", FLOAT, unit=True),
        Dataset("relative_permittivity", FLOAT, required=False),
    ),
    groups=(
        Group(
            "cutoff",
            required=False,
            datasets=(
                Dataset("potential_modifier", STRING, enum=COULOMB_MODIFIERS),
                Dataset("switch_distance", FLOAT, required=False, unit=True),
            ),
        ),
        Group(
            "reaction_field",
            required=False,
            datasets=(Dataset("relative_permittivity", FLOAT),),
        ),
    ),
)

# ==========================================================================
# S7.6  Van der Waals interactions
# ==========================================================================
VDW_TYPES = frozenset({"pme", "cutoff"})
VDW_MODIFIERS = frozenset(
    {"none", "potential_shift", "potential_switch", "force_switch"}
)
DISPERSION_CORRECTIONS = frozenset({"energy", "pressure", "energy_and_pressure"})

VDW_INTERACTIONS = Group(
    "vdw_interactions",
    required=False,
    datasets=(
        Dataset("type", STRING, enum=VDW_TYPES),
        Dataset("cutoff_distance", FLOAT, unit=True),
        Dataset(
            "dispersion_correction", STRING, required=False, enum=DISPERSION_CORRECTIONS
        ),
    ),
    groups=(
        Group(
            "cutoff",
            required=False,
            datasets=(
                Dataset("potential_modifier", STRING, enum=VDW_MODIFIERS),
                Dataset("switch_distance", FLOAT, required=False, unit=True),
            ),
        ),
    ),
)

# ==========================================================================
# S7.7  Boundary conditions
# ==========================================================================
BOUNDARY_CONDITIONS = Group(
    "boundary_conditions",
    required=False,
    datasets=(Dataset("periodicity", BOOLEAN, shape=(3,)),),
)

# ==========================================================================
# S7.8  Constraints
# ==========================================================================
CONSTRAINT_ALGORITHMS = frozenset({"lincs", "shake"})
CONSTRAINT_MODES = frozenset({"hydrogens", "all_bonds"})

CONSTRAINTS = Group(
    "constraints",
    required=False,
    datasets=(
        Dataset("algorithm", STRING, enum=CONSTRAINT_ALGORITHMS),
        Dataset("mode", STRING, required=False, enum=CONSTRAINT_MODES),
    ),
    groups=(
        Group(
            "lincs",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("iterations", INTEGER), Dataset("order", INTEGER)),
        ),
        Group(
            "shake",
            required=False,
            attrs=(Attr("doi", STRING, required=False),),
            datasets=(Dataset("tolerance", FLOAT),),
        ),
    ),
)

PARAMETERS = Group(
    "parameters",
    required=False,
    groups=(
        RUN_CONTROL,
        OUTPUT_CONTROL,
        TEMPERATURE_COUPLING,
        PRESSURE_COUPLING,
        COULOMB_INTERACTIONS,
        VDW_INTERACTIONS,
        BOUNDARY_CONDITIONS,
        CONSTRAINTS,
    ),
)

#: Parameter groups whose ``algorithm``/``type``/``simulation_type`` dataset
#: selects which sibling subgroup must be instantiated (S7).  Maps the group
#: path to (selector dataset name, set of selector values that own a subgroup).
CONDITIONAL_SELECTORS = {
    "parameters/run_control": ("simulation_type", SIMULATION_TYPES),
    "parameters/temperature_coupling": ("algorithm", THERMOSTATS),
    "parameters/pressure_coupling": ("algorithm", BAROSTATS),
    "parameters/constraints": ("algorithm", CONSTRAINT_ALGORITHMS),
    # S7.5/S7.6: only ``cutoff`` and ``reaction_field`` own a subgroup; the
    # remaining long-range treatments (pme, fmm, ...) do not instantiate one.
    "parameters/coulomb_interactions": (
        "type",
        frozenset({"cutoff", "reaction_field"}),
    ),
    "parameters/vdw_interactions": ("type", frozenset({"cutoff"})),
}

#: The full predefined value set of each selector dataset (S7).  This is a
#: superset of CONDITIONAL_SELECTORS: ``coulomb_interactions/type = pme`` is a
#: predefined value even though ``pme`` instantiates no subgroup.
SELECTOR_ENUMS = {
    "parameters/run_control": ("simulation_type", SIMULATION_TYPES),
    "parameters/temperature_coupling": ("algorithm", THERMOSTATS),
    "parameters/pressure_coupling": ("algorithm", BAROSTATS),
    "parameters/constraints": ("algorithm", CONSTRAINT_ALGORITHMS),
    "parameters/coulomb_interactions": ("type", COULOMB_TYPES),
    "parameters/vdw_interactions": ("type", VDW_TYPES),
}

#: Datasets that must have the same length as the ``groups`` dataset of their
#: enclosing coupling category (S7.3: "one value per group").
PER_GROUP_DATASETS = {
    "parameters/temperature_coupling": [
        "reference_temperature",
        "stochastic_velocity_rescaling/coupling_time",
        "andersen/coupling_time",
        "nose_hoover/period_time",
        "berendsen/coupling_time",
    ],
}

# ==========================================================================
# S6  Observables
# ==========================================================================
#: Predefined observables and the shape of a single record (S6, Table SI).
#: All are optional and may be stored time-dependent or time-independent.
OBSERVABLES = {
    "density": (),
    "potential_energy": (),
    "kinetic_energy": (),
    "total_energy": (),
    "conserved_energy": (),
    "enthalpy": (),
    "temperature": (),
    "pressure": (),
    "pressure_tensor": (3, 3),
    "virial_tensor": (3, 3),
    "volume": (),
}


# ==========================================================================
# S8  Force field terms
# ==========================================================================
#: Interaction types, their compound parameter layout, and whether they are
#: bonded (referenced from ``topology/<mol>/interactions``) or nonbonded.
#:
#: ``n_particles`` is the number of particles an interaction entry relates, used
#: to check the second dimension of ``particle_ids``.  NOTE: the specification's
#: Table SI (S8) does not state these arities; see deviations.SPEC_GAP_ARITY.
@dataclass(frozen=True)
class InteractionType:
    name: str
    params: tuple
    bonded: bool
    n_particles: int | None
    #: Alternative compound layouts the spec permits ("c6/c12 or sigma/epsilon").
    alt_params: tuple | None = None


F = FLOAT

INTERACTION_TYPES = {
    it.name: it
    for it in (
        InteractionType(
            "bond",
            (Field("equilibrium_distance", F), Field("force_constant", F)),
            True,
            2,
        ),
        InteractionType(
            "fourth_power_bond",
            (Field("equilibrium_squared_distance", F), Field("force_constant", F)),
            True,
            2,
        ),
        InteractionType(
            "angle",
            (Field("equilibrium_angle", F), Field("force_constant", F)),
            True,
            3,
        ),
        InteractionType(
            "cosine_harmonic_angle",
            (Field("equilibrium_cosine", F), Field("force_constant", F)),
            True,
            3,
        ),
        InteractionType(
            "proper_dihedral",
            (
                Field("phase", F),
                Field("force_constant", F),
                Field("multiplicity", INTEGER),
            ),
            True,
            4,
        ),
        InteractionType(
            "improper_dihedral",
            (Field("equilibrium_angle", F), Field("force_constant", F)),
            True,
            4,
        ),
        InteractionType(
            "periodic_improper_dihedral",
            (
                Field("phase", F),
                Field("force_constant", F),
                Field("multiplicity", INTEGER),
            ),
            True,
            4,
        ),
        InteractionType(
            "ryckaert_bellemans_dihedral",
            tuple(Field(f"c{i}", F) for i in range(6)),
            True,
            4,
        ),
        InteractionType(
            "urey_bradley",
            (
                Field("equilibrium_angle", F),
                Field("angle_force_constant", F),
                Field("equilibrium_distance", F),
                Field("bond_force_constant", F),
            ),
            True,
            3,
        ),
        InteractionType(
            "morse_bond",
            (
                Field("equilibrium_distance", F),
                Field("well_depth", F),
                Field("steepness", F),
            ),
            True,
            2,
        ),
        InteractionType("constraint", (Field("distance", F),), True, 2),
        InteractionType(
            "settle", (Field("distance_OH", F), Field("distance_HH", F)), True, 3
        ),
        InteractionType(
            "lennard_jones_14",
            (Field("c6", F), Field("c12", F)),
            True,
            2,
            alt_params=(Field("sigma", F), Field("epsilon", F)),
        ),
        InteractionType("virtual_site_3", (Field("a", F), Field("b", F)), True, 4),
        InteractionType(
            "lennard_jones_short_range",
            (Field("c6", F), Field("c12", F)),
            False,
            None,
            alt_params=(Field("sigma", F), Field("epsilon", F)),
        ),
        InteractionType(
            "buckingham",
            (Field("constant_a", F), Field("constant_b", F), Field("constant_c", F)),
            False,
            None,
        ),
    )
}

#: Optional provenance attributes on the ``force_field`` group (S8).
FORCE_FIELD_ATTRS = (
    Attr("source", STRING, required=False),
    Attr("modifications", STRING, required=False),
    Attr("polymer_types", STRING, shape=(VAR,), required=False),
    Attr("polymer_force_field", STRING, required=False),
    Attr("solvent_model", STRING, required=False),
    Attr("solvent_component", STRING, required=False),
    Attr("solvent_force_field", STRING, required=False),
    Attr("ion_type", STRING, required=False),
    Attr("ion_force_field", STRING, required=False),
    Attr("membrane_components", STRING, shape=(VAR,), required=False),
    Attr("membrane_force_field", STRING, required=False),
    Attr("ligand_components", STRING, shape=(VAR,), required=False),
    Attr("ligand_force_field", STRING, required=False),
)

# ==========================================================================
# S10  Topology
# ==========================================================================
PARTICLE_DATA_FIELDS = (
    Field("name", STRING),
    Field("species", INTEGER),
    Field("charge", FLOAT),
    Field("mass", FLOAT),
    Field("component_id", INTEGER),
)

SEQUENCE_FIELDS = (
    Field("component_number", INTEGER),
    Field("component_name", INTEGER),
)

MOLECULE_BLOCK_FIELDS = (
    Field("name", STRING),
    Field("count", INTEGER),
)

# ==========================================================================
# S4.1  Root
# ==========================================================================
#: Root groups the module either defines or extends.  All are optional; the
#: only mandatory element is ``h5md/modules/biophysics`` with its version.
ROOT_GROUPS = (
    "connectivity",
    "observables",
    "parameters",
    "particles",
    "force_field",
    "selections",
    "topology",
)

MODULE_NAME = "biophysics"
MODULE_VERSION = (1, 0)
