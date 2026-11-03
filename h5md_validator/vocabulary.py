"""
Translation between GROMACS vocabulary and the ``biophysics`` schema vocabulary.

This module is the concrete, auditable statement of what the reference
implementation claims when it calls the schema engine-independent: every entry
is a GROMACS-specific spelling on the left and the schema-level concept it
denotes on the right.

Used by:
    Validation 2: round-trip test
    Validation 6: cross-engine mapping
"""

from __future__ import annotations

# MD parameter -> run_control/simulation_type map
SIMULATION_TYPE = {
    "md": "molecular_dynamics",
    "md-vv": "molecular_dynamics",
    "md-vv-avek": "molecular_dynamics",
    "sd": "stochastic_dynamics",
    "bd": "stochastic_dynamics",
    "steep": "energy_minimization",
    "cg": "energy_minimization",
    "l-bfgs": "energy_minimization",
}

# MD parameter -> <simulation_type>/integrator map; sd is a leap-frog scheme and
# bd an Euler scheme, the stochastic equation itself being stochastic_dynamics/dynamics
INTEGRATOR = {
    "md": "leap_frog",
    "md-vv": "velocity_verlet",
    "md-vv-avek": "velocity_verlet",
    "sd": "leap_frog",
    "bd": "euler",
    "steep": "steepest_descent",
    "cg": "conjugate_gradient",
    "l-bfgs": "l_bfgs",
}

# MD parameter -> run_control/stochastic_dynamics/dynamics map
SD_DYNAMICS = {
    "sd": "langevin",
    "bd": "brownian",
}

# MD parameter -> temperature_coupling/algorithm map
THERMOSTAT = {
    "No": None,
    "Berendsen": "berendsen",
    "Nose-Hoover": "nose_hoover",
    "V-rescale": "stochastic_velocity_rescaling",
    "Andersen": "andersen",
    "Andersen-massive": "andersen_massive",
}

# MD parameter -> thermostat time field map
THERMOSTAT_TIME_FIELD = {
    "stochastic_velocity_rescaling": "coupling_time",
    "andersen": "coupling_time",
    "andersen_massive": "coupling_time",
    "berendsen": "coupling_time",
    "nose_hoover": "period_time",
}

# MD parameter -> pressure_coupling/algorithm map
BAROSTAT = {
    "No": None,
    "Berendsen": "berendsen",
    "C-rescale": "stochastic_cell_rescaling",
    "Parrinello-Rahman": "parrinello_rahman",
    "MTTK": "mttk",  # not a predefined value of the schema
}

# MD parameter -> pressure_coupling/algorithm time field map
BAROSTAT_TIME_FIELD = {
    "stochastic_cell_rescaling": "coupling_time",
    "berendsen": "coupling_time",
    "parrinello_rahman": "period_time",
}

# MD parameter -> coulomb_interactions/type map
COULOMB_TYPE = {
    "Cut-off": "cutoff",
    "PME": "pme",
    "P3M-AD": "pppm",
    "Ewald": "ewald",
    "Reaction-Field": "reaction_field",
    "Reaction-Field-zero": "reaction_field",
    "FMM": "fmm",
}

# MD parameter -> vdw_interactions/type map
VDW_TYPE = {
    "Cut-off": "cutoff",
    "PME": "pme",
}

# MD parameter -> cutoff/potential_modifier map
POTENTIAL_MODIFIER = {
    "None": "none",
    "Potential-shift": "potential_shift",
    "Potential-shift-Verlet": "potential_shift",
    "Potential-switch": "potential_switch",
    "Force-switch": "force_switch",
}

# MD parameter -> vdw_interactions/dispersion_correction map
DISPERSION_CORRECTION = {
    "No": None,
    "EnerPres": "energy_and_pressure",
    "Ener": "energy",
    "AllEnerPres": "energy_and_pressure",
    "AllEner": "energy",
}

# MD parameter -> center_of_mass_motion_removal/mode map
COM_MODE = {
    "None": None,
    "Linear": "translation",
    "Angular": "translation_and_rotation",
    "Linear-acceleration-correction": "translation",
}

# MD parameter -> constraints/algorithm map
CONSTRAINT_ALGORITHM = {
    "Lincs": "lincs",
    "Shake": "shake",
}

# Interaction type map from GROMACS functype name to schema interaction type
INTERACTION = {
    "BONDS": "bond",
    "HARMONIC": "bond",
    "MORSE": "morse_bond",
    "ANGLES": "angle",
    "UREY_BRADLEY": "urey_bradley",
    "PDIHS": "proper_dihedral",
    "IDIHS": "improper_dihedral",
    "RBDIHS": "ryckaert_bellemans_dihedral",
    "LJ14": "lennard_jones_14",
    "CONSTR": "constraint",
    "CONSTRNC": "constraint",
    "SETTLE": "settle",
    "VSITE3": "virtual_site_3",
    "VSITE3FD": "virtual_site_3",
    "VSITE3FAD": "virtual_site_3",
    "VSITE3OUT": "virtual_site_3",
    "LJ_SR": "lennard_jones_short_range",
    "BHAM": "buckingham",
    "PIDIHS": "periodic_improper_dihedral",
    "G96BONDS": "fourth_power_bond",
    "G96ANGLES": "cosine_harmonic_angle",
    "CMAP": "cmap",
}

# Parameter field map from GROMACS functype field name to schema compound field name
PARAMETER_FIELDS = {
    "BONDS": {"b0A": "equilibrium_distance", "cbA": "force_constant"},
    "G96BONDS": {"b0A": "equilibrium_squared_distance", "cbA": "force_constant"},
    "HARMONIC": {"b0A": "equilibrium_distance", "cbA": "force_constant"},
    "MORSE": {"b0A": "equilibrium_distance", "cbA": "well_depth", "betaA": "steepness"},
    "ANGLES": {"thA": "equilibrium_angle", "ctA": "force_constant"},
    "G96ANGLES": {"thA": "equilibrium_cosine", "ctA": "force_constant"},
    "UREY_BRADLEY": {
        "thetaA": "equilibrium_angle",
        "kthetaA": "angle_force_constant",
        "r13A": "equilibrium_distance",
        "kUBA": "bond_force_constant",
    },
    "PDIHS": {"phiA": "phase", "cpA": "force_constant", "mult": "multiplicity"},
    "PIDIHS": {"phiA": "phase", "cpA": "force_constant", "mult": "multiplicity"},
    "IDIHS": {"xiA": "equilibrium_angle", "cxA": "force_constant"},
    "RBDIHS": {f"rbcA[{i}]": f"c{i}" for i in range(6)},
    "LJ14": {"c6A": "c6", "c12A": "c12"},
    "CONSTR": {"dA": "distance"},
    "CONSTRNC": {"dA": "distance"},
    "SETTLE": {"doh": "distance_OH", "dhh": "distance_HH"},
    "VSITE3": {"a": "a", "b": "b"},
    "LJ_SR": {"c6": "c6", "c12": "c12"},
    "BHAM": {"a": "constant_a", "b": "constant_b", "c": "constant_c"},
}

# GROMACS functypes whose entries are chemical bonds, for ``connectivity/bonds``.
#
# Taken from GROMACS itself: these are exactly the functypes declared with the
# IF_CHEMBOND flag in src/gromacs/topology/ifunc.cpp, which is that engine's own
# statement of what constitutes a chemical connection. Deriving the set from the
# flag rather than from the ones a test system happens to use keeps a topology
# using an unusual bond type from silently losing its connectivity.
#
# Excluded deliberately, despite the names:
#   HARMONIC      a generic harmonic potential between two atoms (def_bonded);
#                 GROMACS does not flag it as a chemical bond
#   CONSTRNC      "Constr. No Co." -- a constraint that explicitly does NOT
#                 represent a connection (def_shk, no IF_CHEMBOND)
#   TABBONDSNC    the no-connection variant of TABBONDS
#   RESTRAINTPOT  a restraint, not a bond
#
# NOTE: SETTLE carries IF_CHEMBOND too, but its entry is a rigid water of three
# atoms that expands to two O-H bonds, so roundtrip.py handles it separately.
BOND_LIKE = {
    "BONDS",
    "G96BONDS",
    "MORSE",
    "CUBICBONDS",
    "CONNBONDS",  # a connection carrying no interaction
    "FENEBONDS",
    "TABBONDS",
    "POLARIZATION",
    "ANHARM_POL",  # shell particle bonded to its parent atom
    "CONSTR",  # a bond converted to a constraint keeps its chemistry
}


# GROMACS .edr energy term -> observables element (S6).
# Only terms the schema names are mapped; the rest of an .edr (per-interaction
# energy decomposition, Constr. rmsd, lambda bookkeeping) has no schema slot and
# is deliberately not compared.
OBSERVABLE = {
    "Potential": "potential_energy",
    "Kinetic En.": "kinetic_energy",
    "Total Energy": "total_energy",
    "Conserved En.": "conserved_energy",
    "Temperature": "temperature",
    "Pressure": "pressure",
    "Volume": "volume",
    "Density": "density",
    "Enthalpy": "enthalpy",
}

#: Tensor observables, assembled from the nine per-component .edr terms.
OBSERVABLE_TENSOR = {
    "virial_tensor": "Vir-{}{}",
    "pressure_tensor": "Pres-{}{}",
}
TENSOR_AXES = ("X", "Y", "Z")
