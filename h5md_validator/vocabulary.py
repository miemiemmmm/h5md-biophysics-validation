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

# MD parameter -> <simulation_type>/integrator map
INTEGRATOR = {
    "md": "leap_frog",
    "md-vv": "velocity_verlet",
    "md-vv-avek": "velocity_verlet",
    "sd": "langevin",
    "bd": "brownian",
    "steep": "steepest_descent",
    "cg": "conjugate_gradient",
    "l-bfgs": "l_bfgs",
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
    "MTTK": "mttk",          # not a predefined value of the schema
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

# TODO: MD parameter -> pressure_coupling/reference_coordinate_scaling map
#:
#: S7.4 describes this element but, alone among the categorical parameters of
#: S7, defines no set of allowed values. The spellings on the right are the ones
#: the reference implementation writes, recovered by inspection rather than from
#: the specification: see deviations.DEV-REFCOORD-NO-VOCABULARY.
REFCOORD_SCALING = {
    "No": "none",
    "All": "all",
    "COM": "center_of_mass",
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
# NOTE: SETTLE describes a rigid water and expands to its two O-H bonds.
BOND_LIKE = {"BONDS", "G96BONDS", "HARMONIC", "MORSE", "CONSTR", "CONSTRNC"}
