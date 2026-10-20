"""
The validation test matrix: chemical systems x simulation protocols.

Validation 1 asks whether the schema can represent the metadata encountered in
realistic biomolecular MD.  That has two independent axes, and this module keeps
them separate:

* CHEMICAL SYSTEMS -- what is in the box (protein, nucleic acid, membrane,
  ligand, ion, solvent), each maps to a topology on validation systems.
* PROTOCOLS -- how the simulation is run (integrator, thermostat, barostat,
  constraints, electrostatics, boundary conditions), expressed as .mdp overrides.

Every protocol exercises a different branch of the ``parameters`` tree of S7, so
the product of the two axes covers the schema rather than merely repeating one
configuration on many molecules.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field

#: The extracted validation dataset: the systems the harness runs on. Download
#: it, unpack it, and point this at the directory that results.
#:
#:     export H5MD_VALIDATION_DATA=/scratch/$USER/validation_systems
#:
#: The default is the copy that sits beside this repository, so a checkout runs
#: without any environment set.
VALIDATION_SYSTEMS = os.environ.get(
    "H5MD_VALIDATION_DATA",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "validation_systems",
    ),
)

#: Where a run writes. Both output trees default to the dataset directory; when
#: that cannot be written to (a read-only copy), to ``h5md_validation/`` under
#: the system temporary directory instead. DATA_DIR and RESULTS_DIR each
#: override one of them. Every stage imports these names, so a redirected run
#: cannot end up half in one tree and half in another.
OUTPUT_ROOT = (
    VALIDATION_SYSTEMS
    if os.access(VALIDATION_SYSTEMS, os.W_OK)
    else os.path.join(tempfile.gettempdir(), "h5md_validation")
)

#: Everything the harness produces on the way to a result. Not part of the
#: dataset, which the leading underscore is there to say. Several GB per full
#: run; safe to delete.
#:
#:     _data/tprs/             grompp output, the round-trip source
#:     _data/resolved_mdp/     grompp -po: the .mdp with every default filled in
#:     _data/reference_files/  the .h5md corpus mdrun writes
#:     _data/tpr_dumps/        cached `gmx dump -s` text, so it runs once per tpr
#:     _data/scratch/          per-pair scratch, deleted after each pair; a
#:                             leftover directory names a run that was interrupted
DATA_ROOT = os.environ.get("DATA_DIR", os.path.join(OUTPUT_ROOT, "_data"))

#: What a run concludes: the JSON measurements and the figures. Small, and the
#: thing worth keeping when ``_data`` is deleted.
RESULTS_DIR = os.environ.get("RESULTS_DIR", os.path.join(OUTPUT_ROOT, "_results"))


def output_locations() -> str:
    """Where this run reads systems from and writes to, for printing."""
    lines = [
        f"  systems : {VALIDATION_SYSTEMS}",
        f"  data    : {DATA_ROOT}",
        f"  results : {RESULTS_DIR}",
    ]
    if OUTPUT_ROOT != VALIDATION_SYSTEMS:
        lines.append(
            f"  (the dataset directory is not writable; outputs default to {OUTPUT_ROOT})"
        )
    return "\n".join(lines)


#: The GROMACS binary that grompps the systems and dumps the .tpr for the
#: round-trip. GMX_BIN points at a build that is not on PATH.
GMX = os.environ.get("GMX_BIN", "gmx")

#: Bases for the negative tests (Validation 4), chosen for different shapes so
#: that as many corruptions as possible apply: a protein with a non-contiguous
#: selection, a nucleic acid, a membrane protein, an ion/solvent box and a
#: protein-DNA complex. A base missing from the corpus is skipped, not an error.
NEGATIVE_TEST_BASES = (
    "bena_12k__selection-noncontiguous",
    "dna_amber__nvt-vrescale",
    "memprotmd_5bkf__nvt-vrescale",
    "solvent_ions__nvt-vrescale",
    "protein_dna__nvt-vrescale",
)

#: Pairs for the cross-engine comparison (Validation 6): (reference file key,
#: AMBER topology relative to VALIDATION_SYSTEMS). Both sides MUST describe the
#: same molecule -- a mismatched pair compares two different systems and its
#: verdict means nothing. Every hexABC and MISATO system ships the AMBER
#: topology its GROMACS topology was converted from (.prmtop, or the deposited
#: production.top.gz, read as is), so each is paired with its own source.
#: Declaring the pair here is what makes that invariant checkable instead of a
#: matter of passing the right flags.
AMBER_INTEROP_PAIRS = (
    (
        "hexabc_a0017b_dna__nvt-vrescale",
        "hexabc/hexabc_a0017b_dna/A017B_topology.prmtop",
    ),
    (
        "hexabc_a0017p_rna__nvt-vrescale",
        "hexabc/hexabc_a0017p_rna/A017P_topology.prmtop",
    ),
    (
        "hexabc_a0017u_rna__nvt-vrescale",
        "hexabc/hexabc_a0017u_rna/A017U_topology.prmtop",
    ),
    (
        "hexabc_a00186_dna__nvt-vrescale",
        "hexabc/hexabc_a00186_dna/A0186_topology.prmtop",
    ),
    (
        "hexabc_a001av_dna_protein__nvt-vrescale",
        "hexabc/hexabc_a001av_dna_protein/A01AV_topology.prmtop",
    ),
    (
        "hexabc_a001fc_dna_protein__nvt-vrescale",
        "hexabc/hexabc_a001fc_dna_protein/A01FC_topology.prmtop",
    ),
    (
        "hexabc_a001tx_rna_protein__nvt-vrescale",
        "hexabc/hexabc_a001tx_rna_protein/MCV1900433_topology.prmtop",
    ),
    (
        "hexabc_a001u2_rna_protein__nvt-vrescale",
        "hexabc/hexabc_a001u2_rna_protein/MCV1900438_topology.prmtop",
    ),
    (
        "misato_1fl6__nvt-vrescale",
        "misato/misato_1fl6/production.top.gz",
    ),
    (
        "misato_1tfq__nvt-vrescale",
        "misato/misato_1tfq/production.top.gz",
    ),
    (
        "misato_1z3j__nvt-vrescale",
        "misato/misato_1z3j/production.top.gz",
    ),
    (
        "misato_2xde__nvt-vrescale",
        "misato/misato_2xde/production.top.gz",
    ),
    (
        "misato_2zyn__nvt-vrescale",
        "misato/misato_2zyn/production.top.gz",
    ),
    (
        "misato_3f78__nvt-vrescale",
        "misato/misato_3f78/production.top.gz",
    ),
    (
        "misato_3lq5__nvt-vrescale",
        "misato/misato_3lq5/production.top.gz",
    ),
    (
        "misato_3ptg__nvt-vrescale",
        "misato/misato_3ptg/production.top.gz",
    ),
    (
        "misato_3s2v__nvt-vrescale",
        "misato/misato_3s2v/production.top.gz",
    ),
    (
        "misato_3unz__nvt-vrescale",
        "misato/misato_3unz/production.top.gz",
    ),
    (
        "misato_4d0w__nvt-vrescale",
        "misato/misato_4d0w/production.top.gz",
    ),
    (
        "misato_4q1n__nvt-vrescale",
        "misato/misato_4q1n/production.top.gz",
    ),
    (
        "misato_4xua__nvt-vrescale",
        "misato/misato_4xua/production.top.gz",
    ),
    (
        "misato_4zw5__nvt-vrescale",
        "misato/misato_4zw5/production.top.gz",
    ),
    (
        "misato_5ivt__nvt-vrescale",
        "misato/misato_5ivt/production.top.gz",
    ),
    (
        "misato_5nhj__nvt-vrescale",
        "misato/misato_5nhj/production.top.gz",
    ),
    (
        "misato_6dyw__nvt-vrescale",
        "misato/misato_6dyw/production.top.gz",
    ),
    (
        "misato_6i77__nvt-vrescale",
        "misato/misato_6i77/production.top.gz",
    ),
    (
        "misato_6mim__nvt-vrescale",
        "misato/misato_6mim/production.top.gz",
    ),
    (
        "misato_6nwc__nvt-vrescale",
        "misato/misato_6nwc/production.top.gz",
    ),
)


@dataclass
class System:
    """A chemical system, defined by whatever is needed to grompp it."""

    id: str
    description: str
    #: dimensions of Validation 1's coverage table that this system exercises
    dimensions: tuple
    #: the system's directory in the validation dataset
    location: str = ""
    structure: str = ""
    topology: str = ""
    index: str = ""
    force_field: str = ""
    water_model: str = ""
    #: .mdp settings this system requires regardless of protocol (coupling
    #: groups, mostly), applied on top of the protocol's settings
    mdp: dict = field(default_factory=dict)
    #: index groups this system needs passed to grompp
    use_index: bool = False
    #: False for a vacuum system: no periodic boundary, infinite cut-offs
    periodic: bool = True
    #: protocol ids this system is run under; None means every protocol
    protocols: tuple = None


@dataclass
class Protocol:
    """A set of .mdp overrides exercising one branch of the parameters tree."""

    id: str
    description: str
    #: schema elements this protocol is meant to populate
    exercises: tuple
    mdp: dict = field(default_factory=dict)


# ==========================================================================
# Protocols -- each one lights up a different part of S7
# ==========================================================================
_BASE_MD = {
    "integrator": "md",
    "dt": "0.002",
    "nsteps": "500",
    "nstxout": "100",
    "nstvout": "0",
    "nstfout": "0",
    "nstenergy": "100",
    "nstlog": "100",
    "nstlist": "20",
    "cutoff-scheme": "Verlet",
    "coulombtype": "PME",
    "rcoulomb": "1.2",
    "vdwtype": "Cut-off",
    "vdw-modifier": "Potential-shift",
    "rvdw": "1.2",
    "pbc": "xyz",
    "constraints": "h-bonds",
    "constraint-algorithm": "lincs",
    "lincs-order": "4",
    "lincs-iter": "1",
    "tcoupl": "V-rescale",
    "tc-grps": "System",
    "tau-t": "0.1",
    "ref-t": "300",
    "pcoupl": "no",
    "comm-mode": "Linear",
    "nstcomm": "100",
    "gen-vel": "no",
}


def _md(**overrides):
    out = dict(_BASE_MD)
    out.update({k.replace("_", "-"): str(v) for k, v in overrides.items()})
    return out


PROTOCOLS = [
    Protocol(
        "nvt-vrescale",
        "NVT, leap-frog, stochastic velocity rescaling",
        (
            "run_control/molecular_dynamics",
            "temperature_coupling/stochastic_velocity_rescaling",
        ),
        _md(),
    ),
    Protocol(
        "npt-prahman",
        "NPT, Parrinello-Rahman barostat",
        ("pressure_coupling/parrinello_rahman",),
        _md(
            pcoupl="Parrinello-Rahman",
            tau_p="2.0",
            ref_p="1.0",
            compressibility="4.5e-5",
            pcoupltype="isotropic",
            refcoord_scaling="com",
        ),
    ),
    Protocol(
        "npt-crescale",
        "NPT, stochastic cell rescaling barostat",
        ("pressure_coupling/stochastic_cell_rescaling",),
        _md(
            pcoupl="C-rescale",
            tau_p="2.0",
            ref_p="1.0",
            compressibility="4.5e-5",
            pcoupltype="isotropic",
        ),
    ),
    Protocol(
        "npt-berendsen",
        "NPT, Berendsen barostat",
        ("pressure_coupling/berendsen",),
        _md(
            pcoupl="Berendsen",
            tau_p="2.0",
            ref_p="1.0",
            compressibility="4.5e-5",
            pcoupltype="isotropic",
        ),
    ),
    Protocol(
        "npt-semiisotropic",
        "NPT, semi-isotropic coupling (anisotropic tensor)",
        ("pressure_coupling/reference_pressure as tensor",),
        _md(
            pcoupl="Parrinello-Rahman",
            pcoupltype="semiisotropic",
            tau_p="5.0",
            ref_p="1.0 1.0",
            compressibility="4.5e-5 4.5e-5",
        ),
    ),
    Protocol(
        "nvt-nosehoover",
        "NVT, Nose-Hoover thermostat with velocity Verlet",
        ("temperature_coupling/nose_hoover",),
        _md(
            integrator="md-vv", tcoupl="nose-hoover", tau_t="1.0", nh_chain_length="10"
        ),
    ),
    Protocol(
        "nvt-andersen",
        "NVT, Andersen (massive) thermostat with velocity Verlet",
        ("temperature_coupling/andersen_massive",),
        _md(
            integrator="md-vv",
            tcoupl="andersen-massive",
            tau_t="0.1",
            constraints="none",
            dt="0.0005",
            nsttcouple="10",
        ),
    ),
    Protocol(
        "nvt-berendsen",
        "NVT, Berendsen thermostat",
        ("temperature_coupling/berendsen",),
        _md(tcoupl="berendsen", tau_t="0.5"),
    ),
    Protocol(
        "nve",
        "NVE, no temperature or pressure coupling",
        ("run_control without coupling groups",),
        _md(tcoupl="no", pcoupl="no"),
    ),
    Protocol(
        "sd-langevin",
        "Stochastic dynamics (Langevin) integrator",
        ("run_control/stochastic_dynamics",),
        _md(integrator="sd", tau_t="1.0", ld_seed="1993"),
    ),
    Protocol(
        "em-steep",
        "Energy minimisation, steepest descent",
        ("run_control/energy_minimization",),
        {
            "integrator": "steep",
            "nsteps": "50",
            "emtol": "1000.0",
            "emstep": "0.01",
            "nstenergy": "10",
            "nstlog": "10",
            "nstxout": "10",
            "cutoff-scheme": "Verlet",
            "nstlist": "20",
            "coulombtype": "PME",
            "rcoulomb": "1.2",
            "vdwtype": "Cut-off",
            "rvdw": "1.2",
            "pbc": "xyz",
            "constraints": "none",
        },
    ),
    Protocol(
        "em-cg",
        "Energy minimisation, conjugate gradient",
        ("run_control/energy_minimization",),
        {
            "integrator": "cg",
            "nsteps": "50",
            "emtol": "1000.0",
            "emstep": "0.01",
            "nstcgsteep": "100",
            "nstenergy": "10",
            "nstlog": "10",
            "nstxout": "10",
            "cutoff-scheme": "Verlet",
            "nstlist": "20",
            "coulombtype": "PME",
            "rcoulomb": "1.2",
            "vdwtype": "Cut-off",
            "rvdw": "1.2",
            "pbc": "xyz",
            "constraints": "none",
        },
    ),
    Protocol(
        "constraints-none",
        "No bond constraints",
        ("constraints absent; bonds retained as interactions",),
        _md(constraints="none", dt="0.0005"),
    ),
    Protocol(
        "constraints-allbonds",
        "All bonds constrained",
        ("constraints/mode = all_bonds",),
        _md(constraints="all-bonds"),
    ),
    Protocol(
        "constraints-shake",
        "SHAKE constraint algorithm",
        ("constraints/shake",),
        _md(
            constraint_algorithm="shake",
            shake_tol="1e-5",
            integrator="md-vv",
            tcoupl="no",
        ),
    ),
    Protocol(
        "coulomb-rf",
        "Reaction-field electrostatics",
        ("coulomb_interactions/reaction_field",),
        _md(coulombtype="Reaction-Field", epsilon_rf="78", rcoulomb="1.2"),
    ),
    Protocol(
        "coulomb-cutoff",
        "Plain cut-off electrostatics",
        ("coulomb_interactions/cutoff",),
        _md(coulombtype="Cut-off", coulomb_modifier="Potential-shift"),
    ),
    Protocol(
        "vdw-forceswitch",
        "Van der Waals force switch",
        ("vdw_interactions/cutoff with force_switch",),
        _md(vdw_modifier="Force-switch", rvdw_switch="1.0"),
    ),
    Protocol(
        "vdw-potswitch",
        "Van der Waals potential switch",
        ("vdw_interactions/cutoff with potential_switch",),
        _md(vdw_modifier="Potential-switch", rvdw_switch="1.0"),
    ),
    Protocol(
        "dispcorr",
        "Analytical dispersion correction",
        ("vdw_interactions/dispersion_correction",),
        _md(DispCorr="EnerPres"),
    ),
    Protocol(
        "selection-noncontiguous",
        "Non-contiguous energy group, exercising index-based selections",
        ("selections/<group>/indices", "selections/<group>/ranges"),
        _md(energygrps="C-alpha"),
    ),
    Protocol(
        "tcoupl-multigroup",
        "Two temperature coupling groups",
        ("temperature_coupling with per-group arrays",),
        _md(tc_grps="Protein Non-Protein", tau_t="0.1 0.1", ref_t="300 300"),
    ),
    Protocol(
        "mass-repartition",
        "Hydrogen mass repartitioning",
        ("run_control/mass_repartitioning",),
        _md(mass_repartition_factor="3", dt="0.004"),
    ),
]

PROTOCOLS_BY_ID = {p.id: p for p in PROTOCOLS}
