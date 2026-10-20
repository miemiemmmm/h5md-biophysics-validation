"""
Registry of specification-vs-implementation mismatches.

``SI-Specification.pdf`` is normative here. If the reference files write something
that is out of the specification, the validator raises a DEVIATION finding with one
of the codes below, so the mismatch waits for an editorial decision instead of
being reported as a file error.

Each entry cites the specification section concerned rather than quoting its
text, so the registry does not go stale as the text is revised. ``status`` is
OPEN until a decision is taken; ``resolution`` records what it was.

Used by:
    Validation 3: internal consistency
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Deviation:
    code: str
    title: str
    #: Section(s) of the specification concerned, e.g. "S8".
    section: str
    #: What the reference implementation does, and where.
    implementation: str
    #: Why it matters for a reader or a third-party implementer.
    impact: str
    #: Suggested way to close it.  Not a decision.
    suggestion: str
    status: str = "OPEN"
    resolution: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


DEVIATIONS = {
    d.code: d
    for d in (
        Deviation(
            code="DEV-FF-BONDED-LAYER",
            title="force_field has an undocumented bonded/nonbonded layer",
            section="S8",
            implementation=(
                "Files are written as force_field/bonded/<type>/parameters and "
                "force_field/nonbonded/<type>/parameters."
            ),
            impact=(
                "A third-party reader written from the specification alone would "
                "look for force_field/<type>/parameters and find nothing."
            ),
            suggestion=(
                "Document the bonded/ and nonbonded/ subgroups in S8 and replace "
                "the asterisk convention with the explicit path, or drop the layer "
                "from the writer."
            ),
        ),
        Deviation(
            code="DEV-OBS-GROUP-NESTING",
            title="observables may nest per-selection-group subgroups",
            section="S6",
            implementation=(
                "Per-coupling-group observables are written one level deeper, as "
                "observables/<group_name>/<observable> (e.g. observables/solvent/temperature)."
            ),
            impact=(
                "A reader enumerating observables/ gets a mixture of observables and "
                "group names with no documented way to tell them apart."
            ),
            suggestion=(
                "Document the nesting in S6 and state the rule that distinguishes a "
                "group subgroup from an observable (e.g. an observable group always "
                "contains a 'value' element)."
            ),
        ),
        Deviation(
            code="DEV-CONSTRAINTS-MODE-UNWRITTEN",
            title="constraints/mode is specified but never written",
            section="S7.8",
            implementation=(
                "No run writes the element. GROMACS applies the mdp 'constraints' "
                "setting during preprocessing and does not retain it in the run "
                "input, so the writer has nothing to copy. Verified across three "
                "runs of the same system with constraints = none, h-bonds and "
                "all-bonds: all three produce an identical parameters/constraints "
                "group holding only algorithm and the lincs subgroup."
            ),
            impact=(
                "The information is not lost -- it is recoverable by comparing the "
                "bond and constraint interaction counts in the topology, which differ "
                "exactly as expected (604/0, 321/283, 0/604 for the three runs). But "
                "it is not where the specification says it is, and recovering it "
                "requires an inference the specification never describes. A reader "
                "consulting only parameters/constraints, as the specification invites, "
                "cannot tell the three runs apart.\n"
                "The underlying engine difference is documented in the literature: "
                "Abraham et al., 'Sharing Data from Molecular Simulations', J. Chem. "
                "Inf. Model. 59(10):4093-4099 (2019), doi:10.1021/acs.jcim.9b00665, "
                "note that 'different programs store and represent molecular "
                "information in different ways (e.g. SHAKE information for "
                "constraining bonds is represented in the molecular topology in "
                "GROMACS, while it is a simulation command parameter in NAMD and "
                "AMBER)'. So this is a cross-engine semantic difference the schema "
                "has to take a position on, not merely a GROMACS writer gap."
            ),
            suggestion=(
                "Either derive mode in the writer from the topology (the counts make "
                "it unambiguous), or mark the element as not recoverable from a "
                "GROMACS run input and say in S7.8 how to infer it from the "
                "interaction lists."
            ),
        ),
        Deviation(
            code="DEV-CONSTRAINTS-ALGORITHM-UNCONDITIONAL",
            title="constraints/algorithm is written even when nothing is constrained",
            section="S7.8",
            implementation=(
                "The element is written unconditionally from the GROMACS "
                "constraint-algorithm setting, which carries its default value whether "
                "or not any bond is constrained. A run with constraints = none, which "
                "converts zero bonds to constraints, still reports algorithm = 'lincs'."
            ),
            impact=(
                "The element reads as an assertion that LINCS constrained this "
                "simulation, which is false for such a run: the only constrained "
                "degrees of freedom left are rigid waters handled by SETTLE, which is "
                "a different algorithm and is recorded elsewhere. Of the schema's "
                "claims this is the one that can mislead rather than merely omit."
            ),
            suggestion=(
                "Write the constraints group only when the topology contains "
                "constraint interactions, or state in S7.8 that algorithm records the "
                "configured solver rather than evidence that constraints were applied."
            ),
        ),
        Deviation(
            code="DEV-REFCOORD-NO-VOCABULARY",
            title="reference_coordinate_scaling has no defined value set",
            section="S7.4",
            implementation=(
                "The writer emits 'none', 'all' and 'center_of_mass', corresponding to "
                "the GROMACS refcoord-scaling settings No, All and COM."
            ),
            impact=(
                "A reader cannot know what strings to expect, and a second "
                "implementation has no way to choose the same spellings. This is the "
                "one place in the parameters tree where the schema is not "
                "machine-checkable, and it was found only because an independent "
                "mapping written from the specification guessed 'no' where the "
                "implementation writes 'none'."
            ),
            suggestion=(
                "Add the defined values to S7.4 in the same style as the neighbouring "
                "elements: none, all, center_of_mass."
            ),
        ),
        Deviation(
            code="DEV-NO-ATOM-NAME-CONVENTION",
            title="atom and component names have no canonical convention",
            section="S10",
            implementation=(
                "Names are copied from whatever the preparation tool produced. The "
                "same B-DNA dodecamer prepared in GROMACS and in AMBER gives 48 of "
                "758 atoms different names for the same atom -- O1P/OP1, O2P/OP2, "
                "H3T/HO3', H5T/HO5' -- and terminal residues named DC5/DC3 in one and "
                "DC in the other."
            ),
            impact=(
                "Two files describing the same molecule cannot be compared or merged "
                "without an external alias table. Element identity, masses, charges "
                "and the covalent graph all agree once the atoms are aligned, so the "
                "obstacle is purely nomenclature -- but it is an obstacle to exactly "
                "the reuse the schema exists to enable."
            ),
            suggestion=(
                "Name a convention in S10 (the PDB Chemical Component Dictionary is "
                "the usual choice), or add an optional attribute recording which "
                "convention the names follow, so a reader can translate rather than "
                "guess."
            ),
        ),
        Deviation(
            code="DEV-FF-PROVENANCE-UNWRITTEN",
            title="no file records which force field it used",
            section="S8",
            implementation=(
                "None of the thirteen is written in any of the 158 reference files, "
                "which between them use CHARMM36, CHARMM27, AMBER99SB-ILDN, OPLS-AA "
                "and GROMOS54A7 with TIP3P, TIP4P and SPC. force_field carries no "
                "attributes at all. The name lives in the #include lines of the "
                ".top, which grompp resolves and then discards."
            ),
            impact=(
                "They are optional, so their absence is not a conformance error -- "
                "but they answer the first question anyone reusing a trajectory "
                "asks. A file carries every force-field parameter and does not say "
                "which force field they came from.\n"
                "There is a precedent, and it is not encouraging. TNG, the other "
                "structured trajectory format GROMACS writes, has a dedicated "
                "force-field name in the general info block of every file, with "
                "tng_forcefield_name_set()/_get() exposed in C, C++ and Fortran "
                "and covered by the TNG test suite. GROMACS has never called it: "
                "the only occurrences of the setter anywhere in the source tree "
                "are the declaration, the definition, the two language bindings "
                "and the library's own test. Because it is never set, the field is "
                "malloc'd to the empty string and written empty "
                "(tng_io.c:2091, written at tng_io.c:2361), so every TNG file "
                "GROMACS has ever produced carries a present, empty force-field "
                "name. Confirmed here by writing one: the general info block holds "
                "'gmx mdrun 2027.0-dev-...', the user name and the host name, all "
                "of which GROMACS does set, and then an empty force-field "
                "string.\n"
                "That is the sentinel design already run as an experiment for over "
                "a decade, in a format with the slot built in. A declared-but-"
                "uncaptured field stays empty, for the same reason in both "
                "formats: grompp resolves the name and discards it, so no writer "
                "downstream has anything to put there.\n"
                "Two measurements bound the design space "
                "(scripts/provenance_recoverability.py):\n"
                "  Inferring it from the parameters does not work. Across 78 system "
                "pairs, the lowest bonded-parameter overlap between two systems "
                "sharing a force field is 23.8% (AMBER protein against AMBER RNA) "
                "while the highest between two that do not is 90.3% (CHARMM36 "
                "against CHARMM27; AMBER against OPLS reaches 85%). Overlap tracks "
                "what is in the box, not what parametrised it, so the signature is "
                "not separable and no fingerprint recovers the name.\n"
                "  Reading it from the topology does work, for most files. The "
                "force field and the water model are both named in an #include "
                "path for 11 of 15 systems -- every system prepared through "
                "pdb2gmx, which is the reference implementation's own documented "
                "route. The four exceptions were assembled externally "
                "(CHARMM-GUI), where the includes name local files only."
            ),
            suggestion=(
                "Three changes, in decreasing order of return. The TNG precedent "
                "argues for doing the first one rather than only the third: a slot "
                "without a capture path is what TNG already has.\n"
                "  1. Have grompp keep the resolved #include paths and carry them "
                "into the .tpr, and the writer emit them as source and "
                "polymer_force_field / solvent_model. This is free for the user "
                "and exact, and it covers the 11 of 15 systems above.\n"
                "  2. Add .mdp fields for the same attributes, for systems "
                "assembled elsewhere and for anyone who wants to override or "
                "enrich what the paths say. This is the manual annotation route, "
                "and it should be the fallback rather than the only mechanism -- "
                "as an optional attribute with no capture path it produced 0 of "
                "158.\n"
                "  3. Distinguish 'not declared' from 'absent'. An omitted "
                "attribute is silent: a reader cannot tell an undeclared force "
                "field from one the writer had no element for. Requiring source to "
                "be present, with a defined value meaning undeclared, makes the "
                "gap visible in the file instead of only in a validator that knows "
                "to look."
            ),
        ),
        Deviation(
            code="DEV-FF-NO-EXTENSION-CONVENTION",
            title="S8 gives no convention for interaction types outside the table",
            section="S7, S8",
            implementation=(
                "The writer now emits a dedicated group for every term outside the "
                "table, so nothing is silently lost (defects.IMP-SILENT-INTERACTION-"
                "LOSS, FIXED). What it cannot do is name the convention it follows: "
                "the groups carry no doi or url, because S8 asks for none."
            ),
            impact=(
                "The set of bonded functional forms in use across MD engines is open, "
                "and enumerating all of them was never the intent -- so the "
                "convention is what carries the extensibility claim for force-field "
                "terms, exactly as S7's convention carries it for simulation "
                "parameters. Right now that half of the claim rests on a convention "
                "the document does not state, which is why the reference "
                "implementation does not follow one."
            ),
            suggestion=(
                "Add to S8 the clause S7 already has, adapted to interaction types: "
                "a term outside Table SI is written as a group under force_field "
                "named after the term, holding its own parameters with a units "
                "attribute and a doi or url for the functional form, with a matching "
                "interactions group under the molecule type. Give cmap a row of its "
                "own stating the grid layout, since it is already named in the "
                "preamble and is the only type the convention currently has to carry."
            ),
            resolution=(
                "Partly taken. The revision of 2026-10-01 gives fourth_power_bond, "
                "cosine_harmonic_angle and periodic_improper_dihedral rows of their "
                "own in Table SI, so three of the four terms this validation found "
                "are now predefined rather than extensions. The clause itself is "
                "still absent and cmap still has no row, so the deviation stands for "
                "those two points."
            ),
        ),
        Deviation(
            code="DEV-FF-NO-FUNCTIONAL-FORM",
            title="interaction types do not state their functional form",
            section="S8",
            implementation=(
                "The names imply the harmonic forms, and the writer treats them that "
                "way, but nothing in the document says so."
            ),
            impact=(
                "Table SI now names the pairs that used to be ambiguous -- bond "
                "against fourth_power_bond, angle against cosine_harmonic_angle, "
                "improper_dihedral against periodic_improper_dihedral -- so the "
                "common cases are decidable from the names. What is still not "
                "stated is the potential each one denotes, which is what an "
                "implementer needs for a form that is in neither column. A quartic "
                "bond and a harmonic bond both carry an equilibrium distance and a "
                "force constant and are not interconvertible; the distinction is "
                "now drawn by the row names alone, and a reader has to infer from "
                "those names which potential is meant."
            ),
            suggestion=(
                "Give Table SI a column stating each potential explicitly. It is "
                "the smallest change that makes the boundary of the predefined set "
                "decidable, and it is a precondition for the convention in "
                "DEV-FF-NO-EXTENSION-CONVENTION to be applied consistently."
            ),
        ),
        Deviation(
            code="DEV-SPEC-GAP-ARITY",
            title="interaction particle arity is not specified",
            section="S8, S10",
            implementation=(
                "Arity follows the physics of each term (angle 3, dihedrals 4, "
                "settle 3, ...) and is implicit in the writer."
            ),
            impact=(
                "A validator or reader cannot check or allocate particle_ids from the "
                "specification alone. The validator supplies the arities itself; they "
                "are an assumption, not a specified value."
            ),
            suggestion="Add an N_particles column to Table SI.",
        ),
    )
}


def get(code: str) -> Deviation:
    return DEVIATIONS[code]


def open_deviations():
    return [d for d in DEVIATIONS.values() if d.status == "OPEN"]
