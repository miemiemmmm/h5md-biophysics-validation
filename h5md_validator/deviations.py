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
    )
}


def get(code: str) -> Deviation:
    return DEVIATIONS[code]


def open_deviations():
    return [d for d in DEVIATIONS.values() if d.status == "OPEN"]
