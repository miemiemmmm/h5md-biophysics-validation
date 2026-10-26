#!/usr/bin/env python
"""
Validation 5 figure: what the metadata costs to store and to write.

    validation5_cost.png  A  metadata stored against system size
                          B  metadata bytes per unique atom (connectivity excluded)
                          C  each phase's share of the write time, from the timing
                             data in data/decompose_run*.json (see data/README.md)

Written to <results>/figures/.

Used by:
    Validation 5: storage overhead
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics as st
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from config.matrix import RESULTS_DIR


def _r(name: str) -> str:
    return str(pathlib.Path(RESULTS_DIR) / name)


OUT = Path(RESULTS_DIR) / "figures"

COLORS = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#5DD3B6",
    "#ff7f00",
    "#999999",
    "#a65628",
    "#f781bf",
    "#A72703",
]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#999999"
RULE = "#d9d8d3"

plt.rcParams["font.size"] = 16

# --------------------------------------------------------------------------
# How a system is classified. At module level so that the listing printed by
# main() is produced by the same function that colours the plot, and cannot
# drift from it.
# --------------------------------------------------------------------------
#: Membrane systems are not named "membrane": they are the MemProtMD deposits
#: and two of the GROMACS benchmarks. Matching on the word alone silently put
#: every one of them in the protein family.
MEMBRANE = ("memprotmd", "membrane", "aquaporin", "ion_channel")

#: The hexABC deposits arrive solvent-stripped, so every particle is solute.
#: They belong in their own category: metadata scales with solute, not with
#: particle count, and a solvated box is roughly 90% water that costs almost
#: nothing to describe -- the topology stores one water type and a count.
#: Mixed with the rest they read as outliers; separated, they are the only
#: systems showing what describing a biomolecule actually costs.
SOLVENT_FREE = ("hexabc",)

#: Coarse-grained systems, where one MARTINI bead stands for about four heavy
#: atoms. Their per-particle cost is low for a reason of resolution, not of
#: chemistry: there are simply fewer particles to describe, each carrying the
#: same kind of entry. Grouped with the atomistic membranes they pulled that
#: family's lower bound down and made the comparison a mix of two things.
COARSE_GRAINED = ("_cg",)

#: Systems carrying a nucleic acid whose name does not say so. Satellite tobacco
#: mosaic virus is a protein capsid around an RNA genome -- its .tpr holds 619
#: adenine and 330 uracil residues alongside the capsid protein -- so matching
#: on "rna" alone filed it as protein only.
NUCLEIC_ACID = ("dna", "rna", "stmv")

#: Satellite tobacco mosaic virus is kept apart from the other nucleic-acid
#: systems. Its topology declares 190 molecule types of one copy each, so unlike
#: its groupmates (repetition factor 21-41) nothing in it repeats as declared,
#: and it writes 163 ms of metadata against their 3.7-5.5 ms. Grouped with them
#: it is invisible behind the median while dragging it; alone it is the case
#: that shows what a topology declaring no reuse costs. Checked before
#: NUCLEIC_ACID, which also matches "stmv".
CAPSID = ("stmv",)
CAPSID_FAMILY = "virus capsid"

#: Names that contain a nucleic-acid token but no nucleic acid. Ribonuclease is
#: an enzyme that degrades RNA; bench_rnase holds moltypes Protein, SOL and CL
#: and not a single nucleotide. Checked before NUCLEIC_ACID.
NOT_NUCLEIC_ACID = ("rnase",)

#: Families drawn with a categorical hue. Solvation and resolution are marked
#: by shape instead -- they are properties of the representation rather than of
#: the chemistry -- so the scatter needs only three of the palette's slots.
ORDER = ["protein in solvent", "nucleic acid", "membrane (all-atoms)"]
SOLVENT_ONLY = "solvent-only"
ALL_SOLUTE = "solute-only"
COARSE = "membrane (coarse-grained)"

#: Legend order, which is not draw order. Drawing goes hue-coded families first
#: and shape-coded ones after, so the two membrane entries would otherwise be
#: separated by the solvation extremes. Related families read together here.
LEGEND_ORDER = [*ORDER, CAPSID_FAMILY, COARSE, ALL_SOLUTE, SOLVENT_ONLY]


def family(system: str) -> str:
    if any(token in system for token in SOLVENT_FREE):
        return ALL_SOLUTE
    if any(system.endswith(token) or token in system for token in COARSE_GRAINED):
        return COARSE
    if any(token in system for token in MEMBRANE):
        return "membrane (all-atoms)"
    if any(token in system for token in CAPSID):
        return CAPSID_FAMILY
    if any(token in system for token in NUCLEIC_ACID) and not any(
        token in system for token in NOT_NUCLEIC_ACID
    ):
        return "nucleic acid"
    if system == "solvent_ions":
        return SOLVENT_ONLY
    return "protein in solvent"


# --------------------------------------------------------------------------
# Writing-time decomposition, measured by timers inside the GROMACS writer.
# Process wall time and the 'Write traj.' cycle counter both bracket the file
# initialisation together with per-frame output, so neither isolates the
# metadata; these timers do.
# --------------------------------------------------------------------------
#: Phase as the writer reports it -> label, in the stacking order of the bars
#: (first at the bottom).
PHASES = (
    ("topology", "Topology"),
    ("connectivity", "Bond connectivity"),
    ("parameters", "MD parameters"),
    ("force_field", "Force field"),
    ("selections", "Selection groups"),
    ("trajectory_scaffold", "Trajectory scaffold"),
    ("h5md_metadata", "H5MD metadata"),
)

#: Colour per phase, keyed so that reordering the stack does not repaint it.
PHASE_COLORS = {
    "h5md_metadata": COLORS[0],
    "trajectory_scaffold": COLORS[1],
    "topology": COLORS[2],
    "connectivity": COLORS[3],
    "force_field": COLORS[4],
    "parameters": COLORS[5],
    "selections": COLORS[6],
}


#: Benchmark property name in the gtest JSON -> phase name used here.
GTEST_PHASES = {
    "metadata": "h5md_metadata",
    "particle_block": "trajectory_scaffold",
    "topology": "topology",
    "bond_connectivity": "connectivity",
    "force_field": "force_field",
    "md_parameters": "parameters",
    "selection_groups": "selections",
}


def load_timing(path):
    """Per-system median milliseconds for each phase, over the repeats.

    ``path`` is either a glob of the gtest JSON files written by
    ``h5md-test --gtest_filter=*BenchmarkDecomposeWritingProcess*
    --gtest_output=json:...`` (one file per repeat), or a times.tsv table.
    Reading the JSON directly keeps the panel tied to the benchmark output
    itself rather than to a derived table that can go missing.

    Returns (systems sorted by size, {phase: [ms per system]}, n_repeats) or
    None when nothing is found -- the storage panels stand without it.
    """
    import csv
    import glob
    from collections import defaultdict

    per = defaultdict(lambda: defaultdict(list))
    natoms, repeats = {}, defaultdict(set)
    files = sorted(glob.glob(str(path)))
    if not files:
        return None
    for repeat, name in enumerate(files):
        if name.endswith(".json"):
            for suite in json.loads(Path(name).read_text())["testsuites"]:
                for test in suite["testsuite"]:
                    if "metadata_elapsed_us" not in test:
                        continue
                    # "bena_12k__npt-prahman.tpr" is one protocol of "bena_12k"
                    system = test["value_param"].strip('"').split("__")[0]
                    natoms[system] = int(test["natoms"])
                    repeats[system].add(repeat)
                    for key, phase in GTEST_PHASES.items():
                        per[system][phase].append(
                            int(test[f"{key}_elapsed_us"]) / 1000.0
                        )
        else:
            with open(name) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    per[row["system"]][row["phase"]].append(float(row["ms"]))
                    natoms[row["system"]] = int(row["natoms"])
                    repeats[row["system"]].add(row["repeat"])
    if not per:
        return None
    systems = sorted(per, key=lambda s: natoms[s])
    series = {
        phase: [st.median(per[s].get(phase, [0.0])) for s in systems]
        for phase, _ in PHASES
    }
    # The usual repeat count, not the minimum: a run still in progress leaves
    # its last system with one repeat and would otherwise retitle the panel.
    from collections import Counter

    counts = Counter(len(r) for r in repeats.values())
    return systems, natoms, series, counts.most_common(1)[0][0]


#: The eight systems drawn in the writing-time panel. Seven are one per family,
#: each the member whose phase profile is closest to its family's median, so no
#: bar is a composite of systems that do not resemble each other. The eighth,
#: the largest solvated system, is added because the text contrasts it with the
#: capsid. Between them they span the corpus: repetition factor 1 to 4,165 and
#: 2 to 190 molecule types. Per-family distributions, with every system named,
#: are in the supplementary figure.
EXEMPLARS = (
    ("solvent_ions", "solvent\nand ions"),
    ("hexabc_a00186_dna", "small\nsolute"),
    ("memprotmd_7e6u_cg", "membrane\ncoarse-gr."),
    ("dna_amber", "nucleic\nacid"),
    ("misato_6dyw", "protein\nin solvent"),
    ("memprotmd_5nc4", "membrane\nall-atom"),
    ("bena_826k", "largest\nsolvated"),
    ("bench_stmv", "virus\ncapsid"),
)


def timing_panel(ax, timing, exemplars=EXEMPLARS):
    """Phase shares of the write time for a few named systems.

    Each bar is one measured system rather than an average over a family. A
    family median had to be renormalised to sum to 100, which produces a
    composition no system has, and it concealed the range it was drawn from:
    within the solvated proteins the connectivity share runs from 2.8% to 44%.
    """
    systems, natoms, series, _ = timing
    index = {name: i for i, name in enumerate(systems)}

    drawn = [(name, label) for name, label in exemplars if name in index]
    missing = [name for name, _ in exemplars if name not in index]
    if missing:
        # Silence here would show a figure with fewer bars than intended.
        print(
            f"  note: timing panel is missing {len(missing)} exemplar(s): "
            f"{', '.join(missing)}"
        )

    #: Ordered by repetition factor would be ideal, but the panel only receives
    #: the timings; EXEMPLARS is therefore already written in the intended order.
    x = np.arange(len(drawn))
    bottom = np.zeros(len(drawn))
    for phase, label in PHASES:
        pct = []
        for name, _ in drawn:
            i = index[name]
            total = sum(series[q][i] for q, _ in PHASES)
            pct.append(100.0 * series[phase][i] / total if total else 0.0)
        pct = np.array(pct)
        ax.bar(x, pct, 0.72, bottom=bottom, label=label, color=PHASE_COLORS[phase])
        bottom += pct

    ax.set_ylabel("Share of write time (%)")
    ax.set_xticks(x)
    # Eight bars in a narrow panel leave no room for horizontal two-line names;
    # angled single-line labels stay legible and cannot run into each other.
    ax.set_xticklabels(
        [
            f"{label.replace(chr(10), ' ')} ({natoms[name] // 1000}K)"
            for name, label in drawn
        ],
        fontsize=9,
        rotation=35,
        ha="right",
        rotation_mode="anchor",
    )
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0, fontsize=12)


# ==========================================================================
def cost_figure(path_stem: Path, rows, timing=None):
    """What the metadata costs: in bytes stored, and in time to write.

    One marker per chemical system, not per reference file: across the
    protocols of a system the stored metadata varies by under 1%, so plotting
    every file would show tight clusters and imply a precision the measurement
    does not have.
    """
    import statistics as st
    from collections import defaultdict

    per_system = defaultdict(list)
    for r in rows:
        per_system[r["system"]].append(r)

    points = []
    for name, rs in per_system.items():
        n_ua = rs[0]["n_unique_atoms"]
        points.append(
            {
                "system": name,
                "family": family(name),
                "n_particles": rs[0]["n_particles"],
                "n_unique_atoms": n_ua,
                "metadata_kb": st.median(r["metadata_bytes"] for r in rs) / 1000,
                # Connectivity is excluded: it is a global bond list, so it grows
                # with every bond in the system, waters included, and not with the
                # unique atoms this panel normalises by.
                "bytes_per_unique_atom": st.median(
                    r["metadata_bytes"] - r["connectivity_bytes"] for r in rs
                )
                / n_ua,
                "n_files": len(rs),
            }
        )
    points.sort(key=lambda p: p["n_particles"])

    order = ORDER
    #: One colour per family, shared by its markers and its annotations.
    colors = {
        order[0]: COLORS[0],
        order[1]: COLORS[1],
        order[2]: COLORS[2],
        CAPSID_FAMILY: COLORS[3],
        SOLVENT_ONLY: "black",
        ALL_SOLUTE: COLORS[5],
        COARSE: COLORS[6],
    }
    #: human-readable names for the points worth calling out
    LABELS = {
        "bench_stmv": "STMV",
        "solvent_ions": "water and ions",
        "bena_826k": "largest\nprotein-in-solvent",
        "hexabc_a00186_dna": "nucleic-acid\nsolvent-stripped",
        "prot_gromos": "GROMOS54A7\n(united atom)",
    }

    n_panels = 3 if timing else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    #: (x, y) of each storage panel. A answers "what will my system cost" on the
    #: axis a user knows; B answers "what drives it" on the axis that explains it.
    PANELS = (
        (axes[0], "n_particles", "metadata_kb"),
        (axes[1], "n_unique_atoms", "bytes_per_unique_atom"),
    )
    for ax, xkey, _ in PANELS:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        ax.set_xscale("log")
        # Derived, not fixed: a hard lower bound silently clipped the six
        # smallest solvent-free deposits (711-6,747 particles) off the plot.
        ax.set_xlim(
            min(p[xkey] for p in points) * 0.6, max(p[xkey] for p in points) * 1.4
        )
    axes[0].set_xlabel("particles in the system")
    axes[1].set_xlabel(r"unique atoms across molecule types, $N_{ua}$")

    for name in order:
        ps = [p for p in points if p["family"] == name]
        if not ps:
            continue
        for ax, xkey, ykey in PANELS:
            ax.scatter(
                [p[xkey] for p in ps],
                [p[ykey] for p in ps],
                s=26,
                color=colors[name],
                label=name,
                edgecolor="white",
                linewidth=0.7,
                zorder=3,
            )

    # The two extremes of solvation are marked by shape rather than by a fourth
    # and fifth hue: what separates them from the rest is how much of the box is
    # solute, not what chemistry it contains.
    shaped = (
        (
            CAPSID_FAMILY,
            {
                "marker": "*",
                "s": 90,
                "facecolor": colors[CAPSID_FAMILY],
                "edgecolor": "white",
                "linewidth": 0.6,
            },
        ),
        (
            SOLVENT_ONLY,
            {
                "marker": "o",
                "s": 30,
                "facecolor": "none",
                "edgecolor": colors[SOLVENT_ONLY],
                "linewidth": 0.6,
            },
        ),
        (
            ALL_SOLUTE,
            {
                "marker": "s",
                "s": 30,
                "facecolor": colors[ALL_SOLUTE],
                "edgecolor": "white",
                "linewidth": 0.6,
            },
        ),
        (
            COARSE,
            {
                "marker": "^",
                "s": 38,
                "facecolor": colors[COARSE],
                "edgecolor": "white",
                "linewidth": 0.6,
            },
        ),
    )
    for fam, style in shaped:
        ps = [p for p in points if p["family"] == fam]
        if not ps:
            continue
        for ax, xkey, ykey in PANELS:
            ax.scatter(
                [p[xkey] for p in ps],
                [p[ykey] for p in ps],
                label=fam,
                zorder=4,
                **style,
            )

    axes[0].set_yscale("log")
    axes[0].set_ylabel("metadata stored (kB)")
    # axes[0].set_title("Metadata growth with the system size", loc="center", color=INK)

    # Per unique atom, everything but connectivity is near-constant for the
    # atomistic systems: solvated proteins, nucleic acids, all-atom membranes and
    # solute-only deposits all sit on one line. The median of those families is
    # drawn so the flatness is read against a reference, not inferred.
    lo = min(p["bytes_per_unique_atom"] for p in points)
    hi = max(p["bytes_per_unique_atom"] for p in points)
    axes[1].set_ylabel("metadata bytes per unique atom\n(excluding connectivity)")
    axes[1].set_yscale("log")
    axes[1].set_ylim(lo * 0.7, hi * 1.6)
    atomistic = [
        p["bytes_per_unique_atom"]
        for p in points
        if p["family"] in [*ORDER, ALL_SOLUTE]
    ]
    if atomistic:
        median = st.median(atomistic)
        axes[1].axhline(
            median, color="black", linewidth=0.8, linestyle="--", alpha=0.6, zorder=2
        )
        axes[1].annotate(
            f"atomistic median {median:.0f} B",
            xy=(1.0, median),
            xycoords=("axes fraction", "data"),
            xytext=(-4, 5),
            textcoords="offset points",
            fontsize=10,
            alpha=0.75,
            ha="right",
            va="bottom",
        )

    # Every system named in LABELS is annotated, on both storage panels unless
    # ONLY_IN says otherwise, so a point can be followed from one axis to the
    # other. Previously only the two extremes were labelled, so an entry added
    # for any other system had no effect and looked like the point was missing
    # from the plot.
    #: system -> the one panel index it is labelled on. GROMOS is unremarkable
    #: in A and stands out only per unique atom, where united atoms show.
    ONLY_IN = {"prot_gromos": 1}
    #: Panel B collapses the atomistic systems onto one band, so a label set
    #: directly under a band member lands on its neighbours. Those are pulled
    #: into the empty space below the band, with a leader line back to the point.
    #: (panel index, system) -> (dx, dy, ha) in points.
    MOVED = {
        (1, "bena_826k"): (30, -35, "left"),
        (1, "hexabc_a00186_dna"): (-40, -55, "right"),
        (1, "prot_gromos"): (-25, 15, "right"),
    }
    for p in points:
        if p["system"] not in LABELS:
            continue
        color = colors.get(p["family"], INK2)
        for i, (ax, xkey, ykey) in enumerate(PANELS):
            if ONLY_IN.get(p["system"], i) != i:
                continue
            dx, dy, ha = MOVED.get((i, p["system"]), (0, -5, "center"))
            arrow = (
                {
                    "arrowstyle": "-",
                    "color": color,
                    "alpha": 0.6,
                    "linewidth": 0.7,
                    "shrinkA": 0,
                    "shrinkB": 3,
                }
                if (i, p["system"]) in MOVED
                else None
            )
            ax.annotate(
                LABELS[p["system"]],
                xy=(p[xkey], p[ykey]),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=10,
                color=color,
                alpha=0.75,
                ha=ha,
                va="top",
                linespacing=1.3,
                arrowprops=arrow,
            )

    # A name in LABELS that matches no system is a typo or a renamed system; say
    # so rather than silently drawing nothing.
    unplaced = set(LABELS) - {p["system"] for p in points}
    if unplaced:
        print(
            f"  note: LABELS names {len(unplaced)} system(s) not in the data: "
            f"{', '.join(sorted(unplaced))}"
        )

    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered = [(by_label[name], name) for name in LEGEND_ORDER if name in by_label]
    # Anything not named in LEGEND_ORDER still appears, rather than vanishing
    # silently if a family is added and the list is not updated.
    ordered += [
        (h, label) for label, h in by_label.items() if label not in LEGEND_ORDER
    ]
    axes[0].legend(
        [h for h, _ in ordered],
        [label for _, label in ordered],
        loc="upper left",
        fontsize=10,
        labelcolor=INK2,
        handletextpad=0.35,
        borderaxespad=0.3,
        borderpad=0.2,
    )

    if timing:
        timing_panel(axes[2], timing)
        axes[2].set_title(f"Median over {timing[3]} runs", color=INK)

    for ax, letter in zip(axes, "ABC"):
        ax.text(
            -0.1, 1.05, letter, transform=ax.transAxes, fontsize=20, fontweight="bold"
        )

    fig.tight_layout()
    fig.savefig(f"{path_stem}.pdf", dpi=350, bbox_inches="tight")
    fig.savefig(f"{path_stem}.png", dpi=350, bbox_inches="tight")
    plt.close(fig)


def unique_atoms(h5md_path) -> int:
    """N_ua: atoms summed over the molecule types, each type counted once.

    Read from the file's own topology -- one ``particle_count`` per molecule
    type group -- so the panel needs nothing beyond the reference file. It
    agrees with the count taken from the .tpr for all 58 systems.
    """
    import h5py

    with h5py.File(h5md_path, "r") as f:
        topology = f["topology"]
        return sum(
            int(np.atleast_1d(topology[name].attrs["particle_count"])[0])
            for name in topology
            if isinstance(topology[name], h5py.Group)
        )


# ==========================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--storage", default=_r("storage.json"))
    ap.add_argument("--out", default=str(OUT / "validation5_cost"))
    ap.add_argument(
        "--timing",
        default=str(
            Path(__file__).resolve().parents[1] / "data" / "decompose_run*.json"
        ),
        help="gtest JSON of BenchmarkDecomposeWritingProcess, one file "
        "per repeat (a glob), or a times.tsv table",
    )
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    storage = json.loads(Path(args.storage).read_text())["files"]

    rows = []
    n_unique = {}
    for rec in storage.values():
        if rec.get("n_particles") and rec.get("metadata_bytes"):
            # N_ua depends only on the system, so read it from one of its files.
            if rec["system"] not in n_unique:
                n_unique[rec["system"]] = unique_atoms(rec["file"])
            rows.append(
                {
                    "system": rec["system"],
                    "protocol": rec["protocol"],
                    "n_particles": rec["n_particles"],
                    "n_unique_atoms": n_unique[rec["system"]],
                    "metadata_bytes": rec["metadata_bytes"],
                    "connectivity_bytes": rec["per_group_bytes"].get("connectivity", 0),
                    "bytes_per_particle": rec["metadata_bytes_per_particle"],
                }
            )

    if not rows:
        print("  no storage measurements found; nothing to plot")
        return 0

    print_families(rows)

    timing = load_timing(args.timing)
    if timing is None:
        print(
            f"  note: no writing-time decomposition at {args.timing}; panel C omitted"
        )
    else:
        systems, _, _, n_repeats = timing
        missing = sorted({r["system"] for r in rows} - set(systems))
        print(
            f"  writing-time decomposition: {len(systems)} systems, {n_repeats} repeats"
        )
        if missing:
            print(
                f"  note: {len(missing)} system(s) have storage but no timing: "
                f"{', '.join(missing)}"
            )

    cost_figure(Path(args.out), rows, timing=timing)
    print(f"  wrote {args.out}.png  ({len(rows)} storage points)")
    return 0


def print_families(rows):
    """List every system under the family it is plotted as.

    Classification is by name, so a system that is renamed or added can fall
    into the wrong family without the figure looking wrong. Printing it makes
    the grouping checkable instead of implicit.
    """
    import statistics as st
    from collections import defaultdict

    per_system = defaultdict(list)
    for r in rows:
        per_system[r["system"]].append(r)

    grouped = defaultdict(list)
    for name, rs in per_system.items():
        grouped[family(name)].append(
            (name, rs[0]["n_particles"], st.median(r["bytes_per_particle"] for r in rs))
        )

    print(f"\n  {len(per_system)} systems in {len(grouped)} families\n")
    for fam in [*ORDER, COARSE, SOLVENT_ONLY, ALL_SOLUTE]:
        members = sorted(grouped.get(fam, []), key=lambda m: m[1])
        if not members:
            continue
        bpp = [m[2] for m in members]
        print(
            f"  {fam}  ({len(members)} systems, "
            f"{min(bpp):.0f}-{max(bpp):.0f} B/particle)"
        )
        for name, n, b in members:
            print(f"      {name:30s} {n:>9,} particles   {b:6.1f} B/particle")
        print()


if __name__ == "__main__":
    raise SystemExit(main())
