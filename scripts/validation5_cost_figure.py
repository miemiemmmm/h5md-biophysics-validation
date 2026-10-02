#!/usr/bin/env python
"""
Validation 5 figure: what the metadata costs to store.

    validation5_cost.png  storage against system size

Written to <results>/figures/. The validation_framework schematic is a diagram
rather than a measurement and lives in renderer/render_framework_figure.py.

Used by:
    Validation 5: storage overhead
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import pathlib
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

# Palette and rcParams mirror gromacs_MR2/plot_benchmark_metadata.py, so the
# two figures of this work read as one set.
COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#5DD3B6',
          '#ff7f00', '#999999', '#a65628', '#f781bf', '#A72703']
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#999999"
RULE = "#d9d8d3"

plt.rcParams['font.size'] = 16

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
LEGEND_ORDER = ORDER + [COARSE, ALL_SOLUTE, SOLVENT_ONLY]


def family(system: str) -> str:
    if any(token in system for token in SOLVENT_FREE):
        return ALL_SOLUTE
    if any(system.endswith(token) or token in system for token in COARSE_GRAINED):
        return COARSE
    if any(token in system for token in MEMBRANE):
        return "membrane (all-atoms)"
    if (any(token in system for token in NUCLEIC_ACID)
            and not any(token in system for token in NOT_NUCLEIC_ACID)):
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
#: Phase as the writer reports it -> label, in the stacking order of the bars.
PHASES = (
    ("h5md_metadata", "H5MD metadata"),
    ("trajectory_scaffold", "Trajectory scaffold"),
    ("topology", "Topology"),
    ("connectivity", "Bond connectivity"),
    ("force_field", "Force field"),
    ("parameters", "MD parameters"),
    ("selections", "Selection groups"),
)


def load_timing(path):
    """Per-system median milliseconds for each phase, over the repeats.

    Returns (systems sorted by size, {phase: [ms per system]}, n_repeats) or
    None when the file is absent -- the storage panels stand without it.
    """
    import csv
    import statistics as st
    from collections import defaultdict

    path = Path(path)
    if not path.is_file():
        return None
    per = defaultdict(lambda: defaultdict(list))
    natoms, repeats = {}, defaultdict(set)
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            per[row["system"]][row["phase"]].append(float(row["ms"]))
            natoms[row["system"]] = int(row["natoms"])
            repeats[row["system"]].add(row["repeat"])
    systems = sorted(per, key=lambda s: natoms[s])
    series = {phase: [st.median(per[s].get(phase, [0.0])) for s in systems]
              for phase, _ in PHASES}
    # The usual repeat count, not the minimum: a run still in progress leaves
    # its last system with one repeat and would otherwise retitle the panel.
    from collections import Counter
    counts = Counter(len(r) for r in repeats.values())
    return systems, natoms, series, counts.most_common(1)[0][0]


def timing_panel(ax, timing):
    """Stacked bars: each phase as a percentage of that system's write time."""
    systems, natoms, series, _ = timing
    x = np.arange(len(systems))
    totals = [sum(series[phase][i] for phase, _ in PHASES)
              for i in range(len(systems))]
    bottom = np.zeros(len(systems))
    for (phase, label), color in zip(PHASES, COLORS):
        pct = np.array([100.0 * series[phase][i] / totals[i] if totals[i] else 0.0
                        for i in range(len(systems))])
        ax.bar(x, pct, 0.9, bottom=bottom, label=label, color=color)
        bottom += pct
    ax.set_xlabel("System size (Katoms)")
    ax.set_ylabel("Percentage of total write time (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{natoms[s] // 1000}" for s in systems],
                       rotation=75, ha="right", fontsize=6)
    ax.set_ylim(0, 100)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0,
              fontsize=14)


# ==========================================================================
def cost_figure(path_stem: Path, rows, timing=None):
    """What the metadata costs: in bytes stored, and in time to write.

    One marker per chemical system, not per reference file: across the
    protocols of a system the stored metadata varies by under 1%, so plotting
    every file would show tight clusters and imply a precision the measurement
    does not have.
    """
    import statistics as st  # noqa: F401  (used by both halves)
    from collections import defaultdict

    per_system = defaultdict(list)
    for r in rows:
        per_system[r["system"]].append(r)

    points = []
    for name, rs in per_system.items():
        points.append({
            "system": name,
            "family": family(name),
            "n_particles": rs[0]["n_particles"],
            "metadata_kib": st.median(r["metadata_bytes"] for r in rs) / 1024,
            "bytes_per_particle": st.median(r["bytes_per_particle"] for r in rs),
            "n_files": len(rs),
        })
    points.sort(key=lambda p: p["n_particles"])

    order = ORDER
    colors = {order[0]: COLORS[0], order[1]: COLORS[1], order[2]: COLORS[2]}
    #: human-readable names for the points worth calling out
    LABELS = {
        "bench_stmv": "STMV",
        "solvent_ions": "water and ions",
        "bena_826k": "largest\nprotein-in-solvent",
        "hexabc_a00186_dna": "nucleic-acid\nsolvent-stripped",
    }

    n_panels = 3 if timing else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    for ax in axes[:2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        ax.set_xscale("log")
        # Derived, not fixed: a hard lower bound silently clipped the six
        # smallest solvent-free deposits (711-6,747 particles) off the plot.
        ax.set_xlim(min(p["n_particles"] for p in points) * 0.6,
                    max(p["n_particles"] for p in points) * 1.4)
        ax.set_xlabel("particles in the system")

    for name in order:
        ps = [p for p in points if p["family"] == name]
        if not ps:
            continue
        for ax, key in ((axes[0], "metadata_kib"), (axes[1], "bytes_per_particle")):
            ax.scatter([p["n_particles"] for p in ps], [p[key] for p in ps],
                       s=26, color=colors[name], label=name,
                       edgecolor="white", linewidth=0.7, zorder=3)

    # The two extremes of solvation are marked by shape rather than by a fourth
    # and fifth hue: what separates them from the rest is how much of the box is
    # solute, not what chemistry it contains.
    shaped = ((SOLVENT_ONLY, {"marker": "o", "s": 30, "facecolor": "none",
                              "edgecolor": "black", "linewidth": 0.6}),
              (ALL_SOLUTE, {"marker": "s", "s": 30, "facecolor": COLORS[5],
                            "edgecolor": "white", "linewidth": 0.6}),
              (COARSE, {"marker": "^", "s": 38, "facecolor": COLORS[6],
                        "edgecolor": "white", "linewidth": 0.6}))
    for fam, style in shaped:
        ps = [p for p in points if p["family"] == fam]
        if not ps:
            continue
        for ax, key in ((axes[0], "metadata_kib"), (axes[1], "bytes_per_particle")):
            ax.scatter([p["n_particles"] for p in ps], [p[key] for p in ps],
                       label=fam, zorder=4, **style)

    axes[0].set_yscale("log")
    axes[0].set_ylabel("metadata stored (KiB)")
    # axes[0].set_title("Metadata growth with the system size", loc="center", color=INK)

    lo = min(p["bytes_per_particle"] for p in points)
    hi = max(p["bytes_per_particle"] for p in points)
    axes[1].set_ylabel("metadata bytes per particle")
    # axes[1].set_title(f"Per particle: {lo:.0f}-{hi:.0f} bytes, set by solute fraction",
    #                   loc="left", color=INK)
    # Log: most systems fall in a decade near the bottom of the range, and on a
    # linear axis they compress into a line while the solvent-free band at the
    # top takes three quarters of the height.
    axes[1].set_yscale("log")
    axes[1].set_ylim(lo * 0.7, hi * 1.6)

    # Shade the solvated band only. The solvent-free deposits sit above it, and
    # including them would make the band span the whole axis and say nothing.
    solvated = [p for p in points if p["family"] != ALL_SOLUTE]
    if solvated:
        axes[1].axhspan(min(p["bytes_per_particle"] for p in solvated),
                        max(p["bytes_per_particle"] for p in solvated),
                        color=MUTED, alpha=0.07, zorder=0)

    # Name the extremes, so the spread is read as solvation rather than noise:
    # water is most of a solvated box and costs almost nothing to describe --
    # the topology stores one water type and a count -- so per-particle cost
    # tracks how much of the system is solute, not how large it is.
    # Every system named in LABELS is annotated. Previously only the two
    # extremes were, so an entry added for any other system had no effect and
    # looked like the point was missing from the plot.
    # Set to the left of the marker and vertically centred on it, so a label
    # never sits on the point above or below it in the crowded low-cost band.
    for p in points:
        if p["system"] not in LABELS:
            continue
        axes[1].annotate(LABELS[p["system"]],
                         xy=(p["n_particles"], p["bytes_per_particle"]),
                         xytext=(0, -5), textcoords="offset points",
                         fontsize=10, color=colors.get(p["family"], INK2),
                         alpha=0.75, ha="center", va="top", linespacing=1.3)

    # A name in LABELS that matches no system is a typo or a renamed system; say
    # so rather than silently drawing nothing.
    unplaced = set(LABELS) - {p["system"] for p in points}
    if unplaced:
        print(f"  note: LABELS names {len(unplaced)} system(s) not in the data: "
              f"{', '.join(sorted(unplaced))}")

    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ordered = [(by_label[name], name) for name in LEGEND_ORDER if name in by_label]
    # Anything not named in LEGEND_ORDER still appears, rather than vanishing
    # silently if a family is added and the list is not updated.
    ordered += [(h, l) for l, h in by_label.items()
                if l not in LEGEND_ORDER]
    axes[0].legend([h for h, _ in ordered], [l for _, l in ordered],
                   loc="upper left", fontsize=10, labelcolor=INK2,
                   handletextpad=0.35, borderaxespad=0.3, borderpad=0.2)

    if timing:
        timing_panel(axes[2], timing)
        axes[2].set_title(f"Writing time breakdown\n(median over "
                          f"{timing[3]} runs)", color=INK)

    for ax, letter in zip(axes, "ABC"):
        ax.text(-0.1, 1.05, letter, transform=ax.transAxes,
                fontsize=20, fontweight="bold")

    fig.tight_layout()
    fig.savefig(f"{path_stem}.pdf", dpi=350, bbox_inches="tight")
    fig.savefig(f"{path_stem}.png", dpi=350, bbox_inches="tight")
    plt.close(fig)


# ==========================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validation", default=_r("validation.json"))
    ap.add_argument("--out", default=str(OUT / "validation5_cost"))
    ap.add_argument("--timing", default="/tmp/claude-1000/validate/timing/times.tsv",
                    help="per-phase write timings from the instrumented writer")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    results = json.loads(Path(args.validation).read_text())
    ok = {k: v for k, v in results.items() if v.get("status") == "ok"}

    rows = []
    for key, v in ok.items():
        st = v.get("storage") or {}
        if st.get("n_particles") and st.get("metadata_bytes"):
            rows.append({"system": v["system"], "protocol": v["protocol"],
                         "n_particles": st["n_particles"],
                         "metadata_bytes": st["metadata_bytes"],
                         "bytes_per_particle": st["metadata_bytes_per_particle"]})

    if not rows:
        print("  no storage measurements found; nothing to plot")
        return 0

    print_families(rows)

    timing = load_timing(args.timing)
    if timing is None:
        print(f"  note: no writing-time decomposition at {args.timing}; "
              "panel C omitted")
    else:
        systems, _, _, n_repeats = timing
        missing = sorted({r["system"] for r in rows} - set(systems))
        print(f"  writing-time decomposition: {len(systems)} systems, "
              f"{n_repeats} repeats")
        if missing:
            print(f"  note: {len(missing)} system(s) have storage but no timing: "
                  f"{', '.join(missing)}")

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
        grouped[family(name)].append((
            name, rs[0]["n_particles"],
            st.median(r["bytes_per_particle"] for r in rs)))

    print(f"\n  {len(per_system)} systems in {len(grouped)} families\n")
    for fam in ORDER + [COARSE, SOLVENT_ONLY, ALL_SOLUTE]:
        members = sorted(grouped.get(fam, []), key=lambda m: m[1])
        if not members:
            continue
        bpp = [m[2] for m in members]
        print(f"  {fam}  ({len(members)} systems, "
              f"{min(bpp):.0f}-{max(bpp):.0f} B/particle)")
        for name, n, b in members:
            print(f"      {name:30s} {n:>9,} particles   {b:6.1f} B/particle")
        print()


if __name__ == "__main__":
    raise SystemExit(main())
