#!/usr/bin/env python
"""
Molecule-block-based against system-based topology.

    python scripts/render_layout_comparison.py
    python scripts/render_layout_comparison.py --out somewhere/figure

Two ways to store the same topology. The system-based layout walks every atom of
the system once; the molecule-block-based layout stores each molecule type once
and lets the molecule blocks say how often it repeats. Which is smaller is not a
property of the schema but of the system: the saving is whatever the repetition
factor is, so a solvated protein gains a great deal and a crystal of distinct
chains gains nothing.

Plotted against particle count because that is the axis a reader expects, but
particle count is not what drives the ratio -- the annotation on the parity line
is there to say so, and the correlation printed by main() quantifies it.

Reads the log of the two layout cases of the benchmark in
data/h5md_benchmark.patch (gtest, two variants per system); the repetition
factor comes from each system's nvt-vrescale reference file.

Drawn in the style of panel A of scripts/validation5_cost_figure.py (family
colours, markers, legend and labels are imported from it), so the two figures
read alike.
"""

from __future__ import annotations

import argparse
import collections
import re
import statistics as st
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from validation5_cost_figure import (
    FAMILY_COLORS,
    INK2,
    LABEL_STYLE,
    REFERENCE_LINE,
    family,
    family_legend,
    leader,
    scatter_families,
    style_axes,
)

from config.matrix import DATA_ROOT, RESULTS_DIR

#: The benchmark prints one line per variant per system. "compound" is the
#: system-based variant (writeCompoundTopologyData over a whole-system
#: AtomRange), "multids" the molecule-block-based one
#: (writeMoleculeTypesAsCompound over the molecule types).
LINE = re.compile(
    r"\[(compound|multids)\s*\]\s+(\S+)\s+"
    r"atoms=\s*(\d+)\s+time=(\d+) us\s+size=(\d+) B"
)

VARIANT = {"compound": "system", "multids": "block"}


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def repetition_factor(system: str) -> float | None:
    """Particles over the particles of the unique molecule types, from the
    system's nvt-vrescale reference file; None when that file is missing."""
    ref = Path(DATA_ROOT) / "reference_files" / f"{system}__nvt-vrescale.h5md"
    if not ref.exists():
        return None
    with h5py.File(ref, "r") as f:
        topo = f["topology"]
        blocks = [(_s(r["name"]), int(r["count"])) for r in topo["molecule_blocks"][()]]
        total = sum(len(topo[n]["particle_data"]) * c for n, c in blocks)
        unique = sum(len(topo[n]["particle_data"]) for n in {n for n, _ in blocks})
    return total / unique if unique else None


def load(log_path: Path):
    """One record per system: particles, both sizes, and their ratio.

    Each gtest case writes its own temporary file holding only that variant's
    data, so the two sizes are directly comparable. Repeats are medians.
    """
    times = collections.defaultdict(list)
    sizes: dict[tuple[str, str], int] = {}
    particles: dict[str, int] = {}

    with open(log_path) as handle:
        for line in handle:
            match = LINE.search(line)
            if not match:
                continue
            variant = VARIANT[match.group(1)]
            # "bena_12k__nve.tpr" is one protocol of the system "bena_12k"
            system = match.group(2).split("__")[0]
            times[(variant, system)].append(int(match.group(4)))
            sizes[(variant, system)] = int(match.group(5))
            particles[system] = int(match.group(3))

    #: Repetition is what actually drives the ratio, and it is not in the
    #: benchmark log; it comes from the reference files.
    repetition = {system: repetition_factor(system) for system in particles}

    points = []
    for system, n_particles in particles.items():
        if ("system", system) not in sizes or ("block", system) not in sizes:
            continue
        points.append(
            {
                "system": system,
                "family": family(system),
                "n_particles": n_particles,
                "system_based_mb": sizes[("system", system)] / 1e6,
                "block_based_mb": sizes[("block", system)] / 1e6,
                "size_ratio": sizes[("system", system)] / sizes[("block", system)],
                "time_ratio": (
                    st.median(times[("system", system)])
                    / st.median(times[("block", system)])
                ),
                "repetition": repetition.get(system),
            }
        )
    points.sort(key=lambda p: p["n_particles"])
    return points


#: Call-outs: the systems the text names, in their family's colour, with the
#: ratio filled in for {}. system -> (text, (dx, dy) in points, ha, va, leader).
#: The top of the plot is the free space: bena_457k and the water box are set
#: there, above their points; the largest protein, whose surroundings are
#: crowded, is pulled into the empty space below with a leader line.
CALLOUTS = {
    "bena_457k": ("bena_457k {}", (0, 5), "center", "bottom", False),
    "solvent_ions": ("water and ions {}", (30, 5), "center", "bottom", False),
    "bena_826k": ("largest protein\nin solvent {}", (-15, -58), "center", "top", True),
    "bench_stmv": ("STMV {}", (0, -5), "center", "top", False),
}

#: Wide enough that, at \linewidth of a two-column page, text and markers print
#: at the size of panel A of validation5_cost.pdf at \textwidth.
FIGSIZE = (8, 4.5)

#: Text of the call-outs and of the parity line: that of panel A, a size up.
ANNOTATION = {**LABEL_STYLE, "fontsize": 13.5}


def figure(path_stem: Path, points):
    """Style of panel A of scripts/validation5_cost_figure.py, from the same
    shared definitions."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    style_axes(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(
        min(p["n_particles"] for p in points) * 0.6,
        max(p["n_particles"] for p in points) * 1.4,
    )

    # Parity: above it the molecule-block layout is the smaller of the two,
    # below it the system-based one is.
    ax.axhline(1.0, **REFERENCE_LINE)
    ax.annotate(
        "equal size",
        xy=(1.0, 1.0),
        xycoords=("axes fraction", "data"),
        xytext=(-4, 5),
        textcoords="offset points",
        ha="right",
        va="bottom",
        **ANNOTATION,
    )

    scatter_families(ax, points, "n_particles", "size_ratio")

    # Room below parity for the labels, and above the data so that the legend
    # in the upper left clears every point.
    ax.set_ylim(0.4, 100)
    ax.set_xlabel("Number of particles in the system")
    ax.set_ylabel("Storage ratio (flat / block)")
    # Behind the marks: a point under the legend stays visible.
    family_legend(ax).set_zorder(1)

    # The systems the text names, and the count of those below parity; a number
    # on every point would say nothing the axis does not already say.
    for p in points:
        if p["system"] not in CALLOUTS:
            continue
        text, offset, ha, va, lead = CALLOUTS[p["system"]]
        ratio = p["size_ratio"]
        color = FAMILY_COLORS.get(p["family"], INK2)
        ax.annotate(
            text.format(f"{ratio:.2f}×" if ratio < 1 else f"{ratio:.0f}×"),
            xy=(p["n_particles"], ratio),
            xytext=offset,
            textcoords="offset points",
            color=color,
            ha=ha,
            va=va,
            linespacing=1.3,
            arrowprops=leader(color) if lead else None,
            **ANNOTATION,
        )
    below = [p for p in points if p["size_ratio"] < 1.0]
    if below:
        ax.annotate(
            f"{len(below)} systems below parity",
            xy=(min(p["n_particles"] for p in below), min(p["size_ratio"] for p in below)),
            xytext=(0, -8),
            textcoords="offset points",
            ha="left",
            va="top",
            **ANNOTATION,
        )
    unplaced = set(CALLOUTS) - {p["system"] for p in points}
    if unplaced:
        print(f"  note: CALLOUTS names systems not in the data: {sorted(unplaced)}")

    fig.tight_layout()
    fig.savefig(f"{path_stem}.pdf", dpi=350, bbox_inches="tight")
    fig.savefig(f"{path_stem}.png", dpi=350, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--layout-log",
        default=str(Path(__file__).resolve().parents[1] / "data" / "layout_benchmark.log"),
        help="gtest output of the BenchmarkAsFlatTopology and "
        "BenchmarkAsMoleculeBlockTopology cases (data/h5md_benchmark.patch)",
    )
    ap.add_argument(
        "--out", default=str(Path(RESULTS_DIR) / "figures" / "layout_comparison")
    )
    args = ap.parse_args()

    log = Path(args.layout_log)
    if not log.exists():
        print(f"  no benchmark log at {log}")
        return 1
    points = load(log)
    if not points:
        print(f"  no paired measurements in {log}")
        return 1

    ratios = [p["size_ratio"] for p in points]
    print(
        f"  {len(points)} systems, "
        f"{min(p['n_particles'] for p in points)}-"
        f"{max(p['n_particles'] for p in points)} particles"
    )
    print(
        f"  size ratio: median {st.median(ratios):.1f}x, "
        f"range {min(ratios):.1f}-{max(ratios):.1f}x"
    )
    print(f"  below parity: {sum(r < 1 for r in ratios)} systems")

    log_p = np.log10([p["n_particles"] for p in points])
    log_r = np.log10(ratios)
    print(
        f"  log-log correlation of ratio with particles : "
        f"r = {np.corrcoef(log_p, log_r)[0, 1]:.3f}"
    )
    with_rep = [p for p in points if p["repetition"]]
    if with_rep:
        log_rep = np.log10([p["repetition"] for p in with_rep])
        log_rr = np.log10([p["size_ratio"] for p in with_rep])
        print(
            f"  log-log correlation of ratio with repetition: "
            f"r = {np.corrcoef(log_rep, log_rr)[0, 1]:.3f}  "
            f"({len(with_rep)} systems)"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    figure(out, points)
    print(f"  wrote {out}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
