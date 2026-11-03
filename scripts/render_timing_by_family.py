#!/usr/bin/env python
"""
Writing time decomposed, every system named.

    python scripts/render_timing_by_family.py
    python scripts/render_timing_by_family.py --out somewhere/figure

The writing-time panel of the cost figure shows eight systems chosen to span
the corpus. This is the same measurement for all of them, one panel per family,
so that the range each exemplar stands for can be read rather than taken on
trust -- within the solvated proteins alone the connectivity share runs from
under 3% to 44%.

Reads the shipped benchmark output (data/decompose_run*.json), the same data
as panel C of the cost figure, and each system's nvt-vrescale reference file
for its size.

Style mirrors scripts/validation5_cost_figure.py.
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
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
    GTEST_PHASES,
    INK,
    INK2,
    LEGEND_ORDER,
    PHASE_COLORS,
    PHASES,
    family,
)

from config.matrix import DATA_ROOT, RESULTS_DIR

plt.rcParams["font.size"] = 11


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def load(pattern: str):
    """One record per system: its phase shares and the sizes that drive them.

    Timings are the shipped benchmark output (data/decompose_run*.json), the
    same data as panel C of the cost figure; the sizes come from each system's
    nvt-vrescale reference file.
    """
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    natoms = {}
    for name in sorted(glob.glob(pattern)):
        for suite in json.loads(Path(name).read_text())["testsuites"]:
            for test in suite["testsuite"]:
                if "metadata_elapsed_us" not in test:
                    continue
                system = test["value_param"].strip('"').split("__")[0]
                natoms[system] = int(test["natoms"])
                for key, phase in GTEST_PHASES.items():
                    per[system][phase].append(int(test[f"{key}_elapsed_us"]) / 1000.0)
    times = {(s, ph): v for s, d in per.items() for ph, v in d.items()}
    dims = {}
    for system, n in natoms.items():
        ref = Path(DATA_ROOT) / "reference_files" / f"{system}__nvt-vrescale.h5md"
        with h5py.File(ref, "r") as f:
            topo = f["topology"]
            blocks = [_s(r["name"]) for r in topo["molecule_blocks"][()]]
            unique = sum(len(topo[b]["particle_data"]) for b in set(blocks))
        dims[system] = {
            "n_particles": n,
            "n_moltypes": len(set(blocks)),
            "repetition": n / max(unique, 1),
        }

    points = []
    for system, d in dims.items():
        med = {
            phase: st.median(times[(system, phase)])
            for phase, _ in PHASES
            if times.get((system, phase))
        }
        total = sum(med.values())
        if not total:
            continue
        points.append(
            {
                "system": system,
                "family": family(system),
                "shares": {
                    phase: 100.0 * med.get(phase, 0.0) / total for phase, _ in PHASES
                },
                "total_ms": total,
                **d,
            }
        )
    return points


#: Which collection a system was taken from. The solvated proteins are two
#: thirds of the corpus and come from four sources whose preparation differs;
#: blocking them by source rather than by size shows whether the decomposition
#: is a property of the schema or of how one group built its inputs.
SOURCES = (
    ("bena_", "Bena"),
    ("misato_", "MISATO"),
    ("bench_", "GROMACS benchmark"),
    ("hexabc_", "hexABC"),
    ("memprotmd_", "MemProtMD"),
)
SELF_BUILT = "self-generated"

#: Headers are often wider than the block they name, so the long family names
#: are shortened and consecutive headers sit at two heights.
SHORT = {
    "membrane (all-atoms)": "membrane",
    "membrane (coarse-grained)": "membrane CG",
    "virus capsid": "capsid",
}


def source(system: str) -> str:
    for prefix, name in SOURCES:
        if system.startswith(prefix):
            return name
    return SELF_BUILT


def _draw_row(ax, blocks, phases_legend=False):
    """One row of bars, in labelled blocks. Returns the x extent used."""
    GAP = 1.2
    positions, spans, cursor, ordered = [], [], 0.0, []
    for label, ps in blocks:
        first = cursor
        for point in ps:
            positions.append(cursor)
            ordered.append(point)
            cursor += 1.0
        spans.append((label, first, cursor - 1.0, len(ps)))
        cursor += GAP

    x = np.array(positions)
    bottom = np.zeros(len(ordered))
    for phase, label in PHASES:
        pct = np.array([p["shares"][phase] for p in ordered])
        ax.bar(x, pct, 0.82, bottom=bottom, label=label, color=PHASE_COLORS[phase])
        bottom += pct

    for row, (label, first, last, n) in enumerate(spans):
        rule, text = (102.5, 104.0) if row % 2 == 0 else (108.0, 109.5)
        ax.plot(
            [first - 0.45, last + 0.45],
            [rule, rule],
            color=INK2,
            linewidth=1.0,
            clip_on=False,
        )
        ax.text(
            (first + last) / 2,
            text,
            f"{SHORT.get(label, label)} ($n$={n})",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK,
            clip_on=False,
        )

    ax.set_ylim(0, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("Share of write time (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([p["system"] for p in ordered], rotation=90, fontsize=7)
    ax.tick_params(axis="x", length=0)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return cursor - GAP


def figure(path_stem: Path, points):
    """Two rows of bars of one width, the solvated proteins apart from the rest.

    Small multiples were tried first and abandoned: matplotlib spreads a panel's
    bars over its full width, so a family of one drew a bar thirty-seven times
    wider than a family of thirty-seven, and the eye read the width as meaning.
    Both rows here share one scale of data units per inch, so every bar in the
    figure is the same width whichever row it sits in.
    """
    solvated = [p for p in points if p["family"] == "protein in solvent"]
    rest = [p for p in points if p["family"] != "protein in solvent"]

    # Row 1: by the collection the system came from, each block ordered by
    # repetition, which is the gradient the shares follow.
    by_source = collections.defaultdict(list)
    for p in solvated:
        by_source[source(p["system"])].append(p)
    for ps in by_source.values():
        ps.sort(key=lambda p: p["repetition"])
    row1 = [
        (name, by_source[name])
        for name in sorted(by_source, key=lambda k: -len(by_source[k]))
    ]

    # Row 2: the other families, in the order the rest of the paper uses.
    by_family = collections.defaultdict(list)
    for p in rest:
        by_family[p["family"]].append(p)
    for ps in by_family.values():
        ps.sort(key=lambda p: p["repetition"])
    row2 = [(f, by_family[f]) for f in LEGEND_ORDER if f in by_family]

    slots = max(len(solvated) + 1.2 * len(row1), len(rest) + 1.2 * len(row2))
    fig, axes = plt.subplots(2, 1, figsize=(0.30 * slots + 4.5, 11.0))
    # The rows need a wide gap: the upper one carries rotated system names below
    # it and the lower one carries its block headers above, and the two met.
    fig.subplots_adjust(hspace=0.62)

    _draw_row(axes[0], row1)
    _draw_row(axes[1], row2)
    # The shorter row keeps the longer row's scale instead of stretching to fill,
    # which is what would make its bars wider than those above them.
    for ax in axes:
        ax.set_xlim(-0.8, slots - 0.2)

    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.005, 0.98),
        borderaxespad=0,
        fontsize=11,
        labelcolor=INK2,
    )

    # "tight" trims the outer margins, which is what keeps the legend on the
    # canvas; it does not disturb the hspace set above.
    fig.savefig(f"{path_stem}.pdf", dpi=350, bbox_inches="tight")
    fig.savefig(f"{path_stem}.png", dpi=350, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--timing",
        default=str(
            Path(__file__).resolve().parents[1] / "data" / "decompose_run*.json"
        ),
        help="glob of the benchmark JSON files",
    )
    ap.add_argument(
        "--out", default=str(Path(RESULTS_DIR) / "figures" / "timing_by_family")
    )
    args = ap.parse_args()

    if not glob.glob(args.timing):
        print(f"  no benchmark output at {args.timing}")
        return 1
    points = load(args.timing)
    if not points:
        print(f"  no usable rows in {args.timing}")
        return 1

    print(f"  {len(points)} systems in {len({p['family'] for p in points})} families")
    for phase, label in PHASES:
        vals = [p["shares"][phase] for p in points]
        print(
            f"    {label:<20} median {st.median(vals):5.1f}%   "
            f"range {min(vals):5.1f}-{max(vals):5.1f}%"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    figure(out, points)
    print(f"  wrote {out}.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
