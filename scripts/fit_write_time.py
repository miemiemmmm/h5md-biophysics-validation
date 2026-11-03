#!/usr/bin/env python
"""
Correlation and model fitting for the metadata writing time.

    python scripts/fit_write_time.py
    python scripts/fit_write_time.py --json fit.json
    python scripts/fit_write_time.py --timing 'other/decompose_run*.json'

T is the total time to write the metadata of one system: the seven phases timed
by the BenchmarkDecomposeWritingProcess gtest in h5md_benchmark.cpp, each taken
as its median over the repeats and then summed. It is set against four
structural counts of the same system:

    N_ua     unique atoms: atoms summed over molecule types, each type once
    N_b      bonds in the whole system (the global connectivity list)
    N_um     molecule types
    natoms   particles in the system

The counts are read from the reference H5MD file written from the same .tpr,
so the script needs the validation dataset but not GROMACS. Taken from the
.tpr instead, all four agree for all 58 systems.

What is printed, and why each part is there:

* Pearson r with and without the excluded system (STMV by default). STMV has
  32x the unique atoms of any other system and ~97% of the variance in both T
  and N_ua, so it alone drives r towards 1; the value without it describes the
  ordinary systems. Log-log Pearson r and Spearman's rho are given alongside
  because no single point can dominate either.
* Linear models compared by leave-one-out error, not only R^2: R^2 never
  decreases when a term is added, while leave-one-out error shows whether the
  term predicts systems the fit did not see.
* Per-phase correlations, which tie each model term to the part of the writer
  that produces it: topology to N_ua, connectivity to N_b.
* The extrapolation check: what a fit to the ordinary systems predicts for the
  excluded one, so a reader can see how far the model can be trusted outside
  the range it was fitted on.

Not part of the validation suite: it reads a recorded benchmark, which the
suite does not produce.
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
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config.matrix import DATA_ROOT

#: Property name in the gtest JSON -> phase name used in the output.
PHASES = {
    "metadata": "h5md_metadata",
    "particle_block": "trajectory_scaffold",
    "topology": "topology",
    "bond_connectivity": "connectivity",
    "force_field": "force_field",
    "md_parameters": "parameters",
    "selection_groups": "selections",
}

PREDICTORS = ("N_ua", "N_b", "N_um", "natoms")

#: Models compared, as lists of predictors; an intercept is always included.
MODELS = (["N_ua"], ["N_ua", "N_b"], ["N_um", "N_ua", "N_b"])

#: Which phase each scaling term is expected to drive.
MECHANISM = (
    ("topology", "N_ua"),
    ("connectivity", "N_b"),
    ("connectivity", "N_ua"),
    ("parameters", "N_ua"),
)


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def load_timing(pattern: str):
    """{tpr name: {phase: [ms per repeat]}} from the gtest JSON files."""
    runs = collections.defaultdict(lambda: collections.defaultdict(list))
    files = sorted(glob.glob(pattern))
    for name in files:
        for suite in json.loads(Path(name).read_text())["testsuites"]:
            for test in suite["testsuite"]:
                if "metadata_elapsed_us" not in test:
                    continue
                tpr = test["value_param"].strip('"')
                for key, phase in PHASES.items():
                    runs[tpr][phase].append(int(test[f"{key}_elapsed_us"]) / 1000.0)
    return runs, len(files)


def structural_counts(h5md_path: Path) -> dict:
    """N_ua, N_b, N_um and natoms from one reference H5MD file."""
    with h5py.File(h5md_path, "r") as f:
        topology = f["topology"]
        molecule_types = [m for m in topology if isinstance(topology[m], h5py.Group)]
        return {
            "N_ua": sum(
                int(np.atleast_1d(topology[m].attrs["particle_count"])[0])
                for m in molecule_types
            ),
            "N_b": f["connectivity/bonds"].shape[0] if "connectivity/bonds" in f else 0,
            "N_um": len(molecule_types),
            "natoms": int(np.atleast_1d(topology.attrs["total_particle_count"])[0]),
        }


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def pearson(x, y) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def fit(columns, y):
    """Ordinary least squares with an intercept.

    Returns coefficients, R^2, adjusted R^2, leave-one-out RMSE and the median
    leave-one-out relative error. Leave-one-out residuals come from the hat
    matrix, e_i / (1 - h_ii), which equals refitting without each point.
    """
    A = np.column_stack([np.ones(len(y)), *columns])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    n, p = A.shape
    ss_tot = float((y - y.mean()) @ (y - y.mean()))
    r2 = 1.0 - float(resid @ resid) / ss_tot
    adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p)
    loo = resid / (1.0 - np.diag(A @ np.linalg.pinv(A)))
    return {
        "coefficients": coef.tolist(),
        "r2": r2,
        "adj_r2": adj,
        "loo_rmse_ms": float(np.sqrt(np.mean(loo**2))),
        "loo_median_rel_err": float(np.median(np.abs(loo) / y)),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--timing",
        default=str(
            Path(__file__).resolve().parents[1] / "data" / "decompose_run*.json"
        ),
        help="gtest JSON files, one per repeat (a glob)",
    )
    ap.add_argument(
        "--reference-dir",
        default=str(Path(DATA_ROOT) / "reference_files"),
        help="reference H5MD files, named <tpr stem>.h5md",
    )
    ap.add_argument(
        "--exclude",
        default="bench_stmv",
        help="system reported separately as a high-leverage point",
    )
    ap.add_argument("--json", help="also write every number printed here")
    args = ap.parse_args()

    runs, n_files = load_timing(args.timing)
    if not runs:
        print(f"  no benchmark output matches {args.timing}")
        return 1

    reference = Path(args.reference_dir)
    keys, missing = [], []
    for tpr in sorted(runs):
        (
            keys if (reference / tpr.replace(".tpr", ".h5md")).is_file() else missing
        ).append(tpr)
    if missing:
        print(
            f"  note: {len(missing)} timed system(s) have no reference file and are "
            f"left out: {', '.join(missing)}"
        )

    counts = {
        k: structural_counts(reference / k.replace(".tpr", ".h5md")) for k in keys
    }
    system = [k.split("__")[0] for k in keys]
    medians = {k: {p: st.median(v) for p, v in runs[k].items()} for k in keys}
    T = np.array([sum(medians[k].values()) for k in keys])
    X = {name: np.array([counts[k][name] for k in keys], float) for name in PREDICTORS}
    keep = np.array([s != args.exclude for s in system])
    if keep.all():
        print(
            f"  note: --exclude {args.exclude!r} matches no system; "
            "the 'without' columns repeat the full set"
        )

    n_repeats = collections.Counter(len(runs[k]["topology"]) for k in keys).most_common(
        1
    )[0][0]
    totals = [
        [
            sum(runs[k][p][r] for p in PHASES.values())
            for r in range(len(runs[k]["topology"]))
        ]
        for k in keys
    ]
    cv = np.median([np.std(t) / np.mean(t) for t in totals])
    out = {
        "n_systems": len(keys),
        "n_files": n_files,
        "n_repeats": n_repeats,
        "excluded": args.exclude,
        "T_ms": {
            "median": float(np.median(T)),
            "min": float(T.min()),
            "max": float(T.max()),
            "run_to_run_cv": float(cv),
        },
    }

    print(f"\n  {len(keys)} systems, {n_repeats} repeats ({n_files} files)")
    print(
        f"  T = sum of the seven phase medians: median {np.median(T):.2f} ms, "
        f"range {T.min():.2f}-{T.max():.2f} ms, run-to-run CV {100 * cv:.1f}%"
    )
    if (~keep).any():
        print(
            f"  {args.exclude}: T = {T[~keep][0]:.1f} ms, "
            f"N_ua = {X['N_ua'][~keep][0]:,.0f}"
        )

    # ---- correlations ------------------------------------------------------
    print("\n  correlation of T with each count")
    print(
        f"  {'':8s} {'Pearson r':>10s} {'r w/o ' + args.exclude[:10]:>17s} "
        f"{'r log-log':>10s} {'Spearman':>9s}"
    )
    out["correlation"] = {}
    for name in PREDICTORS:
        x = X[name]
        row = {
            "pearson": pearson(x, T),
            "pearson_without_excluded": pearson(x[keep], T[keep]),
            "pearson_log_log": pearson(np.log10(x), np.log10(T)),
            "spearman": float(stats.spearmanr(x, T)[0]),
        }
        out["correlation"][name] = row
        print(
            f"  {name:8s} {row['pearson']:10.3f} {row['pearson_without_excluded']:17.3f} "
            f"{row['pearson_log_log']:10.3f} {row['spearman']:9.3f}"
        )

    # ---- models --------------------------------------------------------------
    out["models"] = {}
    for label, mask in (("all", np.ones(len(T), bool)), ("without_excluded", keep)):
        print(
            f"\n  models, {'all systems' if label == 'all' else 'without ' + args.exclude}"
        )
        print(
            f"  {'T = a + ...':22s} {'R2':>6s} {'adj R2':>7s} {'LOO RMSE':>9s} "
            f"{'LOO err':>8s}   coefficients"
        )
        for cols in MODELS:
            name = " + ".join(cols)
            res = fit([X[c][mask] for c in cols], T[mask])
            out["models"].setdefault(label, {})[name] = res
            coef = "  ".join(f"{c:.3g}" for c in res["coefficients"])
            print(
                f"  {name:22s} {res['r2']:6.3f} {res['adj_r2']:7.3f} "
                f"{res['loo_rmse_ms']:7.2f} ms {100 * res['loo_median_rel_err']:7.1f}%   [{coef}]"
            )

    # ---- which phase each term drives ------------------------------------------
    print(f"\n  per-phase Pearson r (without {args.exclude})")
    out["mechanism"] = {}
    for phase, name in MECHANISM:
        y = np.array([medians[k][phase] for k in keys])
        r = pearson(X[name][keep], y[keep])
        out["mechanism"][f"{phase}~{name}"] = r
        print(f"  {phase:14s} vs {name:6s}  r = {r:.3f}")

    # ---- extrapolation -----------------------------------------------------
    if (~keep).any():
        res = fit([X[c][keep] for c in ("N_ua", "N_b")], T[keep])
        a, b, c = res["coefficients"]
        predicted = a + b * X["N_ua"][~keep][0] + c * X["N_b"][~keep][0]
        out["extrapolation"] = {
            "predicted_ms": float(predicted),
            "measured_ms": float(T[~keep][0]),
        }
        print(
            f"\n  a + N_ua + N_b fitted without {args.exclude} predicts it at "
            f"{predicted:.0f} ms; measured {T[~keep][0]:.0f} ms"
        )

    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f"\n  wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
