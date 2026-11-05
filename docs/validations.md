# Validation steps

Each step of `scripts/run_all.sh`, one validation at a time. Setup and the quick
start are in the [README](../README.md). Paths below are relative to the output
root, by default `$H5MD_VALIDATION_DATA`; `DATA_DIR` and `RESULTS_DIR` move
`_data` and `_results`.

## Build the reference files

Every validation reads this corpus, so reference files are built only once.
Each validation then runs on its own, in any order: none reads another's
results. 
This produces one H5MD file per (system, protocol) pair, together with
the `.tpr` it came from. Both files are kept: the `.tpr` is the source of the
round-trip comparison.

```bash
$PY_EXE scripts/build_reference_files.py --list   # show the test matrix
$PY_EXE scripts/build_reference_files.py --all    # grompp + mdrun every pair
```

Existing `.tpr` and `.h5md` files are reused; add `--force` flag to rebuild them.
`--ntomp N` sets the OpenMP threads of each `mdrun` (default 8).

To build only selected pairs, name the systems and protocols; both options can
be repeated, and each listed system is built under every listed protocol that
applies to it. For example, the four example files published with the dataset:

```bash
$PY_EXE scripts/build_reference_files.py --system hexabc_a001u2_rna_protein --protocol nvt-vrescale
$PY_EXE scripts/build_reference_files.py --system bench_villin \
    --protocol tcoupl-multigroup --protocol sd-langevin --protocol selection-noncontiguous
```

With only some pairs built, `run_all.sh` stops at Validation 4 when none of the
negative-test bases declared in `config/matrix.py` is among them; run the
validations one by one instead, and Validation 4 with `--all`.

**Result artifacts**

| Path | What it is |
|:--|:--|
| `_results/reference_files.json` | the manifest: one entry per pair, with its dimensions and protocol targets |
| `_data/tprs/` | grompp output, the round-trip source side |
| `_data/reference_files/` | the H5MD corpus |
| `_data/resolved_mdp/` | grompp `-po`: the complete `.mdp` with every default filled in |
| `_data/tpr_dumps/` | cached `gmx dump -s` text, written on first read |
| `_data/scratch/` | per-pair scratch, deleted after each pair |

A leftover directory under `_data/scratch/` names a pair whose run was
interrupted.

---

## Validation 1 — schema coverage

**Objective.** Establish that the schema can hold the metadata real simulations
produce, and say exactly which parts of it no run exercised. An element that no
file ever wrote is an untested branch, not a working one, and the difference
between those two is what separates a coverage claim from an assumption.

It compares the specification model with the reference files alone. The first
script asks which specified elements the files write, and which chemical and
protocol dimensions they cover; the second asks how much of a flat GROMACS
parameter file is actually in force.

```bash
$PY_EXE scripts/validation1_schema_coverage.py
$PY_EXE scripts/validation1_parameter_surface.py     # every .mdp in _data/resolved_mdp
```

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/coverage.json` | `element_coverage` (specified / observed / never observed / constant), `dimensions`, `protocol_targets` |
| `_results/parameter_surface.json` | per-option classification of the resolved `.mdp` files |

Every never-observed element is listed by name under
`element_coverage.never_observed`, and each needs an explanation (not written
by the reference implementation, or not exercised by the test set) before the
coverage figure is quoted.

---

## Validation 2 — round-trip fidelity

**Objective.** Show that metadata written to H5MD and read back is unchanged.
The comparison runs against `gmx dump` output rather than the writer's own view,
so the two sides of the test share no code path.

The first script compares every metadata item; the second counts the bonded
interaction entries of each type in the run input against the file. Dumps are
cached after the first run.

```bash
$PY_EXE scripts/validation2_roundtrip.py
$PY_EXE scripts/validation2_interaction_coverage.py
```

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/roundtrip.json` | `files`: per file, every compared item with a status of `match`, `mismatch`, `missing_in_file` or `missing_in_source`, and a `faithful` verdict; `item_fidelity`: the same items aggregated over all files (represented / missing / incorrect) |
| `_results/interaction_coverage.json` | `per_system`, `per_file`, `unsupported` |

Numerical comparison uses a per-value tolerance bounded by whichever is coarser
— `gmx dump`'s text precision or the file's float width — rather than one global
epsilon. The `cmap` grids cannot be compared: `gmx dump` prints the CMAP type
index and never the grid values, so there is nothing on the source side.

---

## Validation 3 — internal consistency

**Objective.** Show that each file is conformant on its own terms: declared
types and shapes, enumerated values, unit attributes, time alignment, and the
referential integrity of every index it contains.

Each file is checked against the specification model alone.

```bash
$PY_EXE scripts/validation3_conformance.py
```

**Artifacts.** `_results/conformance.json`: under `files`, one entry per
reference file with a list of findings, each with a stable code, a severity and
the HDF5 path it is anchored to, plus a `conformant` verdict. The codes are
listed in [check-codes.md](check-codes.md).

---

## Validation 4 — negative tests

**Objective.** Show the validator has teeth. Proving that valid files pass is
only half of a schema claim; the other half is that invalid metadata is caught.
Each corruption injects exactly one defect into a copy of a valid file and
declares the finding code the validator is required to raise.

```bash
$PY_EXE  scripts/validation4_negative_tests.py FILE.h5md    # explicitly define a H5MD file
$PY_EXE  scripts/validation4_negative_tests.py --quick      # go through the seven reference files in config/matrix.py
$PY_EXE  scripts/validation4_negative_tests.py --all        # go through all reference files, takes a while
```

`run_all.sh` passes the bases declared in `config/matrix.py`, chosen for
different shapes so that as many corruptions as possible apply.

**Artifacts.** `_results/negative_tests.json` — `by_category`, the per-file
outcomes, and the detection counts.

The full list of corruptions, with the finding each must trigger, is in
[corruptions.md](corruptions.md).

A corruption that cannot apply to a particular base (no barostat group, no
selections) is skipped and reported as such rather than counted as a pass.

---

## Validation 5 — what the metadata costs to store and to write

**Objective.** Show that carrying the metadata is cheap enough to be
unobjectionable. A schema nobody adopts because it is expensive has not solved
the problem.

Storage is measured at the HDF5 layer: the stored size of every dataset of each
reference file, attributed to its root group, and so to the metadata or to the
trajectory.

```bash
$PY_EXE scripts/validation5_storage.py
```

Write *time* is measured inside the GROMACS writer instead: whole-process wall
time and GROMACS's `Write traj.` counter both mix the one-off metadata write with
start-up and per-frame output. A GROMACS benchmark times each phase of the
metadata write (H5MD metadata, trajectory scaffold, topology, bond connectivity,
force field, MD parameters, selection groups) for the `nvt-vrescale` run input
of every system, ten times; its output ships in `data/` (see
[data/README.md](../data/README.md)). The cost of writing the metadata is read
as each phase's share of that time, median over the repeats.

`scripts/validation5_cost_figure.py` (see [Figures](#figures)) plots both:
metadata stored against system size (panel A), metadata per unique atom without
connectivity (panel B) and each phase's share of the write time (panel C).

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/storage.json` | per file: stored bytes per root group, metadata and trajectory bytes, file size, particle and frame counts, and two derived ratios (metadata bytes per particle; metadata over trajectory bytes, which reflects the few frames of a reference run rather than a production trajectory) |
| `data/decompose_run*.json` | shipped benchmark output: per system and repeat, the elapsed time of each phase of the metadata write |

---

## Validation 6 — cross-engine interoperability

**Objective.** Show that the schema's concepts are not GROMACS concepts under
another name. Every system that was set up in AMBER (the hexABC and MISATO
groups of the dataset) is compared item by item: its deposited AMBER topology
against the H5MD file GROMACS writes after ParmEd converted that topology.

```bash
$PY_EXE scripts/validation6_amber_interop.py --run-all     # every pair in AMBER_INTEROP_PAIRS
$PY_EXE scripts/validation6_amber_interop.py --h5md <file.h5md> --prmtop <file.prmtop>
```

Both sides must describe the **same molecule** — a mismatched pair compares two
different systems and its verdict means nothing. `--run-all` takes the pairs
from `AMBER_INTEROP_PAIRS` in `config/matrix.py`, where each H5MD file is
declared alongside the AMBER topology it corresponds to, so the invariant is
stated in one place rather than left to the person passing the flags.

**Artifacts.** `_results/amber_interop.json` — one entry per pair under
`pairs`, each with per-item rows under `chemistry` and `parametrisation` and an
`interoperable` verdict, plus `n_pairs` and `n_interoperable`.

The chemistry items are particle and residue counts, residue sequence,
molecular formula, atom names, atomic numbers, bond count and bond topology. The
parametrisation items are total charge, partial charges and masses, and every
force-field term: bonds, angles, proper and improper dihedrals (term by term),
1-4 pairs and the Lennard-Jones type matrix. The H5MD side is rebuilt to the
full system from `topology/molecule_blocks`. The H–H "bond" AMBER gives rigid
water keeps it rigid under SHAKE and is not chemistry, so it is removed from
the AMBER side of the bond-topology comparison and the number removed is
reported; among the force-field terms it is matched to the H–H distance of the
water's SETTLE.

For the force-field terms the prmtop is read directly, not through ParmEd,
which did the conversion under test. AMBER units and conventions are converted
before comparing (kcal → kJ, Å → nm, rad → degree; bond and angle force
constants × 2, since GROMACS writes E = k/2 (x − x₀)²; 1-4 Lennard-Jones divided
by SCNB and the Coulomb scale is 1/SCEE), at a relative tolerance of 1e-4. A
bond may be written as a bond, a constraint or a SETTLE distance. Three AMBER
devices are reported rather than counted as differences: dihedral terms with
zero force constant (no energy; `grompp` removes them from the run input), a
10-12 pair with zero coefficients (tleap's way of switching off water O–H
Lennard-Jones), and an improper written with its last two atoms swapped when its phase is 0° or 180°,
which leaves the energy unchanged.

---

## Figures

```bash
$PY_EXE scripts/validation5_cost_figure.py
```

| Path | What it is |
|:--|:--|
| `_results/figures/` | `validation5_cost.png`: metadata storage per system (panels A, B) and the write-time decomposition (panel C) |

The figures are derived from the JSON beside them. Re-running the validation
on different hardware, or with a different GROMACS build, means regenerating
them rather than assuming the shipped numbers still describe that run.

`RESULTS_DIR` (and `DATA_DIR` for `_data`) redirects every stage at once, which
is how to work on a set of results kept elsewhere:

```bash
RESULTS_DIR=/some/other/results $PY_EXE scripts/validation5_cost_figure.py
```
