# Validation steps

Each step of `scripts/run_all.sh`, one validation at a time. Setup and the quick
start are in the [README](../README.md). Paths below are relative to the output
root, by default `$H5MD_VALIDATION_DATA`; `DATA_DIR` and `RESULTS_DIR` move
`_data` and `_results`.

## Build the reference files

Every validation reads this corpus, so reference files are built only once. 
This produces one H5MD file per (system, protocol) pair, together with
the `.tpr` it came from. Both files are kept: the `.tpr` is the source of the
round-trip comparison.

```bash
$PY_EXE scripts/build_reference_files.py --list   # show the test matrix
$PY_EXE scripts/build_reference_files.py --all    # grompp + mdrun every pair
```

Existing `.tpr` and `.h5md` files are reused; add `--force` flag to rebuild them.

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

Three stages, in increasing sharpness. The first asks which specified elements
appear anywhere; the second counts individual bonded interaction entries in the
run input against the file; the third asks how much of a flat GROMACS parameter
file is actually in force.

```bash
$PY_EXE scripts/validation1_schema_coverage.py
$PY_EXE scripts/validation1_interaction_coverage.py
$PY_EXE scripts/validation1_parameter_surface.py     # every .mdp in _data/resolved_mdp
```

The first reads `validation.json`, so run Validation 2/3 first. The second reads
the `.tpr` dumps directly and needs no other stage.

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/coverage.json` | `element_coverage` (specified / observed / never observed / constant), `item_fidelity`, `dimensions`, `protocol_targets` |
| `_results/interaction_coverage.json` | `per_system`, `per_file`, `unsupported` |
| `_results/parameter_surface.json` | per-option classification of the resolved `.mdp` files |

Two results need reading rather than quoting. Every never-observed element is
listed by name under `element_coverage.never_observed`, and each needs an
explanation (a registered deviation or a gap in the test set) before the
coverage figure is quoted. Separately,
`cmap.parameters` cannot be verified at all: `gmx dump` emits the CMAP type
index and never the grid values, so there is nothing on the source side to
compare against.

---

## Validation 2 — round-trip fidelity

**Objective.** Show that metadata written to H5MD and read back is unchanged.
The comparison runs against `gmx dump` output rather than the writer's own view,
so the two sides of the test share no code path.

```bash
$PY_EXE scripts/check_reference_files.py
```

One pass per file does Validations 2, 3 and 5 together, because all three need
the file opened and the `.tpr` dumped. Dumps are cached after the first run.

**Artifacts.** `_results/validation.json`, keyed by reference file. The
`roundtrip` section of each entry holds every compared item with a status of
`match`, `mismatch`, `missing_in_file` or `missing_in_source`, plus the
per-category totals and a `faithful` verdict.

Numerical comparison uses a per-value tolerance bounded by whichever is coarser
— `gmx dump`'s text precision or the file's float width — rather than one global
epsilon.

---

## Validation 3 — internal consistency

**Objective.** Show that each file is conformant on its own terms: declared
types and shapes, enumerated values, unit attributes, time alignment, and the
referential integrity of every index it contains.

Same command as Validation 2; the two are produced together.

```bash
$PY_EXE scripts/check_reference_files.py
```

**Artifacts.** `_results/validation.json` → the `conformance` section of each
entry: a list of findings with a stable code, a severity and the HDF5 path it is
anchored to, plus a `conformant` verdict.

Deviations are counted separately and do not make a file non-conformant. Each
one is a registered disagreement between the specification text and the writer,
awaiting an editorial decision rather than describing a bad file; they are
listed in `h5md_validator/deviations.py`.

---

## Validation 4 — negative tests

**Objective.** Show the validator has teeth. Proving that valid files pass is
only half of a schema claim; the other half is that invalid metadata is caught.
Each corruption injects exactly one defect into a copy of a valid file and
declares the finding code the validator is required to raise.

```bash
$PY_EXE  scripts/validation4_negative_tests.py FILE.h5md    # explicitly define a H5MD file
$PY_EXE  scripts/validation4_negative_tests.py --quick      # go through five reference files in the config/matrix.py
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

## Validation 5 — what the metadata costs to store

**Objective.** Show that carrying the metadata is cheap enough to be
unobjectionable. A schema nobody adopts because it is expensive has not solved
the problem.

Measured at the HDF5 layer: the stored size of every dataset, attributed to the
metadata groups or to the trajectory. It comes from the same pass as
Validations 2 and 3.

```bash
$PY_EXE scripts/check_reference_files.py
```

Write *time* is measured inside the GROMACS writer instead: whole-process wall
time and GROMACS's `Write traj.` counter both mix the one-off metadata write with
start-up and per-frame output. A GROMACS benchmark times each phase of the
metadata write for every system; its output ships in `data/` (see
[data/README.md](../data/README.md)) and is panel C of the figure.

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/validation.json` → `storage` | per file: metadata bytes, trajectory bytes, bytes per particle, ratio |

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
molecular formula, atom names, atomic numbers, bond count and bond topology; the
parametrisation items are total charge, partial charges and masses. The H5MD
side is rebuilt to the full system from `topology/molecule_blocks`. The H–H
"bond" AMBER gives rigid water keeps it rigid under SHAKE and is not chemistry,
so it is removed from the AMBER side and the number removed is reported.

---

## Supporting analysis

Not one of the six, but run by `run_all.sh`:

```bash
$PY_EXE scripts/provenance_recoverability.py
```

Measures whether force-field provenance could be recovered automatically from a
file's contents, which is the evidence behind `DEV-FF-PROVENANCE-UNWRITTEN`.
Writes `_results/provenance_recoverability.json`.

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
