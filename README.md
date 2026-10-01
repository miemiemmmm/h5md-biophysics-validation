# Technical validation of the H5MD `biophysics` metadata schema

This directory holds everything needed to reproduce the Technical Validation of the H5MD `biophysics` metadata schema. 
The specification (`si-specification.tex` in the supplementary material) is the normative source of the validation suite.
The schema model in `h5md_validator/spec.py` is transcribed from the specification rather than the GROMACS reference implementation.
By generating the H5MD outputs from the GROMACS reference implementation, the validation suite can compare them against the schema model to identify any discrepancies.
Where the two disagree, the disagreement is reported as a finding rather than absorbed into the model; the list is in `$H5MD_VALIDATION_DATA/_results/specification_findings.md`.

## The six validations

| | Question | Artifact |
|:--|:--|:--|
| 1 | Can the schema represent the metadata of realistic MD systems? | `coverage.json`, `interaction_coverage.json`, `parameter_surface.json` |
| 2 | Does serialising and parsing preserve it? | `validation.json` → `roundtrip` |
| 3 | Is a file internally consistent? | `validation.json` → `validation` |
| 4 | Are invalid files detected? | `negative_tests.json` |
| 5 | What does the metadata cost? | `validation.json` → `storage`, `write_overhead.json` |
| 6 | Do other engines' parameters map onto the same concepts? | `amber_interop.json` |

Fidelity and coverage are kept apart throughout. A file can carry back everything
it holds without change and still have come from a source the schema could not
fully represent; conflating the two would let a coverage gap read as a success.

## Before doing validation

### Setup

Requires `h5py` and `numpy`, and a GROMACS build with the `biophysics` H5MD
writer. Three variables are needed by every step below:

```bash
export PY_EXE=/path/to/python/env/bin/python   # needs h5py and numpy
export GMX_BIN=/path/to/gromacs/bin/gmx        # a build with the biophysics writer
```


To run everything in order, including the report:

```bash
bash scripts/run_all.sh              # analysis only, reference files must exist
bash scripts/run_all.sh --generate   # build the reference files first (hours)
```

The rest of this file is that chain, one validation at a time.

### Download the validation dataset
Download the validation dataset and unpack it. Then setup the `H5MD_VALIDATION_DATA` environment variable to point at the extracted directory.

```bash
export H5MD_VALIDATION_DATA=/path/to/validation_systems
```

Main results and artifacts will be run put into:

```text
$H5MD_VALIDATION_DATA/
    _data/      Reference files, tprs and output H5MD files. Transient, several GB, safe to delete.
    _results/   Measurements, reports, the tables and the figures.
```


### Build the reference files

Every validation reads this corpus, so it is built once.

**Objective.** Produce one H5MD file per (system, protocol) pair, together with
the `.tpr` it came from. Both sides are kept: the `.tpr` is the source of the
round-trip comparison, so the test never compares the writer against itself.

```bash
$PY_EXE scripts/build_reference_files.py --list   # show the test matrix
$PY_EXE scripts/build_reference_files.py --all    # grompp + mdrun every pair
```

Existing `.tpr` and `.h5md` files are reused; `--force` rebuilds them.

**Artifacts**

| Path | What it is |
|:--|:--|
| `_results/reference_files.json` | the manifest: one entry per pair, with its dimensions and protocol targets |
| `_data/tprs/` | grompp output, the round-trip source side |
| `_data/resolved_mdp/` | grompp `-po`: the complete `.mdp` with every default filled in |
| `_data/reference_files/` | the H5MD corpus |
| `_data/tpr_dumps/` | cached `gmx dump -s` text, written on first read |
| `_data/scratch/` | per-pair scratch, deleted after each pair |

A clean run reports **379/379 reference files available**. A leftover directory
under `_data/scratch/` names a pair whose run was interrupted.

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
$PY_EXE scripts/validation1_parameter_surface.py "$H5MD_VALIDATION_DATA"/_data/resolved_mdp/*.mdp
```

The first reads `validation.json`, so run Validation 2/3 first. The second reads
the `.tpr` dumps directly and needs no other stage.

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/coverage.json` | `element_coverage` (specified / observed / never observed / constant), `item_fidelity`, `dimensions`, `protocol_targets` |
| `_results/interaction_coverage.json` | `per_system`, `per_file`, `unsupported` |
| `_results/parameter_surface.json` | per-option classification of the resolved `.mdp` files |

**A clean run.** 123 of 145 specified elements observed, 103 of them with more
than one distinct value; 92 of 93 source items represented; **0** interaction
types unsupported and **0** entries dropped.

Two numbers deserve reading rather than quoting. The 22 never-observed elements
are each explained in the report by a registered deviation or limitation, and
any element with no recorded explanation is flagged there by name. The one
unrepresented item is `cmap.parameters`, which cannot be verified at all:
`gmx dump` emits the CMAP type index and never the grid values, so there is
nothing on the source side to compare against.

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

**A clean run.** 379 files compared, **379 faithful**, 0 mismatches, 0 items
missing from a file. Numerical comparison uses a per-value tolerance bounded by
whichever is coarser — `gmx dump`'s text precision or the file's float width —
rather than one global epsilon.

---

## Validation 3 — internal consistency

**Objective.** Show that each file is conformant on its own terms: declared
types and shapes, enumerated values, unit attributes, time alignment, and the
referential integrity of every index it contains.

Same command as Validation 2; the two are produced together.

```bash
$PY_EXE scripts/check_reference_files.py
```

**Artifacts.** `_results/validation.json` → the `validation` section of each
entry: a list of findings with a stable code, a severity and the HDF5 path it is
anchored to, plus a `conformant` verdict.

**A clean run.** **379/379 conformant, 0 ERROR-level findings.**

Deviations are counted separately and do not make a file non-conformant. Each
one is a registered disagreement between the specification text and the writer,
awaiting an editorial decision rather than describing a bad file; they are
listed in `specification_findings.md`.

---

## Validation 4 — negative tests

**Objective.** Show the validator has teeth. Proving that valid files pass is
only half of a schema claim; the other half is that invalid metadata is caught.
Each corruption injects exactly one defect into a copy of a valid file and
declares the finding code the validator is required to raise.

```bash
$PY_EXE scripts/validation4_negative_tests.py <a valid .h5md> [more .h5md ...]
```

`run_all.sh` passes the bases declared in `config/matrix.py`, chosen for
different shapes so that as many corruptions as possible apply.

**Artifacts.** `_results/negative_tests.json` — `by_category`, the per-file
outcomes, and the counts below.

**A clean run.** 51 corruption kinds over 5 base files, **249 tested, 249
detected, 0 missed**, detection rate 100%. A corruption that cannot apply to a
particular base (no barostat group, no selections) is skipped and reported as
such rather than counted as a pass.

---

## Validation 5 — what the metadata costs

**Objective.** Show that carrying the metadata is cheap enough to be
unobjectionable, in bytes and in write time. Both halves matter: a schema
nobody adopts because it is expensive has not solved the problem.

Storage comes from the same pass as Validations 2 and 3. Write time is a
separate measurement, because it means running each system twice — once writing
H5MD with metadata and once writing a trajectory only, from the same `.tpr`.

```bash
$PY_EXE scripts/check_reference_files.py                       # storage
$PY_EXE scripts/validation5_write_overhead.py --repeats 5      # time
```

**Artifacts**

| Path | Key contents |
|:--|:--|
| `_results/validation.json` → `storage` | per file: metadata bytes, trajectory bytes, bytes per particle, ratio |
| `_results/write_overhead.json` | per system: both wall times and their difference, over `--repeats` runs |

**A clean run.** Metadata costs **11–363 bytes per particle**, and writing it
adds a **median of +8 ms, at most +152 ms** over trajectory-only output.

---

## Validation 6 — cross-engine interoperability

**Objective.** Show that the schema's concepts are not GROMACS concepts under
another name. The same molecule is prepared independently in two engines with
different force fields, and the chemistry each one records is compared through
the schema.

```bash
$PY_EXE scripts/validation6_amber_interop.py --h5md <file.h5md> --prmtop <file.prmtop>
```

Both sides must describe the **same molecule** — a mismatched pair compares two
different systems and its verdict means nothing. `run_all.sh` takes the pairs
from `AMBER_INTEROP_PAIRS` in `config/matrix.py`, where each H5MD file is
declared alongside the AMBER topology it corresponds to, so the invariant is
stated in one place rather than left to the person passing the flags.

**Artifacts.** `_results/amber_interop.json` — per-item rows under `TOPOLOGY`
and `PARAMETRISATION`, the alias table used to align atom names, and an
`interoperable` verdict.

**A clean run.** **9 of 9 chemistry items agree**: particle and residue counts,
molecular formula, bond count and bond topology, total charge, masses and
atomic numbers. Naming differences are resolved by an alias table and reported
rather than hidden — `O1P`/`OP1`, `H3T`/`HO3*`, and terminal residues named
`DA5`/`DA3` against plain `DA`.

---

## Supporting analysis

Not one of the six, but run by `run_all.sh` and cited by the report:

```bash
$PY_EXE scripts/provenance_recoverability.py
```

Measures whether force-field provenance could be recovered automatically from a
file's contents, which is the evidence behind `DEV-FF-PROVENANCE-UNWRITTEN`.
Writes `_results/provenance_recoverability.json`.

## Report and figures

```bash
$PY_EXE scripts/make_figures.py
$PY_EXE scripts/report.py
```

| Path | What it is |
|:--|:--|
| `_results/technical_validation.md` | the assembled report; every number comes from a result file |
| `_results/specification_findings.md` | the deviation registry, rendered |
| `_results/tables/*.tex` | LaTeX tables for the paper |
| `_results/figures/` | `validation_framework` and `metadata_scaling`, as PDF and PNG |

Both are derived from the JSON beside them. Re-running the validation on
different hardware, or with a different GROMACS build, means regenerating the
report rather than assuming the shipped numbers still describe that run.

`RESULTS_DIR` redirects every stage at once, which is how to report on a set of
results kept elsewhere:

```bash
RESULTS_DIR=/some/other/results $PY_EXE scripts/report.py
```

## Validating a single file

```python
from h5md_validator import validate
report = validate("trajectory.h5md")
print(report.summary())
for finding in report.errors:
    print(finding)
```

`report.conformant` is true when no ERROR-level finding was raised. Findings
come at four levels: `ERROR` (violates a MUST of the specification), `DEVIATION`
(a registered specification-versus-implementation mismatch, pending an editorial
decision), `WARNING` (violates a SHOULD), and `INFO`.

## Layout

```
h5md_validator/      the library: schema model and the checks, no knowledge of the dataset
    spec.py          schema model, transcribed from the specification
    validate.py      conformance and internal-consistency checks
    findings.py      finding and report types
    units.py         H5MD unit-notation grammar
    gmxdump.py       parser for `gmx dump -s`, the round-trip source side
    vocabulary.py    GROMACS <-> schema translation tables
    roundtrip.py     item-by-item comparison of source against file
    corrupt.py       the corruption catalogue for the negative tests
    deviations.py    registry of specification-versus-implementation mismatches
config/              what this test set is, and where it lives
    matrix.py        protocols, dataset paths, negative-test bases, interop pairs
    systems.py       the 58 chemical systems of the test matrix
records/             authored claims the report quotes; nothing here is measured
    limitations.py   registry of configurations the test set could not produce
scripts/             the stages, each writing one JSON artifact
    build_reference_files.py          grompp + mdrun every (system, protocol) pair
    check_reference_files.py          Validations 2, 3 and 5, one pass per file
    validation1_schema_coverage.py    element and item coverage over the corpus
    validation1_interaction_coverage.py  bonded interactions: source against file
    validation1_parameter_surface.py  how much of a flat .mdp is actually in force
    validation4_negative_tests.py     corruption tests
    validation5_write_overhead.py     H5MD against trajectory-only output, same .tpr
    validation6_amber_interop.py      the same molecule in GROMACS and in AMBER
    provenance_recoverability.py      where force-field provenance could come from
    make_figures.py                   the two figures
    report.py                         assemble reports and LaTeX tables
    run_all.sh                        chain all of the above
$H5MD_VALIDATION_DATA/
    _data/           transient: tprs, resolved_mdp, reference_files, tpr_dumps, scratch
    _results/        the JSON measurements, the rendered report, tables and figures
```

The prepared systems live outside this directory, in `../validation_systems/`:
one directory per system, holding the grompp inputs the validation was computed
from. They are kept separate because they are large and are archived rather
than versioned.

## Two registries

Findings are kept in two separate registries rather than one list, because they
need different actions from different people:

| Registry | What it holds | Who resolves it |
|:--|:--|:--|
| `h5md_validator/deviations.py` | the specification text and the implementation disagree | an editorial decision on the specification |
| `records/limitations.py` | the test set could not produce a configuration at all | more test systems, or a stated caveat |

Defects found in the reference implementation during validation were fixed in
the writer rather than tracked here; the corpus is regenerated after each fix,
so the current results describe the fixed implementation.

## Check codes

Findings carry a stable code so they can be cited and tracked.

| Prefix | Family |
|:--|:--|
| `MOD` | module declaration and version |
| `STRUCT` | presence, type and shape of specified elements |
| `ENUM` | values of enumerated elements |
| `UNIT`, `UNITS` | unit attributes and their notation |
| `COND` | conditional subgroup instantiation |
| `TIME` | step/time/value alignment |
| `TOPO` | topology invariants |
| `REF` | referential integrity of indices |
| `SEL`, `CONN` | selection and connectivity bounds |
| `FILE` | the file cannot be read at all |
| `DEV-*` | registered specification deviation |
