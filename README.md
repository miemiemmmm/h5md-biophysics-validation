# Technical validation of the H5MD `biophysics` metadata schema

This directory holds everything needed to reproduce the Technical Validation
section: a validator, a round-trip test, a corruption test suite, and the
drivers that generate the reference files they run against.

The specification (`../H5MD_metadata_paper_15Sep26/si-specification.tex`) is the
normative source. 
The schema model in `h5md_validator/spec.py` is transcribed
from it, section by section, and nothing in it is derived from the GROMACS
reference implementation. 
Where the two disagree, the disagreement is reported
as a finding rather than absorbed into the model; the list is in
`reports/specification_findings.md`.

## The five validations

| | Question | Implemented by |
|:--|:--|:--|
| 1 | Can the schema represent the metadata of realistic MD systems? | `scripts/validation1_schema_coverage.py`, `scripts/validation1_interaction_coverage.py` |
| 2 | Does serialising and parsing preserve it? | `h5md_validator/roundtrip.py` |
| 3 | Is a file internally consistent? | `h5md_validator/validate.py` |
| 4 | Are invalid files detected? | `h5md_validator/corrupt.py` |
| 5 | What does the metadata cost? | `scripts/run_validation.py` (storage), `scripts/validation5_write_overhead.py` (time) |
| 6 | Do other engines' parameters map onto the same concepts? | `h5md_validator/cross_engine.py`, `scripts/validation6_amber_interop.py` |

Fidelity and coverage are kept apart throughout. A file can carry back everything
it holds without change and still have come from a source the schema could not
fully represent; conflating the two would let a coverage gap read as a success.
`scripts/validation1_interaction_coverage.py` is the sharpest form of the coverage question:
it counts bonded interaction entries in the run input against the file.

## What is in here, and what is not

This directory is the supplementary artifact: everything needed to re-run the
validation and read what it found, and nothing that the run itself produces.

| | |
|:--|:--|
| `h5md_validator/` | the validator, the schema model and the finding registries |
| `scripts/` | the drivers that build the reference files and analyse them |
| `reports/` | the written findings, the LaTeX tables and the figures |
| `README.md` | this file |

Two things live outside it:

- **`../validation_systems/`** — the prepared systems, one directory per system,
  holding the grompp inputs the validation was computed from. Kept separate
  because they are large and are archived rather than versioned.
- **`../validation_local/`** — the JSON results of a particular run, and working
  notes. These are outputs, not inputs: `scripts/run_all.sh` regenerates them.

A consequence worth knowing: the contents of `reports/` were generated *from*
those results. Someone re-running the validation on different hardware, or with
a different GROMACS build, will regenerate the results and should regenerate the
reports with `scripts/report.py` rather than assume the shipped numbers still
describe their run.

## Running it

Requires `h5py` and `numpy`, and a GROMACS build with the `biophysics` H5MD
writer. The environment used here:

```bash
PY_EXE=/path/to/python/env/python
export GMX_BIN=/path/to/gromacs/bin/gmx
```

<!-- ```bash
PY_EXE=/home/yzhang/mamba/envs/nearl-dev-env/bin/python
export GMX_BIN=/home/yzhang/gromacs_dev/gromacs/build/bin/gmx
``` -->

Everything after the reference files is chained by one script:

```bash
bash scripts/run_all.sh              # analysis only
bash scripts/run_all.sh --generate   # build the reference files first (hours)
```

or step by step:

```bash
# 1. build the reference files: systems x protocols -> .tpr + .h5md
$PY_EXE scripts/generate.py --list          # show the test matrix
$PY_EXE scripts/generate.py --all

# 2. conformance + round-trip + storage, per reference file
$PY_EXE scripts/run_validation.py

# 3. coverage, over the corpus and over bonded interactions
$PY_EXE scripts/validation1_schema_coverage.py
$PY_EXE scripts/validation1_interaction_coverage.py
$PY_EXE scripts/provenance_recoverability.py

# 4. negative tests
$PY_EXE scripts/validation4_negative_tests.py <a valid .h5md> [more .h5md ...]

# 5. what the metadata costs to write
$PY_EXE scripts/validation5_write_overhead.py --repeats 5

# 6. cross-engine: the same molecule prepared in AMBER
$PY_EXE scripts/validation6_amber_interop.py --h5md <dna .h5md> --prmtop <dna.prmtop>

# 7. figures and report
$PY_EXE scripts/make_figures.py
$PY_EXE scripts/report.py
```

Generated data — .tpr files, H5MD files and `gmx dump` caches — goes to
`$DATA_ROOT` (`/home/yzhang/Data/h5md_validation` by default), outside this
directory, because it runs to several GB.

The analysis scripts read and write `results/` relative to the working
directory. `RESULTS_DIR` redirects all of them at once, which is how to report
on a set of results kept elsewhere:

```sh
export RESULTS_DIR=../validation_local/results
$PY scripts/report.py                       # or: --results <dir>
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
h5md_validator/
    spec.py          schema model, transcribed from the specification
    validate.py      conformance and internal-consistency checks
    findings.py      finding and report types
    units.py         H5MD unit-notation grammar
    gmxdump.py       parser for `gmx dump -s`, the round-trip source side
    vocabulary.py    GROMACS <-> schema translation tables
    roundtrip.py     item-by-item comparison of source against file
    corrupt.py       the corruption catalogue for the negative tests
    deviations.py    static report: registry of specification-versus-implementation mismatches
    defects.py       static report: registry of reference-implementation defects
    limitations.py   static report: registry of configurations the test set could not produce
    cross_engine.py  static report: concept mapping across GROMACS/AMBER/NAMD/OpenMM/LAMMPS
scripts/
    matrix.py        simulation protocols of the test matrix
    systems.py       chemical systems of the test matrix, and how to build them
    generate.py      produce the reference files
    run_validation.py      validate + round-trip + measure every reference file
    coverage.py            element and item coverage over the corpus
    interaction_coverage.py  bonded interactions: source counts against file counts
    provenance_recoverability.py  where force-field provenance could come from
    run_negative_tests.py  corruption tests
    measure_write_overhead.py  H5MD against trajectory-only output, same .tpr
    amber_interop.py       the same molecule prepared in GROMACS and in AMBER
    make_figures.py        the two figures
    report.py              assemble reports and LaTeX tables
    run_all.sh             chain all of the above
results/             machine-readable results (JSON)
reports/             the written report, LaTeX tables and figures
```

## Three registries

Findings are kept in three separate registries rather than one list, because
they need different actions from different people:

| Registry | What it holds | Who resolves it |
|:--|:--|:--|
| `deviations.py` | the specification text and the implementation disagree | an editorial decision on the specification |
| `defects.py` | the implementation does not do what both it and the specification say | a code change in the writer |
| `limitations.py` | the test set could not produce a configuration at all | more test systems, or a stated caveat |

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
