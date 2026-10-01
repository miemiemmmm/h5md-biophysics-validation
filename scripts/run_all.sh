#!/bin/bash -eu
# Reproduce the whole technical validation from the reference files onward.
#
#   bash scripts/run_all.sh            # analysis only (reference files must exist)
#   bash scripts/run_all.sh --generate # build the reference files first (hours)

set -eu

# The extracted validation dataset. Unset means the copy beside this checkout;
# set it to wherever you unpacked the download. Everything the run produces goes
# into <dataset>/_data, which is several GB and safe to delete afterwards.
if [ -n "${H5MD_VALIDATION_DATA:-}" ]; then
    export H5MD_VALIDATION_DATA
    DATASET="$H5MD_VALIDATION_DATA"
else
    DATASET="$(cd "$(dirname "$0")/../.." && pwd)/validation_systems"
fi
REF="$DATASET/_data/reference_files"

# Both PY_EXE and GMX_BIN are required as environment variables.
for var in PY_EXE GMX_BIN; do
    if [ -z "${!var:-}" ]; then
        cat >&2 <<USAGE
$var is not set. Both of these are required:
  export PY_EXE=/path/to/python      # needs h5py and numpy
  export GMX_BIN=/path/to/gmx        # a build with the biophysics H5MD writer
USAGE
        exit 1
    fi
done

# Set but wrong is as common as unset, and fails much later and less clearly.
command -v "$PY_EXE" >/dev/null 2>&1 || {
    echo "PY_EXE=$PY_EXE is not executable" >&2; exit 1; }
command -v "$GMX_BIN" >/dev/null 2>&1 || {
    echo "GMX_BIN=$GMX_BIN is not executable" >&2; exit 1; }
"$PY_EXE" -c "import h5py, numpy" 2>/dev/null || {
    echo "PY_EXE=$PY_EXE cannot import h5py and numpy" >&2; exit 1; }
export GMX_BIN
cd "$(dirname "$0")/.."

if [ "${1:-}" = "--generate" ]; then
    echo "==> generating reference files"
    $PY_EXE -u scripts/generate.py --all --ntomp 16
fi

echo "==> conformance, round-trip and storage, per reference file"
$PY_EXE scripts/run_validation.py

echo "==> coverage analysis"
$PY_EXE scripts/coverage.py
$PY_EXE scripts/interaction_coverage.py
$PY_EXE scripts/provenance_recoverability.py
$PY_EXE scripts/parameter_surface.py "$DATASET"/_data/mdout/*.mdp || echo "(skipped)"

echo "==> negative tests"
# Three bases with different shapes, so that as many corruptions as possible apply:
# a small protein, a nucleic acid, and a multicomponent membrane system.
BASES=""
for name in bena_12k__selection-noncontiguous dna_amber__nvt-vrescale \
            memprotmd_5bkf__nvt-vrescale solvent_ions__nvt-vrescale \
            protein_dna__nvt-vrescale; do
    f="$REF/$name.h5md"
    [ -f "$f" ] && BASES="$BASES $f"
done
$PY_EXE scripts/run_negative_tests.py $BASES

echo "==> write overhead"
$PY_EXE scripts/measure_write_overhead.py --repeats 5 || echo "(skipped)"

echo "==> cross-engine interoperability"
PRMTOP=${PRMTOP:-}
if [ -n "$PRMTOP" ] && [ -f "$PRMTOP" ]; then
    $PY_EXE scripts/amber_interop.py \
        --h5md "$REF/dna_amber__nvt-vrescale.h5md" \
        --prmtop "$PRMTOP"
else
    echo "(set PRMTOP=... to include the AMBER comparison)"
fi

echo "==> figures"
$PY_EXE scripts/make_figures.py

echo "==> report"
$PY_EXE scripts/report.py
