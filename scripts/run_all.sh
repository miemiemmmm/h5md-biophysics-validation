#!/bin/bash -eu
# Reproduce the whole technical validation from the reference files onward.
#
#   bash scripts/run_all.sh            # analysis only (reference files must exist)
#   bash scripts/run_all.sh --generate # build the reference files first (hours)

set -eu

# The extracted validation dataset. Unset means the copy beside this checkout;
# set it to wherever you unpacked the download. Where the run writes is decided
# in config/matrix.py (DATA_DIR and RESULTS_DIR override it) and printed below.
if [ -n "${H5MD_VALIDATION_DATA:-}" ]; then
    export H5MD_VALIDATION_DATA
fi

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

REF="$($PY_EXE -c "from config.matrix import DATA_ROOT; print(DATA_ROOT)")/reference_files"
echo "==> locations"
$PY_EXE -c "from config.matrix import output_locations; print(output_locations())"

if [ "${1:-}" = "--generate" ]; then
    echo "==> generating reference files"
    $PY_EXE -u scripts/build_reference_files.py --all --ntomp 16
fi

# Each validation reads the reference files (and, for Validation 2, the GROMACS
# run inputs) and writes its own result; none depends on another's output.
echo "==> Validation 1: schema coverage"
$PY_EXE scripts/validation1_schema_coverage.py
$PY_EXE scripts/validation1_parameter_surface.py || echo "(skipped)"

echo "==> Validation 2: round trip"
$PY_EXE scripts/validation2_roundtrip.py
$PY_EXE scripts/validation2_interaction_coverage.py

echo "==> Validation 3: conformance"
$PY_EXE scripts/validation3_conformance.py

echo "==> Validation 4: negative tests"
# The bases are declared in config/matrix.py so that a renamed system carries
# through here instead of silently dropping out of the test.
BASES=""
for name in $($PY_EXE -c "from config.matrix import NEGATIVE_TEST_BASES; print(' '.join(NEGATIVE_TEST_BASES))"); do
    f="$REF/$name.h5md"
    [ -f "$f" ] && BASES="$BASES $f"
done
$PY_EXE scripts/validation4_negative_tests.py $BASES

echo "==> Validation 5: storage"
$PY_EXE scripts/validation5_storage.py

echo "==> Validation 6: cross-engine interoperability"
# Pairs come from AMBER_INTEROP_PAIRS in config/matrix.py, where each reference
# file is declared with the AMBER topology of the SAME molecule.
$PY_EXE scripts/validation6_amber_interop.py --run-all

echo "==> supporting analysis: force-field provenance"
$PY_EXE scripts/provenance_recoverability.py

echo "==> figures"
$PY_EXE scripts/validation5_cost_figure.py
