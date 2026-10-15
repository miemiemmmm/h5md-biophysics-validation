#!/bin/bash
# format_python.sh — lint and format the validation suite with Ruff
# Usage: scripts/format_python.sh [update|check]   (default: check)

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"   # Ruff resolves first-party imports from the working directory
MODE="${1:-check}"

if [[ "${MODE}" != "update" && "${MODE}" != "check" ]]; then
    echo "Usage: $0 [update|check]" >&2
    exit 2
fi

if ! command -v ruff >/dev/null 2>&1; then
    echo "Ruff is not installed." >&2
    exit 2
fi

COMMON=(--target-version py39 --line-length 88)
LINT=(--select E,F,I,B,UP,SIM,RUF --ignore E501,B005,B008,RUF022,E402,BLE001,RUF001)

mapfile -d '' PYTHON_FILES < <(
    find "${ROOT}" \
        \( -name ".git" -o -name "__pycache__" -o -name ".ruff_cache" -o -name ".pytest_cache" \) -prune -o \
        -type f -name "*.py" -print0
)

if (( ${#PYTHON_FILES[@]} == 0 )); then
    echo "No Python files found in ${ROOT}."
    exit 0
fi

if [[ "${MODE}" == "check" ]]; then
    status=0

    echo "Checking lint rules..."
    ruff check "${COMMON[@]}" "${LINT[@]}" "${PYTHON_FILES[@]}" || status=1

    echo "Checking formatting..."
    ruff format --check "${COMMON[@]}" "${PYTHON_FILES[@]}" || status=1

    exit "${status}"
fi

echo "Applying lint fixes..."
ruff check --fix "${COMMON[@]}" "${LINT[@]}" "${PYTHON_FILES[@]}"

echo "Formatting..."
ruff format "${COMMON[@]}" "${PYTHON_FILES[@]}"
