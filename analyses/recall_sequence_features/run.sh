#!/usr/bin/env bash
# Builds and executes the recall / sequence-feature notebooks from their jupytext sources.
# Usage: ./run.sh [python] [notebook ...]   (default: all notebooks in order)
set -euo pipefail
cd "$(dirname "$0")"
PY="${1:-python}"
shift || true
NOTEBOOKS=("$@")
if [ ${#NOTEBOOKS[@]} -eq 0 ]; then
    NOTEBOOKS=(01_read_accounting 02_feature_table 03_groups_and_tests 04_models 05_simulation_recall 06_error_profiles 07_relaxed_decoding)
fi
for nb in "${NOTEBOOKS[@]}"; do
    echo "== $nb"
    "$PY" -m jupytext --to ipynb --output "$nb.ipynb" "$nb.py"
    "$PY" -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 "$nb.ipynb"
done
