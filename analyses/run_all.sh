#!/usr/bin/env bash
# Re-run every manuscript analysis from layouts, decoys, and QUIK outputs.
# Order matters: real_arrays reads the phi_eff values written by precision_calibration.
# Usage: ./run_all.sh [python]   (default: python3; needs the packages in requirements.txt)
set -euo pipefail
cd "$(dirname "$0")"
PY=${1:-python3}
ROOT="$(cd .. && pwd)"

run_notebook() {
    (cd "$(dirname "$1")" && "$PY" -m jupyter nbconvert --to notebook --execute --inplace "$(basename "$1")")
}

# Regenerate structure-matched decoys (repo-root decoys/, not data/decoys/)
(cd "$ROOT/decoys/scripts" && ./run_all.sh)

run_notebook precision_calibration/precision_calibration.ipynb
run_notebook real_arrays/real_arrays.ipynb
run_notebook decoy_quality/decoy_distances.ipynb

(cd recall_sequence_features && "$PY" read_accounting.py && "$PY" features.py && "$PY" decoding_rescue.py && "$PY" error_profiles.py && "$PY" analyse.py && "$PY" g_position.py)
run_notebook recall_sequence_features/read_accounting.ipynb
run_notebook recall_sequence_features/recall_sequence_features.ipynb
run_notebook recall_sequence_features/decoding_rescue.ipynb
run_notebook recall_sequence_features/error_profiles.ipynb

(cd dual_library_diversity && "$PY" compare_dual_libraries.py \
    --printed "$ROOT/data/1in1/layout_235x363_dual_4x4.csv" \
    --lib1 "$ROOT/inhouse_barcodes/bar_1in1_column_major_first18.txt" \
    --lib2 "$ROOT/inhouse_barcodes/bar_1in1_column_major_last18.txt" \
    --out dual_library_diversity.json)
