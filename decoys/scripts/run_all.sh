#!/usr/bin/env bash
# Generate the structure-matched decoys used in the manuscript and their QC.
#   dual/                     full-density Dual (1:2, 8 decoys per complete 4x4 block)
#   dual_1in2/                half-density Dual checkerboard (same per-block rule)
#   zipcodes_balanced_khalf/  Zipcodes V2 (1:1, Bresenham, 181 virtual columns x 117 virtual rows)
#   greedy/                   Greedy + 2-OPT (1:2, drawn from an independent random pool)
#   zipcodes/                 earlier Zipcodes set (1:2, linspace); kept for comparison only
# Usage: bash run_all.sh [OUT_ROOT]   (default ../output; pass another folder to avoid overwriting)
set -euo pipefail
cd "$(dirname "$0")"

ROOT=../..
REAL=$ROOT/data/1in1
HALF=$ROOT/data/1in2
OLD=$ROOT/data/decoy_layouts/barcodes/1in1
LIB=$ROOT/inhouse_barcodes/bar_1in1_column_major_first18.txt
OUT=${1:-../output}
QC=$OUT/../qc
mkdir -p "$QC"

python3 make_dual_decoys.py --layout "$REAL/layout_235x363_dual_4x4.csv" \
    --library "$LIB" --fraction 0.5 --seed 2026 --out-dir "$OUT/dual"
python3 make_dual_decoys.py --layout "$HALF/layout_1in2_dual_4x4.csv" \
    --library "$LIB" --fraction 0.5 --seed 2026 --out-dir "$OUT/dual_1in2"
python3 make_zip_decoys.py --layout "$REAL/layout_235x363_zipcodes_nosep.csv" \
    --spacing bresenham --n-virtual-columns 181 --n-virtual-rows 117 \
    --candidates 4000 --seed 2026 --out-dir "$OUT/zipcodes_balanced_khalf"
python3 make_zip_decoys.py --layout "$REAL/layout_235x363_zipcodes_nosep.csv" \
    --fraction 0.5 --seed 2026 --out-dir "$OUT/zipcodes"
python3 make_greedy_decoys.py --real "$REAL/layout_235x363_randombarcodes_opt.txt" \
    --decoys "$OLD/greedy_decoybar_1in1_column_major.txt" \
    --fraction 0.5 --seed 2026 --out-dir "$OUT/greedy"

python3 check_decoys.py --real "$REAL/layout_235x363_dual_4x4.txt" \
    --decoys "$OUT/dual/dual_decoys.txt" --label structure_matched \
    --decoys "$OLD/layout_235x363_dual_4x4.txt" --label old_independent \
    --match-size --out "$QC/dual_qc.json"
python3 check_decoys.py --real "$REAL/layout_235x363_zipcodes_nosep.txt" \
    --decoys "$OUT/zipcodes_balanced_khalf/zip_decoys.txt" --label structure_matched \
    --decoys "$OLD/layout_235x363_zipcodes_nosep.txt" --label old_independent \
    --match-size --out "$QC/zipcodes_balanced_khalf_qc.json"
python3 check_decoys.py --real "$REAL/layout_235x363_randombarcodes_opt.txt" \
    --decoys "$OUT/greedy/greedy_decoys.txt" --label pool_subset \
    --out "$QC/greedy_qc.json"
