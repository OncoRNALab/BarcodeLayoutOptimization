#!/usr/bin/env bash
# Generate structure-matched decoys for the four new 235×363 layouts.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DECOY_SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
OUT="$ROOT/decoys/output"
PY="${PYTHON:-python3}"
mkdir -p "$OUT"

echo "=== 1) Dual pseudo-replicate (shared_pool G10, 4×4) ==="
$PY "$DECOY_SCRIPTS/make_dual_decoys.py" \
  --layout "$ROOT/Dual_barcodes/layout_235x363_dual_shared_pool_g10.csv" \
  --library "$ROOT/inhouse_barcodes/bar_1in1_column_major_first18.txt" \
  --block-size 4 --bar1-len 18 --bar2-len 18 \
  --bar2-mode library \
  --g-max 0.10 --max-homopolymer 5 \
  --seed 2026 \
  --out-dir "$OUT/dual_shared_pool_g10"

echo "=== 2) Dual fixed_bar2 8×8 (mutate_bar1) ==="
$PY "$DECOY_SCRIPTS/make_dual_decoys.py" \
  --layout "$ROOT/Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10.csv" \
  --block-size 8 --bar1-len 18 --bar2-len 18 \
  --bar2-mode mutate_bar1 --mutate-flips-min 3 --mutate-flips-max 6 \
  --g-max 0.10 --max-homopolymer 5 \
  --seed 2026 \
  --out-dir "$OUT/dual_fixed_bar2_8x8_g10"

echo "=== 3) Zipcodes flip G10 ==="
$PY "$DECOY_SCRIPTS/make_zip_decoys.py" \
  --layout "$ROOT/Zipcodes/output_235x363_flip_g10_h5_d8-12_far10/layout_235x363_zipcodes_nosep.csv" \
  --segment-length 18 --fraction 1.0 --spacing bresenham \
  --g-max 0.10 --max-homopolymer 5 \
  --candidates 12000 --seed 2026 \
  --out-dir "$OUT/zipcodes_flip_g10"

echo "=== 4) HSDB spot-level (24+12, mutate bar2) ==="
$PY "$DECOY_SCRIPTS/make_dual_decoys.py" \
  --layout "$ROOT/HSDB/outputs/235x363_g10_h5_d8-12_far8_24+12/layout_235x363_hsdb.csv" \
  --block-size 4 --bar1-len 24 --bar2-len 12 \
  --bar2-mode mutate --mutate-flips-min 2 --mutate-flips-max 4 \
  --g-max 0.10 --max-homopolymer 5 \
  --seed 2026 \
  --out-dir "$OUT/hsdb_g10_24+12"

echo "Done. Outputs under $OUT"
