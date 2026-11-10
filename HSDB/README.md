# HSDB — Hierarchical Spatial Dual Barcodes

Implementation of the combined Idea 1 + Idea 2 + Idea 3 from `IMPROVEMENT_PLAN.md`.

## Supported array sizes

| Array | Spots | Blocks (4×4) | bar_1 config |
|---|---|---|---|
| 118 × 182 | 21,476 | 30 × 46 = **1,380** | `bar1_zipcode/block_zipcode_config.json` |
| 235 × 363 | 85,305 | 59 × 91 = **5,369** | `bar1_zipcode/block_zipcode_config_235x363.json` |
| 253 × 363 | 91,839 | 64 × 91 = **5,824** | `bar1_zipcode/block_zipcode_config_253x363.json` |

All sizes use the **same 36 nt barcode format** — only the bar_1 library differs.

## Barcode structure

```
[  28 nt bar_1 (block ID, Zipcode spatial gradient)  ] + [ 8 nt bar_2 (spot ID within block) ]
```

- **bar_1**: unique per block. Axis codes are generated in **flip mode** (each candidate is the preceding code with a few bit flips). Soft distance targets and the far-code minimum are enforced with `require_tolerance`. G content and junction homopolymers are rejected **inside** the candidate loop. There is no post-hoc G→T rewrite.
- **bar_2**: 16 fixed 8-nt sequences, reused in every block, min Hamming distance = 5, G ≤ 20%, homopolymer ≤ 3 (enforced during selection).
- Every `(bar_1, bar_2)` pair is globally unique by construction.
- Combined barcodes keep G ≤ 25%. Homopolymers at the bar_1/bar_2 junction can still exceed 3 nt.

## Decoding strategy

**Hierarchical decoder** (instead of flat 1-of-N):

1. Decode first 28 nt against bar_1 library → block ID.
2. Decode last 8 nt against bar_2 library (1-of-16) → spot within block.

## Directory layout

```
HSDB/
├── bar2_library/
│   ├── generate_spot_barcodes.py
│   └── spot_barcodes_8nt.txt
├── bar1_zipcode/
│   ├── block_zipcode_config*.json
│   ├── generate_block_barcodes.py
│   ├── 118x182/   235x363/   …   # per-size bar_1 outputs
├── layout_generator/
│   └── dual_barcode_hierarchical.py
├── simulation/
│   └── benchmark.py
└── outputs/
```

## Running the pipeline

### 235 × 363 array

```bash
# Phase 1 — bar_2 library (one-time, shared)
python bar2_library/generate_spot_barcodes.py

# Phase 2 — bar_1 library
python bar1_zipcode/generate_block_barcodes.py \
    --config bar1_zipcode/block_zipcode_config_235x363.json \
    --out-dir bar1_zipcode/235x363/

# Phase 3 — assemble layout
python layout_generator/dual_barcode_hierarchical.py \
    --array-rows 235 --array-cols 363 \
    --bar1 bar1_zipcode/235x363/bar1_barcodes_28nt.txt \
    --out-dir outputs/235x363/
```

### 118 × 182 array

```bash
python bar1_zipcode/generate_block_barcodes.py \
    --config bar1_zipcode/block_zipcode_config.json \
    --out-dir bar1_zipcode/118x182/

python layout_generator/dual_barcode_hierarchical.py \
    --array-rows 118 --array-cols 182 \
    --bar1 bar1_zipcode/118x182/bar1_barcodes_28nt.txt \
    --out-dir outputs/118x182/
```

## What was fixed (2026)

Previously bar_1 used the Zipcodes **random** candidate mode and then rewrote barcodes that broke the G or homopolymer limits. That destroyed the column/row structure and the distance gradient. Generation now uses the same flip-mode + in-loop composition checks as the updated Zipcodes designs. See `../AGENTS.md`.
