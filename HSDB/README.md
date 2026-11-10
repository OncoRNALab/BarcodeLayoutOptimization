# HSDB — Hierarchical Spatial Dual Barcodes

Implementation of the combined Idea 1 + Idea 2 + Idea 3 from `IMPROVEMENT_PLAN.md`.

## Supported array sizes

| Array | Spots | Blocks (4×4) | bar_1 library | bar_2 library |
|---|---|---|---|---|
| 118 × 182 | 21,476 | 30 × 46 = **1,380** | `bar1_barcodes_28nt.txt` | (shared) |
| 253 × 363 | 91,839 | 64 × 91 = **5,824** | `bar1_barcodes_28nt_253x363.txt` | (shared) |

Both sizes use the **same 36 nt barcode format** — only the bar_1 library differs.

## Barcode structure

```
[  28 nt bar_1 (block ID, Zipcode spatial gradient)  ] + [ 8 nt bar_2 (spot ID within block) ]
```

- **bar_1**: unique per block, generated with Zipcode spatial constraints so that
  adjacent blocks differ by ~4 bases and far blocks differ by ≥ 16 bases.
- **bar_2**: 16 fixed sequences, reused in every block, min Hamming distance = 5.
- Every `(bar_1, bar_2)` pair is globally unique by construction.
- G content ≤ 25% per barcode (Idea 3 constraint enforced at generation time).

## Decoding strategy

**Hierarchical decoder** (instead of flat 1-of-N):

1. Decode first 28 nt against bar_1 library → block ID.
2. Decode last 8 nt against bar_2 library (1-of-16) → spot within block.

## Directory layout

```
HSDB/
├── bar2_library/
│   ├── generate_spot_barcodes.py       # Phase 1: generate 16 spot barcodes
│   └── spot_barcodes_8nt.txt           # output (shared by all array sizes)
│
├── bar1_zipcode/
│   ├── block_zipcode_config.json            # config for 30×46 block grid (118×182)
│   ├── block_zipcode_config_253x363.json    # config for 64×91 block grid (253×363)
│   └── generate_block_barcodes.py      # Phase 2: run Zipcode on block grid
│
├── layout_generator/
│   └── dual_barcode_hierarchical.py    # Phase 3: assemble full 36 nt layout
│
├── simulation/
│   └── benchmark.py                    # Phase 4: benchmark vs current 18+18
│
└── outputs/                            # generated layout files and QC reports
    ├── layout_118x182_hsdb.csv
    ├── barcodes_hsdb_118x182_36nt.txt
    ├── layout_253x363_hsdb.csv
    └── barcodes_hsdb_253x363_36nt.txt
```

## Running the pipeline

### 118 × 182 array

```bash
# Phase 1 — bar_2 library (one-time, shared)
python bar2_library/generate_spot_barcodes.py

# Phase 2 — bar_1 library for 118×182
python bar1_zipcode/generate_block_barcodes.py \
    --config bar1_zipcode/block_zipcode_config.json \
    --out-dir bar1_zipcode/

# Phase 3 — assemble layout
python layout_generator/dual_barcode_hierarchical.py \
    --array-rows 118 --array-cols 182 \
    --bar1 bar1_zipcode/bar1_barcodes_28nt.txt \
    --out-dir outputs/

# Phase 4 — benchmark
python simulation/benchmark.py \
    --hsdb-barcodes outputs/barcodes_hsdb_118x182_36nt.txt
```

### 253 × 363 array

```bash
# Phase 2 — bar_1 library for 253×363  (bar_2 already generated)
python bar1_zipcode/generate_block_barcodes.py \
    --config bar1_zipcode/block_zipcode_config_253x363.json \
    --out-dir bar1_zipcode/253x363/

# Phase 3 — assemble layout
python layout_generator/dual_barcode_hierarchical.py \
    --array-rows 253 --array-cols 363 \
    --bar1 bar1_zipcode/253x363/bar1_barcodes_28nt.txt \
    --out-dir outputs/

# Phase 4 — benchmark
python simulation/benchmark.py \
    --hsdb-barcodes outputs/barcodes_hsdb_253x363_36nt.txt \
    --hsdb-bar1    bar1_zipcode/253x363/bar1_barcodes_28nt.txt \
    --skip-dual
```

## Implementation order (from IMPROVEMENT_PLAN.md)

1. **Idea 3 first** — G ≤ 25% enforced at generation time in all scripts.
2. Phase 1 — `bar2_library/generate_spot_barcodes.py` (quick, done once)
3. Phase 2 — run for each array size with its config
4. Phase 3 — run for each array size
5. Phase 4 — benchmark each design
