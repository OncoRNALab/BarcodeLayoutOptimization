# HSDB — Hierarchical Spatial Dual Barcodes

Implementation of the combined Idea 1 + Idea 2 + Idea 3 from `IMPROVEMENT_PLAN.md`.

## Supported array sizes

| Array | Spots | Blocks (4×4) | bar_1 config |
|---|---|---|---|
| 118 × 182 | 21,476 | 30 × 46 = **1,380** | `bar1_zipcode/block_zipcode_config.json` (legacy 28 nt axes) |
| 235 × 363 | 85,305 | 59 × 91 = **5,369** | `bar1_zipcode/block_zipcode_config_235x363_g0_d4-10_far8_24nt.json` |
| 253 × 363 | 91,839 | 64 × 91 = **5,824** | `bar1_zipcode/block_zipcode_config_253x363.json` (legacy 28 nt axes) |

All sizes use a **36 nt** barcode. The current format is **24 nt bar_1 + 12 nt bar_2**.

## Barcode structure

```
[  24 nt bar_1 (block ID, Zipcode spatial gradient)  ] + [ 12 nt bar_2 (spot ID within block) ]
```

- **bar_1**: unique per block. Axis codes are 24-bit (= 12 nt) column + row Zipcodes generated in **flip mode**. Soft distance targets and the far-code minimum are enforced with `require_tolerance`. G content and junction homopolymers are rejected **inside** the candidate loop. There is no post-hoc G→T rewrite.
- **bar_2**: 16 fixed 12-nt sequences, reused in every block. Default zero-G set: `spot_barcodes_12nt_g0.txt` (min Hamming 5, leading run ≤ 1, homopolymer ≤ 5).
- Every `(bar_1, bar_2)` pair is globally unique by construction.
- Homopolymers are controlled across the bar_1/bar_2 junction (bar_2 leading runs ≤ 1; bar_1 endings leave room; the layout only assigns safe pairs).

## Decoding strategy

**Hierarchical decoder** (instead of flat 1-of-N):

1. Decode first 24 nt against bar_1 library → block ID.
2. Decode last 12 nt against bar_2 library (1-of-16) → spot within block.

## Directory layout

```
HSDB/
├── bar2_library/
│   ├── generate_spot_barcodes.py
│   ├── spot_barcodes_12nt_g0.txt          # current zero-G bar_2
│   └── spot_barcodes_8nt*.txt             # legacy 8 nt sets
├── bar1_zipcode/
│   ├── block_zipcode_config*_24nt.json
│   ├── generate_block_barcodes.py
│   ├── 235x363_g0_d4-10_far8_24+12/       # current zero-G bar_1
│   └── 118x182/  235x363/  …              # other / legacy outputs
├── layout_generator/
│   └── dual_barcode_hierarchical.py
├── simulation/
│   └── benchmark.py
└── outputs/
    └── 235x363_g0_d4-10_far8_24+12/
```

## Running the pipeline

### 235 × 363 array (current 24+12, zero G)

```bash
# Phase 1 — bar_2 library (12 nt, zero G)
python bar2_library/generate_spot_barcodes.py \
    --length 12 --n 16 --g-max 0 --max-homopolymer 5 --max-leading 1 \
    --out bar2_library/spot_barcodes_12nt_g0.txt

# Phase 2 — bar_1 library (24 nt = 12+12)
python bar1_zipcode/generate_block_barcodes.py \
    --config bar1_zipcode/block_zipcode_config_235x363_g0_d4-10_far8_24nt.json \
    --out-dir bar1_zipcode/235x363_g0_d4-10_far8_24+12 \
    --bar2 bar2_library/spot_barcodes_12nt_g0.txt

# Phase 3 — assemble layout
python layout_generator/dual_barcode_hierarchical.py \
    --array-rows 235 --array-cols 363 \
    --bar1-len 24 --bar2-len 12 \
    --bar1 bar1_zipcode/235x363_g0_d4-10_far8_24+12/bar1_barcodes.txt \
    --bar2 bar2_library/spot_barcodes_12nt_g0.txt \
    --out-dir outputs/235x363_g0_d4-10_far8_24+12 \
    --max-homopolymer 5
```

Far=10 fails under zero G at 24 bits; far=8 works.

## What was fixed (2026)

Previously bar_1 used the Zipcodes **random** candidate mode and then rewrote barcodes that broke the G or homopolymer limits. That destroyed the column/row structure and the distance gradient. Generation now uses the same flip-mode + in-loop composition checks as the updated Zipcodes designs. See `../AGENTS.md`.
