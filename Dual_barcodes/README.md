# Dual Barcode Layout Generator

Block-wise Dual barcodes for photolithographic spatial arrays: each \(B \times B\) block shares one `bar_1`; spots inside the block get distinct `bar_2`s. The combined barcode is `bar_1 + bar_2`. Global uniqueness is the **pair**, so `bar_2` can be reused across blocks.

## Modes

### 1. `shared_pool` (legacy / printed)

- One predefined half-barcode library for both `bar_1` and `bar_2`.
- No distance design on `bar_2`; no combined G/HP filter unless `sequence_qc` is set.
- Default block size **4×4**.
- Config: `config_235x363.json`

```bash
python generate_layout_235x363.py config_235x363.json
```

### 2. `fixed_bar2` (redesigned)

- `n_bar2 = B²` distance-designed repertoire (min Hamming), reused in every block.
- `bar_1` from a file pool or newly generated; each accepted `bar_1` must be QC-compatible with the `bar_2` set.
- Optional combined QC (G, homopolymer, junction via leading-run limit on `bar_2`).
- Config: `config_235x363_fixed_bar2_g10.json`

```bash
python generate_layout_235x363.py config_235x363_fixed_bar2_g10.json
```

Generate only the `bar_2` set:

```bash
python generate_bar2_repertoire.py --block-size 4 --length 18 --min-hamming 8 \
  --g-max 0.10 --max-homopolymer 5 --out bar2_repertoire_4x4_18nt_g10.txt
```

## Configurable block size

`grid_configuration.block_size` is passed into `DualBarcodeLayout` (not documentation-only). Edge blocks are truncated.

## Key files

| File | Role |
|---|---|
| `dual_barcode_layout.py` | Layout engine (`shared_pool` / `fixed_bar2`) |
| `sequence_qc.py` | G / HP / GC / leading-run helpers |
| `generate_bar2_repertoire.py` | Distance-designed `bar_2` set |
| `generate_layout_235x363.py` | Config-driven runner |

## Python API

```python
from dual_barcode_layout import DualBarcodeLayout, load_barcodes_from_file
from sequence_qc import SequenceQC

# Legacy
pool = load_barcodes_from_file("../inhouse_barcodes/bar_1in1_column_major_first18.txt")
layout = DualBarcodeLayout(235, 363, pool, block_size=4)
layout.generate_layout(seed=42)

# Fixed bar_2 + QC
qc = SequenceQC(g_content_max=0.10, max_homopolymer=5, max_leading=1)
layout = DualBarcodeLayout(
    235, 363, bar1_pool,
    block_size=4,
    bar2_barcodes=bar2_set,
    sequence_qc=qc,
)
layout.generate_layout(seed=42)
```
