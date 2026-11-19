# barcode_design

Design and QC of spatial barcode layouts for photolithographic spatial transcriptomics arrays (MAS 2.0 maskless synthesizer, 235 × 363 = 85,305 features, 36 nt barcodes, A > C > G > T deposition order).

## Branches

| Branch | Purpose |
|--------|---------|
| **`main`** (this branch) | Living toolkit — four NN≈8 redesigned arrays, generators, decoys, read simulation. Use for new work. |
| **`manuscript`** | Frozen snapshot — printed layouts (Dual, Zipcodes, Greedy) + structure-matched decoys to reproduce the BMC Bioinformatics paper. See `REPRODUCE.md` on that branch. |

Large files (1M-read FASTQ pairs, oversized intermediates) are deposited on Zenodo. See [`docs/ZENODO.md`](docs/ZENODO.md).

---

## Four current arrays (NN≈8)

All 235 × 363, 36 nt, G < 10%, HP ≤ 5. Full details in [`docs/NewFourArrays_Summary.md`](docs/NewFourArrays_Summary.md) and [`SummaryNewDesigns.md`](SummaryNewDesigns.md).

| Array | Layout file | Real NN μ |
|-------|------------|-----------|
| Dual shared_pool minH6 | `Dual_barcodes/layout_235x363_dual_shared_pool_g10_minH6.csv` | **8.6** |
| Dual fixed_bar2 8×8 | `Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10_nn8.csv` | 7.0 |
| Zipcodes flip | `Zipcodes/output_235x363_flip_g10_h5_d12-16_far10/` | 7.0 |
| HSDB colored | `HSDB/outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored/` | 8.0 |

Structure-matched decoys: `decoys/output/dual_shared_pool_g10_minH6/`, `dual_fixed_bar2_8x8_g10_nn8/`, `zipcodes_flip_g10_nn8/`, `hsdb_g10_nn8_colored/`.

NN histograms and activation masks: `decoys/figures/`.

---

## Quickstart

### Regenerate four layouts + decoys + figures

```bash
# Zipcodes flip
cd Zipcodes
python generate_layout_235x363_zipcodes.py \
  config_235x363_flip_g10_h5_d12-16_far10.json \
  output_235x363_flip_g10_h5_d12-16_far10

# HSDB (bar₁ → color → layout)
cd ../HSDB
python bar2_library/generate_spot_barcodes.py \
  --n 16 --length 12 --min-hamming 8 --g-max 0.10 --max-homopolymer 5 \
  --gc-min 0 --gc-max 1 --out bar2_library/spot_barcodes_12nt_g10_nn8.txt
python bar1_zipcode/generate_block_barcodes.py \
  --config bar1_zipcode/block_zipcode_config_235x363_g10_h5_d12-16_far8_24nt.json \
  --out-dir bar1_zipcode/235x363_g10_h5_d12-16_far8_24+12_nn8 \
  --bar2 bar2_library/spot_barcodes_12nt_g10_nn8.txt
python layout_generator/color_bar2_for_nn.py \
  --bar1 bar1_zipcode/235x363_g10_h5_d12-16_far8_24+12_nn8/bar1_barcodes.txt \
  --target-nn 8 --out-dir outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored

# Dual shared_pool minH6
cd ../Dual_barcodes
python generate_layout_235x363.py config_235x363_shared_pool_g10_minH6.json

# Dual fixed_bar2 8×8 NN≈8
python generate_layout_235x363.py config_235x363_fixed_bar2_8x8_g10_nn8.json

# Decoys + NN histograms + activation masks
cd ../decoys/scripts
# see decoys/README.md for per-array commands
python plot_nn_distance_histograms.py --out-dir ../figures --sample 1000 --seed 1
python plot_activation_masks.py --cycle 1 --out-dir ../figures
```

### Simulated reads (1M pairs, sub=0.04 ins=0.05 del=0.14)

```bash
cd simulation
for p in dual_shared_pool_g10_minH6 dual_fixed_bar2_8x8_g10_nn8 zipcodes_flip_g10_nn8 hsdb_g10_nn8_colored; do
  python make_reads.py --preset $p --out-dir output/$p --n-reads 1000000 --seed 1 \
    --srate 0.04 --irate 0.05 --drate 0.14
done
```

Outputs land in `simulation/output/` (gitignored). Pre-generated reads are on Zenodo (see [`docs/ZENODO.md`](docs/ZENODO.md)).

---

## Repository structure

```
AGENTS.md                     Entry point for agents / collaborators (authoritative)
docs/
  methods.md                  Manuscript Methods excerpt (printed designs, authoritative)
  NewFourArrays_Summary.md    Full methods for four NN≈8 arrays
  ZENODO.md                   Zenodo deposit checklist and DOI
SummaryNewDesigns.md          Quick reference table for four arrays
Dual_barcodes/                Dual barcode generator + four configs/layouts
Zipcodes/                     Zipcode generator + flip-mode configs/outputs
HSDB/                         Hierarchical Spatial Dual Barcode generator + outputs
Quadtree/                     Quadtree method (non-manuscript)
inhouse_barcodes/             Random 18-nt input libraries
decoys/
  scripts/                    Decoy generators, QC, NN histogram and activation mask plots
  output/                     Generated decoy files (four arrays + manuscript)
  qc/                         Real-barcode lists and QC JSONs
  figures/                    NN histograms, activation masks
simulation/
  make_reads.py               Unequal-walk read simulator (R1 + R2=RC(R1), Phred, FASTQ)
  error_model.py              Walk channel (sub/ins/del with configurable rates)
  README.md                   Channel parameters and regenerate commands
  output/                     (gitignored) FASTQ outputs → Zenodo
data/
  1in1/                       Printed full-density layouts (Dual, Zip, Greedy)
  1in2/                       Printed half-density Dual layout
  decoy_layouts/              Original independent-random decoy pools
archive/historical/           Obsolete notes and scripts kept for reference
```

---

## Methods

- **Printed designs (BMC paper):** see `docs/methods.md` and the `manuscript` branch.
- **Four NN≈8 redesigns:** see `docs/NewFourArrays_Summary.md`.
- **Synthesis model:** A > C > G > T deposition order, barcodes start in the same cycle; two-half printing for Dual and Zipcodes (bar₁/column first, then bar₂/row).
- **Decoys:** structure-matched; φ = 1/3 for all four arrays.

## Printed layout files (ground truth for paper)

> These are also on the `manuscript` branch under `data/1in1/`.

| Method | File |
|--------|------|
| Dual | `data/1in1/layout_235x363_dual_4x4.csv` |
| Zipcodes | `data/1in1/layout_235x363_zipcodes_nosep.csv` |
| Greedy | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` |
| Dual half-density | `data/1in2/layout_1in2_dual_4x4.csv` |
