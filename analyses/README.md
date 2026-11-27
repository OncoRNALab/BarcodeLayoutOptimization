# Manuscript analyses

Python scripts and notebooks that turn QUIK decode outputs into the paper’s
rates, calibration checks, and figure inputs.

## Prerequisites

1. Layouts and decoys from this `manuscript` branch (`data/`, `decoys/`).
2. Decoded QUIK outputs from [`docs/decode/README.md`](../docs/decode/README.md)
   placed under `data/quik_out_manuscript/` (same folder layout as in the paper
   working copy: e.g. `dual_combined_decoy_thr7/`, `simreads_thr7_1M/`, …).
3. Optional: raw FASTQs under `data/raw_fastq_files/` where a script reads them
   directly (see individual scripts).

## Run everything

```bash
cd analyses
pip install -r requirements.txt
./run_all.sh          # or: ./run_all.sh python3
```

Order matters: `precision_calibration` writes φ_eff values that `real_arrays`
consumes.

## Main pieces

| Folder | Role |
|--------|------|
| `precision_calibration/` | Simulation-based precision / φ_eff |
| `real_arrays/` | Printed-array assignment and decoy rates |
| `decoy_quality/` | Decoy distance / local-share checks |
| `recall_sequence_features/` | Full-length accounting, sequence features |
| `dual_checkerboard/` | Half-density Dual metrics |
| `zipcodes_v2/` | Zipcodes recompute helpers |
| `library_complexity/` | Library complexity panels |
| `figures_manuscript/` | Results / schematic figure builders |
| `dual_library_diversity/` | Dual half-library comparison |

Figure → script mapping (main-text Figures 1–9): [`REPRODUCE.md`](../REPRODUCE.md) §9.

See the parent [`REPRODUCE.md`](../REPRODUCE.md) workflow overview.
