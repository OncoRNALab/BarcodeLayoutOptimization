# Zenodo deposit

Large files from this repository are deposited on Zenodo because they are too
big for git (FASTQs; full set on the order of 1–2 GB).

**DOI:** [10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834)

---

## What is deposited (manuscript branch)

### Real sequencing FASTQs

Demultiplexed MiSeq paired-end reads (47-8-8-47). Barcode occupies positions
11–46 of R1 (1-based) → QUIK `--barcode_start 10 --barcode_length 36`.

| Layout | Sample ID | Expected R1 / R2 basename pattern |
|--------|-----------|-----------------------------------|
| Dual full density | `26_018_S1` | `26_018_S1_*_R1_001.fastq(.gz)`, `…_R2_…` |
| Greedy + 2-OPT | `26_018_S3` | `26_018_S3_*_R{1,2}_001.fastq(.gz)` |
| Zipcodes | `26_018_S4` | `26_018_S4_*_R{1,2}_001.fastq(.gz)` |
| Dual half density | `26_022_S1` | `26_022_S1_*_R{1,2}_001.fastq(.gz)` |

Exact uploaded names may include a run or project prefix; match by sample ID.

### Simulated FASTQs (1M read pairs each)

Bare 36-nt barcodes (no handles) → QUIK `--barcode_start 0`.

| Layout | Regime | Files |
|--------|--------|-------|
| Dual full | high (3% / 15% / 14%) | `1M_high_dual_R1.fastq`, `1M_high_dual_R2.fastq` (+ optional `1M_high_dual_answers`) |
| Greedy + 2-OPT | high | `1M_high_greedy_R{1,2}.fastq` (+ answers) |
| Zipcodes | high | `1M_high_zip_R{1,2}.fastq` (+ answers) |
| Dual half | medium (4% / 5% / 14%) | `1M_med_dual_checker_R{1,2}.fastq` (+ answers) |

### Machine-readable summary tables (optional)

CSV companions to Additional files 2–9 may also be deposited under names like
`summary_tables_Additional_file_*.csv`.

---

## How to use the FASTQs

Decode with QUIK via
[BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)
(`--tool quik` only). Instructions and params templates:

→ [`docs/decode/README.md`](decode/README.md)

Whitelists are committed in this branch under `decoys/output/`.

To **regenerate** simulated reads instead of downloading them, see
[`REPRODUCE.md`](../REPRODUCE.md) §4.2 (RandomBarcodes / Press `makeerrors`).
The `simulation/` code on branch `main` is for NN≈8 redesigns only, not these
files.

---

## Checklist

- [ ] Real FASTQs for S1 / S3 / S4 / half-density Dual present
- [ ] Simulated `1M_high_*` and `1M_med_dual_checker_*` present
- [ ] MD5 sums recorded in Zenodo metadata
- [ ] DOI linked from manuscript Availability and from `REPRODUCE.md`
