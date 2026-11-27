# Decode FASTQ files with QUIK

Paper results used **QUIK only** (`strategy = 4_mer_gpu_v4`), run through the
same Nextflow wrapper as our companion benchmarking paper:

[OncoRNALab/BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)

That pipeline can call QUIK, RandomBarcodes, or Columba. For this manuscript,
set `--tool quik` (or `"tool": "quik"` in a params JSON) and ignore the other
decoders.

---

## 1. Prerequisites

- **Nextflow** ≥ 23.04
- **GPU** + either:
  - `-profile local,singularity` (Singularity/Apptainer; recommended for
    reproducibility), or
  - `-profile slurm` (HPC modules for CMake / CUDA / PyTorch — see the
    benchmark README)
- This repository checked out on branch **`manuscript`**
- FASTQs from Zenodo ([10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834))

---

## 2. Clone the decoding pipeline

```bash
git clone https://github.com/OncoRNALab/BarcodeCalling_benchmark.git
cd BarcodeCalling_benchmark
```

Optional: pin a known-good commit if you need bit-identical tooling later.

For Singularity cache (avoid filling `$HOME`):

```bash
export APPTAINER_TMPDIR=/path/to/scratch/.apptainer/tmp
export APPTAINER_CACHEDIR=/path/to/scratch/.apptainer/cache
export NXF_SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
```

---

## 3. Inputs for this paper

### 3.1 Whitelists (committed here)

Always decode against **real barcodes + layout-matched decoys** (same files used
in the paper):

| Layout | Whitelist |
|--------|-----------|
| Dual full density | `decoys/output/dual/dual_whitelist_real_plus_decoys.txt` |
| Dual half density | `decoys/output/dual_1in2/dual_whitelist_real_plus_decoys.txt` |
| Zipcodes | `decoys/output/zipcodes_balanced_khalf/zip_whitelist_real_plus_decoysv2.txt` |
| Greedy + 2-OPT | `decoys/output/greedy/greedy_whitelist_real_plus_decoys.txt` |

One barcode sequence per line. Paths are relative to this repo root.

### 3.2 Real sequencing FASTQs (Zenodo)

Demultiplexed MiSeq paired-end files (47-cycle R1/R2 with 10-nt handles).
Paper sample IDs:

| Layout | Sample ID | Typical Illumina basename |
|--------|-----------|---------------------------|
| Dual full | `26_018_S1` | `26_018_S1_L001_R{1,2}_001.fastq(.gz)` |
| Greedy + 2-OPT | `26_018_S3` | `26_018_S3_L001_R{1,2}_001.fastq(.gz)` |
| Zipcodes | `26_018_S4` | `26_018_S4_L001_R{1,2}_001.fastq(.gz)` |
| Dual half (checkerboard) | `26_022_S1` | `26_022_S1_*_R{1,2}_001.fastq(.gz)` |

Exact filenames on Zenodo may include a project prefix; point `--r1_fastq` /
`--r2_fastq` at the files you download.

**QUIK window for real reads:** `--barcode_start 10 --barcode_length 36`
(barcode occupies positions 11–46 of the 47-nt read, 1-based).

### 3.3 Simulated FASTQs (Zenodo)

| Layout | Error regime | Files |
|--------|--------------|-------|
| Dual full | high (3% sub / 15% ins / 14% del) | `1M_high_dual_R1.fastq`, `1M_high_dual_R2.fastq` |
| Greedy + 2-OPT | high | `1M_high_greedy_R1.fastq`, `1M_high_greedy_R2.fastq` |
| Zipcodes | high | `1M_high_zip_R1.fastq`, `1M_high_zip_R2.fastq` |
| Dual half | medium (4% / 5% / 14%) | `1M_med_dual_checker_R1.fastq`, `1M_med_dual_checker_R2.fastq` |

Simulated reads are **bare 36-nt barcodes** (no handles).

**QUIK window for simulated reads:** `--barcode_start 0 --barcode_length 36`.

Optional ground-truth files (`1M_high_*_answers`, etc.) can be passed as
`--ground_truth` with `--data_mode simulated` if you want pipeline precision
metrics; the manuscript also scores assignments offline from filtered FASTQs.

---

## 4. QUIK parameters used in the paper

| Parameter | Value |
|-----------|-------|
| `tool` | `quik` |
| `strategy` | `4_mer_gpu_v4` |
| `distance_measure` | `SEQUENCE_LEVENSHTEIN` |
| `rejection_threshold` (`r`) | full density: `{5, 6, 7}`; half density: `{6, 7, 8}` |
| `barcode_length` | `36` |
| `barcode_start` | `10` (real) / `0` (simulated) |
| `data_mode` | `real` (no ground truth) or `simulated` (with answers) |

Primary Results figures use **`r = 7`** unless noted otherwise.

---

## 5. Run QUIK only (example)

Replace paths with your local clone of this repo and your Zenodo download
directory.

### Real Dual full density, `r = 7`

```bash
LAYOUT=/path/to/BarcodeLayoutOptimization   # manuscript branch
FASTQ=/path/to/zenodo_fastqs
OUT=results/dual_full_r7

nextflow run main.nf \
  --tool quik \
  --barcode_file "$LAYOUT/decoys/output/dual/dual_whitelist_real_plus_decoys.txt" \
  --r1_fastq "$FASTQ/26_018_S1_L001_R1_001.fastq.gz" \
  --r2_fastq "$FASTQ/26_018_S1_L001_R2_001.fastq.gz" \
  --sample_id 26_018_S1_quik \
  --barcode_start 10 \
  --barcode_length 36 \
  --strategy 4_mer_gpu_v4 \
  --distance_measure SEQUENCE_LEVENSHTEIN \
  --rejection_threshold 7 \
  --data_mode real \
  --outdir "$OUT" \
  -profile local,singularity
```

### Simulated Dual full density, `r = 7`

```bash
nextflow run main.nf \
  --tool quik \
  --barcode_file "$LAYOUT/decoys/output/dual/dual_whitelist_real_plus_decoys.txt" \
  --r1_fastq "$FASTQ/1M_high_dual_R1.fastq" \
  --r2_fastq "$FASTQ/1M_high_dual_R2.fastq" \
  --sample_id quik_dual_combined_1M \
  --barcode_start 0 \
  --barcode_length 36 \
  --strategy 4_mer_gpu_v4 \
  --distance_measure SEQUENCE_LEVENSHTEIN \
  --rejection_threshold 7 \
  --data_mode real \
  --outdir results/sim_dual_r7 \
  -profile local,singularity
```

(`data_mode real` here means “no ground-truth scoring inside the pipeline”;
the simulated answers file is optional.)

### Params JSON

Example JSON files are in [`params/`](params/). Edit absolute paths, then:

```bash
nextflow run main.nf \
  -params-file /path/to/BarcodeLayoutOptimization/docs/decode/params/dual_full_real_r7.json \
  -profile local,singularity
```

On SLURM:

```bash
nextflow run main.nf \
  -params-file .../dual_full_real_r7.json \
  -profile slurm
```

---

## 6. Full matrix (what the paper ran)

For each layout × threshold below, run QUIK once against the matching
whitelist.

| Layout | FASTQ | Whitelist | `r` values | `barcode_start` |
|--------|-------|-----------|------------|-----------------|
| Dual full | `26_018_S1` | `decoys/output/dual/...` | 5, 6, 7 | 10 |
| Greedy | `26_018_S3` | `decoys/output/greedy/...` | 5, 6, 7 | 10 |
| Zipcodes | `26_018_S4` | `decoys/output/zipcodes_balanced_khalf/...` | 5, 6, 7 | 10 |
| Dual half | `26_022_S1` | `decoys/output/dual_1in2/...` | 6, 7, 8 | 10 |
| Dual full sim | `1M_high_dual_*` | dual whitelist | 5, 6, 7 | 0 |
| Greedy sim | `1M_high_greedy_*` | greedy whitelist | 5, 6, 7 | 0 |
| Zipcodes sim | `1M_high_zip_*` | zipcodes_balanced_khalf whitelist | 5, 6, 7 | 0 |
| Dual half sim | `1M_med_dual_checker_*` | dual_1in2 whitelist | 6, 7, 8 | 0 |

---

## 7. Outputs

Under `--outdir`, QUIK publishes (names follow `--sample_id`):

- `{sample_id}_R1_filtered.fastq` / `_R2_filtered.fastq` — assigned pairs
- `{sample_id}_barcode_calling_stats.txt` — QUIK log (strategy, `r`, rates)
- `{sample_id}_per_barcode.csv` — counts per whitelist entry
- `{sample_id}_stats_summary.csv` / `_stats_report.txt` — assignment summary

Downstream paper tables (assignment / decoy rates, spatial maps) are built from
these filtered FASTQs and per-barcode counts. See the parent [`REPRODUCE.md`](../../REPRODUCE.md).

---

## 8. Notes

- Do **not** pass `--tool randombarcodes` or `--tool columba` for manuscript
  reproduction.
- Do **not** use real barcodes alone if you need decoy-based precision; always
  use the `*_whitelist_real_plus_decoys*.txt` files above.
- Gzipped FASTQs are fine if Nextflow/QUIK on your profile accepts them; if a
  run fails on `.gz`, gunzip first.
- Full pipeline docs: [BarcodeCalling_benchmark README](https://github.com/OncoRNALab/BarcodeCalling_benchmark#readme)
  and [`docs/REAL_DATA_GUIDE.md`](https://github.com/OncoRNALab/BarcodeCalling_benchmark/blob/main/docs/REAL_DATA_GUIDE.md).
