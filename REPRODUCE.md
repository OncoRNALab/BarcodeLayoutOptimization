# Reproduce paper results

This file describes how to regenerate every file needed to reproduce the
BMC Bioinformatics paper results from this `manuscript` branch.

---

## Workflow overview

End-to-end path from raw inputs to paper metrics. For most reproduction work,
follow the four steps below using the **deposited** layouts (this GitHub
branch) and FASTQs (Zenodo). Optional from-scratch generation is documented
afterward.

```text
┌─────────────────────────────────────────────────────────────────────────┐
│  1. Layouts + decoys                                                     │
│     → this repo (committed)  OR  regenerate with scripts here (§1–3)     │
├─────────────────────────────────────────────────────────────────────────┤
│  2. FASTQ data                                                           │
│     Real:     Zenodo                                                     │
│     Simulated: Zenodo  OR  generate with BarcodeCalling_benchmark (§4)   │
├─────────────────────────────────────────────────────────────────────────┤
│  3. Decode                                                               │
│     QUIK only, via BarcodeCalling_benchmark (--tool quik) (§5)           │
├─────────────────────────────────────────────────────────────────────────┤
│  4. Analyse                                                              │
│     Python scripts in analyses/ (run_all.sh) (§6)                        │
└─────────────────────────────────────────────────────────────────────────┘
```

### 1. Get the layouts (and decoys)

| Option | What to do |
|--------|------------|
| **Use paper files (recommended)** | Clone this `manuscript` branch. Printed layouts are in `data/1in1/` and `data/1in2/`; structure-matched decoy whitelists are in `decoys/output/`. |
| **Generate from scratch** | Run the Dual / Zipcodes generators and decoy scripts (§2–3). Greedy + 2-OPT needs a 150k barcode pool then [BarLay](https://github.com/Jatzkowski/BarLay) (§2). |

Layouts are versioned in **this GitHub repository**, not on Zenodo.

### 2. Get the FASTQ data

| Data | Option A (recommended) | Option B |
|------|------------------------|----------|
| **Real sequencing** | Download from Zenodo ([10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834)) | — |
| **Simulated reads** | Download the same Zenodo deposit | Generate from the printed barcode `.txt` lists using [BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark) (`bin/press/`: simulate → paired FASTQ; §4.2) |

Filenames: [`docs/ZENODO.md`](docs/ZENODO.md).

### 3. Decode FASTQs with QUIK

Use the Nextflow pipeline from the benchmarking paper:

[OncoRNALab/BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)

- Set **`--tool quik`** only (do not run the other decoders).
- Whitelist = real barcodes + layout-matched decoys (`decoys/output/…`).
- Real reads: `barcode_start=10`; simulated: `barcode_start=0`.
- Full density: `r ∈ {5,6,7}`; Dual half density: `r ∈ {6,7,8}`.

Step-by-step commands and example params: §5 and [`docs/decode/README.md`](docs/decode/README.md).

### 4. Run the analysis scripts

After decoding, compute paper tables and figures with the Python analyses in
this branch:

```bash
cd analyses
pip install -r requirements.txt   # once
./run_all.sh                      # or: ./run_all.sh python3
```

`run_all.sh` regenerates decoy QC where needed, then runs the notebooks and
scripts for precision calibration, real-array rates, decoy quality, recall /
sequence features, and related summaries. Individual folders under `analyses/`
can also be run on their own. Details: §6.

---

## Generating data from scratch (optional)

The numbered subsections below explain how to **rebuild** layouts, decoys,
simulated FASTQs, and related inputs from generator scripts. They are for
parameter changes, verification, or full regeneration.

**We recommend skipping this section for standard reproduction.** Use the
files already deposited instead:

| Input | Where to get it |
|-------|-----------------|
| Printed layouts + decoy whitelists | This GitHub repository, branch `manuscript` (`data/`, `decoys/output/`) |
| Real and simulated FASTQs | Zenodo ([10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834)); see [`docs/ZENODO.md`](docs/ZENODO.md) |

Then decode with QUIK (§5) and run `analyses/run_all.sh` (§6).

---

### 1. Printed layout files (ground truth — already committed)

| Method | File | n barcodes |
|--------|------|------------|
| Dual 4×4 (full density) | `data/1in1/layout_235x363_dual_4x4.csv` | 85,305 |
| Dual 4×4 (half density) | `data/1in2/layout_1in2_dual_4x4.csv` | 42,471 |
| Zipcodes (full density) | `data/1in1/layout_235x363_zipcodes_nosep.csv` | 85,305 |
| Greedy + 2-OPT | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` | 85,305 |

`.txt` versions (one barcode per line, column-major): same directory, same base name.

---

### 2. Regenerate printed layouts from configs

#### Dual (full density)

```bash
cd Dual_barcodes
python generate_layout_235x363.py config_235x363.json
# → layout_235x363_dual_4x4.csv (seed inside config; identical to data/1in1/ copy)
```

Config: `Dual_barcodes/config_235x363.json`

#### Dual (half density)

```bash
cd Dual_barcodes
python generate_checkerboard_235x363.py
# → data/1in2/layout_1in2_dual_4x4.csv
```

#### Zipcodes (full density)

```bash
cd Zipcodes
python generate_layout_235x363_zipcodes.py config_235x363_nosep.json output_235x363_nosep
# → output_235x363_nosep/layout_235x363_zipcodes_nosep_colmajor.csv
#   (identical to data/1in1/layout_235x363_zipcodes_nosep.csv)
```

Config: `Zipcodes/config_235x363_nosep.json`  
Seeds: row + column axes both use `random_seed: 42`.

#### Greedy + 2-OPT

##### Step A — candidate pool (input to BarLay)

Before placement, a **candidate pool of 150,000 unique 36-nt barcodes** was
built with the barcode-generation step in
[BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)
(`bin/press/scripts/1_BarcodeGen_split.py` /
[`BENCHMARK_README.md`](https://github.com/OncoRNALab/BarcodeCalling_benchmark/blob/main/bin/press/BENCHMARK_README.md));
set `N: 150000`, `M: 36`, `homomax: 3`, and a G limit of 25% to match the
manuscript wording.

A **second, independent** pool of 150,000 barcodes was generated with the same
procedure for Greedy decoys (not used as BarLay input). Decoys were subsampled
from that pool; the committed source used by `decoys/scripts/make_greedy_decoys.py`
is `data/decoy_layouts/barcodes/1in1/greedy_decoybar_1in1_column_major.txt`
(85,305 sequences from that second pool).

##### Step B — placement with BarLay

The first 150k pool was the input library for
[BarLay](https://github.com/Jatzkowski/BarLay) (external repository; commit
`e22a877603` as cited in the paper). BarLay ran greedy placement + 2-OPT local
search with default settings on NVIDIA A100 GPUs at the VSC and placed
**85,305** barcodes on the 235 × 363 grid. The printed layout is committed in
`data/1in1/`:

- `layout_235x363_randombarcodes_opt_coords.csv`
- `layout_235x363_randombarcodes_opt.txt`

---

### 3. Regenerate structure-matched decoys

All decoys are regenerated by:

```bash
cd decoys/scripts
bash run_all.sh
# writes to decoys/output/ (overwrites committed copies)
# QC reports to decoys/qc/
```

The regenerated files are **byte-identical** to the committed copies (verified Oct 2026).

#### Individual scripts

| Set | Script | φ used |
|-----|--------|--------|
| `dual/` | `make_dual_decoys.py --layout data/1in1/layout_235x363_dual_4x4.csv …` | 1/3 (global) |
| `dual_1in2/` | Same rule on half-density layout | 1/3 |
| `zipcodes_balanced_khalf/` | `make_zip_decoys.py … --spacing bresenham` | 1/3 (local, per axis) |
| `greedy/` | `make_greedy_decoys.py …` | 1/3 (global) |

See `decoys/scripts/run_all.sh` for exact arguments.

---

### 4. Simulated and real FASTQs

#### 4.1 Pre-generated FASTQs (recommended)

**Simulated** and **real sequencing** FASTQs used in the paper are on Zenodo
([10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834)).
Filenames: [`docs/ZENODO.md`](docs/ZENODO.md).

| | Full-density Dual / Greedy / Zipcodes | Dual half density |
|--|---------------------------------------|-------------------|
| Error level (benchmark naming) | **high**: sub 0.03, ins 0.15, del 0.14 | **medium-like**: sub 0.04, ins 0.05, del 0.14 |
| Reads per layout | 1,000,000 paired-end | 1,000,000 paired-end |
| Source barcodes | printed layouts only (no decoys) | same |
| Form | bare 36-nt sequences (no handles) | same |

Real sequencing: demultiplexed MiSeq FASTQs for Dual, Greedy + 2-OPT,
Zipcodes (full density) and Dual half density (checkerboard).

**This `manuscript` branch does not ship simulation code.** For paper-identical
inputs, download Zenodo. To regenerate yourself, use the companion benchmarking
repository (next subsection).

> The `simulation/` directory on branch `main` is for separate NN≈8 redesign
> experiments only. Do **not** use it for these manuscript reads.

#### 4.2 Regenerate simulated FASTQs
([BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark))

All read simulation and FASTQ conversion for this paper uses the same tooling
documented in the benchmarking repository under `bin/press/`:

| Step | Script | Docs in that repo |
|------|--------|-------------------|
| 1. Simulate reads (plain text + answer indices) | `bin/press/scripts/2_SimGen_split.py` | [`bin/press/BENCHMARK_README.md`](https://github.com/OncoRNALab/BarcodeCalling_benchmark/blob/main/bin/press/BENCHMARK_README.md) |
| 2. Build paired FASTQ (R1 + R2) | `bin/press/scripts/3_Generate_FASTQ_pairs.py` | [`bin/press/FASTQ_GENERATION_README.md`](https://github.com/OncoRNALab/BarcodeCalling_benchmark/blob/main/bin/press/FASTQ_GENERATION_README.md) |

Follow those two READMEs for dependencies, YAML configs, and job scripts. Below
is only what differs for **this** paper’s layouts.

#### Clone the benchmarking pipeline

```bash
git clone https://github.com/OncoRNALab/BarcodeCalling_benchmark.git
cd BarcodeCalling_benchmark
```

#### Barcode libraries (inputs from this branch)

Use the printed barcode lists (one sequence per line), **not** decoy whitelists:

| Layout | File in this repo |
|--------|-------------------|
| Dual full | `data/1in1/layout_235x363_dual_4x4.txt` |
| Zipcodes | `data/1in1/layout_235x363_zipcodes_nosep.txt` |
| Greedy + 2-OPT | `data/1in1/layout_235x363_randombarcodes_opt.txt` |
| Dual half | `data/1in2/layout_1in2_dual_4x4.txt` |

Point the simulation YAML `filename_1` (barcode library) at the matching
`.txt` path, with `N` = number of lines and `M: 36`. Simulation loads the
companion pickle (`filename_1 + ".pkl"`) produced when the library is prepared
in the benchmarking workflow — see `BENCHMARK_README.md` (barcode generation /
library pickle step).

#### Error rates and read count

Use the same levels documented in
[`bin/press/BENCHMARK_README.md`](https://github.com/OncoRNALab/BarcodeCalling_benchmark/blob/main/bin/press/BENCHMARK_README.md)
(`Q: 1000000`, `nave: 4`):

| Paper set | Config rates (`srate` / `irate` / `drate`) |
|-----------|---------------------------------------------|
| Dual / Greedy / Zipcodes (full density) | **high** — `0.03` / `0.15` / `0.14` |
| Dual half density | `0.04` / `0.05` / `0.14` (edit the medium template if needed) |

Example high-error YAML fields (paths are illustrative):

```yaml
N: 85305
M: 36
filename_1: "/path/to/BarcodeLayoutOptimization/data/1in1/layout_235x363_dual_4x4.txt"
Q: 1000000
nave: 4
srate: 0.03
irate: 0.15
drate: 0.14
readsfilename: "/path/to/out/1M_high_dual"
answersfilename: "/path/to/out/1M_high_dual_answers"
```

#### Step 1 — simulate reads

From the benchmarking repo (after the library pickle exists):

```bash
python bin/press/scripts/2_SimGen_split.py \
  --config /path/to/your_sim_config.yaml
```

Writes a plain-text read list (`readsfilename`) and answer indices
(`answersfilename`). Details: `bin/press/BENCHMARK_README.md`.

#### Step 2 — convert to paired FASTQ

```bash
python bin/press/scripts/3_Generate_FASTQ_pairs.py \
  --reads /path/to/out/1M_high_dual \
  --answers /path/to/out/1M_high_dual_answers \
  --output-r1 /path/to/out/1M_high_dual_R1.fastq \
  --output-r2 /path/to/out/1M_high_dual_R2.fastq
```

Or pass the same YAML with `--config` (see
`bin/press/FASTQ_GENERATION_README.md`). That step produces:

- **R1** — simulated read
- **R2** — reverse complement of R1  
Headers follow `@read_<index>_barcode_<barcode_index>`.

These FASTQs are what QUIK consumes in §5 (`barcode_start: 0`).

Regenerated files match the documented rates and workflow but are not
guaranteed byte-identical to Zenodo (seeds / sampling stream not pinned here).
Prefer Zenodo for manuscript-identical inputs.

---

## 5. Decode FASTQs with QUIK

Paper decoding used **QUIK only**, via the companion Nextflow pipeline
[BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)
(`--tool quik`; do not run RandomBarcodes or Columba for manuscript
reproduction).

Full step-by-step instructions, the layout × threshold matrix, and example
params JSON files:

→ **[`docs/decode/README.md`](docs/decode/README.md)**

Summary:

```bash
git clone https://github.com/OncoRNALab/BarcodeCalling_benchmark.git
cd BarcodeCalling_benchmark

# Real Dual full density, r = 7 (edit paths)
nextflow run main.nf \
  --tool quik \
  --barcode_file /path/to/BarcodeLayoutOptimization/decoys/output/dual/dual_whitelist_real_plus_decoys.txt \
  --r1_fastq /path/to/zenodo/26_018_S1_L001_R1_001.fastq.gz \
  --r2_fastq /path/to/zenodo/26_018_S1_L001_R2_001.fastq.gz \
  --sample_id 26_018_S1_quik \
  --barcode_start 10 \
  --barcode_length 36 \
  --strategy 4_mer_gpu_v4 \
  --distance_measure SEQUENCE_LEVENSHTEIN \
  --rejection_threshold 7 \
  --data_mode real \
  --outdir results/dual_full_r7 \
  -profile local,singularity
```

| | Real sequencing | Simulated |
|--|-----------------|-----------|
| `barcode_start` | `10` | `0` |
| Whitelist | `decoys/output/*/…whitelist_real_plus_decoys*.txt` | same |
| `r` full density | `{5, 6, 7}` | `{5, 6, 7}` |
| `r` Dual half | `{6, 7, 8}` | `{6, 7, 8}` |

Or use a params file from `docs/decode/params/` (replace `/ABS/PATH/...`):

```bash
nextflow run main.nf \
  -params-file /path/to/BarcodeLayoutOptimization/docs/decode/params/dual_full_real_r7.json \
  -profile local,singularity
```

---

## 6. Analysis (Python scripts)

After QUIK decoding (§5), compute paper rates and figure inputs with the
scripts under [`analyses/`](analyses/):

```bash
cd analyses
pip install -r requirements.txt
./run_all.sh                 # or: ./run_all.sh python3
```

**Expected inputs**

| Input | Location |
|-------|----------|
| Layouts | `data/1in1/`, `data/1in2/` |
| Decoy whitelists | `decoys/output/` |
| QUIK decode outputs | `data/quik_out_manuscript/` (place §5 `outdir` results here, using the folder names expected by the notebooks — see [`analyses/README.md`](analyses/README.md)) |

`run_all.sh` runs precision calibration, real-array metrics, decoy-quality
checks, recall / sequence-feature scripts, and related notebooks (in that
order). Individual folders under `analyses/` can be run alone. Figure builders
live in `analyses/figures_manuscript/`.

---

## 7. File checklist

Tick before submitting or sharing:

- [ ] `data/1in1/layout_235x363_dual_4x4.csv` present
- [ ] `data/1in1/layout_235x363_zipcodes_nosep.csv` present
- [ ] `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` present
- [ ] `data/1in2/layout_1in2_dual_4x4.csv` present
- [ ] `decoys/output/dual/dual_decoys.txt` and `dual_whitelist_real_plus_decoys.txt` present
- [ ] `decoys/output/dual_1in2/dual_decoys.txt` present
- [ ] `decoys/output/zipcodes_balanced_khalf/zip_decoysv2.txt` and `zip_whitelist_real_plus_decoysv2.txt` present
- [ ] `decoys/output/greedy/greedy_decoys.txt` present
- [ ] `decoys/scripts/run_all.sh` runs successfully on a fresh clone
- [ ] Zenodo FASTQs downloaded; QUIK decode steps in `docs/decode/` runnable
- [ ] `analyses/run_all.sh` present; QUIK outputs available under `data/quik_out_manuscript/`
- [ ] Zenodo DOI listed in `docs/ZENODO.md`

---

## 8. Key parameters

| Parameter | Value |
|-----------|-------|
| Array size | 235 × 363 = 85,305 spots |
| Barcode length | 36 nt |
| Synthesis order | A > C > G > T |
| Bit mapping | 00→A, 01→C, 10→T, 11→G |
| Dual block size | 4 × 4 |
| Dual random seed | 42 |
| Greedy input pool | 150,000 unique 36-nt barcodes (homopolymer ≤ 3, G ≤ 25%) |
| Zipcodes random seed | 42 |
| Zipcodes soft targets (bits) | 2, 4, 6, 8 at offsets 1–4 |
| Zipcodes hard min (bits) | ≥ 8 at offset ≥ 5 |
| Decoy φ (Dual, Greedy) | 1/3 (global) |
| Decoy φ (Zipcodes) | 1/3 (local, per axis) |
| QUIK strategy | `4_mer_gpu_v4` |
| QUIK threshold r | {5, 6, 7} full density; {6, 7, 8} half density |
| QUIK `barcode_start` | 10 (real); 0 (simulated) |

---

## 9. Manuscript figures → scripts

Main-text figures are written as PDFs under `figures/` (submission package
renames them to `Figure1.pdf` … `Figure9.pdf`). Run the analysis pipeline
first (§6) so QUIK-derived tables exist, then:

```bash
# Figures 1–4, 6–7
python analyses/figures_manuscript/make_results_figures.py

# Figure 5
python analyses/library_complexity/library_complexity.py \
  --data-dir data/1in1 --out-dir figures

# Figures 8–9
python analyses/figures_manuscript/make_checkerboard_figures.py
```

| Figure | Label | Output PDF | Script |
|--------|-------|------------|--------|
| **1** | `fig:masks` | `figures/masks.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`masks()`) |
| **2** | `fig:decoy_distances` | `figures/decoy_distances.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`decoy_distances()`) |
| **3** | `fig:precision_calibration` | `figures/precision_calibration.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`calibration()`) |
| **4** | `fig:assignment` | `figures/assignment_vs_precision.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`assignment()`) |
| **5** | `fig:complexity` | `figures/library_complexity.pdf` | `analyses/library_complexity/library_complexity.py` |
| **6** | `fig:spatial` | `figures/spatial_reads.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`spatial()`) |
| **7** | `fig:gcontent` | `figures/g_content.pdf` | `analyses/figures_manuscript/make_results_figures.py` (`g_content()`) |
| **8** | `fig:checkerboard_spatial` | `figures/checkerboard_spatial.pdf` | `analyses/figures_manuscript/make_checkerboard_figures.py` |
| **9** | `fig:checkerboard_assignment` | `figures/checkerboard_assignment.pdf` | `analyses/figures_manuscript/make_checkerboard_figures.py` |

Related scripts (not main-text figures in the current manuscript):

| Output | Script | Notes |
|--------|--------|-------|
| `figures/encoding_schematic.pdf` | `analyses/figures_manuscript/make_encoding_schematic.py` | Dual / Zipcodes encoding schematic |
| `figures/greedy_principle.pdf` | `analyses/figures_manuscript/make_greedy_schematic.py` | Greedy 8-neighbour schematic |
| `figures/dual_full_vs_half_all_reads_assignment_vs_precision.pdf` | `analyses/figures_manuscript/make_dual_density_all_reads.py` | All-reads Dual full vs half comparison |
