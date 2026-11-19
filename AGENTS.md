# AGENTS.md — barcode_design

Read this first. It is the entry point for any agent or new collaborator working in this repository.

## Branch policy

| Branch | Purpose |
|--------|---------|
| **`main`** | Living toolkit. Four NN≈8 redesigned arrays are the current designs. Use for all new work. |
| **`manuscript`** | **Frozen.** Printed layouts (Dual, Zipcodes, Greedy) + structure-matched decoys to reproduce the BMC Bioinformatics paper. See `REPRODUCE.md` on that branch. Do not add new design work here. |

Large files (simulated FASTQs, oversized intermediates) are **not in git**. See [`docs/ZENODO.md`](docs/ZENODO.md) for what is deposited on Zenodo and how to obtain it.

---

## What this project is

This repository designs **spatial barcode layouts** for photolithographic (maskless, MAS 2.0) spatial-transcriptomics arrays. It also generates the **decoy whitelists** used to estimate decoding precision on printed arrays.

- **Array:** full density is 235 rows × 363 columns, giving 85,305 features of 13.7 µm with 1 µm gaps. Half density is a checkerboard of 42,471 printed features.
- **Barcode:** every feature carries a 36-nt barcode, read from Read 1 at positions 11–46, right after the 10-nt UMI.
- **Synthesis model:** probes are synthesized 5′→3′ with the fixed printing (deposition) order **A > C > G > T**, repeated every 4 cycles, and each barcode is embedded left-most. Features that share a barcode prefix are illuminated (co-activated) together, which can reduce stray-light insertions.
- **Decoding:** QUIK, using sequence-Levenshtein distance and the 4-mer single-stage filter, with rejection threshold r ∈ {5, 6, 7}, or {6, 7, 8} at half density.

The work is published in a BMC Bioinformatics manuscript, which lives outside this repo at `~/Documents/phd_thesis_revision/layoutopt_manuscript/`. `docs/methods.md` is an excerpt of its Methods and is the **authoritative description** of the printed designs and decoys. If code or older READMEs disagree with `docs/methods.md`, trust `docs/methods.md`, then the data, then the code.

## Current four NN≈8 arrays (`main`)

The active designs on `main`. Full documentation in [`docs/NewFourArrays_Summary.md`](docs/NewFourArrays_Summary.md).

| Array | Layout | Real NN μ |
|-------|--------|-----------|
| Dual shared_pool minH6 | `Dual_barcodes/layout_235x363_dual_shared_pool_g10_minH6.csv` | **8.6** |
| Dual fixed_bar2 8×8 | `Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10_nn8.csv` | 7.0 |
| Zipcodes flip | `Zipcodes/output_235x363_flip_g10_h5_d12-16_far10/` | 7.0 |
| HSDB colored | `HSDB/outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored/` | 8.0 |

## Printed layouts (`manuscript` branch / `data/1in1/`)

Three layouts were printed and described in the BMC paper. Their ground-truth files are in `data/1in1/` on **both** branches. The `manuscript` branch also carries their structure-matched decoys. See `docs/methods.md` for the authoritative description.

| Method | File | Real NN μ |
|--------|------|-----------|
| Dual 4×4 | `data/1in1/layout_235x363_dual_4x4.csv` | 10.0 |
| Zipcodes | `data/1in1/layout_235x363_zipcodes_nosep.csv` | 6.6 |
| Greedy | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` | 14.5 |
| Dual half-density | `data/1in2/layout_1in2_dual_4x4.csv` | — |

## Methods in this repository

All methods are first-class, and none are archived. Three were printed and appear in the manuscript. The others are kept for future work.

| Folder | Method | Status |
|---|---|---|
| `Dual_barcodes/` | Dual barcodes: 4×4 blocks share bar_1 | Printed at full and half density (manuscript); redesigned with minH6 and fixed_bar2 on `main` |
| `Zipcodes/` | Zipcodes: [column code][row code] | Printed at full density (manuscript); flip-mode NN≈8 redesign on `main` |
| *(not in this repo)* | Greedy + 2-OPT (BarLay, GPU; third-party baseline) | Printed at full density (manuscript). The algorithm belongs to a **separate repository** that was not developed in this project; only its printed layout and decoy pool are kept here |
| `HSDB/` | Hierarchical Spatial Dual Barcodes (24-nt Zipcode bar_1 + 12-nt bar_2) | NN≈8 colored design on `main`; not printed |
| `Quadtree/` | Hierarchical quadtree + RLL, hybrid 24 nt + 12 nt | Layouts exist; printed variant `data/1in1/layout_235x363_quadtree_8x3_12nt.*`; **excluded from the manuscript** |
| `HSDB/` | Hierarchical Spatial Dual Barcodes (24-nt Zipcode bar_1 + 12-nt bar_2 = 36 nt) | Designed, not printed. bar_1 uses flip-mode Zipcodes with G/junction enforced in the candidate loop (see below) |
| `inhouse_barcodes/` | Random barcode libraries (filtered random 36-nt pool and its 12/18-nt prefixes and suffixes) | Input for Dual, Quadtree and decoys |
| `decoys/` | Structure-matched decoy generation + QC | Used for every manuscript precision estimate |
| `figures/` | Early 118×182 schematic figures | Historical |

### Dual barcodes (printed + redesigned)
- **Layout:** blocks of size \(B \times B\) (printed: **4×4**; 59 × 91 = 5,369 blocks on 235×363, edge blocks truncated). Barcode = bar_1 + bar_2.
- **Modes** (`Dual_barcodes/generate_layout_235x363.py`):
  - **`shared_pool` (legacy / printed):** one predefined half library for both halves (`inhouse_barcodes/bar_1in1_column_major_first18.txt`). No combined G/HP filter. Config `Dual_barcodes/config_235x363.json`.
  - **`fixed_bar2` (new):** `n_bar2 = B²` distance-designed repertoire (min Hamming), reused in every block; bar_1 from file or generated and filtered for compatibility with the bar_2 set; optional `sequence_qc` on combined barcodes. Config `Dual_barcodes/config_235x363_fixed_bar2_g10.json` (G≤10%, HP≤5, 16×18-nt bar_2 with min Hamming 8).
- **Assignment:** unique bar_1 per block; distinct bar_2 within a block; global uniqueness is the pair `(bar_1, bar_2)`.
- **Block size:** `grid_configuration.block_size` is passed into `DualBarcodeLayout` (no longer hardcoded-only).
- **Half density (printed):** first 117 rows checkerboard; see `generate_checkerboard_235x363.py`.

### Zipcodes (printed)
- **Barcode:** [column code 18 nt][row code 18 nt], with no separator and no random tail.
- **Axis codes:** each is a **36-bit codeword**, converted at 2 bits per nt with **00→A, 01→C, 10→T, 11→G**.
- **Generation:** a randomized search with backtracking (`Zipcodes/spatial_zipcodes/axis_code_generator.py`, seed 42), testing up to 10,000 candidates per position.
  - **Soft targets:** bit distance to the codes 1, 2, 3 and 4 positions back should be 2, 4, 6 and 8 bits. A candidate is accepted if the mean squared deviation is ≤ 16.
  - **Hard constraint:** ≥ 8 bits from every code 5 or more positions back.
  - **Backtracking:** up to 5 positions.
- **Measured distances:** adjacent column codes differ by 8.8 bits on average (range 4–12), 10.5 bits at offset 4 and 17.9 bits beyond. The targets are **not reached**.
  - **Why:** every candidate is a fresh random 36-bit code, about 18 ± 3 bits from any other, so a 2-bit neighbour is practically never drawn.
  - **Hidden failure:** when nothing met the tolerance, the generator silently kept the best of the 10,000 draws. Only 19 of 362 printed positions (5%) are within tolerance.
- **Code:** `Zipcodes/generate_layout_235x363_zipcodes.py Zipcodes/config_235x363_nosep.json` (a copy of the config is in `data/1in1/config_files/zip_config_235x363_nosep.json`).

### Zipcodes candidate modes (`axis_code` section of the config)
- `candidate_mode`:
  - `"random"` (the default, used for the printed array) keeps `random_binary_code`.
  - `"flip"` builds each candidate from the preceding code with `flip_neighbour_code`, flipping a number of bits drawn from `flip_counts` (default `[1, 2, 3]`). Position 0 is still random.
  - In flip mode the first candidate that hits all four targets exactly is taken; otherwise the best-scoring one.
- `require_tolerance: true` backtracks instead of keeping an out-of-tolerance code.
- `positions_outside_tolerance` is now reported by `verify_constraints` and printed by the layout scripts.
- **Composition is enforced inside the candidate loop**, in both modes, before distances are scored:
  - the G cap per 18-nt axis code (`sequence_qc.g_content_max`);
  - the homopolymer limit (`max_homopolymer_length`), including the run across the column/row junction (`reject_sequences_with_long_homopolymers: true`).
  - The GC range (`gc_content_min` / `gc_content_max`) is still only reported, not enforced. About half the flip-design axis codes fall outside 0.4–0.6.
- **Run:** `cd Zipcodes && python generate_layout_235x363_zipcodes.py <config.json> <output_dir>` (about 1 minute in flip mode).

### Zipcodes flip-mode designs (235×363, not printed yet)

Both designs use seed 42, the printed bit mapping (00→A, 01→C, 10→T, 11→G), `require_tolerance: true` and homopolymers ≤ 3 nt including the junction. All 85,305 barcodes are unique in each.

| | `flip` | `flip_g15_d6-12` |
|---|---|---|
| Config | `Zipcodes/config_235x363_flip.json` | `Zipcodes/config_235x363_flip_g15_d6-12.json` |
| Output | `Zipcodes/output_235x363_flip/` | `Zipcodes/output_235x363_flip_g15_d6-12/` |
| Targets, offsets 1–4 (bits) | 2, 4, 6, 8 | 6, 8, 10, 12 |
| Min distance ≥ 5 apart (bits) | 8 | 12 (raised so distant codes are never closer than offset 4) |
| `flip_counts` | 1, 2, 3 | 4–8 |
| G cap per axis code | 25% (max reached 22%) | 15%, i.e. ≤ 2 G per 18 nt (max reached 11%, mean 9%) |
| Achieved, offsets 1–4 (bits) | exactly 2, 4, 6, 8 everywhere | 5.9, 8.0, 10.0, 12.0 on average; 93–100% exactly on target, rest within ±2 |
| Achieved, offset 5 / offset 20 (bits) | 9.7 / 17.6 | 13.5 / 17.6 |
| Bases that differ, offset 1 | — | mean 5.4 (3–8), mostly 5–6 of 18 |
| Positions outside tolerance / backtracks | 0 / 0 | 0 / 0 |
| Row codes equal to a column code | 12 | 0 |

Distances are in **bits**, and each base holds 2 bits. A bit distance d changes between d/2 bases (both bits of the same bases) and d bases (one bit per base), usually close to d. To control the difference in bases directly, the soft-constraint score would have to be computed on the DNA sequences instead.

### HSDB (Hierarchical Spatial Dual Barcodes)
- **Structure:** bar_1 (24 nt = 12 nt column + 12 nt row Zipcodes on the **block** grid) + bar_2 (12 nt spot ID, 16 codes reused in every 4×4 block) = 36 nt. (Older 28+8 outputs remain in the repo for comparison.)
- **bar_1 generator:** `HSDB/bar1_zipcode/generate_block_barcodes.py`. It calls `generate_column_and_row_codes` from `Zipcodes/`.
  - **Was broken:** candidates were random bit strings; G ≤ 25% and homopolymers were applied afterwards by rewriting finished barcodes (G→T, forced substitutions). That broke the X+Y structure and the distance gradient (adjacent DNA mean 6.6, far min 1).
  - **Fixed:** `candidate_mode: "flip"`, `require_tolerance: true`, G and junction homopolymers rejected in the candidate loop. No post-hoc rewrite. QC asserts uniqueness and intact X+Y structure.
- **bar_2:** `bar2_library/generate_spot_barcodes.py` (default 12 nt). G, GC, homopolymer and leading-run ≤ 1 are enforced during selection.
- **Configs (all flip mode, printed bit mapping):**

| Array | Config | Block grid | Axis bits | Targets (bits) | Far min (bits) |
|---|---|---|---|---|---|
| 118×182 | `block_zipcode_config.json` | 30×46 | 28 (legacy) | 4, 8, 12, 16 | 12 (was 16; 16 failed under the G cap) |
| 235×363 | `block_zipcode_config_235x363.json` | 59×91 | 28 (legacy) | 2, 4, 6, 8 | 8 |
| 235×363 zero-G 24+12 | `block_zipcode_config_235x363_g0_d4-10_far8_24nt.json` | 59×91 | **24** | 4, 6, 8, 10 | **8** (10 failed under zero G) |
| 253×363 | `block_zipcode_config_253x363.json` | 64×91 | 28 (legacy) | 2, 4, 6, 8 | 8 |
| 295×363 extra | `block_zipcode_config_253x363_extra.json` | 85×91 | 28 (legacy) | 2, 4, 6, 8 | 8 |

- **bar_1/bar_2 junction:** bar_2 leading runs are ≤ 1 nt. Row (Y) codes are generated with `downstream_start_runs` so their trailing run plus a bar_2 start stays ≤ the homopolymer limit. The layout assembler only assigns junction-safe bar_2 codes and checks the combined homopolymer.
- **Current zero-G 235×363 (24+12, targets 4/6/8/10, far 8):** `HSDB/bar1_zipcode/235x363_g0_d4-10_far8_24+12/` and `HSDB/outputs/235x363_g0_d4-10_far8_24+12/`. bar_2: `spot_barcodes_12nt_g0.txt` (min Hamming 5). **No G** anywhere; homopolymers ≤ 5; adjacent DNA mean 3.6 bases; far DNA min 5.
- **Legacy 28+8 zero-G 235×363 (targets 4/6/8/10, far 10):** `HSDB/bar1_zipcode/235x363_g0_d4-10_far10/` and `HSDB/outputs/235x363_g0_d4-10_far10/` with `spot_barcodes_8nt_g0.txt`.
- **Legacy 28+8 G≤25% 235×363 (targets 2/4/6/8):** `HSDB/bar1_zipcode/235x363/` and `HSDB/outputs/235x363/`.
- **118×182 bar_1 (legacy 28 nt):** `HSDB/bar1_zipcode/118x182/`. Soft targets 4/8/12/16 hit exactly; adjacent DNA mean 3.7 bases.

### Greedy + 2-OPT (printed baseline; separate repository)
Greedy + 2-OPT is a third-party method (BarLay, Jatzkowski et al. 2024) used as the baseline. Its code lives in its own repository (local clones: `~/Documents/Barlay_new/BarLay`, `~/Documents/layout_opt/BarLay`) and should not be added to or modified from this repo. This repo only consumes its output: the printed layout in `data/1in1/` and the second pool for decoys.
- **Pool and placement:** a random pool of unique 36-nt barcodes (homopolymers ≤ 3 nt, G ≤ 25%) was placed by BarLay (Jatzkowski et al. 2024) with default settings.
- **Cost:** each barcode is a left-most schedule binary vector, and the cost is the summed schedule Hamming distance over all 8-neighbour pairs.
- **Greedy step:** fills one position at a time with the unused barcode closest to its already placed neighbours.
- **2-OPT step:** batched best-improvement swaps between placed barcodes.
- **Pools:** the pool of 150,000 barcodes was generated twice. One pool fed Greedy; 85,305 barcodes of the second pool (`data/decoy_layouts/barcodes/1in1/greedy_decoybar_1in1_column_major.txt`) are the source of the Greedy decoys.

### Synthesis / co-activation facts (verified from printing masks)
- All barcodes start in the same cycle, right after the UMI.
- **Both co-designed layouts are printed in two synchronized halves.** The printer finishes the first 18-nt half on **every** feature before it starts the second half:
  - Dual barcodes: all bar_1 codes, then all bar_2 codes.
  - Zipcodes: all column codes, then all row codes.
- **Why it matters:** the second half starts in the same cycle everywhere. Features sharing a first half (a Dual block, a Zipcodes column) stay in phase during the first half. Features sharing a row code co-activate during the second half. Without the wait, the second half would start at different cycles and its co-activation would drift towards random.
- **Models must include the wait:** any simulation of masks, cycles or co-activation should model it. Some older analysis code assumes no wait, e.g. `phd_thesis_revision/layoutopt_manuscript/analyses/recall_sequence_features/02_feature_table.py` and the Dual co-activation animation.

## Decoys (`decoys/`)

Precision on printed arrays is estimated as P̂ = 1 − (D/A)/φ:

- **D:** reads assigned to decoys.
- **A:** all assigned reads.
- **φ:** the expected share of misassigned reads that land on decoys.

Decoys must reproduce each layout's local structure; independent random decoys lie too far away and underestimate errors.

| Set (`decoys/output/…`) | Rule | n decoys | φ used |
|---|---|---|---|
| `dual/` | Keep each block's bar_1; add round(0.5 × n_real) extra bar_2 from the same 18-nt library (8 per complete block) | 42,652 | 1/3 (global) |
| `dual_1in2/` | Same per-block rule on the half-density layout | 21,236 | 1/3 |
| `zipcodes_balanced_khalf/` (**the manuscript set**, files `*v2*`) | Virtual axis codes in the gaps (Bresenham): 181 columns, 117 rows. Each is made by bit flips of a flanking real code, distance-matched at offsets 1–4 (4,000 candidates). Pairings are virtual column × all real rows and virtual row × all real columns; **no virtual × virtual** | 85,006 | **1/3 (local, per axis)**, not the global ≈ 1/2 |
| `greedy/` | Random subset of the independent second pool | 42,652 | 1/3 |
| `zipcodes/` | Earlier 1:2 linspace Zipcodes set; **not** used in the manuscript | 42,802 | — |

Run everything (≈ 2 min, needs only `numpy`):

```bash
cd decoys/scripts
bash run_all.sh                 # overwrites decoys/output and decoys/qc
bash run_all.sh /tmp/x/output   # writes elsewhere (QC goes to /tmp/x/qc)
```

The regenerated outputs are byte-identical to the committed ones (verified Oct 2026). The script writes `zip_decoys.txt`; the committed manuscript copy was renamed `zip_decoysv2.txt`. `check_decoys.py` reports uniqueness, overlap with real barcodes, decoy→nearest-real distance, local decoy share among the 10 nearest neighbours, G content and homopolymers.

## Repository map

```
AGENTS.md                 <- this file (read first)
docs/
  methods.md              manuscript Methods excerpt (authoritative for printed designs)
  NewFourArrays_Summary.md  full methods for four NN≈8 designs (main only)
  ZENODO.md               Zenodo deposit checklist and DOI
SummaryNewDesigns.md      quick reference for four NN≈8 arrays (main only)
data/
  1in1/                   printed full-density layouts (.csv with coords, .txt barcodes) + config_files/
  1in2/                   half-density Dual layout (printer x,y)
  decoy_layouts/barcodes/1in1/   old independent decoys + second Greedy pool
decoys/
  scripts/                make_dual_decoys.py, make_zip_decoys.py, make_greedy_decoys.py,
                          check_decoys.py, run_all.sh (manuscript),
                          plot_nn_distance_histograms.py, plot_activation_masks.py (main)
  output/                 generated decoy files
  qc/                     real-barcode lists + QC JSONs
  figures/                NN histograms, activation masks (main)
  README.md               decoy documentation
Dual_barcodes/            Dual barcode generator + configs/layouts (printed + NN≈8 redesigns)
Zipcodes/                 Zipcode generator + configs/outputs (printed + flip-mode redesigns)
HSDB/                     HSDB generator + configs/outputs (NN≈8 colored design)
Quadtree/                 Quadtree method (non-manuscript)
inhouse_barcodes/         random 18-nt input libraries
simulation/               read simulator (main); output/ gitignored → Zenodo
archive/historical/       obsolete notes/scripts kept for reference
figures/                  early 118×182 schematic figures (historical)
```

**Printed layout files** (ground truth; on both `main` and `manuscript`):

| Method | File | Format |
|---|---|---|
| Dual | `data/1in1/layout_235x363_dual_4x4.csv` | `row,col,combined` |
| Zipcodes | `data/1in1/layout_235x363_zipcodes_nosep.csv` | `row,col,barcode` |
| Greedy | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` | `barcode,col,row` (no header) |
| Dual half density | `data/1in2/layout_1in2_dual_4x4.csv` | `row,col,combined,x,y` |

The `.txt` versions list one barcode per line, in column-major order.

## Known pitfalls

1. **Zipcodes are not Gray codes.** `README.md`, `figures/figure_3a_gray_code_principle.py` and parts of `Zipcodes/README.md` describe Gray codes and distances in bases. Both are outdated: the distances are in bits, and the codes come from randomized search. Banners mark the READMEs as historical.
2. **Bit-to-base mapping.**
   - The printed design and `make_zip_decoys.py` use 00→A, 01→C, 10→T, 11→G.
   - `spatial_zipcodes/bit_to_dna.py` defaults to 10→G, 11→T.
   - The printed config (`data/1in1/config_files/zip_config_235x363_nosep.json`) labels its mapping "4-bits-per-base" one-hot with an explicit `1000→A, 0100→C, 0010→T, 0001→G`, which yields 10→T, 11→G.
   - `Zipcodes/config_235x363_nosep.json` sets the same scheme **without** a mapping, so it falls back to the code default (0010→G), i.e. 10→G, 11→T. That is **not** the printed mapping.
   - The flip configs set the printed mapping explicitly. Check before reusing any other config.
3. **Mixed Zipcodes runs.**
   - In `Zipcodes/output_235x363_nosep/`, only `layout_235x363_zipcodes_nosep_colmajor.csv` is the printed layout.
   - The other files there come from a **different run** and share no axis codes with the printed array: `layout_235x363_zipcodes_nosep.csv`, `codebook.csv`, `barcodes_235x363_nosep.txt` and `intermediate/axis_codes.json`.
   - No codebook (bit strings) from the printed run is in the repo. The first 182 printed column codes equal the X codes of `Zipcodes/output_l36_nospacer/codebook.csv` (118×182, seed 42). Recover bits from the printed sequences with the mapping in pitfall 2.
   - Always take printed barcodes from `data/1in1/`.
4. **Zipcodes row codes equal column codes.** The 235 row codes are identical to the first 235 column codes (same seed and settings). As a result, 235 barcodes are tandem repeats, and 27,495 pairs are half-swapped (X_aX_b vs X_bX_a).
5. **No composition limits on the printed co-designed layouts.**
   - Zipcodes: up to 53% G, homopolymers up to 11 nt. This is despite `sequence_qc` in the config, which was reported, not enforced.
   - Dual: up to 42% G.
   - High-G Zipcodes axis codes correlate with low-read rows and columns.
   - HSDB adds G ≤ 25%. The current Zipcodes generator enforces the G cap and the junction homopolymer limit during generation (see *Zipcodes candidate modes*).
6. **2-OPT pool exchanges (open question, concerns the external BarLay code).** The manuscript says 2-OPT also swaps placed barcodes with unplaced pool barcodes. The BarLay paper describes swaps between placed barcodes only. Unresolved.
7. **φ for Zipcodes.** Zipcodes misassignments stay on one axis, so use φ = 1/3 (local), not the global ≈ 1/2.
8. **Absolute cluster paths.** These point at `/user/gent/446/vsc44685/ScratchVO_dir/layout_opt/…`. They appear in the `*.sh` scripts and in `barcode_file` / `library_file` in the JSON configs. Adapt them before running locally.
9. **Not tracked in git (see `.gitignore`).** Simulated reads (`*.fastq`, `simulation/output/`, ~1.9 GB) and SLURM logs are gitignored. Pre-generated reads are on Zenodo (see `docs/ZENODO.md`).
10. **`generate_layout_235x363_zipcodes.py` exits with code 1 even on success.**
    - It prints `Warning: Analysis error: 'distance_stats'`, then ends with `name 'analysis' is not defined`.
    - Both are bugs in the summary printout, present since the initial import. Every output file is written before they occur.
    - Check the files, not the exit code.
11. **Bits vs bases.** Distance targets and `flip_counts` are in bits, not nucleotides. See the note under *Zipcodes flip-mode designs*.

## Conventions

- Python 3, plus `numpy` and `pandas`. The decoy scripts need only `numpy`.
- Seeds: layouts use 42; decoys use 2026.
- Figures for publication: editable vector PDF with `matplotlib.rcParams['pdf.fonttype'] = 42` and Arial.
- Methods are named **Dual barcodes**, **Zipcodes** and **Greedy + 2-OPT**; the last is the external baseline. Quadtree and HSDB are kept for future designs.
