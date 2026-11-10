# AGENTS.md — barcode_design

Read this first. It is the entry point for any agent or new collaborator working in this repository.

## What this project is

This repository designs **spatial barcode layouts** for photolithographic (maskless, MAS 2.0) spatial-transcriptomics arrays. It also generates the **decoy whitelists** used to estimate decoding precision on printed arrays.

- **Array:** full density is 235 rows × 363 columns, giving 85,305 features of 13.7 µm with 1 µm gaps. Half density is a checkerboard of 42,471 printed features.
- **Barcode:** every feature carries a 36-nt barcode, read from Read 1 at positions 11–46, right after the 10-nt UMI.
- **Synthesis model:** probes are synthesized 5′→3′ with the deposition order A, C, G, T, and each barcode is embedded left-most. Features that share a barcode prefix are illuminated (co-activated) together, which can reduce stray-light insertions.
- **Decoding:** QUIK, using sequence-Levenshtein distance and the 4-mer single-stage filter, with rejection threshold r ∈ {5, 6, 7}, or {6, 7, 8} at half density.

The work is published in a BMC Bioinformatics manuscript, which lives outside this repo at `~/Documents/phd_thesis_revision/layoutopt_manuscript/`. `docs/methods.md` is an excerpt of its Methods and is the **authoritative description** of the printed designs and decoys. If code or older READMEs disagree with `docs/methods.md`, trust `docs/methods.md`, then the data, then the code.

## Methods in this repository

All methods are first-class, and none are archived. Three were printed and appear in the manuscript. The others are kept for future work.

| Folder | Method | Status |
|---|---|---|
| `Dual_barcodes/` | Dual barcodes: 4×4 blocks share bar_1 | Printed at full and half density (manuscript) |
| `Zipcodes/` | Zipcodes: [column code][row code] | Printed at full density (manuscript) |
| *(external)* | Greedy + 2-OPT (BarLay, GPU) | Printed at full density (manuscript); **generation code is not here** |
| `Quadtree/` | Hierarchical quadtree + RLL, hybrid 24 nt + 12 nt | Layouts exist; printed variant `data/1in1/layout_235x363_quadtree_8x3_12nt.*`; **excluded from the manuscript** |
| `HSDB/` | Hierarchical Spatial Dual Barcodes (28-nt Zipcode bar_1 + 8-nt bar_2, G ≤ 25%) | Designed, not printed; see `IMPROVEMENT_PLAN.md`, `SummaryNewDesigns.md` |
| `inhouse_barcodes/` | Random barcode libraries (filtered random 36-nt pool and its 12/18-nt prefixes and suffixes) | Input for Dual, Quadtree and decoys |
| `decoys/` | Structure-matched decoy generation + QC | Used for every manuscript precision estimate |
| `figures/` | Early 118×182 schematic figures | Historical |

### Dual barcodes (printed)
- **Layout:** 59 × 91 = 5,369 blocks of 4×4, with truncated blocks at the bottom and right edges. Barcode = bar_1 (18 nt, shared within a block) + bar_2 (18 nt).
- **Library:** both halves come from one library, `inhouse_barcodes/bar_1in1_column_major_first18.txt` (85,305 sequences; the 18-nt prefixes of the filtered random 36-nt pool).
- **Assignment:** each block's bar_1 is the least-used sequence, so all 5,369 block IDs are distinct. bar_2 is drawn at random, excluding codes already paired with that block's bar_1. bar_2 codes are reused across blocks.
- **Composition:** no filter was applied to the concatenated barcode. It has up to 42% G, and homopolymers of up to 6 nt at the junction.
- **Code:** `Dual_barcodes/generate_layout_235x363.py` (config `data/1in1/config_files/dual_config_235x363.json`, seed 42). Half-density checkerboard: `Dual_barcodes/generate_checkerboard_235x363.py` / `dual_checkerboard_layout.py`.
- **Half density:** the first 117 rows of the full-density layout (117 × 363 = 42,471), placed where x + y is even. Printer coordinates are y = col + 4, and x = 4 + 2·row for even columns or 5 + 2·row for odd columns. A 4×4 block occupies an 8×4 printer window.

### Zipcodes (printed)
- **Barcode:** [column code 18 nt][row code 18 nt], with no separator and no random tail.
- **Axis codes:** each is a **36-bit codeword**, converted at 2 bits per nt with **00→A, 01→C, 10→T, 11→G**.
- **Generation:** a randomized search with backtracking (`Zipcodes/spatial_zipcodes/axis_code_generator.py`, seed 42), testing up to 10,000 candidates per position.
  - **Soft targets:** bit distance to the codes 1, 2, 3 and 4 positions back should be 2, 4, 6 and 8 bits. A candidate is accepted if the mean squared deviation is ≤ 16.
  - **Hard constraint:** ≥ 8 bits from every code 5 or more positions back.
  - **Backtracking:** up to 5 positions.
- **Measured distances:** adjacent column codes differ by 8.8 bits on average (range 4–12), 10.5 bits at offset 4 and 17.9 bits beyond. The targets are not reached.
- **Code:** `Zipcodes/generate_layout_235x363_zipcodes.py Zipcodes/config_235x363_nosep.json` (a copy of the config is in `data/1in1/config_files/zip_config_235x363_nosep.json`).

### Greedy + 2-OPT (printed, external code)
- **Pool and placement:** a random pool of unique 36-nt barcodes (homopolymers ≤ 3 nt, G ≤ 25%) was placed by BarLay (Jatzkowski et al. 2024) with default settings.
- **Cost:** each barcode is a left-most schedule binary vector, and the cost is the summed schedule Hamming distance over all 8-neighbour pairs.
- **Greedy step:** fills one position at a time with the unused barcode closest to its already placed neighbours.
- **2-OPT step:** batched best-improvement swaps between placed barcodes.
- **Pools:** the pool of 150,000 barcodes was generated twice. One pool fed Greedy; 85,305 barcodes of the second pool (`data/decoy_layouts/barcodes/1in1/greedy_decoybar_1in1_column_major.txt`) are the source of the Greedy decoys.

### Synthesis / co-activation facts (verified from printing masks)
- All barcodes start in the same cycle, right after the UMI.
- **Zipcodes:** the printer **waits until every column code is finished** before starting the row codes, so whole rows co-activate during the second half. Without the wait, phases would drift between features.
- **Dual:** it is unknown whether the printer also waits between bar_1 and bar_2.

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
AGENTS.md                 <- this file
docs/methods.md           <- manuscript Methods excerpt (authoritative)
data/
  1in1/                   printed full-density layouts (.csv with coords, .txt one barcode per line) + config_files/
  1in2/                   half-density layouts (Dual: layout_1in2_dual_4x4.csv with printer x,y)
  decoy_layouts/barcodes/1in1/   old independent decoys + the second Greedy pool
decoys/scripts/           make_dual_decoys.py, make_zip_decoys.py, make_greedy_decoys.py, check_decoys.py, run_all.sh
decoys/output/, decoys/qc/       generated decoys and QC JSONs
Dual_barcodes/ Zipcodes/ Quadtree/ HSDB/   method code + their outputs (incl. 118×182 sets and *_decoy variants)
inhouse_barcodes/         random input libraries
figures/, plot_*.py, generate_presentation.py, combine_full_checkerboard.py, get_layout_bounds.py   plotting/helpers
README.md, README_118x182_generation.md, IMPROVEMENT_PLAN.md, SummaryNewDesigns.md, manuscript.txt   older notes
*.sh                      HPC (VSC, SLURM) job scripts with absolute cluster paths
```

**Printed layout files** (the ground truth):

| Method | File | Format |
|---|---|---|
| Dual | `data/1in1/layout_235x363_dual_4x4.csv` | `row,col,combined` (identical to `Dual_barcodes/layout_235x363_dual_4x4.csv`, which also has `bar_1,bar_2`) |
| Zipcodes | `data/1in1/layout_235x363_zipcodes_nosep.csv` | `row,col,barcode` (identical to `Zipcodes/output_235x363_nosep/layout_235x363_zipcodes_nosep_colmajor.csv`) |
| Greedy | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` | no header: `barcode,col,row` |
| Dual half density | `data/1in2/layout_1in2_dual_4x4.csv` | `row,col,combined,x,y` |

The `.txt` versions list one barcode per line, in column-major order.

## Known pitfalls

1. **Zipcodes are not Gray codes.** `README.md`, `figures/figure_3a_gray_code_principle.py` and parts of `Zipcodes/README.md` describe Gray codes and distances in bases. Both are outdated: the distances are in bits, and the codes come from randomized search. Banners mark the READMEs as historical.
2. **Bit-to-base mapping.**
   - The printed design and `make_zip_decoys.py` use 00→A, 01→C, 10→T, 11→G.
   - `spatial_zipcodes/bit_to_dna.py` defaults to 10→G, 11→T.
   - The config labels its mapping "4-bits-per-base" one-hot (1000→A, 0100→C, 0010→T, 0001→G), which in practice yields the 10→T, 11→G mapping.
   - Check before reusing.
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
   - HSDB and the newer Zipcodes ideas add G ≤ 25%.
6. **2-OPT pool exchanges (open question).** The manuscript says 2-OPT also swaps placed barcodes with unplaced pool barcodes. The BarLay paper describes swaps between placed barcodes only. Unresolved.
7. **φ for Zipcodes.** Zipcodes misassignments stay on one axis, so use φ = 1/3 (local), not the global ≈ 1/2.
8. **Absolute cluster paths.** These point at `/user/gent/446/vsc44685/ScratchVO_dir/layout_opt/…`. They appear in the `*.sh` scripts and in `barcode_file` / `library_file` in the JSON configs. Adapt them before running locally.
9. **Not tracked in git (see `.gitignore`).** Simulated reads (`*.fastq`, `*_read_L36zip*.txt`, ~230 MB) and SLURM logs are local only.

## Conventions

- Python 3, plus `numpy` and `pandas`. The decoy scripts need only `numpy`.
- Seeds: layouts use 42; decoys use 2026.
- Figures for publication: editable vector PDF with `matplotlib.rcParams['pdf.fonttype'] = 42` and Arial.
- Methods are named **Dual barcodes**, **Zipcodes** and **Greedy + 2-OPT**. Quadtree and HSDB are kept for future designs.
