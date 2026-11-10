# Dual Barcodes — Improvement Ideas: Evaluation & Implementation Plan

**Context:** Current state: 36nt barcodes split as 18nt bar_1 + 18nt bar_2.  
Array: 118×182 = 21,476 spots, divided into 1,380 blocks of 4×4.  
Library: 5,369 barcodes of 18nt.  
Results: Precision = 0.9731, Recall = 0.5703 (107,427 erasures out of 250,000 reads).

---

## Idea 1 — Asymmetric split: longer bar_1 (block ID), shorter bar_2 (spot ID)

### The proposal
Change the 18+18 split to 28nt bar_1 + 8nt bar_2 (or similar asymmetric ratio).  
Rationale: bar_2 only needs to distinguish 16 spots within a block. In theory, 2nt suffices (4² = 16); practically, 8nt is generous and allows high minimum Hamming distance among the 16 spot barcodes.

### Evaluation

**Bar_2 feasibility (8nt for 16 unique spot IDs):**
- Sequence space: 4⁸ = 65,536 possible 8-mers.
- Selecting 16 barcodes with minimum Hamming distance ≥ 4 is trivial from this space.
- In practice a greedy or exhaustive search can find 16 barcodes with min-distance ≥ 4 or even ≥ 5 easily.
- This means within-block spot decoding is a 1-of-16 problem with generous error tolerance — easy to decode even under high error rates.

**Bar_1 feasibility (28nt for 1,380 block IDs):**
- Need a new library of 1,380 barcodes of 28nt with high minimum Hamming distance.
- The existing 5,369 barcodes are 18nt; they cannot be used directly, but can be extended (e.g., append a fixed 10nt suffix, or concatenate two 14nt sub-libraries).
- Alternatively, generate a new random library with minimum Hamming distance ≥ 8 from scratch (straightforward with 28nt sequences).

**Error tolerance argument:**
- The current 36nt combined barcode is decoded by the RandomBarcodes algorithm as a 1-of-21,476 problem. A single read accumulates on average ~(0.03+0.04+0.14)×36 ≈ 7.6 expected error events. With 21,476 candidates and that error rate, many reads fall into the erasure zone (the 107,427 erasures correspond to a 43% erasure rate).
- With an asymmetric split and **hierarchical decoding**:
  - Step 1: Decode the first 28nt against 1,380 candidates → 1-of-1,380 problem.
  - Step 2: Decode the last 8nt against 16 candidates → 1-of-16 problem.
  - Both sub-problems are far easier than the current 1-of-21,476 problem, so recall should increase substantially.
- Even with the same RandomBarcodes decoder on the full 36nt barcode, the 28nt bar_1 carries more redundancy: if bar_1 minimum Hamming distance is ≥ 8, the decoder can tolerate up to 4 substitutions in bar_1 before misidentifying the block.

**Potential downside:**
- A shorter 8nt bar_2 means fewer bases for spot disambiguation. At error rates sub=0.03, ins=0.04, del=0.14, an 8nt sequence accumulates ~1.7 expected errors. With min-distance ≥ 4 among 16 candidates, this should still be decodable; but with indels it is tighter. Verifying this via simulation is essential.
- If bar_2 is corrupted beyond recovery, the spot is lost but the block is still known — a graceful degradation (block-level assignment could still be useful).

**Verdict: ✅ Sound idea, straightforward to implement. Gains are primarily through enabling hierarchical decoding and allocating bits more efficiently (block ID is harder to decode than within-block spot ID, so it deserves more bits).**

---

## Idea 2 — Zipcode-style inter-block distance constraint on bar_1

### The proposal
Instead of randomly assigning bar_1 values to blocks, apply the Zipcode (Gray-code spatial gradient) algorithm at the **block level**. Treat the block grid (30×46 = 1,380 blocks) as a chip, run the Zipcode algorithm to generate bar_1 sequences, then append a short bar_2 for within-block uniqueness.

### Evaluation

**What this fixes:**
- The main weakness of current Dual Barcodes is **sharp discontinuities at block boundaries**. Because bar_1 is randomly assigned, adjacent blocks get unrelated (high-Hamming-distance) barcodes from a high-diversity library. A synthesis error that activates a spot in the adjacent block produces a totally different barcode → likely a wrong decode or erasure.
- With Zipcode-style bar_1, adjacent blocks differ by ~2 bases. A spot that gets the adjacent block's bar_1 due to a synthesis error would differ from the correct full barcode by only 2 bases (bar_1 part) → likely still decoded correctly or at least assigned close by.

**Spatial profile of the combined barcode:**
| Spot pair | Hamming distance contribution |
|---|---|
| Within block (same bar_1, different bar_2) | Only bar_2 differs → low (≈ min Hamming of bar_2 set) |
| Adjacent block (different bar_1, different bar_2) | bar_1 differs by ~2 bases + bar_2 differs → moderate |
| Far blocks | bar_1 differs by many bases + bar_2 differs → high |

This is a strictly better spatial profile than the current random assignment.

**Interaction with Idea 1:**
Applying Zipcode-style constraints is orthogonal to the bar split ratio. Both ideas can and should be combined:
- **28nt Zipcode-style bar_1** (spatial gradient between blocks)
- **8nt bar_2** (spot ID within block)
- Total: 36nt

**Implementation complexity:**
- The Zipcode package already exists in `Layout_opt/Zipcodes/`.
- Running it on a 30×46 grid instead of 118×182 is a minor configuration change.
- The output would be 1,380 bar_1 sequences of appropriate length with the spatial gradient property.
- The distance targets for the block-level Zipcode need to be re-calibrated (a block is physically 4× larger than a spot, so the "1-step" distance in block space corresponds to 4 spots in physical space).

**Key design decision — bar_1 length vs Zipcode bit budget:**
- The Zipcode algorithm for a 30×46 grid uses ⌈log₂(30)⌉ = 5 bits for Y and ⌈log₂(46)⌉ = 6 bits for X.
- With 2-bits-per-base encoding: 5 bits → 3 bases for Y, 6 bits → 3 bases for X → total ≈ 6 bases for the positional encoding.
- This is far less than 28nt. The remaining bases can be filled with random library barcodes (as in the Quadtree hybrid approach) or by running multiple Zipcode rounds.
- Alternatively, increase Zipcode distances [4, 8, 12, 16] at the block level to use more bits.

**Verdict: ✅ Excellent idea. It directly addresses the core gap in the current algorithm (no inter-block spatial gradient) and is achievable by reusing the existing Zipcode infrastructure.**

---

## Combined Proposal: Hierarchical Spatial Dual Barcodes (HSDB)

Combine both ideas into a single design:

```
[28nt Zipcode bar_1 (block ID with spatial gradient)] + [8nt bar_2 (spot ID within block)]
Total: 36nt
```

**Properties:**
1. **Within block:** All 16 spots share bar_1 → 0 Hamming distance in bar_1; only bar_2 differs.
2. **Adjacent blocks:** bar_1 differs by ~2 bases (Zipcode property); bar_2 also differs → graceful transition.
3. **Far blocks:** bar_1 differs by many bases; bar_2 also differs → easy to distinguish.
4. **Hierarchical decoding:** Decode 28nt bar_1 first (1-of-1,380), then 8nt bar_2 (1-of-16) → much higher recall expected.

---

## Implementation Plan

### Phase 1 — Design the bar_2 mini-library (within-block spot IDs)

**Goal:** Generate 16 barcodes of 8nt with maximum minimum Hamming distance.

**Steps:**
1. Write `generate_spot_barcodes.py`:
   - Enumerate or greedy-search 8nt sequences.
   - Select 16 with minimum Hamming distance ≥ 4 (target ≥ 5).
   - Apply QC: GC content 40–60%, no homopolymers > 3.
   - Export to `spot_barcodes_8nt.txt`.
2. These 16 barcodes are **fixed and reused** across every block (the bar_1 provides block-level uniqueness, so bar_2 can repeat across blocks as long as the (bar_1, bar_2) pair is globally unique).

**Estimated effort:** 1–2 hours.

---

### Phase 2 — Generate bar_1 with Zipcode spatial constraints at block level

**Goal:** Generate 1,380 barcodes of 28nt for the 30×46 block grid with spatial gradient.

**Option A — Run Zipcode on the block grid (recommended):**
1. Create a new config file `block_zipcode_config.json` pointing to a 30×46 "chip".
2. Set distance targets appropriate for block-level spacing.
   - Since 1 block step = 4 spot steps physically, the inter-block Hamming distances can be set higher than spot-level Zipcodes.
   - Suggested: `distances: [4, 8, 12, 16]`, `far_features_min_distance: 16`.
3. Run `python -m spatial_zipcodes.cli full-pipeline --config block_zipcode_config.json --output-dir output_block_zip/`
4. This produces 1,380 barcodes. Convert to 28nt by using 4-bit photolithography encoding or by padding with the random library.

**Option B — Extend existing Zipcode barcodes:**
- Take the existing Zipcode output (`barcodes_l36_zip.txt` or `barcodes_len26.csv`) at the block-corner positions (rows 0,4,8,…; cols 0,4,8,…) and use those as bar_1.
- This reuses already-generated and evaluated barcodes.
- Less clean but quicker.

**Estimated effort:** 2–4 hours for Option A, 1 hour for Option B.

---

### Phase 3 — Build the combined layout generator

**Goal:** New script `dual_barcode_hierarchical.py` that assembles the full 36nt layout.

**Steps:**
1. Load 1,380 block bar_1 sequences (from Phase 2) — already spatially ordered (block_row, block_col) → bar_1.
2. Load 16 spot bar_2 sequences (from Phase 1).
3. For each block (br, bc):
   - Assign `bar_1 = block_bar1[br, bc]` to all 16 spots in the block.
   - Randomly assign each of the 16 `bar_2` sequences to one spot in the block (this is fine because bar_1 already makes all combinations globally unique — no two blocks share the same bar_1).
4. Combined barcode = `bar_1 + bar_2`.
5. Validate global uniqueness.
6. Export CSV.

**Key simplification vs current code:** Because bar_1 values are all unique (one per block) and bar_2 values are the same fixed set of 16, every `(bar_1, bar_2)` combination is automatically unique. No need for the complex `used_combinations` tracking in the current code.

**Estimated effort:** 3–5 hours.

---

### Phase 4 — Simulation and benchmarking

**Goal:** Compare HSDB against current Dual Barcodes and other methods.

**Steps:**
1. Run the same 250,000-read simulation (error rates: sub=0.03, ins=0.04, del=0.14) using the RandomBarcodes decoder.
2. Additionally, implement a **hierarchical decoder**:
   - Decode first 28nt of a read against bar_1 library (1,380 candidates).
   - Decode last 8nt against bar_2 library (16 candidates).
   - Report block and spot ID separately.
3. Compare metrics:

| Metric | Current Dual (18+18) | HSDB flat (28+8, random bar_1) | HSDB spatial (28+8, Zipcode bar_1) |
|---|---|---|---|
| Precision | 0.9731 | TBD | TBD |
| Recall | 0.5703 | TBD | TBD |
| Erasure rate | 0.43 | TBD | TBD |

4. Also compute:
   - Hamming-distance vs. physical-distance correlation (spatial metric).
   - Neighbor Hamming distance distribution (within-block vs. across-block).

**Estimated effort:** 4–6 hours.

---

### Phase 5 — Analysis and iteration

**Goal:** Interpret results and tune parameters.

**Key questions to answer:**
1. Does the asymmetric split (28+8) improve recall vs. 18+18 when using the same flat decoder?
2. Does the hierarchical decoder outperform the flat decoder on the asymmetric layout?
3. Does the Zipcode-style bar_1 improve precision/recall at block boundaries specifically?
4. Is 8nt sufficient for bar_2 under indel errors (del=0.14)?
   - If not, try 10nt bar_2 + 26nt bar_1.
5. What is the optimal split ratio? (Sweep: 20+16, 24+12, 28+8, 30+6)

---

---

## Idea 3 — G-content control (photolithography-specific constraint)

### The problem

In photolithography-based synthesis, coupling efficiency is **not equal across all four bases**. G has measurably lower coupling efficiency than A, C, or T. This means that every G base in a barcode is a position where synthesis is more likely to fail — producing a missing base (effective deletion), a truncated strand, or a carry-over from the previous cycle.

The **current GC-content filter** in both pipelines does not address this, because it treats G and C identically. A barcode with 40% GC content could have all of that as G (40% G, 0% C) and still pass QC.

There are **two distinct levels** at which G content must be controlled:

### Level 1 — Per-barcode G content

Each individual barcode should have G frequency below a threshold. This is the straightforward analog of the GC-content filter.

**Current state (measured):**

| Library / Layout | Mean G% | Barcodes with G > 25% | Barcodes with G > 30% |
|---|---|---|---|
| Dual Barcodes library (18nt, n=5,369) | 21.2% | 30.3% | 9.1% |
| Zipcodes barcodes (36nt, n=21,476) | **31.6%** | **79.5%** | **64.2%** |

Zipcodes has a severe G-content problem: nearly 80% of its barcodes exceed 25% G. Per-position analysis revealed the root cause (see below).

**Target:** G content per barcode ≤ 25% (i.e., at most 6 G's in a 24nt barcode, at most 9 in a 36nt barcode).

### Level 2 — Per-cycle G frequency across the array

For photolithography, all spots on the array are synthesized **one cycle at a time**. At each cycle *i*, every spot receives whichever base its barcode specifies at position *i*. If many spots require G at the same cycle, the G coupling problem is concentrated in that cycle.

**Measured per-cycle G frequency:**

| Layout | Mean G/cycle | Max G/cycle | Cycles with G > 30% |
|---|---|---|---|
| Dual Barcodes (36nt, 118×182) | 21.6% | 23.9% (pos 2) | 0 |
| Zipcodes with GGGG separator (36nt) | 31.6% overall | **100% (pos 16–19)** | **4 cycles at 100%** |

**Target:** No synthesis cycle should have G frequency > 25% across the array.

### Root cause in Zipcodes — the GGGG separator, not the encoding

Per-position analysis of `barcodes_l36_zip.txt` showed:

```
Position:  0────────────────15 | 16─19  | 20────────────────35
Component: [  X-code (16 nt) ] | [GGGG] | [  Y-code (16 nt) ]
G content:       ~25%          |  100%  |       ~22%
```

Positions 16–19 are **100% G** across all 21,476 barcodes. The X and Y code portions themselves have G ≈ 22–25%, which is perfectly acceptable. The problem is entirely caused by the 4-nt `GGGG` spacer inserted between them.

**Arithmetic:** 4 forced G's out of 36 nt adds 11% G on top of the natural ~22–25% → total ≈ 31.6%.

The bit-to-DNA mapping and the Gray code structure are **not the primary problem**. The fix is simply to eliminate or replace the GGGG separator.

**Fix options (in order of preference):**
1. **Remove the separator entirely** — already tested: Precision = 0.9014, Recall = 0.6039. G content drops to ~23%.
2. **Replace with a G-free separator** — if a boundary landmark is still needed for decoding, use `ATAT`, `CACA`, or similar. Same structural benefit, zero G cost.
3. **Use a shorter separator** — e.g., `AT` (2 nt) as a minimal boundary marker with no G.

The bit-to-DNA mapping is a secondary consideration. The X and Y codes themselves have G ~22–25%, which is already within target. No mapping change is needed unless the separator is kept.

#### Zipcodes — per-cycle G frequency monitor

Add a function to `distance_analysis.py` or a new `synthesis_qc.py` module:

```python
def per_cycle_base_frequencies(codebook_df, barcode_col='zip_seq'):
    """
    Returns a DataFrame of per-position base frequencies.
    Each row is a cycle (position), columns are A, C, G, T fractions.
    """
```

Report: max G frequency per cycle, cycles with G > 25%, heatmap of per-cycle G content.

#### Dual Barcodes — library generation and filtering

- Filter the existing library (`barcodes_1in4_half.txt`) to G ≤ 25%: of 5,369 barcodes, 30.3% (1,625) would be removed, leaving ~3,744. This is still sufficient for 1,380 blocks.
- For new libraries generated under Idea 1 (28nt bar_1, 8nt bar_2), enforce G ≤ 25% at generation time.
- Add a per-cycle G frequency check to the layout validation step.

#### Quadtree — RLL permutation selection

The RLL algorithm already selects among 24 permutations to minimize homopolymers. Add a **secondary criterion**: among permutations with equal or near-equal homopolymer violations, prefer the one that minimizes G usage at each cycle. This is a 1-line change in the permutation scoring function.

### Implementation steps

1. **Zipcodes — replace the GGGG separator** in all configs (`spacer_sequence: ""` or `"ATAT"`). This alone resolves most of the G inflation. (5 min)
2. **Add `g_content_max` to `sequence_qc.py`** and update `SequenceQC` and `validate_axis_sequences()`. (1–2 hours)
3. **Update all config files** to include `g_content_max: 0.25`. (15 min)
4. **Add per-cycle G frequency reporting** to the Zipcodes pipeline output. (1–2 hours)
5. **Filter Dual Barcodes library** to G ≤ 25% and update generation scripts. (30 min)
6. **Add G-minimization secondary criterion** to Quadtree RLL permutation selection. (1 hour)
7. **Re-run simulations** for all methods with corrected barcodes.

**Estimated effort:** 3–5 hours across both pipelines.

---

## Summary

| | Idea 1 (Asymmetric split) | Idea 2 (Zipcode bar_1) | Idea 3 (G-content control) | Combined (HSDB + G-control) |
|---|---|---|---|---|
| Addresses sharp block boundaries | ❌ | ✅ | ❌ | ✅ |
| Enables hierarchical decoding | ✅ | ❌ | ❌ | ✅ |
| Reduces synthesis errors (G coupling) | ❌ | ❌ | ✅ | ✅ |
| Reuses existing infrastructure | Partially | ✅ (Zipcode package) | ✅ (add to QC module) | ✅ |
| Applies to Dual Barcodes | ✅ | ✅ | ✅ | ✅ |
| Applies to Zipcodes | ❌ | ✅ | ✅ (urgent) | ✅ |
| Applies to Quadtree | ❌ | ❌ | ✅ (RLL tweak) | — |
| Implementation complexity | Low | Medium | Low–Medium | Medium |
| Expected recall improvement | High | Moderate | Moderate | High |
| Expected precision improvement | Moderate | Moderate | Moderate | High |

**Recommended implementation order:**
1. **Idea 3 first** — G-content control is a correctness issue (especially for Zipcodes) and cuts across all methods. It should be fixed before generating any new barcode sets.
   - Add `g_content_max` to `sequence_qc.py`
   - Update configs and regenerate Zipcodes barcodes with G-friendly mapping
   - Filter Dual Barcodes library
2. Phase 1 (bar_2 mini-library with G ≤ 25%) — quick, independent
3. Phase 3 (layout generator with flat Idea 1 only) — test asymmetric split
4. Phase 4 (simulation, flat decoder) — validate Idea 1 in isolation
5. Phase 2 (Zipcode bar_1 with G control) — add spatial constraint to block-level bar_1
6. Phase 3 update + Phase 4 repeat — validate combined HSDB design
7. Phase 4 hierarchical decoder — test the full HSDB concept
