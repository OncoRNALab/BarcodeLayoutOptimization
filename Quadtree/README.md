# Quadtree Barcode Generation with Run-Length Limiting (RLL)

## Overview

This directory implements quadtree-based spatial barcode generation for photolithography arrays with Run-Length Limiting (RLL) to prevent homopolymer sequences.

## Method Description

### Quadtree Spatial Encoding

The quadtree method encodes spatial positions hierarchically through recursive subdivision:

1. **Coordinate Scaling**: Map (row, col) coordinates to binary space [0, 2^L)
2. **Hierarchical Encoding**: At each level, determine which quadrant contains the position
3. **Base Assignment**: Map quadrant index (0-3) to DNA base (A, C, G, T)

This creates barcodes where spatial proximity correlates with sequence similarity (low Hamming distance).

### Run-Length Limiting (RLL) Algorithm

**Problem**: Without constraints, quadtree encoding can generate homopolymers (e.g., AAA, GGGG) which are problematic for DNA synthesis and sequencing.

**Solution**: At each cycle, intelligently select a permutation that minimizes homopolymer formation while maintaining spatial encoding.

#### RLL Process (Per Cycle)

**Step 1: Build Contingency Tables**

For all 21,476 spots, create a 4×4 table tracking:
- Rows: Current quadrant (0, 1, 2, 3)
- Columns: Previous base (A, C, G, T)
- Values: Count of spots in each combination

**Step 2: Evaluate All 24 Permutations**

Test every possible mapping of quadrants → bases:
```
Example permutations:
- ACGT: Quad0→A, Quad1→C, Quad2→G, Quad3→T
- ATCG: Quad0→A, Quad1→T, Quad2→C, Quad3→G
- ... (24 total)
```

For each permutation P, calculate violations:
```python
violations = sum(contingency[q, P[q]] for q in [0,1,2,3])
```

This counts how many spots would get the same base as their previous base (creating homopolymers).

**Step 3: Select Best Permutation**

Choose P* with minimum violations:
```
P* = argmin(violations)
```

**Step 4: Apply Permutation**

Assign bases to all spots:
```python
for spot in spots:
    quadrant = calculate_quadrant(spot.row, spot.col, cycle)
    base = P*[quadrant]
```

**Step 5: Forced Flips (Fallback)**

If a spot *still* violates max_run after best permutation:
```python
if (base == prev_base) and (run_length >= max_run):
    offset = (row + col + cycle) % 3 + 1  # Deterministic: 1, 2, or 3
    base = (prev_base + offset) % 4        # Force different base
```

**Step 6: Update State**

Track run length for next cycle:
```python
if base == prev_base:
    run_length += 1
else:
    run_length = 1
prev_base = base
```

### Parameters

- **max_run**: Maximum allowed consecutive identical bases
  - `max_run = 1`: No homopolymers allowed (AA forbidden)
  - `max_run = 3`: Up to 3 consecutive allowed (AAA ok, AAAA forbidden)
- **barcode_length**: Number of cycles (bases per barcode)
  - 24nt: 3 rounds of 8 cycles
  - 36nt: 3 rounds of 8 cycles + 12nt random library

### Key Insights

**Why Global Optimization Works:**
- Evaluating all 21,476 spots simultaneously finds permutations that minimize violations across the entire array
- The 24 permutation choices provide sufficient flexibility
- Different permutations per cycle adapt to changing base distributions

**Why It's Efficient:**
- Only 24 permutations to test per cycle
- O(N) contingency table building
- O(1) permutation evaluation (4 lookups)
- Total: O(N × L) for N spots and L cycles

**Forced Flips:**
- Rare: typically <1% of bases need forced flips
- Deterministic: same (row, col, cycle) always gets same flip
- Minimal impact on spatial similarity

## Files in This Directory

- **`quadtree_rll.py`**: Core RLL implementation
  - `generate_barcode_matrix_rll()`: Main generation function
  - `_best_perm_avoid_violations()`: Permutation selection
  - `verify_no_homopolymers()`: Validation

- **`output_utils.py`**: Export utilities
  - `write_barcodes_column_major()`: Column-major TXT export for photolithography
  - `write_statistics_json()`: Quality metrics
  - `validate_grid_coverage()`: Completeness check

- **`generate_repeated_quadtree.py`**: 24nt barcode generator
  - Generates continuous 24nt barcodes with RLL
  - Outputs: CSV, FASTA, column-major TXT, statistics

- **`generate_hybrid_quadtree_library.py`**: 36nt barcode generator
  - 24nt quadtree + 12nt random library
  - High diversity for sequencing, spatial structure for synthesis

## Usage

### Generate 24nt Barcodes

```bash
cd /path/to/Quadtree
python generate_repeated_quadtree.py
```

Output files in `outputs/no_homopylmers/`:
- `layout_118x182_repeated_quadtree_24nt.csv`
- `layout_118x182_repeated_quadtree_24nt.fasta`
- `layout_118x182_repeated_quadtree_24nt_colmajor.txt`
- `layout_118x182_repeated_quadtree_24nt_stats.json`

### Generate 36nt Hybrid Barcodes

```bash
cd /path/to/Quadtree
python generate_hybrid_quadtree_library.py
```

Requires: `bar_1in_4_len12.txt` (12nt barcode library)

## Results

**24nt Repeated Quadtree (118×182 array):**
- Total barcodes: 21,476
- Homopolymers: 0 (100% clean)
- Unique barcodes: 21,474 (99.99%)
- Forced flip rate: 24.1%
- Generation time: ~1 minute

**Performance (250K simulated reads, error rates: sub=0.03, ins=0.04, del=0.14):**
- 36nt hybrid: Precision = 0.8994, Recall = 0.5931
- 24nt repeated: Evaluation in progress

## References

- Press et al. (2022): RandomBarcodes method for comparison
- Quadtree spatial indexing: Classical computer graphics technique
- RLL principle: Adapted from magnetic recording and data storage
