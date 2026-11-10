# Dual Barcode Layout Generator - Checkerboard Pattern

## Overview

This directory contains scripts for generating dual barcode layouts for **checkerboard pattern arrays** used in DNA in-situ transcriptomics. The checkerboard pattern creates two interleaved sub-arrays with spatially adjacent spots that share barcode components.

## Checkerboard Pattern Description

The checkerboard pattern consists of two interleaved arrays offset by one position:

```
Row 0: O . O . O . O .   (starts at column 0, even columns only)
Row 1: . O . O . O . O   (starts at column 1, odd columns only)
Row 2: O . O . O . O .   (starts at column 0, even columns only)
Row 3: . O . O . O . O   (starts at column 1, odd columns only)
...
```

Where `O` represents a spot and `.` represents empty space.

### Key Properties

- **Array 1 (Even Rows):** Rows 0,2,4,... starting at column 0,2,4,...
  - For 235×363 grid: 118×182 = 21,476 spots
  
- **Array 2 (Odd Rows):** Rows 1,3,5,... starting at column 1,3,5,...
  - For 235×363 grid: 117×181 = 21,177 spots

- **Total spots:** 21,476 + 21,177 = 42,653 spots

## Algorithm

### Dual Barcode Approach for Checkerboard

1. **8×8 Checkerboard Blocks**
   - The array is divided into 8×8 blocks in checkerboard space
   - Each block contains approximately 32 spots (half of 64 due to checkerboard pattern)
   - Spots within a block share the same `bar_1` barcode

2. **Two Independent Barcode Sets**
   - **Set 1:** Used for Array 1 (even rows)
   - **Set 2:** Used for Array 2 (odd rows)
   - Each set contains 18nt barcodes

3. **Barcode Assignment**
   - Each 8×8 block is assigned one `bar_1` barcode (18nt)
   - Each spot within the block receives a unique `bar_2` barcode (18nt)
   - Final barcode = `bar_1 + bar_2` = 36nt

4. **Global Uniqueness**
   - All `bar_1 + bar_2` combinations are globally unique across the entire array
   - No two spots have the same combined barcode

## Files

- **`dual_checkerboard_layout.py`** - Core library for checkerboard dual barcode generation
- **`generate_checkerboard_235x363.py`** - Script to generate 235×363 checkerboard layout
- **`barcodes_118x182_dual_set1.txt`** - Barcode set 1 (18nt, for Array 1)
- **`barcodes_118x182_dual_set2.txt`** - Barcode set 2 (18nt, for Array 2)

## Usage

### Generate 235×363 Checkerboard Layout

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Dual_barcodes
python generate_checkerboard_235x363.py
```

### Expected Output

1. **`layout_235x363_checkerboard_colmajor.csv`**
   - Column-major ordered layout (for photolithography synthesis)
   - Columns: `row, col, bar_1, bar_2, barcode`
   - ~42,653 rows (one per spot)

2. **`layout_235x363_checkerboard_stats.json`**
   - Generation timestamp and configuration
   - Statistics on barcode usage
   - Validation metrics

### Output Format Example

```csv
row,col,bar_1,bar_2,barcode
0,0,ATGATACCCAGTGTCGAC,TGCGGGCTCATAGTTCCT,ATGATACCCAGTGTCGACTGCGGGCTCATAGTTCCT
0,2,ATGATACCCAGTGTCGAC,ACACTCGTGTGGAGCAAT,ATGATACCCAGTGTCGACACACTCGTGTGGAGCAAT
0,4,ATGATACCCAGTGTCGAC,TGACTGCCGGACGATAT,ATGATACCCAGTGTCGACATGACTGCCGGACGATAT
...
1,1,TGTGCATTAGTACTCGCT,AGTTGAACAAGAGTCGCT,TGTGCATTAGTACTCGCTAGTTGAACAAGAGTCGCT
1,3,TGTGCATTAGTACTCGCT,TCGCCGCGTGCTCTAACT,TGTGCATTAGTACTCGCTTCGCCGCGTGCTCTAACT
...
```

## Requirements

### Barcode Requirements

For a 235×363 checkerboard array:
- **Array 1 spots:** 21,476
- **Array 2 spots:** 21,177
- **8×8 blocks per array:** ~336 blocks each
- **Minimum barcodes needed:** ~1,350 per set (conservative estimate)
- **Available barcodes:** 21,477 per set (more than sufficient)

### Python Dependencies

```bash
pip install numpy
```

## Algorithm Details

### Checkerboard Block Assignment

For an 8×8 checkerboard block starting at position (block_row, block_col):

```python
# For even rows (Array 1)
positions = []
for r in range(block_row*8, (block_row+1)*8, 2):  # Even rows only
    for c in range(block_col*8, (block_col+1)*8, 2):  # Even cols only
        positions.append((r, c))

# For odd rows (Array 2)
positions = []
for r in range(block_row*8+1, (block_row+1)*8, 2):  # Odd rows only
    for c in range(block_col*8+1, (block_col+1)*8, 2):  # Odd cols only
        positions.append((r, c))
```

### Barcode Assignment Strategy

1. **Group spots into 8×8 checkerboard blocks**
2. **For each block:**
   - Assign one `bar_1` from available pool (least-used)
   - Track which `bar_2` values have been used with this `bar_1`
   - Assign unique `bar_2` to each spot in block
   - Ensure no `bar_1 + bar_2` combination repeats globally

## Advantages of Checkerboard Dual Barcode Approach

1. **Spatial Error Tolerance**
   - Synthesis errors affecting `bar_1` will group nearby spots together
   - Spots within an 8×8 block share `bar_1`, providing spatial redundancy

2. **Efficient Barcode Usage**
   - Only need ~2,700 unique 18nt barcodes total (both sets)
   - Much fewer than the 42,653 spots in the array

3. **Simple Pattern**
   - Deterministic checkerboard positioning
   - Easy to verify and visualize

4. **Global Uniqueness**
   - Every spot has a unique 36nt barcode
   - No ambiguity in spot identification

## Validation

The generated layout is automatically validated for:

1. **Coverage:** All checkerboard positions have barcodes
2. **Block consistency:** All spots in a block share the same `bar_1`
3. **Block uniqueness:** All spots in a block have different `bar_2`
4. **Global uniqueness:** No duplicate combined barcodes across entire array

## Visualization

To visualize which spots are activated at each synthesis cycle:

```bash
python plot_cycle_activation.py layout_235x363_checkerboard_colmajor.csv 0
```

This will show which spots are activated during synthesis cycle 0 (base A).

## Performance

- **Layout generation:** ~5-15 seconds
- **CSV export:** ~5-10 seconds
- **Memory usage:** ~100-200 MB
- **Output file size:** ~3-4 MB

## Comparison: Regular vs Checkerboard

| Property | Regular 4×4 | Checkerboard 8×8 |
|----------|-------------|------------------|
| **Pattern** | Continuous grid | Two interleaved arrays |
| **Block size** | 4×4 = 16 spots | 8×8 checkerboard = ~32 spots |
| **Spatial grouping** | Tight (2μm blocks) | Wider (8μm blocks) |
| **Barcode efficiency** | High | Very high |
| **Synthesis tolerance** | Good | Excellent |

## Troubleshooting

### Error: "Insufficient barcodes"
- **Cause:** Not enough barcodes in input files
- **Solution:** Ensure each barcode file has at least 1,350 barcodes
- **Check:** Count lines in barcode files: `wc -l barcodes_118x182_dual_set*.txt`

### Error: "Cannot find enough unique bar_2 barcodes"
- **Cause:** A `bar_1` has been paired with too many `bar_2` values
- **Solution:** Increase total number of barcodes or reduce block size
- **Note:** Should not occur with provided barcode sets

### Error: "Duplicate combined barcodes"
- **Cause:** Algorithm error or insufficient validation
- **Solution:** Report as bug with reproduction details

## Python API

### CheckerboardDualBarcodeLayout Class

```python
from dual_checkerboard_layout import CheckerboardDualBarcodeLayout, load_barcodes_from_file

# Load barcodes
barcodes_set1 = load_barcodes_from_file('barcodes_118x182_dual_set1.txt')
barcodes_set2 = load_barcodes_from_file('barcodes_118x182_dual_set2.txt')

# Create layout
layout = CheckerboardDualBarcodeLayout(
    rows=235,
    cols=363,
    barcodes_set1=barcodes_set1,
    barcodes_set2=barcodes_set2,
    block_size=8
)

# Generate
grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=42)

# Export
layout.export_layout('output.csv', format='csv', column_major=True)

# Statistics
stats = layout.get_statistics()
print(f"Total spots: {stats['total_spots']}")
print(f"Unique barcodes: {stats['unique_combined_barcodes']}")
```

## Contact

For questions or issues regarding the checkerboard dual barcode approach, please refer to the main project documentation.
