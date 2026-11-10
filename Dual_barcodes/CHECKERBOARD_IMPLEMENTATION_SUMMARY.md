# Checkerboard Dual Barcode Implementation Summary

## Overview

I have successfully adapted the dual barcode approach to work with checkerboard pattern arrays. The implementation generates barcode layouts for arrays where spots are arranged in a checkerboard pattern with two interleaved sub-arrays.

## What Was Created

### 1. Core Library: `dual_checkerboard_layout.py`

A new Python module implementing the `CheckerboardDualBarcodeLayout` class with the following features:

- **Checkerboard Pattern Generation**: Automatically generates valid spot positions where:
  - Even rows (0,2,4,...) have spots at even columns (0,2,4,...)
  - Odd rows (1,3,5,...) have spots at odd columns (1,3,5,...)
  
- **8×8 Block Assignment**: Groups spots into 8×8 checkerboard blocks where:
  - Each block contains approximately 32 spots (half of 64 due to checkerboard)
  - All spots in a block share the same `bar_1` barcode (18nt)
  - Each spot in a block has a unique `bar_2` barcode (18nt)
  
- **Two Independent Barcode Sets**:
  - Set 1 used for Array 1 (even rows)
  - Set 2 used for Array 2 (odd rows)
  
- **Global Uniqueness**: Ensures all `bar_1 + bar_2` combinations are globally unique across the entire array

- **No External Dependencies**: Pure Python implementation (no numpy required)

### 2. Generation Script: `generate_checkerboard_235x363.py`

A complete script to generate the 235×363 checkerboard layout:

- Loads two barcode sets from files
- Creates checkerboard layout with 8×8 blocks
- Validates global uniqueness
- Exports in column-major order (for photolithography)
- Generates statistics JSON file

### 3. Test Script: `test_checkerboard_small.py`

A verification script that:
- Tests the algorithm on a small 16×16 grid
- Verifies checkerboard pattern correctness
- Validates block consistency
- Confirms global uniqueness
- **Status**: ✅ Test passed successfully

### 4. Documentation: `README_CHECKERBOARD.md`

Comprehensive documentation including:
- Checkerboard pattern description
- Algorithm details
- Usage instructions
- API reference
- Troubleshooting guide

## Key Algorithm Features

### Checkerboard Pattern

```
Row 0: O . O . O . O .   (even columns only)
Row 1: . O . O . O . O   (odd columns only)
Row 2: O . O . O . O .   (even columns only)
Row 3: . O . O . O . O   (odd columns only)
```

For a 235×363 grid:
- **Array 1** (even rows): 118 rows × 182 cols = 21,476 spots
- **Array 2** (odd rows): 117 rows × 181 cols = 21,177 spots
- **Total spots**: 42,653 spots

### Block Structure

8×8 checkerboard blocks group spatially adjacent spots:
- Block size: 8 rows × 8 columns (in checkerboard space)
- Spots per block: ~32 (due to checkerboard pattern)
- All spots in a block share `bar_1` (spatial redundancy)
- Each spot has unique `bar_2` (identification)

### Barcode Assignment

1. Divide array into two interleaved sub-arrays (even/odd rows)
2. For each sub-array:
   - Group spots into 8×8 checkerboard blocks
   - Assign one `bar_1` (18nt) to entire block
   - Assign unique `bar_2` (18nt) to each spot in block
   - Track globally to ensure no `bar_1 + bar_2` repeats
3. Final barcode: `bar_1 + bar_2` = 36nt

## Usage Instructions

### Step 1: Prepare Barcode Files

Ensure you have two barcode files with 18nt barcodes:
```
Dual_barcodes/barcodes_118x182_dual_set1.txt  (~21,477 barcodes)
Dual_barcodes/barcodes_118x182_dual_set2.txt  (~21,477 barcodes)
```

### Step 2: Run Test (Optional but Recommended)

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Dual_barcodes
python test_checkerboard_small.py
```

Expected output: `TEST PASSED!`

### Step 3: Generate 235×363 Checkerboard Layout

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Dual_barcodes
python generate_checkerboard_235x363.py
```

This will create:
1. `layout_235x363_checkerboard_colmajor.csv` (~3-4 MB)
2. `layout_235x363_checkerboard_stats.json`

### Step 4: Verify Output

```python
import pandas as pd

# Load layout
df = pd.read_csv('layout_235x363_checkerboard_colmajor.csv')

# Check statistics
print(f"Total spots: {len(df)}")  # Should be ~42,653
print(f"Unique barcodes: {df['barcode'].nunique()}")  # Should equal total spots
print(f"Barcode length: {df['barcode'].str.len().iloc[0]}")  # Should be 36

# Verify checkerboard pattern
even_rows = df[df['row'] % 2 == 0]
odd_rows = df[df['row'] % 2 == 1]

print(f"Even rows have even columns: {(even_rows['col'] % 2 == 0).all()}")  # Should be True
print(f"Odd rows have odd columns: {(odd_rows['col'] % 2 == 1).all()}")  # Should be True
```

## Output Format

### CSV Layout (Column-Major Order)

```csv
row,col,bar_1,bar_2,barcode
0,0,ATGATACCCAGTGTCGAC,TGCGGGCTCATAGTTCCT,ATGATACCCAGTGTCGACTGCGGGCTCATAGTTCCT
0,2,ATGATACCCAGTGTCGAC,ACACTCGTGTGGAGCAAT,ATGATACCCAGTGTCGACACACTCGTGTGGAGCAAT
...
```

**Columns:**
- `row`: Row index (0-234)
- `col`: Column index (0-362)  
- `bar_1`: First barcode segment (18nt), shared within block
- `bar_2`: Second barcode segment (18nt), unique within block
- `barcode`: Combined barcode (36nt) = bar_1 + bar_2

**Order:** Column-major (iterate columns first, then rows) for photolithography synthesis

### Statistics JSON

Contains:
- Generation timestamp
- Array configuration (dimensions, pattern, spots)
- Barcode usage statistics
- Validation metrics

## Advantages of This Approach

1. **Spatial Error Tolerance**: Synthesis errors affecting `bar_1` group nearby spots together (within 8×8 block)

2. **Efficient Barcode Usage**: Only need ~2,700 unique 18nt barcodes total (instead of 42,653)

3. **Global Uniqueness**: Every spot has a unique 36nt identifier

4. **Checkerboard Compatible**: Properly handles the two interleaved sub-arrays

5. **No External Dependencies**: Pure Python (no numpy, scipy, etc.)

6. **Well-Tested**: Passed validation on 16×16 test grid

## Differences from Regular Dual Barcode Approach

| Feature | Regular 4×4 | Checkerboard 8×8 |
|---------|-------------|------------------|
| **Pattern** | Continuous grid | Two interleaved arrays |
| **Block size** | 4×4 = 16 spots | 8×8 checkerboard = ~32 spots |
| **Arrays** | Single array | Two independent arrays |
| **Barcode sets** | One set for both arrays | Two sets (one per array) |
| **Spatial grouping** | 4μm blocks | 8μm blocks |

## Files Created

```
Dual_barcodes/
├── dual_checkerboard_layout.py              # Core library
├── generate_checkerboard_235x363.py         # Generation script  
├── test_checkerboard_small.py               # Test script (✅ passed)
├── README_CHECKERBOARD.md                   # Full documentation
└── CHECKERBOARD_IMPLEMENTATION_SUMMARY.md   # This file
```

## Next Steps

1. **Run the generation script**:
   ```bash
   python generate_checkerboard_235x363.py
   ```

2. **Visualize activation patterns**:
   ```bash
   python ../plot_cycle_activation.py layout_235x363_checkerboard_colmajor.csv 0
   ```

3. **Validate with your decoder** (using the 250K decode test you mentioned)

4. **Compare results** with:
   - Dual barcodes 4×4 (Precision: 0.9731, Recall: 0.5703)
   - Zipcodes with separator (Precision: 0.8301, Recall: 0.7467)
   - Zipcodes without separator (Precision: 0.9014, Recall: 0.6039)

## Technical Notes

- **Block assignment**: Based on integer division of coordinates by block_size
- **Bar_1 selection**: Least-used barcode to balance usage
- **Bar_2 selection**: Random from available pool (not yet used with current bar_1)
- **Validation**: Checks block consistency, global uniqueness, and checkerboard pattern
- **Export order**: Column-major for compatibility with photolithography masks

## Questions or Issues?

If you encounter any problems:

1. Check that barcode files exist and have sufficient barcodes
2. Verify barcode lengths are 18nt
3. Run the test script to verify algorithm correctness
4. Check the statistics JSON for detailed metrics

## Summary

✅ **Checkerboard dual barcode approach is fully implemented and tested**
✅ **Ready to generate 235×363 layout with 8×8 blocks**
✅ **Uses two independent barcode sets for the two interleaved arrays**
✅ **Maintains global uniqueness while providing spatial error tolerance**
