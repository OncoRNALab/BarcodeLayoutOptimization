# Dual Barcode Layout Generator - 235×363 Array

## Overview

This directory contains scripts and configuration for generating a dual barcode layout for a 235×363 DNA in-situ array (85,305 spots).

## Files

- **`generate_layout_235x363.py`** - Main generation script
- **`config_235x363.json`** - Configuration file with all parameters
- **`dual_barcode_layout.py`** - Core dual barcode layout algorithm (reusable library)

## Array Specifications

| Parameter | Value |
|-----------|-------|
| Rows | 235 |
| Columns | 363 |
| Total spots | 85,305 |
| Block size | 4×4 |
| Total blocks | 5,369 (59×91) |
| Barcode segments | 18nt + 18nt |
| Full barcode length | 36nt |
| Barcodes available | 85,305 |

## Algorithm

**Dual Barcodes Approach:**
- Array is divided into 4×4 blocks
- Each block shares one `bar_1` barcode (18nt)
- Each position within a block receives a unique `bar_2` barcode (18nt)
- Final barcode = `bar_1 + bar_2` (36nt total)
- **Global uniqueness:** No `bar_1 + bar_2` combination repeats across the entire array

**Key Properties:**
- Spatially adjacent spots within the same block share `bar_1`
- All 85,305 final barcodes are unique
- Optimized for synthesis accuracy (shared sequences within blocks)

## Usage

### Basic execution:
```bash
python generate_layout_235x363.py
```

### With custom config:
```bash
python generate_layout_235x363.py path/to/custom_config.json
```

### Expected runtime:
- Layout generation: ~5-30 seconds
- CSV export: ~10-20 seconds
- Visualization: ~30-60 seconds
- **Total: ~1-2 minutes**

## Output Files

After successful execution, the following files will be created:

1. **`layout_235x363_dual_4x4.csv`** - Full layout in CSV format
   - Columns: `row, col, bar_1, bar_2, combined_barcode`
   - 85,305 rows (one per spot)
   - ~8-12 MB file size

2. **`layout_235x363_dual_4x4_stats.json`** - Statistics and metadata
   - Generation timestamp
   - Configuration used
   - Usage statistics for bar_1 and bar_2
   - Validation metrics

3. **`layout_235x363_dual_4x4_visualization.png`** - Visual representation
   - Heatmap showing spatial barcode patterns
   - Block structure visualization
   - Optional (requires matplotlib)

## Configuration

Edit `config_235x363.json` to customize:
- Array dimensions
- Block size
- Barcode file path
- Output file names
- Random seed for reproducibility

## Example Output Structure

CSV layout format:
```csv
row,col,bar_1,bar_2,combined_barcode
0,0,CTTCTTCTCTGGTCTCAT,ACTCTATTGTTATGCCAC,CTTCTTCTCTGGTCTCATACTCTATTGTTATGCCAC
0,1,CTTCTTCTCTGGTCTCAT,ATCACTGTGGTCTCGCCA,CTTCTTCTCTGGTCTCATATCACTGTGGTCTCGCCA
0,2,CTTCTTCTCTGGTCTCAT,AGCATCGTTTACGTCCCA,CTTCTTCTCTGGTCTCATAGCATCGTTTACGTCCCA
...
```

## Verification

To verify the generated layout:
```python
import pandas as pd

# Load layout
df = pd.read_csv('layout_235x363_dual_4x4.csv')

# Check dimensions
print(f"Total barcodes: {len(df)}")  # Should be 85,305

# Check uniqueness
print(f"Unique combined: {df['combined_barcode'].nunique()}")  # Should be 85,305

# Check barcode length
print(f"Barcode length: {df['combined_barcode'].str.len().unique()}")  # Should be [36]

# Verify block structure (e.g., first block at 0,0)
block_0_0 = df[(df['row'] < 4) & (df['col'] < 4)]
print(f"Unique bar_1 in block (0,0): {block_0_0['bar_1'].nunique()}")  # Should be 1
print(f"Unique bar_2 in block (0,0): {block_0_0['bar_2'].nunique()}")  # Should be 16
```

## Troubleshooting

**Issue:** "Insufficient barcodes" error
- **Solution:** Check that barcode file has at least 5,369 unique barcodes

**Issue:** "Cannot find enough unique bar_2 barcodes"
- **Solution:** Increase the total number of available barcodes

**Issue:** Visualization fails
- **Solution:** Install matplotlib: `pip install matplotlib`
- **Note:** Visualization is optional; the CSV layout is still valid

## Requirements

- Python 3.7+
- NumPy
- Matplotlib (optional, for visualization only)

## Contact

For questions or issues, refer to the main project README or contact the development team.
