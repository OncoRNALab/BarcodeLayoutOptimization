# Dual Barcode Layout Generator

A Python implementation for generating optimal barcode layouts for spatial transcriptomics in-situ arrays using photolithography synthesis.

## Problem Statement

Arrays synthesized using photolithography require careful barcode design to:
- **Minimize cross-talk** during sequencing (barcodes should be physically distant)
- **Reduce synthesis errors** (neighboring spots should share sequence similarity)

## Solution: Dual Barcode Approach

This implementation uses a dual barcode strategy:

1. **bar_1**: Assigned to each 2×2 block (all 4 spots share the same bar_1)
   - Reduces synthesis errors by grouping similar sequences
   
2. **bar_2**: Assigned uniquely within each 2×2 block (all 4 spots have different bar_2)
   - Prevents cross-talk by ensuring unique identifiers

3. **Final barcode**: Concatenation of bar_1 + bar_2
   - Creates unique identifiers while maintaining spatial structure

## Key Features

- ✓ Uses the **same barcode set** for both bar_1 and bar_2
- ✓ Guarantees unique barcodes within each 2×2 block
- ✓ Handles edge cases with odd-dimension grids
- ✓ Validates all constraints automatically
- ✓ Exports to CSV/TSV formats
- ✓ Includes visualization capabilities

## Installation

```bash
# Required
pip install numpy

# Optional (for visualization)
pip install matplotlib
```

## Quick Start

```python
from dual_barcode_layout import DualBarcodeLayout, generate_random_barcodes

# Generate barcodes
barcodes = generate_random_barcodes(n=6, length=8, seed=42)

# Create layout
layout = DualBarcodeLayout(rows=8, cols=8, barcodes=barcodes)
grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=42)

# Get statistics
stats = layout.get_statistics()
print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")

# Export
layout.export_layout('output.csv', format='csv')

# Visualize
layout.visualize_layout('output.png')
```

## Algorithm Overview

### Step 1: bar_1 Assignment (Block-wise)
- Divide grid into 2×2 blocks
- Assign one bar_1 barcode to all 4 positions in each block
- Use cycling through available barcodes for multiple blocks

### Step 2: bar_2 Assignment (Position-wise)
- For each 2×2 block, randomly select 4 different bar_2 barcodes
- Assign one to each position within the block
- Ensures all positions in a block have unique bar_2

### Step 3: Combination
- Concatenate bar_1 + bar_2 for final barcode
- Guarantees uniqueness within each block

### Step 4: Validation
- Verify all 2×2 blocks have same bar_1
- Verify all 2×2 blocks have different bar_2
- Verify all combined barcodes are unique within blocks

## Minimum Requirements

- **Minimum 4 barcodes** required to fill any grid
- With N barcodes: can create N² unique combinations
- Optimal capacity: 4N spots (N blocks × 4 spots/block)

## Examples

See `example_usage.py` for comprehensive examples including:
1. Small 4×4 grid with minimum barcodes
2. Rectangular grids
3. Large arrays (20×20)
4. Edge cases with odd dimensions
5. Constraint validation
6. Multiple runs with different seeds

## File Structure

```
Dual_barcodes/
├── dual_barcode_layout.py      # Main implementation
├── example_usage.py             # Usage examples
├── README.md                    # This file
└── (generated outputs)
    ├── barcode_layout.csv
    └── barcode_layout_visualization.png
```

## API Reference

### `DualBarcodeLayout`

Main class for generating layouts.

**Constructor:**
```python
DualBarcodeLayout(rows: int, cols: int, barcodes: List[str])
```

**Methods:**
- `generate_layout(seed=None)`: Generate the complete layout
- `get_statistics()`: Get layout statistics
- `export_layout(filepath, format='csv')`: Export to file
- `visualize_layout(output_file=None)`: Create visualization

### `generate_random_barcodes`

Utility function to generate random DNA barcodes.

```python
generate_random_barcodes(n: int, length: int = 8, seed: Optional[int] = None) -> List[str]
```

## Usage Tips

1. **Barcode count**: Use at least 4 barcodes, more provides better diversity
2. **Grid size**: Works with any dimensions, automatically handles edges
3. **Reproducibility**: Set seed for reproducible layouts
4. **Validation**: Automatic validation ensures all constraints are met
5. **Same barcode set**: Using same set for bar_1 and bar_2 simplifies synthesis

## Performance

- Time complexity: O(rows × cols)
- Space complexity: O(rows × cols)
- Very fast even for large grids (e.g., 100×100 completes in milliseconds)

## Advantages Over Other Approaches

1. **Block-wise assignment**: More efficient than backtracking algorithms
2. **Same barcode set**: Reduces synthesis complexity
3. **Guaranteed solution**: Always finds valid layout (with ≥4 barcodes)
4. **Edge handling**: Automatically handles incomplete blocks at boundaries
5. **Validation**: Built-in verification of all constraints

## License

This implementation is provided for research purposes in spatial transcriptomics array design.

## Contact

For questions or issues, please refer to the project documentation or create an issue in the repository.
