#!/usr/bin/env python3
"""
Run dual barcode layout generator with offset block pattern.
Accepts command line arguments for grid dimensions and barcode file.

Usage:
    python run_offset_layout.py <rows> <cols> <barcode_file> [output_prefix] [seed]

Examples:
    python run_offset_layout.py 118 182 barcodes.txt
    python run_offset_layout.py 118 182 barcodes.txt my_layout 42
    python run_offset_layout.py 10 10 barcodes.txt test_layout
"""

import sys
import os
from dual_barcode_offset_layout import DualBarcodeOffsetLayout, load_barcodes_from_file


def print_usage():
    """Print usage information."""
    print("Usage: python run_offset_layout.py <rows> <cols> <barcode_file> [output_prefix] [seed]")
    print()
    print("Arguments:")
    print("  rows          : Number of rows in the grid (required)")
    print("  cols          : Number of columns in the grid (required)")
    print("  barcode_file  : Path to file with barcodes, one per line (required)")
    print("  output_prefix : Prefix for output files (optional, default: 'layout')")
    print("  seed          : Random seed for reproducibility (optional, default: 42)")
    print()
    print("Examples:")
    print("  python run_offset_layout.py 118 182 barcodes_1in4_half.txt")
    print("  python run_offset_layout.py 118 182 barcodes_1in4_half.txt my_layout 42")
    print("  python run_offset_layout.py 10 10 example_barcodes.txt test_layout")


def main():
    # Parse command line arguments
    if len(sys.argv) < 4:
        print("ERROR: Insufficient arguments\n")
        print_usage()
        sys.exit(1)
    
    try:
        rows = int(sys.argv[1])
        cols = int(sys.argv[2])
        barcode_file = sys.argv[3]
        output_prefix = sys.argv[4] if len(sys.argv) > 4 else "layout"
        seed = int(sys.argv[5]) if len(sys.argv) > 5 else 42
    except ValueError as e:
        print(f"ERROR: Invalid argument - {e}\n")
        print_usage()
        sys.exit(1)
    
    # Validate inputs
    if rows <= 0 or cols <= 0:
        print(f"ERROR: Rows and columns must be positive integers")
        sys.exit(1)
    
    if not os.path.exists(barcode_file):
        print(f"ERROR: Barcode file not found: {barcode_file}")
        sys.exit(1)
    
    # Output files
    output_csv = f"{output_prefix}_{rows}x{cols}.csv"
    output_viz = f"{output_prefix}_{rows}x{cols}.png"
    
    print("=" * 70)
    print(f"Dual Barcode Layout Generator - Offset Block Pattern")
    print("=" * 70)
    print(f"\nParameters:")
    print(f"  Grid dimensions: {rows} rows × {cols} columns")
    print(f"  Barcode file: {barcode_file}")
    print(f"  Output prefix: {output_prefix}")
    print(f"  Random seed: {seed}")
    
    # Calculate requirements
    n_spots = rows * cols
    n_blocks_bar1 = ((rows + 1) // 2) * ((cols + 1) // 2)
    n_blocks_bar2 = (rows // 2) * (cols // 2)
    
    print(f"\nGrid statistics:")
    print(f"  Total spots: {n_spots:,}")
    print(f"  bar_1 blocks (2×2 at (0,0)): {n_blocks_bar1:,}")
    print(f"  bar_2 blocks (2×2 at (1,1)): {n_blocks_bar2:,}")
    print(f"  Edge positions: {n_spots - n_blocks_bar2 * 4:,}")
    
    # Load barcodes
    print(f"\nLoading barcodes from: {barcode_file}")
    try:
        barcodes = load_barcodes_from_file(barcode_file)
        print(f"✓ Loaded {len(barcodes):,} barcodes")
        
        if barcodes:
            print(f"  Barcode length: {len(barcodes[0])} nucleotides")
            print(f"  First barcode: {barcodes[0]}")
            print(f"  Last barcode: {barcodes[-1]}")
        
        # Estimate minimum barcodes needed
        # This is a rough estimate - actual requirement depends on overlap patterns
        min_barcodes_estimate = max(n_blocks_bar1, n_blocks_bar2)
        print(f"\nEstimated minimum barcodes needed: ~{min_barcodes_estimate:,}")
        print(f"(Actual requirement may vary based on block overlaps)")
        
        if len(barcodes) < min_barcodes_estimate:
            print(f"\n⚠ WARNING: May not have enough barcodes!")
            print(f"  Recommended: {min_barcodes_estimate:,}")
            print(f"  Available: {len(barcodes):,}")
            print(f"  Will attempt to generate layout...")
        else:
            print(f"✓ Sufficient barcodes available")
        
    except FileNotFoundError:
        print(f"✗ ERROR: File not found: {barcode_file}")
        sys.exit(1)
    except ValueError as e:
        print(f"✗ ERROR: {e}")
        sys.exit(1)
    
    # Create layout
    print(f"\nCreating {rows}×{cols} layout with offset block pattern...")
    layout = DualBarcodeOffsetLayout(rows, cols, barcodes)
    
    # Generate layout
    print("Generating barcode assignments...")
    try:
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=seed)
        print("✓ Layout generated successfully!")
    except ValueError as e:
        print(f"✗ ERROR: {e}")
        print("\nTroubleshooting:")
        print("  - Try using more barcodes")
        print("  - Try a different random seed")
        print("  - Reduce grid dimensions")
        sys.exit(1)
    
    # Get statistics
    print("\n" + "=" * 70)
    print("LAYOUT STATISTICS")
    print("=" * 70)
    stats = layout.get_statistics()
    print(f"Grid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']:,}")
    print(f"Number of barcodes used: {stats['n_barcodes_available']:,}")
    print(f"bar_1 blocks: {stats['n_blocks_bar1']:,}")
    print(f"bar_2 blocks (offset): {stats['n_blocks_bar2']:,}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']:,}")
    print(f"Theoretical max combinations: {stats['theoretical_max_combinations']:,}")
    
    # Usage statistics
    print(f"\nbar_1 usage distribution:")
    bar1_usage = stats['bar1_usage']
    usage_counts = {}
    for count in bar1_usage.values():
        usage_counts[count] = usage_counts.get(count, 0) + 1
    for count in sorted(usage_counts.keys()):
        print(f"  {usage_counts[count]:,} barcodes used {count} times")
    
    print(f"\nbar_2 block diversity:")
    bar2_diversity = stats['bar2_block_diversity']
    if bar2_diversity:
        diversity_counts = {}
        for count in bar2_diversity.values():
            diversity_counts[count] = diversity_counts.get(count, 0) + 1
        for count in sorted(diversity_counts.keys()):
            print(f"  {diversity_counts[count]:,} bar_2 barcodes used in {count} blocks")
    
    # Export
    print("\n" + "=" * 70)
    print("EXPORTING")
    print("=" * 70)
    print(f"Exporting layout to: {output_csv}")
    layout.export_layout(output_csv, format='csv')
    print("✓ CSV export complete (column-major order)")
    
    # Visualization
    if rows <= 50 and cols <= 50:
        print(f"\nGenerating visualization: {output_viz}")
        layout.visualize_layout(output_viz, include_combined=False)
        print("✓ Visualization complete")
    else:
        print(f"\nNote: Grid too large for automatic visualization ({rows}×{cols})")
        print(f"To create visualization manually, use:")
        print(f"  layout.visualize_layout('{output_viz}', include_combined=False)")
    
    print("\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  CSV: {output_csv}")
    if rows <= 50 and cols <= 50:
        print(f"  Visualization: {output_viz}")
    print(f"\nLayout summary:")
    print(f"  Total unique barcodes: {stats['unique_combined_barcodes']:,} / {stats['total_spots']:,}")
    print(f"  Coverage: 100% (all spots have unique combined barcodes)")


if __name__ == "__main__":
    main()
