#!/usr/bin/env python3
"""
Test script for checkerboard dual barcode layout with a small grid.

This verifies the algorithm works correctly before running on the full 235×363 array.
"""

import sys
from dual_checkerboard_layout import CheckerboardDualBarcodeLayout

def generate_test_barcodes(n, length, prefix):
    """Generate unique test barcodes"""
    import random
    bases = ['A', 'C', 'G', 'T']
    barcodes = set()
    
    # Use prefix for first part, then randomize the rest
    prefix_len = min(len(prefix), length // 2)
    
    attempts = 0
    while len(barcodes) < n and attempts < n * 100:
        attempts += 1
        # Start with prefix
        barcode = prefix * (prefix_len // len(prefix) + 1)
        barcode = barcode[:prefix_len]
        
        # Fill rest with random bases
        for _ in range(length - prefix_len):
            barcode += random.choice(bases)
        
        if len(barcode) == length:
            barcodes.add(barcode)
    
    if len(barcodes) < n:
        raise ValueError(f"Could only generate {len(barcodes)} unique barcodes, needed {n}")
    
    return list(barcodes)[:n]


def main():
    print("=" * 80)
    print("CHECKERBOARD DUAL BARCODE LAYOUT - SMALL TEST")
    print("=" * 80)
    
    # Small test: 16×16 grid with 8×8 blocks
    ROWS = 16
    COLS = 16
    BLOCK_SIZE = 8
    BARCODE_LENGTH = 18
    N_BARCODES = 100  # More than enough for testing
    
    print(f"\nTest configuration:")
    print(f"  Grid: {ROWS}×{COLS}")
    print(f"  Block size: {BLOCK_SIZE}×{BLOCK_SIZE} checkerboard")
    print(f"  Barcode length: {BARCODE_LENGTH}nt per segment")
    
    # Generate test barcodes
    print(f"\nGenerating test barcodes...")
    barcodes_set1 = generate_test_barcodes(N_BARCODES, BARCODE_LENGTH, "AT")
    barcodes_set2 = generate_test_barcodes(N_BARCODES, BARCODE_LENGTH, "GC")
    
    print(f"  Set 1: {len(barcodes_set1)} barcodes (prefix: AT)")
    print(f"    Example: {barcodes_set1[0]}")
    print(f"  Set 2: {len(barcodes_set2)} barcodes (prefix: GC)")
    print(f"    Example: {barcodes_set2[0]}")
    
    # Create layout
    print(f"\nCreating layout...")
    try:
        layout = CheckerboardDualBarcodeLayout(
            rows=ROWS,
            cols=COLS,
            barcodes_set1=barcodes_set1,
            barcodes_set2=barcodes_set2,
            block_size=BLOCK_SIZE,
            split_barcodes=False  # Test barcodes are already 18nt
        )
        
        # Generate
        print(f"\nGenerating layout...")
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=42)
        
        # Get statistics
        stats = layout.get_statistics()
        
        print(f"\n" + "=" * 80)
        print("TEST RESULTS")
        print("=" * 80)
        print(f"Total spots: {stats['total_spots']}")
        print(f"  Array 1 (even rows): {stats['array1_spots']}")
        print(f"  Array 2 (odd rows): {stats['array2_spots']}")
        print(f"\nBarcodes:")
        print(f"  Unique combined: {stats['unique_combined_barcodes']}")
        print(f"  Unique bar_1 used: {stats['unique_bar1_used']}")
        print(f"  Unique bar_2 used: {stats['unique_bar2_used']}")
        
        # Verify checkerboard pattern
        print(f"\nCheckerboard pattern verification:")
        print(f"  Checking first few positions...")
        
        expected_positions = [
            (0, 0), (0, 2), (0, 4), (0, 6),  # Row 0: even columns
            (1, 1), (1, 3), (1, 5), (1, 7),  # Row 1: odd columns
            (2, 0), (2, 2), (2, 4), (2, 6),  # Row 2: even columns
        ]
        
        for pos in expected_positions:
            if pos in grid_combined:
                print(f"    ✓ Position {pos}: {grid_combined[pos][:20]}...")
            else:
                print(f"    ✗ Position {pos}: MISSING!")
        
        # Check wrong positions (should not exist)
        wrong_positions = [
            (0, 1), (0, 3),  # Row 0 should not have odd columns
            (1, 0), (1, 2),  # Row 1 should not have even columns
        ]
        
        for pos in wrong_positions:
            if pos not in grid_combined:
                print(f"    ✓ Position {pos}: correctly absent")
            else:
                print(f"    ✗ Position {pos}: should not exist!")
        
        # Verify UNIFIED spatial block consistency
        print(f"\nUnified spatial block consistency check:")
        # Check first 8x8 spatial block (rows 0-7, cols 0-7) - includes BOTH arrays
        block1_array1 = [(r, c) for r in [0,2,4,6] for c in [0,2,4,6]]  # Even positions
        block1_array2 = [(r, c) for r in [1,3,5,7] for c in [1,3,5,7]]  # Odd positions
        block1_all_positions = block1_array1 + block1_array2
        
        block1_bar1 = set(grid_bar1.get(pos) for pos in block1_all_positions if pos in grid_bar1)
        block1_bar2 = set(grid_bar2.get(pos) for pos in block1_all_positions if pos in grid_bar2)
        
        print(f"  Spatial Block 1 (rows 0-7, cols 0-7):")
        print(f"    Array 1 spots: {len([p for p in block1_array1 if p in grid_bar1])}")
        print(f"    Array 2 spots: {len([p for p in block1_array2 if p in grid_bar1])}")
        print(f"    Total spots: {len([p for p in block1_all_positions if p in grid_bar1])}")
        print(f"    Unique bar_1: {len(block1_bar1)} (should be 1)")
        print(f"    Unique bar_2: {len(block1_bar2)} (should equal total spots)")
        
        if len(block1_bar1) == 1:
            print(f"    ✓ All spots (both arrays) share same bar_1")
        else:
            print(f"    ✗ ERROR: Multiple bar_1 in spatial block!")
            print(f"    bar_1 values found: {block1_bar1}")
        
        # Also check that adjacent columns share bar_1
        print(f"\n  Adjacent column check (columns 0 and 1):")
        col0_pos = [(0, 0), (2, 0), (4, 0), (6, 0)]  # Array 1
        col1_pos = [(1, 1), (3, 1), (5, 1), (7, 1)]  # Array 2
        col0_bar1 = set(grid_bar1.get(pos) for pos in col0_pos if pos in grid_bar1)
        col1_bar1 = set(grid_bar1.get(pos) for pos in col1_pos if pos in grid_bar1)
        
        if col0_bar1 == col1_bar1 and len(col0_bar1) == 1:
            print(f"    ✓ Columns 0 and 1 share the same bar_1 (spatial block consistency)")
        else:
            print(f"    ✗ ERROR: Columns 0 and 1 have different bar_1!")
            print(f"    Column 0 bar_1: {col0_bar1}")
            print(f"    Column 1 bar_1: {col1_bar1}")
        
        # Export test layout
        output_file = 'test_checkerboard_16x16.csv'
        print(f"\nExporting test layout to {output_file}...")
        layout.export_layout(output_file, format='csv', column_major=True)
        
        print(f"\n" + "=" * 80)
        print("TEST PASSED!")
        print("=" * 80)
        print(f"\nYou can now safely run: python generate_checkerboard_235x363.py")
        
    except Exception as e:
        print(f"\n" + "=" * 80)
        print("TEST FAILED!")
        print("=" * 80)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
