#!/usr/bin/env python3
"""
Verify the generated dual barcode layout for correctness

Checks:
1. Total number of barcodes = 85,305
2. All combined barcodes are unique
3. All barcodes are 36nt long
4. Each 4×4 block has exactly 1 unique bar_1
5. Each 4×4 block has unique bar_2 values
6. No bar_1+bar_2 combination repeats globally
"""

import sys
import pandas as pd
import numpy as np

def verify_layout(csv_file='layout_235x363_dual_4x4.csv'):
    """Verify the generated layout"""
    
    print("=" * 80)
    print("DUAL BARCODE LAYOUT VERIFICATION")
    print("=" * 80)
    print(f"\nReading layout from: {csv_file}")
    
    # Load layout
    df = pd.read_csv(csv_file)
    print(f"✓ Loaded {len(df)} rows")
    
    # Test 1: Total count
    print("\n[Test 1] Checking total number of barcodes...")
    expected_count = 235 * 363
    actual_count = len(df)
    assert actual_count == expected_count, f"Expected {expected_count}, got {actual_count}"
    print(f"✓ Total barcodes: {actual_count} (235 × 363)")
    
    # Test 2: Uniqueness
    print("\n[Test 2] Checking barcode uniqueness...")
    unique_combined = df['combined'].nunique()
    assert unique_combined == len(df), f"Found duplicates: {len(df) - unique_combined}"
    print(f"✓ All {unique_combined} combined barcodes are unique")
    
    # Test 3: Barcode length
    print("\n[Test 3] Checking barcode lengths...")
    lengths = df['combined'].str.len().unique()
    assert len(lengths) == 1 and lengths[0] == 36, f"Expected 36nt, found: {lengths}"
    bar1_len = df['bar_1'].str.len().unique()[0]
    bar2_len = df['bar_2'].str.len().unique()[0]
    print(f"✓ All barcodes are 36nt ({bar1_len}nt bar_1 + {bar2_len}nt bar_2)")
    
    # Test 4: Block structure - bar_1 consistency
    print("\n[Test 4] Checking bar_1 consistency within blocks...")
    errors = 0
    for block_row in range(59):  # (235 + 3) // 4 = 59
        for block_col in range(91):  # (363 + 3) // 4 = 91
            # Get all positions in this block
            row_start, row_end = block_row * 4, min((block_row + 1) * 4, 235)
            col_start, col_end = block_col * 4, min((block_col + 1) * 4, 363)
            
            block_data = df[
                (df['row'] >= row_start) & (df['row'] < row_end) &
                (df['col'] >= col_start) & (df['col'] < col_end)
            ]
            
            unique_bar1 = block_data['bar_1'].nunique()
            if unique_bar1 != 1:
                print(f"  ✗ Block ({block_row},{block_col}): {unique_bar1} different bar_1 values")
                errors += 1
                if errors >= 5:
                    print(f"  ... (stopping after 5 errors)")
                    break
        if errors >= 5:
            break
    
    if errors == 0:
        print(f"✓ All 5,369 blocks have consistent bar_1 within each block")
    else:
        raise AssertionError(f"Found {errors} blocks with inconsistent bar_1")
    
    # Test 5: Block structure - bar_2 uniqueness
    print("\n[Test 5] Checking bar_2 uniqueness within blocks...")
    errors = 0
    for block_row in range(59):
        for block_col in range(91):
            row_start, row_end = block_row * 4, min((block_row + 1) * 4, 235)
            col_start, col_end = block_col * 4, min((block_col + 1) * 4, 363)
            
            block_data = df[
                (df['row'] >= row_start) & (df['row'] < row_end) &
                (df['col'] >= col_start) & (df['col'] < col_end)
            ]
            
            n_positions = len(block_data)
            unique_bar2 = block_data['bar_2'].nunique()
            
            if unique_bar2 != n_positions:
                print(f"  ✗ Block ({block_row},{block_col}): {unique_bar2} unique bar_2 for {n_positions} positions")
                errors += 1
                if errors >= 5:
                    print(f"  ... (stopping after 5 errors)")
                    break
        if errors >= 5:
            break
    
    if errors == 0:
        print(f"✓ All blocks have unique bar_2 values for each position")
    else:
        raise AssertionError(f"Found {errors} blocks with duplicate bar_2")
    
    # Test 6: Global uniqueness of bar_1+bar_2 combinations
    print("\n[Test 6] Checking global uniqueness of bar_1+bar_2 combinations...")
    df['combination_key'] = df['bar_1'] + '|' + df['bar_2']
    duplicates = df[df['combination_key'].duplicated()]
    if len(duplicates) > 0:
        print(f"  ✗ Found {len(duplicates)} duplicate bar_1+bar_2 combinations")
        print(duplicates.head())
        raise AssertionError("Duplicate combinations found")
    print(f"✓ All 85,305 bar_1+bar_2 combinations are globally unique")
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print(f"Total spots: {len(df):,}")
    print(f"Unique bar_1 values: {df['bar_1'].nunique():,}")
    print(f"Unique bar_2 values: {df['bar_2'].nunique():,}")
    print(f"Unique combined barcodes: {df['combined'].nunique():,}")
    
    print("\nbar_1 usage:")
    bar1_usage = df['bar_1'].value_counts()
    print(f"  Min usage: {bar1_usage.min()} spots")
    print(f"  Max usage: {bar1_usage.max()} spots")
    print(f"  Mean usage: {bar1_usage.mean():.1f} spots")
    
    print("\nbar_2 usage:")
    bar2_usage = df['bar_2'].value_counts()
    print(f"  Min usage: {bar2_usage.min()} spots")
    print(f"  Max usage: {bar2_usage.max()} spots")
    print(f"  Mean usage: {bar2_usage.mean():.1f} spots")
    
    # Example block
    print("\n" + "=" * 80)
    print("EXAMPLE: Block (0,0) - First 4×4 block")
    print("=" * 80)
    block_00 = df[(df['row'] < 4) & (df['col'] < 4)]
    print(f"\nAll 16 positions share bar_1: {block_00['bar_1'].iloc[0]}")
    print(f"Each position has unique bar_2:")
    for idx, row in block_00.iterrows():
        print(f"  ({row['row']},{row['col']}): {row['bar_2']}")
    
    print("\n" + "=" * 80)
    print("✓ ALL TESTS PASSED - Layout is valid!")
    print("=" * 80)

if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else 'layout_235x363_dual_4x4.csv'
    try:
        verify_layout(csv_file)
    except Exception as e:
        print(f"\n✗ VERIFICATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
