#!/usr/bin/env python3
"""
Create Checkerboard Layouts - Alternating Filled and Empty Spots

This script creates TRUE checkerboard layouts where:
- Set 1 occupies positions where (row + col) % 2 == 0 (even squares) → ~50% of array
- Set 2 occupies positions where (row + col) % 2 == 1 (odd squares) → other 50% of array

Each set has EMPTY spots in a checkerboard pattern, creating two independent
interleaved layouts that can be used separately or combined.

Usage:
    python create_checkerboard_layouts.py <set1_layout.csv> <set2_layout.csv>
    
Or run with defaults:
    python create_checkerboard_layouts.py
"""

import sys
import pandas as pd
from pathlib import Path


def load_layout(filepath):
    """Load a barcode layout CSV file"""
    df = pd.read_csv(filepath)
    
    # Verify required columns
    required = ['row', 'col', 'barcode']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Layout file must have columns: {required}")
    
    return df


def create_checkerboard_set(df_source, rows, cols, parity, set_name):
    """
    Create checkerboard layout for one set.
    
    Args:
        df_source: Source DataFrame with barcodes
        rows: Number of rows in array
        cols: Number of columns in array
        parity: 0 for even squares, 1 for odd squares
        set_name: Name of the set (for reporting)
    
    Returns:
        DataFrame with checkerboard layout (only filled positions)
    """
    # Create lookup dictionary
    barcode_dict = {(row, col): barcode 
                    for row, col, barcode in zip(df_source['row'], df_source['col'], df_source['barcode'])}
    
    # Build checkerboard layout
    layout = []
    used_count = 0
    expected_count = 0
    
    for row in range(rows):
        for col in range(cols):
            # Check if this position belongs to this set's checkerboard pattern
            if (row + col) % 2 == parity:
                expected_count += 1
                barcode = barcode_dict.get((row, col))
                
                if barcode is None:
                    print(f"  ⚠ Warning: Missing barcode for position ({row}, {col})")
                    continue
                
                layout.append({
                    'row': row,
                    'col': col,
                    'barcode': barcode
                })
                used_count += 1
    
    df_layout = pd.DataFrame(layout)
    
    print(f"  {set_name}: {used_count:,} positions filled (out of {expected_count:,} in checkerboard)")
    
    return df_layout


def verify_checkerboard(df, parity, set_name):
    """Verify all positions follow the checkerboard pattern"""
    errors = []
    
    for idx, row in df.iterrows():
        r, c = row['row'], row['col']
        if (r + c) % 2 != parity:
            errors.append(f"Position ({r},{c}): violates {set_name} checkerboard pattern")
    
    if errors:
        print(f"  ⚠ Warning: {len(errors)} pattern violations:")
        for err in errors[:5]:
            print(f"    {err}")
        if len(errors) > 5:
            print(f"    ... and {len(errors)-5} more")
        return False
    else:
        print(f"  ✓ {set_name} checkerboard pattern verified")
        return True


def export_checkerboard(df, output_name, rows, cols, parity, set_name):
    """Export checkerboard layout in multiple formats"""
    
    # 1. CSV with filled positions only (column-major for fabrication)
    csv_path = f"{output_name}_colmajor.csv"
    df_sorted = df.sort_values(['col', 'row'])
    df_sorted.to_csv(csv_path, index=False)
    print(f"  ✓ Layout CSV: {csv_path}")
    
    # 2. Barcodes only (column-major order)
    barcodes_path = f"{output_name}_barcodes.txt"
    with open(barcodes_path, 'w') as f:
        for barcode in df_sorted['barcode']:
            f.write(f"{barcode}\n")
    print(f"  ✓ Barcodes: {barcodes_path}")
    
    # 3. Pattern visualization
    pattern_path = f"{output_name}_pattern.txt"
    with open(pattern_path, 'w') as f:
        f.write(f"{set_name} Checkerboard Pattern\n")
        f.write("=" * 82 + "\n")
        f.write("Legend: # = filled spot, . = empty spot\n")
        f.write("=" * 82 + "\n\n")
        
        # Create lookup for quick checks
        filled = set((r, c) for r, c in zip(df['row'], df['col']))
        
        vis_rows = min(30, rows)
        vis_cols = min(60, cols)
        
        # Header
        f.write("    ")
        for c in range(vis_cols):
            f.write(f"{c%10}")
        f.write("\n    " + "-" * vis_cols + "\n")
        
        # Grid
        for r in range(vis_rows):
            f.write(f"{r:3d}|")
            for c in range(vis_cols):
                if (r, c) in filled:
                    f.write('#')
                else:
                    f.write('.')
            f.write("\n")
        
        if rows > vis_rows or cols > vis_cols:
            f.write(f"\n(Showing first {vis_rows} rows × {vis_cols} cols of {rows} × {cols} array)\n")
        
        f.write(f"\nStatistics:\n")
        f.write(f"  Filled spots: {len(df):,}\n")
        f.write(f"  Empty spots: {rows * cols - len(df):,}\n")
        f.write(f"  Fill ratio: {len(df) / (rows * cols) * 100:.1f}%\n")
    
    print(f"  ✓ Pattern: {pattern_path}")


def create_combined_visualization(df1, df2, output_name, rows, cols):
    """Create a visualization showing both sets together"""
    
    pattern_path = f"{output_name}_combined_pattern.txt"
    with open(pattern_path, 'w') as f:
        f.write("COMBINED CHECKERBOARD PATTERN\n")
        f.write("=" * 82 + "\n")
        f.write("Legend: 1 = Set 1, 2 = Set 2, . = error (overlap or gap)\n")
        f.write("=" * 82 + "\n\n")
        
        # Create lookups
        set1_pos = set((r, c) for r, c in zip(df1['row'], df1['col']))
        set2_pos = set((r, c) for r, c in zip(df2['row'], df2['col']))
        
        vis_rows = min(30, rows)
        vis_cols = min(60, cols)
        
        # Header
        f.write("    ")
        for c in range(vis_cols):
            f.write(f"{c%10}")
        f.write("\n    " + "-" * vis_cols + "\n")
        
        # Grid
        errors = 0
        for r in range(vis_rows):
            f.write(f"{r:3d}|")
            for c in range(vis_cols):
                in_set1 = (r, c) in set1_pos
                in_set2 = (r, c) in set2_pos
                
                if in_set1 and in_set2:
                    f.write('X')  # Overlap error
                    errors += 1
                elif in_set1:
                    f.write('1')
                elif in_set2:
                    f.write('2')
                else:
                    f.write('.')  # Gap error
                    errors += 1
            f.write("\n")
        
        if rows > vis_rows or cols > vis_cols:
            f.write(f"\n(Showing first {vis_rows} rows × {vis_cols} cols of {rows} × {cols} array)\n")
        
        # Statistics
        overlap = set1_pos & set2_pos
        total_filled = len(set1_pos | set2_pos)
        
        f.write(f"\nStatistics:\n")
        f.write(f"  Set 1 spots: {len(set1_pos):,}\n")
        f.write(f"  Set 2 spots: {len(set2_pos):,}\n")
        f.write(f"  Total filled: {total_filled:,} / {rows * cols:,}\n")
        f.write(f"  Coverage: {total_filled / (rows * cols) * 100:.1f}%\n")
        
        if overlap:
            f.write(f"  ⚠ Overlaps: {len(overlap)}\n")
        else:
            f.write(f"  ✓ No overlaps\n")
        
        if total_filled == rows * cols:
            f.write(f"  ✓ Complete coverage (no gaps)\n")
        else:
            f.write(f"  ⚠ Gaps: {rows * cols - total_filled}\n")
    
    print(f"  ✓ Combined pattern: {pattern_path}")


def main():
    """Main execution"""
    
    print("=" * 80)
    print("CHECKERBOARD LAYOUT GENERATOR")
    print("Creates true checkerboard patterns with empty spots")
    print("=" * 80)
    print()
    
    # Parse arguments or use defaults
    if len(sys.argv) >= 3:
        set1_file = sys.argv[1]
        set2_file = sys.argv[2]
    else:
        # Default files
        set1_file = "output_118x182_nosep_set1/layout_118x182_nosep_set1_colmajor.csv"
        set2_file = "output_118x182_nosep_set2/layout_118x182_nosep_set2_colmajor.csv"
        
        print(f"Using default files:")
        print(f"  Set 1: {set1_file}")
        print(f"  Set 2: {set2_file}")
        print()
    
    # Load layouts
    print("[1/6] Loading source layouts...")
    try:
        df_set1 = load_layout(set1_file)
        print(f"  ✓ Set 1 loaded: {len(df_set1):,} barcodes")
        
        df_set2 = load_layout(set2_file)
        print(f"  ✓ Set 2 loaded: {len(df_set2):,} barcodes")
    except Exception as e:
        print(f"  ✗ Error: {e}")
        sys.exit(1)
    
    # Determine array dimensions
    rows = max(df_set1['row'].max(), df_set2['row'].max()) + 1
    cols = max(df_set1['col'].max(), df_set2['col'].max()) + 1
    total_spots = rows * cols
    
    print(f"\nArray dimensions: {rows} × {cols} = {total_spots:,} spots")
    expected_per_set = (total_spots + 1) // 2  # Round up for odd dimensions
    print(f"Expected per checkerboard set: ~{expected_per_set:,} spots")
    
    # Create Set 1 checkerboard (even squares)
    print(f"\n[2/6] Creating Set 1 checkerboard (even squares: (0,0), (0,2), (1,1)...)...")
    df_cb1 = create_checkerboard_set(df_set1, rows, cols, parity=0, set_name="Set 1")
    
    # Create Set 2 checkerboard (odd squares)
    print(f"\n[3/6] Creating Set 2 checkerboard (odd squares: (0,1), (1,0), (1,2)...)...")
    df_cb2 = create_checkerboard_set(df_set2, rows, cols, parity=1, set_name="Set 2")
    
    # Verify patterns
    print(f"\n[4/6] Verifying checkerboard patterns...")
    verify_checkerboard(df_cb1, parity=0, set_name="Set 1")
    verify_checkerboard(df_cb2, parity=1, set_name="Set 2")
    
    # Check for overlaps
    set1_pos = set((r, c) for r, c in zip(df_cb1['row'], df_cb1['col']))
    set2_pos = set((r, c) for r, c in zip(df_cb2['row'], df_cb2['col']))
    overlap = set1_pos & set2_pos
    total_coverage = len(set1_pos | set2_pos)
    
    if overlap:
        print(f"  ⚠ Warning: {len(overlap)} overlapping positions!")
    else:
        print(f"  ✓ No overlaps between sets")
    
    print(f"  Total coverage: {total_coverage:,} / {total_spots:,} spots ({total_coverage/total_spots*100:.1f}%)")
    
    # Export Set 1 checkerboard
    print(f"\n[5/6] Exporting Set 1 checkerboard...")
    export_checkerboard(df_cb1, "layout_118x182_checkerboard_set1", rows, cols, 0, "Set 1")
    
    # Export Set 2 checkerboard
    print(f"\n[6/6] Exporting Set 2 checkerboard...")
    export_checkerboard(df_cb2, "layout_118x182_checkerboard_set2", rows, cols, 1, "Set 2")
    
    # Create combined visualization
    print(f"\n[Bonus] Creating combined visualization...")
    create_combined_visualization(df_cb1, df_cb2, "layout_118x182_checkerboard", rows, cols)
    
    # Summary
    print("\n" + "=" * 80)
    print("CHECKERBOARD LAYOUTS CREATED")
    print("=" * 80)
    print(f"\nArray: {rows} × {cols} = {total_spots:,} spots")
    print(f"\nSet 1 (even squares):")
    print(f"  • Filled: {len(df_cb1):,} spots")
    print(f"  • Empty: {total_spots - len(df_cb1):,} spots")
    print(f"  • Files: layout_118x182_checkerboard_set1_*")
    
    print(f"\nSet 2 (odd squares):")
    print(f"  • Filled: {len(df_cb2):,} spots")
    print(f"  • Empty: {total_spots - len(df_cb2):,} spots")
    print(f"  • Files: layout_118x182_checkerboard_set2_*")
    
    print(f"\nUsage:")
    print(f"  • View patterns: Open *_pattern.txt files")
    print(f"  • Fabricate Set 1 only: Use layout_118x182_checkerboard_set1_colmajor.csv")
    print(f"  • Fabricate Set 2 only: Use layout_118x182_checkerboard_set2_colmajor.csv")
    print(f"  • Fabricate both: Run two separate fabrications")
    print("=" * 80)


if __name__ == "__main__":
    main()
