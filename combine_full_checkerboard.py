#!/usr/bin/env python3
"""
Combine Two Full Barcode Sets into Checkerboard Pattern

Takes two complete barcode sets (21,476 each = 42,952 total) and places them
in a larger grid with checkerboard pattern (50% filled).

Pattern:
- Even columns (0,2,4...): Set 1 barcodes at EVEN rows (0,2,4...)
- Odd columns (1,3,5...): Set 2 barcodes at ODD rows (1,3,5...)

This creates a checkerboard where adjacent spots are always empty, and the
two sets never touch each other.

Usage:
    python combine_full_checkerboard.py <set1_layout.csv> <set2_layout.csv> <rows> <cols>
    
Example:
    python combine_full_checkerboard.py set1.csv set2.csv 236 364
"""

import sys
import pandas as pd
import numpy as np


def load_layout_ordered(filepath):
    """Load layout in column-major order"""
    df = pd.read_csv(filepath)
    
    required = ['row', 'col', 'barcode']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Layout file must have columns: {required}")
    
    # Sort by column, then row (column-major)
    df_sorted = df.sort_values(['col', 'row'])
    barcodes = df_sorted['barcode'].tolist()
    
    return barcodes


def create_checkerboard_grid(barcodes_set1, barcodes_set2, target_rows, target_cols):
    """
    Create checkerboard grid with both sets.
    
    Pattern:
    - Even columns: Set 1 at even rows
    - Odd columns: Set 2 at odd rows
    
    Args:
        barcodes_set1: List of barcodes from Set 1 (column-major order)
        barcodes_set2: List of barcodes from Set 2 (column-major order)
        target_rows: Number of rows in output grid
        target_cols: Number of columns in output grid
    
    Returns:
        DataFrame with combined layout
    """
    layout = []
    idx_set1 = 0
    idx_set2 = 0
    
    set1_used = 0
    set2_used = 0
    
    # Process column by column (column-major)
    for col in range(target_cols):
        for row in range(target_rows):
            # Determine if this position should have a barcode
            if col % 2 == 0:  # Even column
                # Set 1 at even rows
                if row % 2 == 0:
                    if idx_set1 >= len(barcodes_set1):
                        raise ValueError(f"Ran out of Set 1 barcodes at position ({row},{col})")
                    
                    layout.append({
                        'row': row,
                        'col': col,
                        'barcode': barcodes_set1[idx_set1],
                        'source': 'set1'
                    })
                    idx_set1 += 1
                    set1_used += 1
                # Odd rows: empty
            else:  # Odd column
                # Set 2 at odd rows
                if row % 2 == 1:
                    if idx_set2 >= len(barcodes_set2):
                        raise ValueError(f"Ran out of Set 2 barcodes at position ({row},{col})")
                    
                    layout.append({
                        'row': row,
                        'col': col,
                        'barcode': barcodes_set2[idx_set2],
                        'source': 'set2'
                    })
                    idx_set2 += 1
                    set2_used += 1
                # Even rows: empty
    
    df = pd.DataFrame(layout)
    
    total_spots = target_rows * target_cols
    filled_spots = len(layout)
    
    print(f"Grid creation complete:")
    print(f"  Grid size: {target_rows} × {target_cols} = {total_spots:,} total spots")
    print(f"  Set 1 used: {set1_used:,} / {len(barcodes_set1):,} barcodes")
    print(f"  Set 2 used: {set2_used:,} / {len(barcodes_set2):,} barcodes")
    print(f"  Total filled: {filled_spots:,} ({filled_spots/total_spots*100:.1f}%)")
    print(f"  Empty spots: {total_spots - filled_spots:,} ({(total_spots-filled_spots)/total_spots*100:.1f}%)")
    
    # Check if we used all barcodes
    if idx_set1 < len(barcodes_set1):
        print(f"  ⚠ Warning: {len(barcodes_set1) - idx_set1} Set 1 barcodes unused")
    if idx_set2 < len(barcodes_set2):
        print(f"  ⚠ Warning: {len(barcodes_set2) - idx_set2} Set 2 barcodes unused")
    
    return df


def calculate_optimal_dimensions(n_set1, n_set2):
    """
    Calculate optimal grid dimensions to fit barcodes in checkerboard.
    
    In checkerboard:
    - Even columns have n_rows/2 spots (even rows)
    - Odd columns have n_rows/2 spots (odd rows)
    - Total filled = (n_cols/2 * n_rows/2) + (n_cols/2 * n_rows/2) = n_rows * n_cols / 2
    
    We need: n_rows * n_cols / 2 >= n_set1 + n_set2
    """
    total_barcodes = n_set1 + n_set2
    total_spots_needed = total_barcodes * 2  # Double for checkerboard
    
    # Try to make it roughly square
    side = int(np.sqrt(total_spots_needed))
    
    # Try dimensions close to square
    candidates = []
    for rows in range(side - 50, side + 50):
        for cols in range(side - 50, side + 50):
            filled = (rows // 2) * cols  # Even columns with even rows + odd columns with odd rows
            if filled >= total_barcodes:
                # Check exact fit per set
                even_cols = (cols + 1) // 2
                odd_cols = cols // 2
                set1_capacity = even_cols * (rows // 2 + rows % 2)  # Even rows in even cols
                set2_capacity = odd_cols * ((rows + 1) // 2)  # Odd rows in odd cols
                
                if set1_capacity >= n_set1 and set2_capacity >= n_set2:
                    candidates.append((rows, cols, filled, abs(rows - cols)))
    
    if not candidates:
        raise ValueError("Could not find suitable dimensions")
    
    # Choose most square option
    candidates.sort(key=lambda x: x[3])
    rows, cols, filled, _ = candidates[0]
    
    return rows, cols


def export_layouts(df, output_base, total_rows, total_cols):
    """Export combined layout in multiple formats"""
    
    # 1. Full layout with source (column-major)
    full_path = f"{output_base}_full_colmajor.csv"
    df_sorted = df.sort_values(['col', 'row'])
    df_sorted.to_csv(full_path, index=False)
    print(f"  ✓ Full layout: {full_path}")
    
    # 2. Simple layout without source (column-major)
    simple_path = f"{output_base}_colmajor.csv"
    df_simple = df[['row', 'col', 'barcode']].copy()
    df_simple = df_simple.sort_values(['col', 'row'])
    df_simple.to_csv(simple_path, index=False)
    print(f"  ✓ Simple layout: {simple_path}")
    
    # 3. Barcodes only (column-major)
    barcodes_path = f"{output_base}_barcodes.txt"
    with open(barcodes_path, 'w') as f:
        for barcode in df_simple['barcode']:
            f.write(f"{barcode}\n")
    print(f"  ✓ Barcodes: {barcodes_path}")
    
    # 4. Pattern visualization
    pattern_path = f"{output_base}_pattern.txt"
    with open(pattern_path, 'w') as f:
        f.write("CHECKERBOARD PATTERN\n")
        f.write("=" * 82 + "\n")
        f.write("Legend: 1=Set1, 2=Set2, .=Empty\n")
        f.write("=" * 82 + "\n\n")
        
        # Create lookup
        pos_to_source = {}
        for _, row in df.iterrows():
            pos_to_source[(row['row'], row['col'])] = row['source']
        
        vis_rows = min(40, total_rows)
        vis_cols = min(80, total_cols)
        
        # Header
        f.write("    ")
        for c in range(vis_cols):
            f.write(f"{c%10}")
        f.write("\n    " + "-" * vis_cols + "\n")
        
        # Grid
        for r in range(vis_rows):
            f.write(f"{r:3d}|")
            for c in range(vis_cols):
                source = pos_to_source.get((r, c))
                if source == 'set1':
                    f.write('1')
                elif source == 'set2':
                    f.write('2')
                else:
                    f.write('.')
            f.write("\n")
        
        if total_rows > vis_rows or total_cols > vis_cols:
            f.write(f"\n(Showing first {vis_rows} rows × {vis_cols} cols of {total_rows} × {total_cols})\n")
        
        f.write(f"\nPattern rules:\n")
        f.write(f"  • Even columns (0,2,4...): Set 1 at even rows (0,2,4...)\n")
        f.write(f"  • Odd columns (1,3,5...): Set 2 at odd rows (1,3,5...)\n")
        f.write(f"  • All other positions: Empty\n")
    
    print(f"  ✓ Pattern: {pattern_path}")
    
    # 5. Statistics
    stats_path = f"{output_base}_stats.txt"
    with open(stats_path, 'w') as f:
        set1_count = len(df[df['source'] == 'set1'])
        set2_count = len(df[df['source'] == 'set2'])
        
        f.write("CHECKERBOARD LAYOUT STATISTICS\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Grid dimensions: {total_rows} × {total_cols} = {total_rows * total_cols:,} spots\n\n")
        f.write(f"Set 1 barcodes: {set1_count:,}\n")
        f.write(f"Set 2 barcodes: {set2_count:,}\n")
        f.write(f"Total filled: {len(df):,} ({len(df)/(total_rows*total_cols)*100:.2f}%)\n")
        f.write(f"Empty spots: {total_rows * total_cols - len(df):,}\n\n")
        
        f.write("Pattern verification:\n")
        errors = []
        for _, row in df.iterrows():
            r, c, source = row['row'], row['col'], row['source']
            if c % 2 == 0 and r % 2 == 0 and source != 'set1':
                errors.append(f"({r},{c}): should be set1, got {source}")
            elif c % 2 == 1 and r % 2 == 1 and source != 'set2':
                errors.append(f"({r},{c}): should be set2, got {source}")
            elif (c % 2 == 0 and r % 2 == 1) or (c % 2 == 1 and r % 2 == 0):
                errors.append(f"({r},{c}): should be empty, got {source}")
        
        if errors:
            f.write(f"  ✗ Found {len(errors)} errors\n")
            for err in errors[:20]:
                f.write(f"    {err}\n")
        else:
            f.write(f"  ✓ Perfect checkerboard pattern\n")
    
    print(f"  ✓ Statistics: {stats_path}")


def main():
    """Main execution"""
    
    print("=" * 80)
    print("FULL CHECKERBOARD LAYOUT COMBINER")
    print("Uses all barcodes from both sets")
    print("=" * 80)
    print()
    
    # Parse arguments
    if len(sys.argv) >= 3:
        set1_file = sys.argv[1]
        set2_file = sys.argv[2]
        
        if len(sys.argv) >= 5:
            target_rows = int(sys.argv[3])
            target_cols = int(sys.argv[4])
            auto_size = False
        else:
            auto_size = True
    else:
        # Defaults
        set1_file = "output_118x182_nosep_set1/layout_118x182_nosep_set1_colmajor.csv"
        set2_file = "output_118x182_nosep_set2/layout_118x182_nosep_set2_colmajor.csv"
        auto_size = True
        
        print(f"Using defaults:")
        print(f"  Set 1: {set1_file}")
        print(f"  Set 2: {set2_file}")
        print()
    
    # Load barcodes
    print("[1/4] Loading barcode sets...")
    try:
        barcodes_set1 = load_layout_ordered(set1_file)
        print(f"  ✓ Set 1: {len(barcodes_set1):,} barcodes")
        
        barcodes_set2 = load_layout_ordered(set2_file)
        print(f"  ✓ Set 2: {len(barcodes_set2):,} barcodes")
        
        total_barcodes = len(barcodes_set1) + len(barcodes_set2)
        print(f"  Total: {total_barcodes:,} barcodes")
    except Exception as e:
        print(f"  ✗ Error: {e}")
        sys.exit(1)
    
    # Calculate dimensions
    print(f"\n[2/4] Calculating grid dimensions...")
    if auto_size:
        target_rows, target_cols = calculate_optimal_dimensions(len(barcodes_set1), len(barcodes_set2))
        print(f"  Auto-calculated: {target_rows} × {target_cols}")
    else:
        print(f"  User-specified: {target_rows} × {target_cols}")
    
    total_spots = target_rows * target_cols
    print(f"  Total spots: {total_spots:,} (need {total_barcodes:,} filled)")
    
    # Create checkerboard
    print(f"\n[3/4] Creating checkerboard layout...")
    try:
        df = create_checkerboard_grid(barcodes_set1, barcodes_set2, target_rows, target_cols)
    except Exception as e:
        print(f"  ✗ Error: {e}")
        sys.exit(1)
    
    # Export
    print(f"\n[4/4] Exporting files...")
    output_base = f"layout_{target_rows}x{target_cols}_checkerboard"
    export_layouts(df, output_base, target_rows, target_cols)
    
    # Summary
    print("\n" + "=" * 80)
    print("CHECKERBOARD LAYOUT COMPLETE")
    print("=" * 80)
    print(f"\nGrid: {target_rows} × {target_cols} = {total_spots:,} spots")
    print(f"Filled: {len(df):,} spots ({len(df)/total_spots*100:.1f}%)")
    print(f"Pattern: Even cols→Set1 at even rows, Odd cols→Set2 at odd rows")
    print(f"\nOutput files: {output_base}_*")
    print("=" * 80)


if __name__ == "__main__":
    main()
