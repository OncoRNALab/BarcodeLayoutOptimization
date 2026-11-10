#!/usr/bin/env python3
"""
Generate Checkerboard Dual Barcode Layout for 235×363 Array

This script generates a dual barcode layout for a checkerboard pattern array where:
- Two interleaved arrays: 118×182 (even rows) and 117×181 (odd rows)
- 8×8 checkerboard blocks (approximately 32 spots per block)
- Each block shares one bar_1 barcode (18nt)
- Each position within a block has a unique bar_2 barcode (18nt)
- Final barcode = bar_1 + bar_2 = 36nt
- All bar_1+bar_2 combinations are globally unique

Usage:
    python generate_checkerboard_235x363.py
"""

import sys
import os
import json
from datetime import datetime
from dual_checkerboard_layout import CheckerboardDualBarcodeLayout, load_barcodes_from_file


def save_statistics(stats, filepath):
    """Save layout statistics to JSON file"""
    output = {
        'generation_timestamp': datetime.now().isoformat(),
        'array_configuration': {
            'pattern': 'checkerboard',
            'total_rows': stats['grid_dimensions'][0],
            'total_cols': stats['grid_dimensions'][1],
            'total_spots': stats['total_spots'],
            'array1_spots': stats['array1_spots'],
            'array2_spots': stats['array2_spots'],
            'block_pattern': stats['block_size']
        },
        'barcode_configuration': {
            'barcodes_set1_available': stats['n_barcodes_set1'],
            'barcodes_set2_available': stats['n_barcodes_set2'],
            'barcode_segment_length': '18nt',
            'full_barcode_length': '36nt'
        },
        'layout_statistics': {
            'unique_combined_barcodes': stats['unique_combined_barcodes'],
            'unique_bar1_used': stats['unique_bar1_used'],
            'unique_bar2_used': stats['unique_bar2_used']
        },
        'usage_statistics': {
            'bar1_usage_range': [
                min(stats['bar1_usage'].values()) if stats['bar1_usage'] else 0,
                max(stats['bar1_usage'].values()) if stats['bar1_usage'] else 0
            ],
            'bar2_usage_range': [
                min(stats['bar2_usage'].values()) if stats['bar2_usage'] else 0,
                max(stats['bar2_usage'].values()) if stats['bar2_usage'] else 0
            ]
        }
    }
    
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Statistics saved to {filepath}")


def main():
    """Main execution function"""
    
    # Configuration
    ROWS = 235
    COLS = 363
    BLOCK_SIZE = 8
    SEED = 42
    
    # Barcode files
    BARCODE_FILE_SET1 = 'barcodes_118x182_dual_set1.txt'
    BARCODE_FILE_SET2 = 'barcodes_118x182_dual_set2.txt'
    
    # Output files
    OUTPUT_LAYOUT = 'layout_235x363_checkerboard_colmajor.csv'
    OUTPUT_STATS = 'layout_235x363_checkerboard_stats.json'
    
    print("=" * 80)
    print("CHECKERBOARD DUAL BARCODE LAYOUT GENERATOR - 235×363 ARRAY")
    print("=" * 80)
    print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print(f"\nArray configuration:")
    print(f"  Grid dimensions: {ROWS} × {COLS}")
    print(f"  Pattern: Checkerboard (two interleaved arrays)")
    print(f"  Block size: {BLOCK_SIZE}×{BLOCK_SIZE} checkerboard pattern")
    print(f"  Expected spots: ~{ROWS * COLS // 2}")
    
    # Load barcodes
    print(f"\n[1/5] Loading barcodes...")
    try:
        if not os.path.exists(BARCODE_FILE_SET1):
            raise FileNotFoundError(f"Barcode file not found: {BARCODE_FILE_SET1}")
        if not os.path.exists(BARCODE_FILE_SET2):
            raise FileNotFoundError(f"Barcode file not found: {BARCODE_FILE_SET2}")
        
        barcodes_set1 = load_barcodes_from_file(BARCODE_FILE_SET1)
        barcodes_set2 = load_barcodes_from_file(BARCODE_FILE_SET2)
        
        print(f"  Set 1 (Array 1): {len(barcodes_set1)} barcodes")
        print(f"    First 3: {barcodes_set1[:3]}")
        print(f"    Length: {len(barcodes_set1[0])}nt")
        
        print(f"  Set 2 (Array 2): {len(barcodes_set2)} barcodes")
        print(f"    First 3: {barcodes_set2[:3]}")
        print(f"    Length: {len(barcodes_set2[0])}nt")
        
    except Exception as e:
        print(f"✗ Error loading barcodes: {e}")
        sys.exit(1)
    
    # Validate barcode lengths
    print(f"\n[2/5] Validating barcodes...")
    actual_len_set1 = len(barcodes_set1[0])
    actual_len_set2 = len(barcodes_set2[0])
    
    if actual_len_set1 == 36 and actual_len_set2 == 36:
        print(f"  Barcodes are 36nt - will be split into 18nt halves")
        split_mode = True
    elif actual_len_set1 == 18 and actual_len_set2 == 18:
        print(f"  Barcodes are already 18nt - will use directly")
        split_mode = False
    else:
        print(f"⚠ Warning: Unexpected barcode lengths: {actual_len_set1}nt and {actual_len_set2}nt")
        split_mode = True  # Default to split mode
    
    # Check uniqueness
    unique_set1 = len(set(barcodes_set1))
    unique_set2 = len(set(barcodes_set2))
    
    if unique_set1 < len(barcodes_set1):
        print(f"⚠ Warning: Set 1 has {len(barcodes_set1) - unique_set1} duplicate barcodes")
    if unique_set2 < len(barcodes_set2):
        print(f"⚠ Warning: Set 2 has {len(barcodes_set2) - unique_set2} duplicate barcodes")
    
    print(f"  ✓ Barcodes validated")
    
    # Create layout generator
    print(f"\n[3/5] Initializing checkerboard layout generator...")
    try:
        layout = CheckerboardDualBarcodeLayout(
            rows=ROWS,
            cols=COLS,
            barcodes_set1=barcodes_set1,
            barcodes_set2=barcodes_set2,
            block_size=BLOCK_SIZE,
            split_barcodes=split_mode
        )
        print(f"  ✓ Generator initialized")
        
    except Exception as e:
        print(f"✗ Initialization error: {e}")
        sys.exit(1)
    
    # Generate layout
    print(f"\n[4/5] Generating dual barcode layout...")
    try:
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=SEED)
        print(f"\n✓ Layout generated successfully")
        
    except Exception as e:
        print(f"✗ Generation error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Get and display statistics
    print(f"\n[5/5] Computing statistics and exporting...")
    try:
        stats = layout.get_statistics()
        
        print(f"\nLayout Statistics:")
        print(f"  Grid: {stats['grid_dimensions'][0]}×{stats['grid_dimensions'][1]}")
        print(f"  Total spots (checkerboard): {stats['total_spots']:,}")
        print(f"    Array 1 (even rows): {stats['array1_spots']:,}")
        print(f"    Array 2 (odd rows): {stats['array2_spots']:,}")
        print(f"  Unique combined barcodes: {stats['unique_combined_barcodes']:,}")
        print(f"  Unique bar_1 used: {stats['unique_bar1_used']}")
        print(f"  Unique bar_2 used: {stats['unique_bar2_used']}")
        
        # Save statistics
        save_statistics(stats, OUTPUT_STATS)
        
        # Export layout in column-major order
        print(f"\n  Exporting layout...")
        layout.export_layout(OUTPUT_LAYOUT, format='csv', column_major=True)
        
        # Check file size
        file_size_mb = os.path.getsize(OUTPUT_LAYOUT) / (1024 * 1024)
        print(f"  File size: {file_size_mb:.2f} MB")
        
    except Exception as e:
        print(f"✗ Statistics/export error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    print(f"\nOutput files:")
    print(f"  1. Layout CSV: {OUTPUT_LAYOUT}")
    print(f"  2. Statistics JSON: {OUTPUT_STATS}")
    
    print(f"\nLayout summary:")
    print(f"  • {stats['unique_combined_barcodes']:,} unique barcodes for {stats['total_spots']:,} spots")
    print(f"  • Checkerboard pattern with two interleaved arrays")
    print(f"  • Each barcode is 36nt (18nt bar_1 + 18nt bar_2)")
    print(f"  • Organized in 8×8 checkerboard blocks with shared bar_1 per block")
    print(f"  • All bar_1+bar_2 combinations are globally unique")
    
    print(f"\nUsage:")
    print(f"  • Read layout: import pandas as pd; df = pd.read_csv('{OUTPUT_LAYOUT}')")
    print(f"  • Visualize: python plot_cycle_activation.py {OUTPUT_LAYOUT} 0")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n✗ Generation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
