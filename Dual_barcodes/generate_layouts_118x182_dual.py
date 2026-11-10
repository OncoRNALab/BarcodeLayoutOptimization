#!/usr/bin/env python3
"""
Generate Dual Barcode Layouts for 118×182 Arrays (21,476 spots per set)

This script generates dual barcode layouts for TWO independent sets using 4×4 blocks where:
- Each block shares one bar_1 barcode (18nt)
- Each position within a block has a unique bar_2 barcode (18nt)
- Final barcode is bar_1 + bar_2 = 36nt
- Set 1 and Set 2 use non-overlapping bar_1 libraries to ensure uniqueness

Usage:
    python generate_layouts_118x182_dual.py

Generates both Set 1 and Set 2 layouts using their respective configurations.
"""

import sys
import os
import json
import numpy as np
from datetime import datetime

# Add parent directory to path to import dual_barcode_layout
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dual_barcode_layout import DualBarcodeLayout, load_barcodes_from_file


def load_config(config_path):
    """Load configuration from JSON file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    return config


def generate_set(set_num, config_file):
    """Generate layout for one set"""
    print(f"\n{'='*80}")
    print(f"GENERATING SET {set_num} - DUAL BARCODES (118×182 ARRAY)")
    print(f"{'='*80}\n")
    print(f"Configuration: {config_file}")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Load configuration
    print(f"[1/6] Loading configuration...")
    config = load_config(config_file)
    
    rows = config['array']['rows']
    cols = config['array']['cols']
    bar_1_file = config['barcodes']['bar_1_file']
    bar_2_file = config['barcodes']['bar_2_file']
    
    print(f"  Array: {rows}×{cols} = {rows*cols:,} spots")
    print(f"  Block size: 4×4 = 16 spots/block")
    print(f"  bar_1 library: {bar_1_file}")
    print(f"  bar_2 library: {bar_2_file}")
    
    # Load barcodes
    print(f"\n[2/6] Loading barcode libraries...")
    bar_1 = load_barcodes_from_file(bar_1_file)
    bar_2 = load_barcodes_from_file(bar_2_file)
    
    print(f"  bar_1: {len(bar_1):,} barcodes loaded")
    print(f"  bar_2: {len(bar_2):,} barcodes loaded")
    
    # Calculate requirements
    n_blocks = (rows * cols) // 16  # 4×4 blocks
    
    print(f"\n[3/6] Calculating requirements...")
    print(f"  Total blocks needed: {n_blocks:,}")
    print(f"  bar_1 needed: {n_blocks:,}")
    print(f"  bar_2 needed per block: 16")
    
    if len(bar_1) < n_blocks:
        raise ValueError(f"Insufficient bar_1: need {n_blocks}, have {len(bar_1)}")
    if len(bar_2) < 16:
        raise ValueError(f"Insufficient bar_2: need at least 16, have {len(bar_2)}")
    
    # Generate layout
    print(f"\n[4/6] Generating dual barcode layout...")
    # DualBarcodeLayout uses a single barcode list and generates 4x4 blocks internally
    # We use bar_1 as the barcode source
    layout = DualBarcodeLayout(rows=rows, cols=cols, barcodes=bar_1)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=1)
    
    # Flatten to get column-major ordered list
    barcodes_list = []
    for col in range(cols):
        for row in range(rows):
            barcodes_list.append(grid_combined[row, col])
    
    print(f"  ✓ Layout generated")
    print(f"  Total barcodes: {len(barcodes_list):,}")
    print(f"  Unique barcodes: {len(set(barcodes_list)):,}")
    
    # Verify
    print(f"\n[5/6] Verifying layout...")
    stats = layout.get_statistics()
    
    n_unique = len(set(barcodes_list))
    n_total = len(barcodes_list)
    
    if n_unique == n_total:
        print(f"  ✓ All {n_total:,} barcodes are unique")
    else:
        print(f"  ⚠ Warning: Only {n_unique:,} unique out of {n_total:,} total")
    
    print(f"  ✓ Unique bar_1 used: {len(stats['bar1_usage'])}")
    print(f"  ✓ Unique bar_2 used: {len(stats['bar2_usage'])}")
    
    # Export
    print(f"\n[6/6] Exporting results...")
    layout_file = config['output']['layout_file']
    barcode_file = config['output']['barcode_file']
    stats_file = config['output']['stats_file']
    
    # Export layout CSV with full details (column-major)
    import csv
    with open(layout_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['row', 'col', 'barcode', 'bar_1', 'bar_2'])
        for col in range(cols):
            for row in range(rows):
                combined = grid_combined[row, col]
                b1 = grid_bar1[row, col]
                b2 = grid_bar2[row, col]
                writer.writerow([row, col, combined, b1, b2])
    print(f"  ✓ Layout exported: {layout_file}")
    
    # Export barcode list (column-major order)
    with open(barcode_file, 'w') as f:
        for barcode in barcodes_list:
            f.write(f"{barcode}\n")
    print(f"  ✓ Barcodes exported: {barcode_file}")
    
    # Export statistics
    stats_output = {
        'generation_timestamp': datetime.now().isoformat(),
        'set_number': set_num,
        'configuration': config,
        'results': {
            'total_spots': rows * cols,
            'total_barcodes': len(barcodes_list),
            'unique_barcodes': len(set(barcodes_list)),
            'barcode_length': len(barcodes_list[0]),
            'n_blocks': stats['n_blocks'],
            'unique_bar_1': len(stats['bar1_usage']),
            'unique_bar_2': len(stats['bar2_usage']),
            'all_unique': len(set(barcodes_list)) == len(barcodes_list),
        }
    }
    
    with open(stats_file, 'w') as f:
        json.dump(stats_output, f, indent=2)
    print(f"  ✓ Statistics saved: {stats_file}")
    
    print(f"\n{'='*80}")
    print(f"SET {set_num} GENERATION COMPLETE")
    print(f"{'='*80}\n")
    
    # Return barcodes list and stats for cross-set verification
    return barcodes_list, stats_output


def verify_uniqueness_between_sets(barcodes1, barcodes2):
    """Verify no barcode overlap between two sets"""
    print(f"\n{'='*80}")
    print(f"CROSS-SET VERIFICATION")
    print(f"{'='*80}\n")
    
    set1 = set(barcodes1)
    set2 = set(barcodes2)
    
    overlap = set1 & set2
    
    print(f"Set 1: {len(set1):,} unique barcodes")
    print(f"Set 2: {len(set2):,} unique barcodes")
    print(f"Overlap: {len(overlap):,} barcodes")
    
    if len(overlap) == 0:
        print(f"\n✓ NO OVERLAP - Both sets are completely unique!")
        return True
    else:
        print(f"\n✗ OVERLAP DETECTED - {len(overlap)} barcodes appear in both sets:")
        for i, bc in enumerate(list(overlap)[:10]):
            print(f"  {i+1}. {bc}")
        if len(overlap) > 10:
            print(f"  ... and {len(overlap)-10} more")
        return False


def main():
    """Main execution"""
    # Generate Set 1
    barcodes1, stats1 = generate_set(1, 'config_118x182_set1_decoy.json')
    
    # Generate Set 2
    barcodes2, stats2 = generate_set(2, 'config_118x182_set2_decoy.json')
    
    # Verify uniqueness between sets
    no_overlap = verify_uniqueness_between_sets(barcodes1, barcodes2)
    
    # Final summary
    print(f"\n{'='*80}")
    print(f"FINAL SUMMARY")
    print(f"{'='*80}")
    print(f"Set 1: {len(barcodes1):,} barcodes ({len(set(barcodes1)):,} unique)")
    print(f"Set 2: {len(barcodes2):,} barcodes ({len(set(barcodes2)):,} unique)")
    print(f"Total: {len(barcodes1) + len(barcodes2):,} barcodes")
    print(f"Cross-set uniqueness: {'✓ PASS' if no_overlap else '✗ FAIL'}")
    print(f"{'='*80}\n")
    
    if not no_overlap:
        sys.exit(1)


if __name__ == '__main__':
    main()
