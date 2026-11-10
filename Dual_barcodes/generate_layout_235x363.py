#!/usr/bin/env python3
"""
Generate Dual Barcode Layout for 235×363 Array (85,305 spots)

This script generates a dual barcode layout using 4×4 blocks where:
- Each block shares one bar_1 barcode (18nt)
- Each position within a block has a unique bar_2 barcode (18nt)
- Final barcode is bar_1 + bar_2 = 36nt
- All bar_1+bar_2 combinations are globally unique

Usage:
    python generate_layout_235x363.py [config_file]
    
    If config_file is not provided, uses config_235x363.json
"""

import sys
import os
import json
import numpy as np
from datetime import datetime
from dual_barcode_layout import DualBarcodeLayout, load_barcodes_from_file


def load_config(config_path='config_235x363.json'):
    """Load configuration from JSON file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    return config


def validate_barcodes(barcodes, expected_length, min_count):
    """Validate barcode properties"""
    if len(barcodes) < min_count:
        raise ValueError(
            f"Insufficient barcodes: need at least {min_count}, "
            f"but only {len(barcodes)} available"
        )
    
    # Check barcode lengths
    barcode_lengths = set(len(bc) for bc in barcodes)
    if len(barcode_lengths) > 1:
        print(f"⚠ Warning: Barcodes have varying lengths: {barcode_lengths}")
    
    actual_length = len(barcodes[0])
    if actual_length != expected_length:
        print(f"⚠ Warning: Expected barcode length {expected_length}nt, "
              f"but got {actual_length}nt")
    
    # Check for duplicates
    unique_barcodes = len(set(barcodes))
    if unique_barcodes < len(barcodes):
        print(f"⚠ Warning: {len(barcodes) - unique_barcodes} duplicate barcodes found")
    
    return True


def save_statistics(stats, config, filepath):
    """Save layout statistics to JSON file"""
    output = {
        'generation_timestamp': datetime.now().isoformat(),
        'configuration': config,
        'layout_statistics': {
            'grid_dimensions': stats['grid_dimensions'],
            'total_spots': stats['total_spots'],
            'n_barcodes_available': stats['n_barcodes_available'],
            'n_blocks': stats['n_blocks'],
            'unique_combined_barcodes': stats['unique_combined_barcodes'],
            'theoretical_max_combinations': stats['theoretical_max_combinations']
        },
        'barcode_usage': {
            'bar_1_unique_count': len(stats['bar1_usage']),
            'bar_2_unique_count': len(stats['bar2_usage']),
            'bar_1_usage_distribution': {
                'min': min(stats['bar1_usage'].values()) if stats['bar1_usage'] else 0,
                'max': max(stats['bar1_usage'].values()) if stats['bar1_usage'] else 0,
                'mean': np.mean(list(stats['bar1_usage'].values())) if stats['bar1_usage'] else 0
            },
            'bar_2_usage_distribution': {
                'min': min(stats['bar2_usage'].values()) if stats['bar2_usage'] else 0,
                'max': max(stats['bar2_usage'].values()) if stats['bar2_usage'] else 0,
                'mean': np.mean(list(stats['bar2_usage'].values())) if stats['bar2_usage'] else 0
            }
        }
    }
    
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Statistics saved to {filepath}")


def main():
    """Main execution function"""
    
    # Parse command line arguments
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_235x363.json'
    
    print("=" * 80)
    print("DUAL BARCODE LAYOUT GENERATOR - 235×363 ARRAY")
    print("=" * 80)
    print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load configuration
    print(f"\n[1/6] Loading configuration from {config_file}...")
    try:
        config = load_config(config_file)
        print("✓ Configuration loaded")
    except Exception as e:
        print(f"✗ Error loading configuration: {e}")
        sys.exit(1)
    
    # Extract parameters
    rows = config['array_dimensions']['rows']
    cols = config['array_dimensions']['columns']
    total_spots = config['array_dimensions']['total_spots']
    barcode_file = config['barcode_configuration']['barcode_file']
    barcode_length = config['barcode_configuration']['barcode_length_per_segment']
    output_layout = config['output']['layout_file']
    output_viz = config['output']['visualization_file']
    output_stats = config['output']['statistics_file']
    seed = config.get('random_seed', 42)
    
    print(f"\nArray configuration:")
    print(f"  Dimensions: {rows} × {cols} = {total_spots} spots")
    print(f"  Block size: 4×4")
    print(f"  Total blocks: {config['grid_configuration']['total_blocks']}")
    print(f"  Barcode length: {barcode_length}nt + {barcode_length}nt = {config['barcode_configuration']['full_barcode_length']}nt")
    
    # Load barcodes
    print(f"\n[2/6] Loading barcodes from {barcode_file}...")
    try:
        if not os.path.exists(barcode_file):
            raise FileNotFoundError(f"Barcode file not found: {barcode_file}")
        
        barcodes = load_barcodes_from_file(barcode_file)
        print(f"✓ Loaded {len(barcodes)} barcodes")
        print(f"  First 3: {barcodes[:3]}")
        print(f"  Last 3: {barcodes[-3:]}")
        
    except Exception as e:
        print(f"✗ Error loading barcodes: {e}")
        sys.exit(1)
    
    # Validate barcodes
    print(f"\n[3/6] Validating barcodes...")
    try:
        min_barcodes_needed = config['grid_configuration']['total_blocks']
        validate_barcodes(barcodes, barcode_length, min_barcodes_needed)
        print(f"✓ Barcodes validated")
        print(f"  Available: {len(barcodes)}")
        print(f"  Required (minimum): {min_barcodes_needed}")
        print(f"  Status: {'OK' if len(barcodes) >= min_barcodes_needed else 'INSUFFICIENT'}")
    except Exception as e:
        print(f"✗ Validation error: {e}")
        sys.exit(1)
    
    # Create layout generator
    print(f"\n[4/6] Generating dual barcode layout...")
    try:
        layout = DualBarcodeLayout(rows, cols, barcodes)
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=seed)
        print(f"✓ Layout generated successfully")
        
        # Verify all spots have barcodes
        n_assigned = np.sum(grid_combined != None)
        print(f"  Assigned barcodes: {n_assigned} / {total_spots}")
        
        if n_assigned != total_spots:
            print(f"⚠ Warning: Not all spots received barcodes!")
        
    except Exception as e:
        print(f"✗ Generation error: {e}")
        sys.exit(1)
    
    # Get and display statistics
    print(f"\n[5/6] Computing statistics...")
    try:
        stats = layout.get_statistics()
        print("✓ Statistics computed")
        print(f"\n  Grid: {stats['grid_dimensions'][0]}×{stats['grid_dimensions'][1]}")
        print(f"  Total spots: {stats['total_spots']:,}")
        print(f"  Unique combined barcodes: {stats['unique_combined_barcodes']:,}")
        print(f"  Number of blocks: {stats['n_blocks']:,}")
        print(f"  Barcodes available: {stats['n_barcodes_available']:,}")
        print(f"  Unique bar_1 used: {len(stats['bar1_usage'])}")
        print(f"  Unique bar_2 used: {len(stats['bar2_usage'])}")
        
        # Save statistics
        save_statistics(stats, config, output_stats)
        
    except Exception as e:
        print(f"✗ Statistics error: {e}")
        sys.exit(1)
    
    # Export layout
    print(f"\n[6/6] Exporting results...")
    try:
        # Export CSV layout
        layout.export_layout(output_layout, format='csv')
        print(f"✓ Layout exported to {output_layout}")
        
        # Check file size
        file_size_mb = os.path.getsize(output_layout) / (1024 * 1024)
        print(f"  File size: {file_size_mb:.2f} MB")
        
        # Generate visualization
        try:
            print(f"\n  Generating visualization...")
            layout.visualize_layout(output_viz)
            print(f"✓ Visualization saved to {output_viz}")
        except Exception as viz_error:
            print(f"⚠ Warning: Could not generate visualization: {viz_error}")
            print(f"  (This is optional - layout data is still valid)")
        
    except Exception as e:
        print(f"✗ Export error: {e}")
        sys.exit(1)
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    print(f"\nOutput files:")
    print(f"  1. Layout CSV: {output_layout}")
    print(f"  2. Statistics JSON: {output_stats}")
    print(f"  3. Visualization: {output_viz} (if generated)")
    
    print(f"\nLayout summary:")
    print(f"  • {stats['unique_combined_barcodes']:,} unique barcodes for {stats['total_spots']:,} spots")
    print(f"  • Each barcode is 36nt (18nt bar_1 + 18nt bar_2)")
    print(f"  • Organized in 4×4 blocks with shared bar_1 per block")
    print(f"  • All bar_1+bar_2 combinations are globally unique")
    
    print(f"\nUsage:")
    print(f"  • Read layout: import pandas as pd; df = pd.read_csv('{output_layout}')")
    print(f"  • Access barcode for position (row, col): df[(df['row']==row) & (df['col']==col)]")
    
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
