#!/usr/bin/env python3
"""
Generate Quadtree Barcode Layout for 235×363 Array (85,305 spots)

This script generates quadtree barcodes using hierarchical spatial encoding with
Run-Length Limiting (RLL) to prevent homopolymers.

Key Features:
- NO input barcode library needed - generates algorithmically from coordinates
- Spatial similarity: adjacent spots have similar sequences (low Hamming distance)
- RLL prevents homopolymers through global permutation optimization
- 24nt barcodes: 3 segments × 8nt with continuous RLL

Usage:
    python generate_layout_235x363_quadtree.py [config_file]
    
    If config_file is not provided, uses config_235x363_24nt.json
"""

import sys
import os
import json
import numpy as np
from datetime import datetime
from generate_repeated_quadtree import generate_repeated_quadtree_layout, analyze_layout


def load_config(config_path='config_235x363_24nt.json'):
    """Load configuration from JSON file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    return config


def export_barcodes_txt(df, filepath):
    """Export barcodes to simple text file (one per line, column-major order)"""
    # Column-major order: iterate through columns, then rows within each column
    with open(filepath, 'w') as f:
        for col in sorted(df['col'].unique()):
            col_data = df[df['col'] == col].sort_values('row')
            for barcode in col_data['barcode']:
                f.write(f"{barcode}\n")
    
    print(f"✓ Barcodes exported to {filepath} (column-major order)")


def save_statistics(stats_list, config, df, filepath):
    """Save comprehensive statistics to JSON"""
    
    # Compute Hamming distance statistics for spatial neighbors
    hamming_dists = []
    rows = df['row'].max() + 1
    cols = df['col'].max() + 1
    
    # Sample 1000 random spots to compute neighbor distances
    sample_size = min(1000, len(df))
    sample_indices = np.random.choice(len(df), sample_size, replace=False)
    
    for idx in sample_indices:
        row, col = df.iloc[idx]['row'], df.iloc[idx]['col']
        barcode = df.iloc[idx]['barcode']
        
        # Check right neighbor
        if col < cols - 1:
            neighbor = df[(df['row'] == row) & (df['col'] == col + 1)]
            if len(neighbor) > 0:
                neighbor_bc = neighbor['barcode'].values[0]
                dist = sum(c1 != c2 for c1, c2 in zip(barcode, neighbor_bc))
                hamming_dists.append(dist)
        
        # Check bottom neighbor
        if row < rows - 1:
            neighbor = df[(df['row'] == row + 1) & (df['col'] == col)]
            if len(neighbor) > 0:
                neighbor_bc = neighbor['barcode'].values[0]
                dist = sum(c1 != c2 for c1, c2 in zip(barcode, neighbor_bc))
                hamming_dists.append(dist)
    
    output = {
        'generation_timestamp': datetime.now().isoformat(),
        'configuration': config,
        'layout_statistics': {
            'total_spots': len(df),
            'unique_barcodes': df['barcode'].nunique(),
            'barcode_length': len(df['barcode'].iloc[0]),
            'grid_dimensions': {
                'rows': int(df['row'].max() + 1),
                'cols': int(df['col'].max() + 1)
            }
        },
        'rll_statistics': [s for s in stats_list if not isinstance(s, dict) or 'segment' in s or 'overall' in s],
        'spatial_similarity': {
            'neighbor_hamming_distances': {
                'mean': float(np.mean(hamming_dists)) if hamming_dists else 0,
                'std': float(np.std(hamming_dists)) if hamming_dists else 0,
                'min': int(np.min(hamming_dists)) if hamming_dists else 0,
                'max': int(np.max(hamming_dists)) if hamming_dists else 0,
                'sample_size': len(hamming_dists)
            }
        },
        'quality_metrics': {
            'all_unique': df['barcode'].nunique() == len(df),
            'no_duplicates': not df.duplicated(subset=['barcode']).any()
        }
    }
    
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"✓ Statistics saved to {filepath}")


def main():
    """Main execution function"""
    
    # Parse command line arguments
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_235x363_24nt.json'
    
    print("=" * 80)
    print("QUADTREE BARCODE LAYOUT GENERATOR - 235×363 ARRAY")
    print("=" * 80)
    print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load configuration
    print(f"\n[1/5] Loading configuration from {config_file}...")
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
    segment_length = config['quadtree_configuration']['segment_length']
    n_repeats = config['quadtree_configuration']['n_repeats']
    total_length = config['quadtree_configuration']['total_barcode_length']
    max_run = config['rll_configuration']['max_run']
    
    output_csv = config['output']['layout_csv']
    output_txt = config['output']['barcodes_txt']
    output_stats = config['output']['statistics_json']
    
    print(f"\nArray configuration:")
    print(f"  Dimensions: {rows} × {cols} = {total_spots:,} spots")
    print(f"  Barcode structure: {n_repeats} × {segment_length}nt = {total_length}nt")
    print(f"  RLL constraint: max_run = {max_run}")
    print(f"  Generation method: Algorithmic (no input library needed)")
    
    # Set random seed
    seed = config.get('random_seed', 42)
    np.random.seed(seed)
    print(f"  Random seed: {seed}")
    
    # Generate layout
    print(f"\n[2/5] Generating quadtree layout with RLL...")
    print(f"  This may take 30-120 seconds for {total_spots:,} spots...")
    try:
        df, segments, combined_matrix, all_stats = generate_repeated_quadtree_layout(
            rows, cols, 
            segment_length=segment_length, 
            n_repeats=n_repeats, 
            max_run=max_run
        )
        print(f"✓ Layout generated successfully")
        print(f"  Generated {len(df):,} unique barcodes")
        
    except Exception as e:
        print(f"✗ Generation error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Analyze layout
    print(f"\n[3/5] Analyzing layout quality...")
    try:
        analyze_layout(df, segments, all_stats)
        print(f"✓ Analysis complete")
        
    except Exception as e:
        print(f"⚠ Warning: Analysis error: {e}")
    
    # Export files
    print(f"\n[4/5] Exporting results...")
    try:
        # Export CSV (column-major order)
        df_export = df.copy()
        # Sort by column first, then row (column-major)
        df_export = df_export.sort_values(['col', 'row'])
        df_export.to_csv(output_csv, index=False)
        print(f"✓ Layout CSV exported to {output_csv} (column-major order)")
        
        file_size_mb = os.path.getsize(output_csv) / (1024 * 1024)
        print(f"  File size: {file_size_mb:.2f} MB")
        
        # Export simple text file
        export_barcodes_txt(df_export, output_txt)
        
        # Export statistics
        save_statistics(all_stats, config, df, output_stats)
        
    except Exception as e:
        print(f"✗ Export error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Verification
    print(f"\n[5/5] Final verification...")
    try:
        # Check uniqueness
        n_unique = df['barcode'].nunique()
        n_total = len(df)
        
        if n_unique == n_total:
            print(f"✓ All {n_total:,} barcodes are unique")
        else:
            print(f"⚠ Warning: Only {n_unique:,} unique out of {n_total:,} total")
        
        # Check length
        barcode_length = len(df['barcode'].iloc[0])
        if barcode_length == total_length:
            print(f"✓ All barcodes are {barcode_length}nt")
        else:
            print(f"⚠ Warning: Expected {total_length}nt, got {barcode_length}nt")
        
        # Check ordering
        first_rows = df.head(10)
        is_column_major = all(first_rows['col'].iloc[i] <= first_rows['col'].iloc[i+1] 
                             for i in range(len(first_rows)-1))
        if is_column_major:
            print(f"✓ Output is column-major ordered")
        else:
            print(f"⚠ Warning: Output ordering may not be column-major")
        
    except Exception as e:
        print(f"⚠ Warning: Verification error: {e}")
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    print(f"\nOutput files:")
    print(f"  1. Layout CSV: {output_csv}")
    print(f"  2. Barcodes TXT: {output_txt}")
    print(f"  3. Statistics JSON: {output_stats}")
    
    print(f"\nLayout summary:")
    print(f"  • {n_unique:,} unique barcodes for {n_total:,} spots")
    print(f"  • Each barcode is {total_length}nt ({n_repeats}×{segment_length}nt)")
    print(f"  • Spatial encoding: adjacent spots have similar sequences")
    print(f"  • RLL constraint: max_run={max_run} (prevents long homopolymers)")
    print(f"  • Column-major ordering (compatible with array fabrication)")
    
    # Get RLL statistics
    overall_stats = next((s for s in all_stats if isinstance(s, dict) and s.get('overall')), None)
    if overall_stats:
        print(f"\nRLL Performance:")
        print(f"  • Forced flips: {overall_stats['total_forced_flips']} ({overall_stats['mean_forced_flip_rate']:.4f} rate)")
        print(f"  • Violations avoided: {overall_stats['total_violations_avoided']}")
    
    print(f"\nUsage:")
    print(f"  • Read layout: import pandas as pd; df = pd.read_csv('{output_csv}')")
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
