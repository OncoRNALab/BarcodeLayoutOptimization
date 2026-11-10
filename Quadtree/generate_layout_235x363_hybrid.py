#!/usr/bin/env python3
"""
Generate Hybrid Quadtree+Library Layout for 235×363 Array (85,305 spots, 36nt barcodes)

This script generates 36nt hybrid barcodes combining:
- 24nt quadtree (3×8nt with RLL, restarted from different permutations) - provides spatial similarity
- 12nt random tail - provides high Hamming distance

Usage:
    python generate_layout_235x363_hybrid.py [config_file]
"""

import sys
import os
import json
import numpy as np
from datetime import datetime
from generate_hybrid_quadtree_library import generate_hybrid_layout, analyze_hybrid_layout


def load_config(config_path='config_235x363_hybrid_36nt.json'):
    """Load configuration from JSON file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    return config


def export_barcodes_txt(df, filepath):
    """Export barcodes to simple text file (column-major order)"""
    with open(filepath, 'w') as f:
        for col in sorted(df['col'].unique()):
            col_data = df[df['col'] == col].sort_values('row')
            for barcode in col_data['barcode']:
                f.write(f"{barcode}\n")
    
    print(f"✓ Barcodes exported to {filepath} (column-major order)")


def save_statistics(all_stats, config, df, filepath):
    """Save comprehensive statistics to JSON"""
    
    # Compute Hamming distance statistics
    hamming_dists = []
    rows = df['row'].max() + 1
    cols = df['col'].max() + 1
    
    sample_size = min(1000, len(df))
    sample_indices = np.random.choice(len(df), sample_size, replace=False)
    
    for idx in sample_indices:
        row, col = df.iloc[idx]['row'], df.iloc[idx]['col']
        quadtree_bc = df.iloc[idx]['quadtree_combined']
        
        if col < cols - 1:
            neighbor = df[(df['row'] == row) & (df['col'] == col + 1)]
            if len(neighbor) > 0:
                neighbor_bc = neighbor['quadtree_combined'].values[0]
                dist = sum(c1 != c2 for c1, c2 in zip(quadtree_bc, neighbor_bc))
                hamming_dists.append(dist)
        
        if row < rows - 1:
            neighbor = df[(df['row'] == row + 1) & (df['col'] == col)]
            if len(neighbor) > 0:
                neighbor_bc = neighbor['quadtree_combined'].values[0]
                dist = sum(c1 != c2 for c1, c2 in zip(quadtree_bc, neighbor_bc))
                hamming_dists.append(dist)
    
    output = {
        'generation_timestamp': datetime.now().isoformat(),
        'configuration': config,
        'layout_statistics': {
            'total_spots': len(df),
            'unique_combined_barcodes': df['barcode'].nunique(),
            'total_barcode_length': len(df['barcode'].iloc[0]),
            'quadtree_length': len(df['quadtree_combined'].iloc[0]),
            'tail_length': len(df['library'].iloc[0]),
            'grid_dimensions': {
                'rows': int(df['row'].max() + 1),
                'cols': int(df['col'].max() + 1)
            }
        },
        'quadtree_statistics': {
            'unique_quadtree_codes': df['quadtree_combined'].nunique(),
            'rll_stats': all_stats
        },
        'tail_statistics': {
            'unique_tail_codes': df['library'].nunique(),
            'average_reuse': len(df) / df['library'].nunique()
        },
        'spatial_similarity': {
            'neighbor_hamming_distances_quadtree_only': {
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
    
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_235x363_hybrid_36nt.json'
    
    print("=" * 80)
    print("HYBRID QUADTREE+RANDOM TAIL LAYOUT GENERATOR - 235×363 ARRAY (36nt)")
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
    quadtree_length = config['quadtree_configuration']['quadtree_length']
    
    tail_length = config['library_configuration']['library_length']
    tail_mode = config['library_configuration'].get('mode', 'random')
    library_file = config['library_configuration'].get('library_file')
    assignment_strategy = config['library_configuration'].get('assignment_strategy', 'sequential')
    
    max_run = config['rll_configuration']['max_run']
    total_length = config['total_barcode_length']
    
    output_csv = config['output']['layout_csv']
    output_txt = config['output']['barcodes_txt']
    output_stats = config['output']['statistics_json']
    
    print(f"\nArray configuration:")
    print(f"  Dimensions: {rows} × {cols} = {total_spots:,} spots")
    print(f"  Barcode structure: {quadtree_length}nt quadtree + {tail_length}nt random tail = {total_length}nt")
    print(f"  Quadtree: {n_repeats} × {segment_length}nt with RLL (max_run={max_run})")
    print(f"  Restart mode: {tail_mode}")
    if tail_mode == 'library' and library_file:
        print(f"  Library: {library_file}")
        print(f"  Assignment: {assignment_strategy}")
    
    # Set random seed
    seed = config.get('random_seed', 42)
    np.random.seed(seed)
    print(f"  Random seed: {seed}")
    
    # Generate layout
    print(f"\n[2/5] Generating hybrid layout...")
    print(f"  This may take 1-2 minutes for {total_spots:,} spots...")
    try:
        df, quadtree_segments, all_stats = generate_hybrid_layout(
            rows, cols,
            library_file=library_file,
            segment_length=segment_length,
            n_quadtree_rounds=n_repeats,
            tail_length=tail_length,
            tail_mode=tail_mode,
            random_seed=seed,
            assignment_strategy=assignment_strategy,
            max_run=max_run
        )
        print(f"✓ Layout generated successfully")
        print(f"  Generated {len(df):,} barcodes")
        
    except Exception as e:
        print(f"✗ Generation error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Analyze layout
    print(f"\n[3/5] Analyzing layout quality...")
    try:
        analyze_hybrid_layout(df, n_quadtree_rounds=n_repeats, all_stats=all_stats)
        print(f"✓ Analysis complete")
        
    except Exception as e:
        print(f"⚠ Warning: Analysis error: {e}")
    
    # Export files
    print(f"\n[4/5] Exporting results...")
    try:
        # Export CSV (column-major order)
        df_export = df.copy()
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
        n_unique = df['barcode'].nunique()
        n_total = len(df)
        
        if n_unique == n_total:
            print(f"✓ All {n_total:,} combined barcodes are unique")
        else:
            print(f"⚠ Warning: Only {n_unique:,} unique out of {n_total:,} total")
        
        barcode_length = len(df['barcode'].iloc[0])
        if barcode_length == total_length:
            print(f"✓ All barcodes are {barcode_length}nt ({quadtree_length}nt + {tail_length}nt)")
        else:
            print(f"⚠ Warning: Expected {total_length}nt, got {barcode_length}nt")
        
        tail_unique = df['library'].nunique()
        print(f"✓ Using {tail_unique:,} unique tail barcodes")
        
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
    print(f"  • {n_unique:,} unique 36nt barcodes for {n_total:,} spots")
    print(f"  • Structure: {quadtree_length}nt quadtree + {tail_length}nt random tail")
    print(f"  • Quadtree: spatial encoding with RLL (max_run={max_run})")
    print(f"  • Tail: {tail_unique:,} unique sequences for sequencing quality")
    print(f"  • Column-major ordering (array fabrication compatible)")
    
    overall_stats = all_stats[-1] if all_stats else None
    if overall_stats:
        print(f"\nRLL Performance (Quadtree portion):")
        print(f"  • Forced flips: {overall_stats['total_forced_flips']} ({overall_stats['mean_forced_flip_rate']:.4f} rate)")
        print(f"  • Violations avoided: {overall_stats['total_violations_avoided']}")
    
    print(f"\nUsage:")
    print(f"  • Read layout: import pandas as pd; df = pd.read_csv('{output_csv}')")
    print(f"  • Quadtree portion: df['quadtree_combined']")
    print(f"  • Tail portion: df['library']")
    print(f"  • Full barcode: df['barcode']")
    
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
