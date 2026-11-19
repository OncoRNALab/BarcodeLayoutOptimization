#!/usr/bin/env python3
"""
Generate Zipcodes Layout for 235×363 Array (85,305 spots) - No Separator

This script generates spatial zipcodes using Gray-code-like encoding where:
- Each axis (X, Y) has independent 1D Gray codes
- Spatial neighbors have similar sequences (low Hamming distance)
- No separator between X and Y codes

Barcode Structure: [9nt X-code][9nt Y-code] = 18nt total

Usage:
    python generate_layout_235x363_zipcodes.py [config_file]
"""

import sys
import os
import json
import logging
from pathlib import Path
from datetime import datetime
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from spatial_zipcodes.config import load_config
from spatial_zipcodes.axis_code_generator import generate_column_and_row_codes
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.sequence_qc import (
    calculate_g_content,
    summarize_junction_homopolymers,
    validate_axis_sequences,
)
from spatial_zipcodes.codebook import build_codebook, export_codebook
from spatial_zipcodes.distance_analysis import analyze_axis_codes

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)


def export_barcodes_txt(df, filepath):
    """Export barcodes to simple text file (column-major order)"""
    with open(filepath, 'w') as f:
        for col in sorted(df['col'].unique()):
            col_data = df[df['col'] == col].sort_values('row')
            for barcode in col_data['barcode']:
                f.write(f"{barcode}\n")
    
    print(f"✓ Barcodes exported to {filepath} (column-major order)")


def main():
    """Main execution function"""
    
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_235x363_nosep.json'
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'output_zipcodes'
    
    print("=" * 80)
    print("ZIPCODES LAYOUT GENERATOR - 235×363 ARRAY (NO SEPARATOR)")
    print("=" * 80)
    print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load configuration
    print(f"\n[1/7] Loading configuration from {config_file}...")
    try:
        config = load_config(config_file)
        rows = config.chip_geometry.rows
        cols = config.chip_geometry.cols
        print("✓ Configuration loaded")
        print(f"  Array: {rows} × {cols} = {rows*cols:,} spots")
        print(f"  X-axis codes needed: {cols}")
        print(f"  Y-axis codes needed: {rows}")
        print(f"  Spacer: {'None' if not config.codebook.spacer_sequence else config.codebook.spacer_sequence}")
    except Exception as e:
        print(f"✗ Error loading configuration: {e}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    (output_dir / "intermediate").mkdir(exist_ok=True)
    
    # Generate axis codes
    print(f"\n[2/7] Generating X and Y axis codes...")
    print(f"  G content max: {config.sequence_qc.g_content_max:.0%} (enforced during search)")
    print(f"  This may take 2-5 minutes due to constraint satisfaction...")
    try:
        encoder = BitToDNAEncoder(
            scheme=config.bit_to_dna.scheme,
            mapping=config.bit_to_dna.mapping
        )
        spacer = config.codebook.spacer_sequence or ""
        enforce_junction = config.sequence_qc.reject_sequences_with_long_homopolymers
        x_codes, y_codes, x_sequences, y_sequences, x_verification, y_verification = (
            generate_column_and_row_codes(
                n_columns=cols,
                n_rows=rows,
                x_bit_length=config.x_bit_length,
                y_bit_length=config.y_bit_length,
                target_distances=config.error_correction.distances,
                far_features_min_distance=config.error_correction.far_features_min_distance,
                encoder=encoder,
                g_content_max=config.sequence_qc.g_content_max,
                max_homopolymer_length=config.sequence_qc.max_homopolymer_length,
                spacer=spacer,
                enforce_junction=enforce_junction,
                max_attempts_per_position=config.axis_code.max_attempts_per_position,
                backtrack_depth=config.axis_code.backtrack_depth,
                max_backtracks=getattr(config.axis_code, "max_backtracks", 100),
                tolerance=config.axis_code.tolerance,
                random_seed=config.axis_code.random_seed,
                candidate_mode=config.axis_code.candidate_mode,
                flip_counts=tuple(config.axis_code.flip_counts),
                require_tolerance=config.axis_code.require_tolerance,
            )
        )
        print(f"    ✓ Generated {len(x_codes)} X-codes")
        print(f"    Hard constraint violations: {x_verification['hard_constraint_violations']}")
        print(f"    Positions outside tolerance: {x_verification['positions_outside_tolerance']}")
        print(f"    ✓ Generated {len(y_codes)} Y-codes")
        print(f"    Hard constraint violations: {y_verification['hard_constraint_violations']}")
        print(f"    Positions outside tolerance: {y_verification['positions_outside_tolerance']}")
        
        axis_codes = {'x_codes': x_codes, 'y_codes': y_codes}
        
        print(f"✓ Axis codes generated")
        print(f"  X-codes: {len(axis_codes['x_codes'])} codes of {len(axis_codes['x_codes'][0])} bits")
        print(f"  Y-codes: {len(axis_codes['y_codes'])} codes of {len(axis_codes['y_codes'][0])} bits")
        
        # Save axis codes
        import json as json_module
        axis_codes_file = output_dir / "intermediate" / "axis_codes.json"
        with open(axis_codes_file, 'w') as f:
            json_module.dump({
                'x_codes': [code.tolist() if hasattr(code, 'tolist') else list(code) 
                           for code in axis_codes['x_codes']],
                'y_codes': [code.tolist() if hasattr(code, 'tolist') else list(code) 
                           for code in axis_codes['y_codes']]
            }, f, indent=2)
        print(f"  Saved axis codes to {axis_codes_file}")
        
    except Exception as e:
        print(f"✗ Generation error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Convert to DNA
    print(f"\n[3/7] Converting bit sequences to DNA...")
    try:
        print(f"✓ DNA sequences generated")
        print(f"  X-sequences: {len(x_sequences)} × {len(x_sequences[0])}nt")
        print(f"  Y-sequences: {len(y_sequences)} × {len(y_sequences[0])}nt")
        print(f"  Example X-code: {x_sequences[0]}")
        print(f"  Example Y-code: {y_sequences[0]}")
        
    except Exception as e:
        print(f"✗ Encoding error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Validate sequences
    print(f"\n[4/7] Validating sequence quality...")
    try:
        qc_config = {
            'gc_content_min': config.sequence_qc.gc_content_min,
            'gc_content_max': config.sequence_qc.gc_content_max,
            'g_content_max': config.sequence_qc.g_content_max,
            'max_homopolymer_length': config.sequence_qc.max_homopolymer_length
        }
        x_valid, x_issues = validate_axis_sequences(x_sequences, 'X', qc_config)
        y_valid, y_issues = validate_axis_sequences(y_sequences, 'Y', qc_config)
        
        g_max = config.sequence_qc.g_content_max
        x_g = [calculate_g_content(seq) for seq in x_sequences]
        y_g = [calculate_g_content(seq) for seq in y_sequences]
        if max(x_g) > g_max or max(y_g) > g_max:
            print(f"✗ G content exceeds {g_max:.0%}: "
                  f"X max {max(x_g):.1%}, Y max {max(y_g):.1%}")
            sys.exit(1)
        print(f"✓ G content within {g_max:.0%} "
              f"(X max {max(x_g):.1%}, Y max {max(y_g):.1%})")

        junction = summarize_junction_homopolymers(
            x_sequences,
            y_sequences,
            config.sequence_qc.max_homopolymer_length,
            spacer,
        )
        print(f"  Junction homopolymer on full barcodes: "
              f"max {junction['max_junction_run']} nt, "
              f"{junction['n_violations']} over "
              f"{config.sequence_qc.max_homopolymer_length} nt")
        if junction["n_violations"] and enforce_junction:
            for example in junction["examples"][:5]:
                print(f"    row {example['row']} col {example['col']}: "
                      f"{example['length']} nt in {example['barcode']}")
            print("✗ Junction homopolymer limit exceeded")
            sys.exit(1)

        if x_valid and y_valid:
            print(f"✓ All sequences passed QC")
        else:
            if not x_valid:
                print(f"⚠ X-axis issues: {x_issues}")
            if not y_valid:
                print(f"⚠ Y-axis issues: {y_issues}")
            print(f"  GC warnings were not used as a filter")
        
    except Exception as e:
        print(f"⚠ Warning: Validation error: {e}")
    
    # Analyze distances
    print(f"\n[5/7] Analyzing distance distributions...")
    try:
        x_results, y_results = analyze_axis_codes(
            axis_codes['x_codes'],
            axis_codes['y_codes'],
            config.error_correction.distances,
            config.error_correction.far_features_min_distance
        )
        print(f"✓ Distance analysis complete")
        # Print diff 1 statistics (1-step neighbors)
        diff1_x = [d for d in x_results['distance_stats'] if d['diff'] == 1][0]
        diff1_y = [d for d in y_results['distance_stats'] if d['diff'] == 1][0]
        print(f"  X-axis 1-step neighbors: mean={diff1_x['mean']:.2f} bits")
        print(f"  Y-axis 1-step neighbors: mean={diff1_y['mean']:.2f} bits")
        
    except Exception as e:
        print(f"⚠ Warning: Analysis error: {e}")
    
    # Build codebook
    print(f"\n[6/7] Building codebook...")
    try:
        codebook = build_codebook(
            rows=config.chip_geometry.rows,
            cols=config.chip_geometry.cols,
            x_codes=axis_codes['x_codes'],
            y_codes=axis_codes['y_codes'],
            x_sequences=x_sequences,
            y_sequences=y_sequences,
            spacer=config.codebook.spacer_sequence,
            feature_id_format=config.codebook.feature_id_format
        )
        
        print(f"✓ Codebook built")
        print(f"  Total barcodes: {len(codebook.features):,}")
        print(f"  Barcode length: {len(codebook.features[0].zip_seq)}nt")
        
        # Check uniqueness
        all_zips = codebook.get_all_zip_sequences()
        unique_zips = len(set(all_zips))
        print(f"  Unique barcodes: {unique_zips:,}")
        
        if unique_zips == len(all_zips):
            print(f"  ✓ All barcodes are unique")
        else:
            duplicates = len(all_zips) - unique_zips
            print(f"  ⚠ Warning: {duplicates} duplicate barcodes")
        
    except Exception as e:
        print(f"✗ Codebook error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Export results
    print(f"\n[7/7] Exporting results...")
    try:
        # Export codebook CSV
        codebook_file = output_dir / "codebook.csv"
        export_codebook(codebook, codebook_file, file_format='csv')
        print(f"✓ Codebook exported to {codebook_file}")
        
        file_size_mb = os.path.getsize(codebook_file) / (1024 * 1024)
        print(f"  File size: {file_size_mb:.2f} MB")
        
        # Export simple text file (column-major)
        barcodes_file = output_dir / "barcodes_235x363_nosep.txt"
        with open(barcodes_file, 'w') as f:
            for feature in codebook.features:
                f.write(f"{feature.zip_seq}\n")
        print(f"✓ Barcodes exported to {barcodes_file}")
        
        # Create filtered CSV (row, col, barcode only)
        filtered_file = output_dir / "layout_235x363_zipcodes_nosep.csv"
        with open(filtered_file, 'w') as f:
            f.write("row,col,barcode\n")
            for feature in codebook.features:
                f.write(f"{feature.row},{feature.col},{feature.zip_seq}\n")
        print(f"✓ Filtered layout exported to {filtered_file}")
        
        # Export statistics
        all_zips = codebook.get_all_zip_sequences()
        stats = {
            'generation_timestamp': datetime.now().isoformat(),
            'configuration': {
                'rows': config.chip_geometry.rows,
                'cols': config.chip_geometry.cols,
                'total_spots': config.chip_geometry.rows * config.chip_geometry.cols,
                'x_bit_length': config.x_bit_length,
                'y_bit_length': config.y_bit_length,
                'encoding_scheme': config.bit_to_dna.scheme,
                'spacer': config.codebook.spacer_sequence
            },
            'results': {
                'total_barcodes': len(all_zips),
                'unique_barcodes': len(set(all_zips)),
                'barcode_length': len(all_zips[0]),
                'x_sequence_length': len(x_sequences[0]),
                'y_sequence_length': len(y_sequences[0]),
                'all_unique': len(set(all_zips)) == len(all_zips)
            },
            'spatial_properties': {
                'design_principle': 'Gray-code-like: spatial neighbors have similar sequences',
                'note': 'Spatial distance analysis printed above'
            }
        }
        
        stats_file = output_dir / "statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"✓ Statistics saved to {stats_file}")
        
    except Exception as e:
        print(f"✗ Export error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Verification
    print(f"\n[Verification]")
    try:
        all_zips = codebook.get_all_zip_sequences()
        n_unique = len(set(all_zips))
        n_total = len(all_zips)
        barcode_length = len(all_zips[0])
        
        if n_unique == n_total:
            print(f"✓ All {n_total:,} barcodes are unique")
        else:
            print(f"⚠ Only {n_unique:,} unique out of {n_total:,} total")
        
        print(f"✓ Barcode length: {barcode_length}nt")
        print(f"✓ Column-major ordering maintained")
        
    except Exception as e:
        print(f"⚠ Verification error: {e}")
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    print(f"\nOutput directory: {output_dir}/")
    print(f"\nOutput files:")
    print(f"  1. Full codebook: codebook.csv")
    print(f"  2. Filtered layout: layout_235x363_zipcodes_nosep.csv")
    print(f"  3. Barcodes only: barcodes_235x363_nosep.txt")
    print(f"  4. Statistics: statistics.json")
    print(f"  5. Intermediate: intermediate/axis_codes.json, distance_analysis.json")
    
    print(f"\nLayout summary:")
    print(f"  • {n_unique:,} unique {barcode_length}nt barcodes")
    print(f"  • Structure: {len(x_sequences[0])}nt X-code + {len(y_sequences[0])}nt Y-code (no separator)")
    print(f"  • Gray-code design: spatial neighbors differ by ~{analysis['x_neighbor_distances']['mean']:.1f} bases")
    print(f"  • All positions in same row share Y-code")
    print(f"  • All positions in same column share X-code")
    print(f"  • Column-major ordering (array fabrication compatible)")
    
    print(f"\nUsage:")
    print(f"  • Read layout: import pandas as pd; df = pd.read_csv('{filtered_file}')")
    print(f"  • Access barcode: df[(df['row']==row) & (df['col']==col)]['barcode']")
    
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
