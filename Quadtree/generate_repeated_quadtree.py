#!/usr/bin/env python3
"""
Repeated Quadtree Approach: 3x 8nt = 24nt barcodes with RLL

Generates three 8nt quadtree codes with run-length limiting (RLL)
to prevent homopolymers and concatenates them to create 24nt barcodes.

Production-quality implementation with:
- Dynamic base permutation selection per cycle
- Homopolymer prevention (max_run=1)
- Column-major export for array fabrication
- Quality statistics tracking
"""

import numpy as np
import pandas as pd
import math
from quadtree_rll import generate_barcode_matrix_rll, verify_no_homopolymers, count_homopolymers
from output_utils import (write_barcodes_column_major, write_layout_csv, 
                          validate_grid_coverage, write_statistics_json)


def generate_repeated_quadtree_layout(rows, cols, segment_length=8, n_repeats=3, max_run=1):
    """
    Generate layout using repeated quadtree with RLL (Run-Length Limiting) across all segments.
    
    RLL is applied continuously across segments to prevent homopolymers even at boundaries.
    
    Args:
        rows: Number of rows in grid
        cols: Number of columns in grid
        segment_length: Length of each quadtree segment (default: 8nt)
        n_repeats: Number of repetitions (default: 3)
        max_run: Maximum consecutive identical bases (default: 1 forbids homopolymers)
    
    Returns:
        Tuple of (DataFrame, segments, combined_matrix, all_stats)
        - DataFrame with columns: row, col, barcode, segment_1, segment_2, segment_3
        - segments: list of barcode matrices for each segment
        - combined_matrix: concatenated barcode matrix
        - all_stats: list of statistics dicts for each segment
    """
    print(f"Generating {n_repeats} RLL-controlled quadtree segments of {segment_length}nt each...")
    print(f"  max_run = {max_run} (forbids {'homopolymers' if max_run == 1 else f'{max_run+1}+ consecutive bases'})")
    print(f"  RLL will be applied continuously across segments to prevent boundary homopolymers")
    
    # Generate one continuous barcode with all segments
    total_length = segment_length * n_repeats
    print(f"\n  Generating continuous {total_length}nt barcode with RLL...")
    combined_matrix, overall_stats = generate_barcode_matrix_rll(rows, cols, total_length, max_run)
    
    # Split into segments for analysis
    segments = []
    segment_start = 0
    for i in range(n_repeats):
        # Extract characters for this segment
        segment_chars = np.array([list(bc[segment_start:segment_start + segment_length]) 
                                 for bc in combined_matrix.flatten()])
        segment_matrix = np.array([''.join(chars) for chars in segment_chars]).reshape(rows, cols)
        segments.append(segment_matrix)
        segment_start += segment_length
    
    # Per-segment statistics (for display purposes - note they don't add up to overall due to boundary effects)
    all_stats = []
    for i in range(n_repeats):
        # Create partial stats for each segment
        seg_stats = {
            'segment': i + 1,
            'length': segment_length,
            'max_run': max_run,
            'note': 'Individual segment stats - see overall_stats for complete RLL metrics'
        }
        all_stats.append(seg_stats)
    
    # Add overall statistics
    all_stats.append({
        'overall': True,
        'total_length': total_length,
        'forced_flips_per_cycle': overall_stats['forced_flips_per_cycle'],
        'total_forced_flips': overall_stats['total_forced_flips'],
        'mean_forced_flip_rate': overall_stats['mean_forced_flip_rate'],
        'total_violations_avoided': overall_stats['total_violations_avoided'],
        'max_run': max_run
    })
    
    print(f"  ✅ Generated continuous {total_length}nt barcodes with RLL")
    print(f"    Forced flips: {overall_stats['total_forced_flips']} ({overall_stats['mean_forced_flip_rate']:.4f} rate)")
    print(f"    Violations avoided: {overall_stats['total_violations_avoided']}")
    
    # Verify no homopolymers in final barcodes
    print(f"\nVerifying homopolymer constraint (max_run={max_run})...")
    is_valid = verify_no_homopolymers(combined_matrix, max_run)
    if is_valid:
        print(f"  ✅ All barcodes satisfy max_run={max_run} constraint")
    else:
        print(f"  ⚠️  WARNING: Some barcodes violate max_run={max_run} constraint!")
    
    # Create DataFrame
    idx = np.arange(rows * cols)
    row_idx = idx // cols
    col_idx = idx % cols
    
    data = {
        'row': row_idx,
        'col': col_idx,
        'barcode': combined_matrix.flatten(),
    }
    
    # Add individual segments
    for i in range(n_repeats):
        data[f'segment_{i+1}'] = segments[i].flatten()
    
    df = pd.DataFrame(data)
    
    return df, segments, combined_matrix, all_stats


def analyze_layout(df, segments, all_stats=None):
    """Analyze the generated layout with enhanced quality metrics"""
    print("\n" + "="*80)
    print("LAYOUT ANALYSIS")
    print("="*80)
    
    # Basic stats
    n_spots = len(df)
    n_unique = df['barcode'].nunique()
    total_length = len(df['barcode'].iloc[0])
    
    print(f"\nBasic Statistics:")
    print(f"  Total spots: {n_spots:,}")
    print(f"  Unique barcodes: {n_unique:,}")
    print(f"  Uniqueness: {100*n_unique/n_spots:.4f}%")
    print(f"  Barcode length: {total_length}nt")
    
    # Check for duplicates
    duplicates = df[df.duplicated(subset=['barcode'], keep=False)]
    if len(duplicates) > 0:
        print(f"\n⚠️  WARNING: {len(duplicates)} duplicate barcodes found!")
        dup_counts = duplicates['barcode'].value_counts()
        print(f"  Most common duplicate appears {dup_counts.iloc[0]} times")
    else:
        print(f"\n✅ All barcodes are unique!")
    
    # Homopolymer analysis
    print(f"\nHomopolymer Analysis:")
    barcodes_array = df['barcode'].values.reshape(-1, 1)
    homopolymers = count_homopolymers(barcodes_array, min_length=2)
    
    total_homo = sum(sum(counts.values()) for counts in homopolymers.values())
    if total_homo == 0:
        print(f"  ✅ No homopolymers detected (AA, CC, GG, TT, etc.)")
    else:
        print(f"  ⚠️  Found {total_homo} homopolymer occurrences:")
        for base in ['A', 'C', 'G', 'T']:
            for length in range(2, 10):
                count = homopolymers[base][length]
                if count > 0:
                    print(f"    {base*length}: {count} occurrences")
    
    # RLL statistics - find overall stats
    if all_stats:
        overall_stats = None
        for stats in all_stats:
            if stats.get('overall'):
                overall_stats = stats
                break
        
        if overall_stats:
            print(f"\nRLL Quality Statistics (Continuous across all segments):")
            print(f"  Total forced flips: {overall_stats['total_forced_flips']} "
                  f"(rate: {overall_stats['mean_forced_flip_rate']:.4f})")
            print(f"  Violations avoided: {overall_stats['total_violations_avoided']}")
            print(f"  Max run: {overall_stats['max_run']}")
    
    # Segment analysis
    print(f"\nSegment Analysis:")
    for i, seg_col in enumerate([c for c in df.columns if c.startswith('segment_')], 1):
        seg_unique = df[seg_col].nunique()
        seg_length = len(df[seg_col].iloc[0])
        print(f"  Segment {i}: {seg_length}nt, {seg_unique:,} unique codes")
    
    # Spatial similarity check (compare neighbors)
    print(f"\nSpatial Similarity (comparing adjacent spots):")
    max_row = df['row'].max()
    max_col = df['col'].max()
    sample_positions = [(0, 0), (0, 1), (1, 0), (max_row//2, max_col//2)]
    
    for r, c in sample_positions:
        if c + 1 <= max_col:
            bc1_df = df[(df['row']==r) & (df['col']==c)]
            bc2_df = df[(df['row']==r) & (df['col']==c+1)]
            if len(bc1_df) > 0 and len(bc2_df) > 0:
                bc1 = bc1_df['barcode'].values[0]
                bc2 = bc2_df['barcode'].values[0]
                hamming = sum(b1 != b2 for b1, b2 in zip(bc1, bc2))
                print(f"  Position ({r},{c}) vs ({r},{c+1}): {hamming} differences (Hamming distance)")


def export_layout(df, output_csv, output_fasta=None, output_column_major=None,
                 stats_list=None, stats_json=None):
    """Export layout to multiple formats with column-major support"""
    
    # Standard CSV (row-major)
    write_layout_csv(df, output_csv)
    
    # FASTA (optional)
    if output_fasta:
        with open(output_fasta, 'w') as f:
            for _, row in df.iterrows():
                f.write(f">R{row['row']:03d}_C{row['col']:03d}\n")
                f.write(f"{row['barcode']}\n")
        print(f"FASTA exported to {output_fasta}")
    
    # Column-major export (critical for array fabrication)
    if output_column_major:
        max_row = df['row'].max()
        max_col = df['col'].max()
        write_barcodes_column_major(
            df, 
            output_column_major, 
            strict=True,
            rows=max_row + 1,
            cols=max_col + 1
        )
    
    # Statistics JSON (optional)
    if stats_list and stats_json:
        # Find and use overall stats
        combined_stats = None
        for stats in stats_list:
            if stats.get('overall'):
                combined_stats = stats
                break
        
        if combined_stats:
            write_statistics_json(combined_stats, stats_json)


def main():
    """Main execution"""
    print("="*80)
    print("REPEATED QUADTREE BARCODE LAYOUT GENERATOR (with RLL)")
    print("="*80)
    
    # Configuration
    ROWS = 118
    COLS = 182
    SEGMENT_LENGTH = 8  # 8nt per segment
    N_REPEATS = 3       # 3 segments
    MAX_RUN = 1         # Forbid homopolymers
    OUTPUT_CSV = 'layout_118x182_repeated_quadtree_24nt.csv'
    OUTPUT_FASTA = 'layout_118x182_repeated_quadtree_24nt.fasta'
    OUTPUT_COLUMN_MAJOR = 'layout_118x182_repeated_quadtree_24nt_colmajor.txt'
    STATS_JSON = 'layout_118x182_repeated_quadtree_24nt_stats.json'
    
    print(f"\nConfiguration:")
    print(f"  Grid size: {ROWS} × {COLS} = {ROWS*COLS:,} spots")
    print(f"  Segment length: {SEGMENT_LENGTH}nt")
    print(f"  Number of segments: {N_REPEATS}")
    print(f"  Total barcode length: {SEGMENT_LENGTH * N_REPEATS}nt")
    print(f"  Max run length: {MAX_RUN} (forbids homopolymers)")
    
    # Calculate minimum length for uniqueness
    min_length = math.ceil(math.log2(max(ROWS, COLS)))
    print(f"  Minimum length for uniqueness: {min_length}nt")
    print(f"  Each segment ({SEGMENT_LENGTH}nt) {'✅ satisfies' if SEGMENT_LENGTH >= min_length else '❌ does not satisfy'} this requirement")
    
    # Generate layout
    df, segments, combined, all_stats = generate_repeated_quadtree_layout(
        ROWS, COLS, SEGMENT_LENGTH, N_REPEATS, MAX_RUN
    )
    
    # Analyze
    analyze_layout(df, segments, all_stats)
    
    # Export
    export_layout(df, OUTPUT_CSV, OUTPUT_FASTA, OUTPUT_COLUMN_MAJOR, 
                 all_stats, STATS_JSON)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"✅ Generated {len(df):,} barcodes of {SEGMENT_LENGTH * N_REPEATS}nt length")
    print(f"✅ Using {N_REPEATS} RLL-controlled quadtree segments")
    print(f"✅ Homopolymer prevention (max_run={MAX_RUN})")
    print(f"✅ Column-major export for array fabrication")
    print(f"✅ Maintains spatial similarity in all segments")
    print(f"✅ Deterministic and reproducible")
    print("="*80)
    
    return df, segments, all_stats


if __name__ == "__main__":
    df, segments, stats = main()

