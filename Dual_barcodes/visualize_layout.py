#!/usr/bin/env python3
"""
Visualize dual barcode layout from CSV file.

This script reads a barcode layout CSV file and generates visualizations
showing bar_1 and bar_2 distributions with block boundaries.

Usage:
    python visualize_layout.py <csv_file> [output_image] [--combined]

Examples:
    python visualize_layout.py my_layout_118x182.csv
    python visualize_layout.py my_layout_118x182.csv layout_viz.png
    python visualize_layout.py layout.csv viz.png --combined
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def load_layout_from_csv(filepath):
    """
    Load barcode layout from CSV file.
    
    Args:
        filepath: Path to CSV file
    
    Returns:
        Tuple of (grid_bar1, grid_bar2, grid_combined, rows, cols)
    """
    # Read CSV
    df = pd.read_csv(filepath)
    
    # Get grid dimensions
    rows = df['row'].max() + 1
    cols = df['col'].max() + 1
    
    # Create grids
    grid_bar1 = np.empty((rows, cols), dtype=object)
    grid_bar2 = np.empty((rows, cols), dtype=object)
    grid_combined = np.empty((rows, cols), dtype=object)
    
    # Fill grids
    for _, row in df.iterrows():
        i, j = int(row['row']), int(row['col'])
        grid_bar1[i, j] = row['bar_1']
        grid_bar2[i, j] = row['bar_2']
        grid_combined[i, j] = row['combined']
    
    return grid_bar1, grid_bar2, grid_combined, rows, cols


def visualize_layout(grid_bar1, grid_bar2, grid_combined, rows, cols, 
                     output_file=None, include_combined=False):
    """
    Create visualization of barcode layout.
    
    Args:
        grid_bar1: Bar_1 grid
        grid_bar2: Bar_2 grid
        grid_combined: Combined barcode grid
        rows: Number of rows
        cols: Number of columns
        output_file: Optional path to save figure
        include_combined: If True, include combined layout panel
    """
    # Create color maps
    unique_bar1 = list(set(grid_bar1.flatten()))
    unique_bar2 = list(set(grid_bar2.flatten()))
    
    bar1_colors = {bc: i for i, bc in enumerate(unique_bar1) if bc is not None}
    bar2_colors = {bc: i for i, bc in enumerate(unique_bar2) if bc is not None}
    
    # Determine figure layout
    if include_combined:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        ax_bar1, ax_bar2, ax_combined = axes[0], axes[1], axes[2]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        ax_bar1, ax_bar2 = axes[0], axes[1]
    
    # Plot bar_1
    bar1_numeric = np.array([[bar1_colors.get(grid_bar1[i, j], 0) 
                              for j in range(cols)] 
                             for i in range(rows)])
    im1 = ax_bar1.imshow(bar1_numeric, cmap='tab20', aspect='auto', interpolation='nearest')
    ax_bar1.set_title(f'bar_1 Layout ({rows}×{cols})\n2×2 blocks at (0,0)', 
                     fontsize=12, fontweight='bold')
    ax_bar1.set_xlabel('Column')
    ax_bar1.set_ylabel('Row')
    
    # Draw bar_1 block boundaries (standard 2x2 blocks starting at 0,0)
    for i in range(0, rows, 2):
        ax_bar1.axhline(i - 0.5, color='black', linewidth=1.5, alpha=0.7)
    for j in range(0, cols, 2):
        ax_bar1.axvline(j - 0.5, color='black', linewidth=1.5, alpha=0.7)
    
    # Add statistics text
    n_unique_bar1 = len(unique_bar1)
    ax_bar1.text(0.02, 0.98, f'Unique bar_1: {n_unique_bar1}', 
                transform=ax_bar1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.8))
    
    # Plot bar_2
    bar2_numeric = np.array([[bar2_colors.get(grid_bar2[i, j], 0) 
                              for j in range(cols)] 
                             for i in range(rows)])
    im2 = ax_bar2.imshow(bar2_numeric, cmap='tab20', aspect='auto', interpolation='nearest')
    ax_bar2.set_title(f'bar_2 Layout ({rows}×{cols})\n2×2 blocks offset at (1,1)', 
                     fontsize=12, fontweight='bold')
    ax_bar2.set_xlabel('Column')
    ax_bar2.set_ylabel('Row')
    
    # Draw bar_2 offset block boundaries (2x2 blocks starting at 1,1)
    for i in range(1, rows, 2):
        ax_bar2.axhline(i - 0.5, color='red', linewidth=1.5, alpha=0.7)
    for j in range(1, cols, 2):
        ax_bar2.axvline(j - 0.5, color='red', linewidth=1.5, alpha=0.7)
    
    # Add statistics text
    n_unique_bar2 = len(unique_bar2)
    ax_bar2.text(0.02, 0.98, f'Unique bar_2: {n_unique_bar2}', 
                transform=ax_bar2.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', 
                facecolor='lightblue', alpha=0.8))
    
    # Plot combined (if requested)
    if include_combined:
        unique_combined = list(set(grid_combined.flatten()))
        combined_colors = {bc: i for i, bc in enumerate(unique_combined) if bc is not None}
        
        combined_numeric = np.array([[combined_colors.get(grid_combined[i, j], 0) 
                                     for j in range(cols)] 
                                    for i in range(rows)])
        im3 = ax_combined.imshow(combined_numeric, cmap='tab20', aspect='auto', interpolation='nearest')
        ax_combined.set_title(f'Combined Layout ({rows}×{cols})\nbar_1 + bar_2', 
                            fontsize=12, fontweight='bold')
        ax_combined.set_xlabel('Column')
        ax_combined.set_ylabel('Row')
        
        # Add statistics text
        n_unique_combined = len(unique_combined)
        total_spots = rows * cols
        ax_combined.text(0.02, 0.98, 
                        f'Unique: {n_unique_combined}\nTotal: {total_spots}\nCoverage: {100*n_unique_combined/total_spots:.1f}%', 
                        transform=ax_combined.transAxes, fontsize=10,
                        verticalalignment='top', bbox=dict(boxstyle='round', 
                        facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Visualization saved to: {output_file}")
    else:
        plt.show()
    
    plt.close()


def print_usage():
    """Print usage information."""
    print("Usage: python visualize_layout.py <csv_file> [output_image] [--combined]")
    print()
    print("Arguments:")
    print("  csv_file     : Path to CSV layout file (required)")
    print("  output_image : Path for output PNG file (optional, default: <csv_basename>.png)")
    print("  --combined   : Include combined barcode panel (optional)")
    print()
    print("Examples:")
    print("  python visualize_layout.py my_layout_118x182.csv")
    print("  python visualize_layout.py my_layout_118x182.csv layout_viz.png")
    print("  python visualize_layout.py layout.csv viz.png --combined")


def main():
    # Check for help flag
    if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help', 'help']:
        print_usage()
        sys.exit(0)
    
    # Parse arguments
    csv_file = sys.argv[1]
    
    # Check if csv file exists
    if not os.path.exists(csv_file):
        print(f"ERROR: File not found: {csv_file}")
        sys.exit(1)
    
    # Determine output file
    if len(sys.argv) > 2 and not sys.argv[2].startswith('--'):
        output_file = sys.argv[2]
        include_combined = '--combined' in sys.argv[3:]
    else:
        # Default output: replace .csv with .png
        base_name = os.path.splitext(csv_file)[0]
        output_file = f"{base_name}_visualization.png"
        include_combined = '--combined' in sys.argv[2:]
    
    print("=" * 70)
    print("Dual Barcode Layout Visualization")
    print("=" * 70)
    print(f"\nInput file: {csv_file}")
    print(f"Output file: {output_file}")
    print(f"Include combined panel: {include_combined}")
    
    # Load layout
    print(f"\nLoading layout from CSV...")
    try:
        grid_bar1, grid_bar2, grid_combined, rows, cols = load_layout_from_csv(csv_file)
        print(f"✓ Loaded {rows}×{cols} grid ({rows*cols:,} spots)")
    except Exception as e:
        print(f"✗ ERROR loading CSV: {e}")
        sys.exit(1)
    
    # Calculate statistics
    unique_bar1 = len(set(grid_bar1.flatten()))
    unique_bar2 = len(set(grid_bar2.flatten()))
    unique_combined = len(set(grid_combined.flatten()))
    total_spots = rows * cols
    
    print(f"\nLayout statistics:")
    print(f"  Unique bar_1 barcodes: {unique_bar1}")
    print(f"  Unique bar_2 barcodes: {unique_bar2}")
    print(f"  Unique combined barcodes: {unique_combined:,}")
    print(f"  Total spots: {total_spots:,}")
    print(f"  Coverage: {100*unique_combined/total_spots:.2f}%")
    
    # Check for duplicates
    if unique_combined < total_spots:
        print(f"  ⚠ WARNING: {total_spots - unique_combined} duplicate combined barcodes found!")
    else:
        print(f"  ✓ All combined barcodes are unique")
    
    # Generate visualization
    print(f"\nGenerating visualization...")
    try:
        visualize_layout(grid_bar1, grid_bar2, grid_combined, rows, cols, 
                        output_file, include_combined)
    except Exception as e:
        print(f"✗ ERROR generating visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)


if __name__ == "__main__":
    # Check if required packages are available
    try:
        import pandas
        import matplotlib
    except ImportError as e:
        print(f"ERROR: Missing required package: {e}")
        print("\nInstall required packages with:")
        print("  pip install pandas matplotlib")
        sys.exit(1)
    
    main()
