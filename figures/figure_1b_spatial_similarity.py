#!/usr/bin/env python3
"""
Figure 1B: Dual Barcode Spatial Similarity Map
Shows Hamming distance between adjacent spots, highlighting block boundaries
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
sys.path.append('../Dual_barcodes')
from dual_barcode_layout import DualBarcodeLayout

def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def create_spatial_similarity_map():
    """Create heatmap showing Hamming distances"""
    
    # Generate a small dual barcode layout
    rows, cols = 32, 32  # 8×8 blocks
    block_size = 4
    
    print("Generating dual barcode layout...")
    barcodes_file = '../inhouse_barcodes/barcodes_1in4_half.txt'
    with open(barcodes_file) as f:
        barcodes = [line.strip() for line in f if line.strip()]
    
    layout = DualBarcodeLayout(rows, cols, barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=42)
    
    # Use the combined grid as barcode matrix
    barcode_matrix = grid_combined
    
    # Calculate Hamming distances to right neighbor
    print("Calculating Hamming distances...")
    hamming_right = np.zeros((rows, cols - 1))
    for i in range(rows):
        for j in range(cols - 1):
            hamming_right[i, j] = hamming_distance(
                barcode_matrix[i, j], 
                barcode_matrix[i, j + 1]
            )
    
    # Calculate Hamming distances to bottom neighbor
    hamming_down = np.zeros((rows - 1, cols))
    for i in range(rows - 1):
        for j in range(cols):
            hamming_down[i, j] = hamming_distance(
                barcode_matrix[i, j], 
                barcode_matrix[i + 1, j]
            )
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Left: Horizontal (right) neighbors
    ax1 = axes[0]
    im1 = ax1.imshow(hamming_right, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=36)
    
    # Draw block boundaries
    num_blocks = rows // block_size
    for i in range(0, rows + 1, block_size):
        ax1.axhline(i - 0.5, color='black', linewidth=2)
    for j in range(0, cols, block_size):
        ax1.axvline(j - 0.5, color='black', linewidth=2)
    
    ax1.set_title('Hamming Distance to Right Neighbor', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Column', fontsize=12)
    ax1.set_ylabel('Row', fontsize=12)
    
    # Add colorbar
    cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.set_label('Hamming Distance (bases)', fontsize=12)
    
    # Add annotations
    ax1.text(cols/2, -3, 
            'Low distance (green) within blocks\nHigh distance (red) across block boundaries',
            ha='center', fontsize=10, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    # Right: Vertical (down) neighbors
    ax2 = axes[1]
    im2 = ax2.imshow(hamming_down, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=36)
    
    # Draw block boundaries
    for i in range(0, rows, block_size):
        ax2.axhline(i - 0.5, color='black', linewidth=2)
    for j in range(0, cols + 1, block_size):
        ax2.axvline(j - 0.5, color='black', linewidth=2)
    
    ax2.set_title('Hamming Distance to Bottom Neighbor', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Column', fontsize=12)
    ax2.set_ylabel('Row', fontsize=12)
    
    # Add colorbar
    cbar2 = plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label('Hamming Distance (bases)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('figure_1b_spatial_similarity.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_1b_spatial_similarity.pdf', bbox_inches='tight')
    print("✅ Saved figure_1b_spatial_similarity.png and .pdf")
    
    # Print statistics
    print("\nStatistics:")
    print(f"Within-block horizontal: mean={np.mean(hamming_right[hamming_right < 20]):.1f}, "
          f"std={np.std(hamming_right[hamming_right < 20]):.1f}")
    print(f"Across-block horizontal: mean={np.mean(hamming_right[hamming_right > 20]):.1f}, "
          f"std={np.std(hamming_right[hamming_right > 20]):.1f}")
    
    plt.close()

if __name__ == "__main__":
    create_spatial_similarity_map()
