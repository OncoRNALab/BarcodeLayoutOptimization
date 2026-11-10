#!/usr/bin/env python3
"""
Figure 3B: Spatial Comparison of All Three Methods
Side-by-side comparison of spatial similarity patterns
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
sys.path.append('../Dual_barcodes')
sys.path.append('../Quadtree')

def hamming_distance(s1, s2):
    """Calculate Hamming distance"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_dual_barcodes(rows, cols, block_size=4):
    """Generate simple dual barcodes (without actual library)"""
    barcodes = np.empty((rows, cols), dtype=object)
    
    num_blocks_row = (rows + block_size - 1) // block_size
    num_blocks_col = (cols + block_size - 1) // block_size
    
    # Simple mock barcodes
    bar_1_length = 8
    bar_2_length = 8
    bases = ['A', 'C', 'G', 'T']
    
    for row in range(rows):
        for col in range(cols):
            block_row = row // block_size
            block_col = col // block_size
            block_idx = block_row * num_blocks_col + block_col
            
            local_row = row % block_size
            local_col = col % block_size
            position_idx = local_row * block_size + local_col
            
            # Generate mock barcodes
            np.random.seed(block_idx)
            bar_1 = ''.join(np.random.choice(bases, bar_1_length))
            
            np.random.seed(position_idx + 1000)
            bar_2 = ''.join(np.random.choice(bases, bar_2_length))
            
            barcodes[row, col] = bar_1 + bar_2
    
    return barcodes

def generate_quadtree_barcodes(rows, cols, length=8):
    """Generate simple quadtree barcodes (no RLL for speed)"""
    barcodes = np.empty((rows, cols), dtype=object)
    bases_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    
    scale = 1 << length
    u = np.floor(np.arange(rows)[:, None] * scale / rows).astype(np.int64)
    v = np.floor(np.arange(cols)[None, :] * scale / cols).astype(np.int64)
    
    for i in range(rows):
        for j in range(cols):
            barcode = ''
            for level in range(length):
                shift = length - 1 - level
                row_bit = (u[i, 0] >> shift) & 1
                col_bit = (v[0, j] >> shift) & 1
                quad_idx = (row_bit << 1) | col_bit
                barcode += bases_map[quad_idx]
            barcodes[i, j] = barcode
    
    return barcodes

def binary_to_gray(n):
    """Convert binary to Gray code"""
    return n ^ (n >> 1)

def generate_zipcode_barcodes(rows, cols):
    """Generate simple zipcode barcodes using Gray codes"""
    barcodes = np.empty((rows, cols), dtype=object)
    bases_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    
    # Calculate bits needed
    row_bits = int(np.ceil(np.log2(rows)))
    col_bits = int(np.ceil(np.log2(cols)))
    
    for i in range(rows):
        for j in range(cols):
            # Convert to Gray code
            row_gray = binary_to_gray(i)
            col_gray = binary_to_gray(j)
            
            # Convert to DNA (2 bits per base)
            row_barcode = ''
            for bit_pos in range(0, row_bits, 2):
                two_bits = (row_gray >> bit_pos) & 0b11
                row_barcode += bases_map[two_bits]
            
            col_barcode = ''
            for bit_pos in range(0, col_bits, 2):
                two_bits = (col_gray >> bit_pos) & 0b11
                col_barcode += bases_map[two_bits]
            
            barcodes[i, j] = row_barcode + col_barcode
    
    return barcodes

def create_spatial_comparison():
    """Create comparison figure of all three methods"""
    
    # Generate barcodes for all methods
    rows, cols = 32, 32
    ref_row, ref_col = 16, 16
    
    print("Generating dual barcodes...")
    dual_barcodes = generate_dual_barcodes(rows, cols, block_size=4)
    
    print("Generating quadtree barcodes...")
    quadtree_barcodes = generate_quadtree_barcodes(rows, cols, length=8)
    
    print("Generating zipcode barcodes...")
    zipcode_barcodes = generate_zipcode_barcodes(rows, cols)
    
    # Calculate Hamming distances to reference point
    methods = [
        ('Dual Barcodes (4×4)', dual_barcodes),
        ('Quadtree', quadtree_barcodes),
        ('Zipcodes (Gray)', zipcode_barcodes),
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    for col_idx, (method_name, barcodes) in enumerate(methods):
        print(f"Processing {method_name}...")
        
        ref_barcode = barcodes[ref_row, ref_col]
        
        # Calculate distances
        distances = np.zeros((rows, cols))
        for i in range(rows):
            for j in range(cols):
                distances[i, j] = hamming_distance(ref_barcode, barcodes[i, j])
        
        # Top row: Heatmap
        ax = axes[0, col_idx]
        im = ax.imshow(distances, cmap='RdYlGn_r', aspect='auto')
        
        # Mark reference point
        ax.plot(ref_col, ref_row, 'k*', markersize=20, markeredgewidth=2)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Hamming Distance', fontsize=10)
        
        ax.set_title(f'{method_name}\nDistance from ({ref_row},{ref_col})', 
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Column', fontsize=10)
        ax.set_ylabel('Row', fontsize=10)
        
        # For dual barcodes, overlay block boundaries
        if 'Dual' in method_name:
            block_size = 4
            for i in range(0, rows + 1, block_size):
                ax.axhline(i - 0.5, color='white', linewidth=2, alpha=0.7)
                ax.axvline(i - 0.5, color='white', linewidth=2, alpha=0.7)
        
        # Bottom row: Profile along horizontal line through reference
        ax2 = axes[1, col_idx]
        
        # Horizontal profile (same row as reference)
        horizontal_profile = distances[ref_row, :]
        ax2.plot(horizontal_profile, 'b-', linewidth=2, label='Horizontal')
        
        # Vertical profile (same col as reference)
        vertical_profile = distances[:, ref_col]
        ax2.plot(vertical_profile, 'r-', linewidth=2, label='Vertical')
        
        ax2.axvline(ref_col, color='k', linestyle='--', alpha=0.5, label='Reference')
        ax2.axvline(ref_row, color='k', linestyle='--', alpha=0.5)
        
        ax2.set_xlabel('Position', fontsize=10)
        ax2.set_ylabel('Hamming Distance', fontsize=10)
        ax2.set_title(f'Distance Profiles', fontsize=11, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(alpha=0.3)
        
        # Print statistics
        neighbor_distances = []
        for i in range(rows):
            for j in range(cols):
                # Right neighbor
                if j < cols - 1:
                    neighbor_distances.append(
                        hamming_distance(barcodes[i, j], barcodes[i, j+1])
                    )
                # Bottom neighbor
                if i < rows - 1:
                    neighbor_distances.append(
                        hamming_distance(barcodes[i, j], barcodes[i+1, j])
                    )
        
        mean_neighbor = np.mean(neighbor_distances)
        std_neighbor = np.std(neighbor_distances)
        
        ax2.text(0.5, 0.95, 
                f'Adjacent neighbors: {mean_neighbor:.1f} ± {std_neighbor:.1f}',
                transform=ax2.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('figure_3b_spatial_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_3b_spatial_comparison.pdf', bbox_inches='tight')
    print("✅ Saved figure_3b_spatial_comparison.png and .pdf")
    plt.close()

if __name__ == "__main__":
    np.random.seed(42)
    create_spatial_comparison()
