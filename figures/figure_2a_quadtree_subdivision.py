#!/usr/bin/env python3
"""
Figure 2A: Quadtree Hierarchical Subdivision
Shows how the quadtree recursively subdivides space
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

def create_quadtree_subdivision():
    """Create visualization of quadtree spatial subdivision"""
    
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    
    # Colors for quadrants: NW, NE, SW, SE
    colors = ['#FFB6C1', '#87CEEB', '#90EE90', '#FFD700']  # Pink, Blue, Green, Gold
    base_labels = ['A', 'C', 'G', 'T']
    quadrant_names = ['NW (00)', 'NE (01)', 'SW (10)', 'SE (11)']
    
    size = 8  # Grid size
    
    # Level 0: Whole array
    ax = axes[0, 0]
    rect = Rectangle((0, 0), size, size, facecolor='lightgray', alpha=0.5)
    ax.add_patch(rect)
    ax.text(size/2, size/2, 'Level 0\nWhole Array', 
           ha='center', va='center', fontsize=12, fontweight='bold')
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Cycle 0: Root', fontsize=12, fontweight='bold')
    
    # Level 1: 4 quadrants
    ax = axes[0, 1]
    for qi, (qx, qy, qname, color, base) in enumerate([
        (0, size/2, 'NW', colors[0], 'A'),  # Top-left
        (size/2, size/2, 'NE', colors[1], 'C'),  # Top-right
        (0, 0, 'SW', colors[2], 'G'),  # Bottom-left
        (size/2, 0, 'SE', colors[3], 'T'),  # Bottom-right
    ]):
        rect = Rectangle((qx, qy), size/2, size/2, 
                        facecolor=color, alpha=0.6, edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(qx + size/4, qy + size/4, f'{qname}\n→ {base}', 
               ha='center', va='center', fontsize=10, fontweight='bold')
    
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Cycle 1: 4 Quadrants', fontsize=12, fontweight='bold')
    
    # Level 2: 16 sub-quadrants (show NW quadrant subdivided)
    ax = axes[0, 2]
    # Draw other quadrants lightly
    for qx, qy, color in [(size/2, size/2, colors[1]), (0, 0, colors[2]), (size/2, 0, colors[3])]:
        rect = Rectangle((qx, qy), size/2, size/2, 
                        facecolor=color, alpha=0.3, edgecolor='gray', linewidth=1)
        ax.add_patch(rect)
    
    # Subdivide NW quadrant
    for i, (sqx, sqy, color, label) in enumerate([
        (0, size/2 + size/4, colors[0], 'AA'),  # NW of NW
        (size/4, size/2 + size/4, colors[1], 'AC'),  # NE of NW
        (0, size/2, colors[2], 'AG'),  # SW of NW
        (size/4, size/2, colors[3], 'AT'),  # SE of NW
    ]):
        rect = Rectangle((sqx, sqy), size/4, size/4, 
                        facecolor=color, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(sqx + size/8, sqy + size/8, label, 
               ha='center', va='center', fontsize=9, fontweight='bold')
    
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Cycle 2: NW → 4 Sub-quadrants', fontsize=12, fontweight='bold')
    
    # Level 3: 64 sub-sub-quadrants (show NW-NW subdivided)
    ax = axes[0, 3]
    # Draw background lightly
    rect = Rectangle((0, 0), size, size, facecolor='lightgray', alpha=0.2)
    ax.add_patch(rect)
    
    # Draw NW quadrant subdivisions
    for sqx, sqy, color in [(size/4, size/2 + size/4, colors[1]), 
                            (0, size/2, colors[2]), 
                            (size/4, size/2, colors[3])]:
        rect = Rectangle((sqx, sqy), size/4, size/4, 
                        facecolor=color, alpha=0.3, edgecolor='gray', linewidth=1)
        ax.add_patch(rect)
    
    # Subdivide NW-NW
    for i, (ssx, ssy, color, label) in enumerate([
        (0, size/2 + 3*size/8, colors[0], 'AAA'),
        (size/8, size/2 + 3*size/8, colors[1], 'AAC'),
        (0, size/2 + size/4, colors[2], 'AAG'),
        (size/8, size/2 + size/4, colors[3], 'AAT'),
    ]):
        rect = Rectangle((ssx, ssy), size/8, size/8, 
                        facecolor=color, alpha=0.8, edgecolor='black', linewidth=1)
        ax.add_patch(rect)
        ax.text(ssx + size/16, ssy + size/16, label, 
               ha='center', va='center', fontsize=7, fontweight='bold')
    
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Cycle 3: NW-NW → 4 Regions', fontsize=12, fontweight='bold')
    
    # Bottom row: Show binary encoding process
    ax = axes[1, 0]
    ax.text(0.5, 0.9, 'Binary Encoding Process', ha='center', fontsize=13, fontweight='bold',
           transform=ax.transAxes)
    
    encoding_text = """
Position (x=3, y=6) in 8×8 grid:

Binary representation:
  x = 3 = 011₂
  y = 6 = 110₂

Cycle 1: Extract MSB
  x[2]=0, y[2]=1 → quad=01 → C

Cycle 2: Extract bit 1
  x[1]=1, y[1]=1 → quad=11 → T

Cycle 3: Extract LSB
  x[0]=1, y[0]=0 → quad=10 → G

Result: CTG (3 cycles)
    """
    ax.text(0.1, 0.5, encoding_text, ha='left', va='center', 
           fontsize=9, family='monospace', transform=ax.transAxes)
    ax.axis('off')
    
    # Show legend
    ax = axes[1, 1]
    legend_elements = []
    for color, base, qname in zip(colors, base_labels, quadrant_names):
        patch = mpatches.Patch(facecolor=color, alpha=0.6, edgecolor='black', 
                              label=f'{qname} → {base}')
        legend_elements.append(patch)
    
    ax.legend(handles=legend_elements, loc='center', fontsize=11, frameon=True)
    ax.text(0.5, 0.9, 'Quadrant Mapping', ha='center', fontsize=13, fontweight='bold',
           transform=ax.transAxes)
    ax.axis('off')
    
    # Show full grid at final level
    ax = axes[1, 2]
    grid_size = 8
    barcode_length = 3
    
    # Generate simple quadtree barcodes (no RLL for visualization)
    bases_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    barcodes = np.empty((grid_size, grid_size), dtype=object)
    
    for i in range(grid_size):
        for j in range(grid_size):
            barcode = ''
            for level in range(barcode_length):
                shift = barcode_length - 1 - level
                row_bit = (i >> shift) & 1
                col_bit = (j >> shift) & 1
                quad_idx = (row_bit << 1) | col_bit
                barcode += bases_map[quad_idx]
            barcodes[i, j] = barcode
    
    # Display as text grid
    for i in range(grid_size):
        for j in range(grid_size):
            ax.text(j, grid_size - 1 - i, barcodes[i, j], 
                   ha='center', va='center', fontsize=7, family='monospace')
    
    ax.set_xlim(-0.5, grid_size - 0.5)
    ax.set_ylim(-0.5, grid_size - 0.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title('Final 8×8 Barcode Grid (3-cycle)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Column', fontsize=10)
    ax.set_ylabel('Row', fontsize=10)
    
    # Spatial similarity demonstration
    ax = axes[1, 3]
    ref_i, ref_j = 3, 3
    ref_barcode = barcodes[ref_i, ref_j]
    
    hamming_grid = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            hamming_grid[i, j] = sum(c1 != c2 for c1, c2 in zip(ref_barcode, barcodes[i, j]))
    
    im = ax.imshow(hamming_grid, cmap='RdYlGn_r', vmin=0, vmax=3)
    ax.plot(ref_j, ref_i, 'k*', markersize=15, markeredgewidth=2, markerfacecolor='none')
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Hamming Distance', fontsize=10)
    
    ax.set_title(f'Distance from Position ({ref_i},{ref_j})', fontsize=12, fontweight='bold')
    ax.set_xlabel('Column', fontsize=10)
    ax.set_ylabel('Row', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('figure_2a_quadtree_subdivision.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_2a_quadtree_subdivision.pdf', bbox_inches='tight')
    print("✅ Saved figure_2a_quadtree_subdivision.png and .pdf")
    plt.close()

if __name__ == "__main__":
    create_quadtree_subdivision()
