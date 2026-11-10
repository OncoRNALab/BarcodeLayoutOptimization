#!/usr/bin/env python3
"""
Figure 1A: Dual Barcode Block Structure
Shows how 4×4 blocks are structured with shared bar_1 and unique bar_2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import sys
sys.path.append('../Dual_barcodes')

def create_block_structure_figure():
    """Create visual showing dual barcode block structure"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left panel: Single block with bar_1 + bar_2 structure
    ax1 = axes[0]
    block_size = 4
    
    # Draw grid
    for i in range(block_size + 1):
        ax1.axhline(i, color='black', linewidth=2)
        ax1.axvline(i, color='black', linewidth=2)
    
    # Color entire block (shared bar_1)
    block_color = '#FFE5B4'  # Peach
    rect = Rectangle((0, 0), block_size, block_size, 
                     facecolor=block_color, alpha=0.5, edgecolor='none')
    ax1.add_patch(rect)
    
    # Add position labels (bar_2 identifiers)
    for i in range(block_size):
        for j in range(block_size):
            position_id = i * block_size + j
            # Center of cell
            x = j + 0.5
            y = block_size - i - 0.5
            
            ax1.text(x, y, f'{position_id:02d}', 
                    ha='center', va='center', fontsize=14, fontweight='bold')
    
    # Add annotations
    ax1.text(2, -0.8, 'bar_1 = ACGTACGTACGTACGT', 
            ha='center', fontsize=12, fontweight='bold', 
            bbox=dict(boxstyle='round', facecolor=block_color, alpha=0.7))
    ax1.text(2, 4.5, 'All 16 spots share bar_1', 
            ha='center', fontsize=11, style='italic')
    
    # Highlight example spot
    example_i, example_j = 1, 2
    rect_highlight = Rectangle((example_j, block_size - example_i - 1), 1, 1,
                               facecolor='none', edgecolor='red', linewidth=4)
    ax1.add_patch(rect_highlight)
    
    ax1.text(5, block_size - example_i - 0.5, 
            f'Spot (1,2):\nbar_2 = position 06\nFinal = bar_1 + bar_2', 
            fontsize=10, va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='red'))
    
    ax1.set_xlim(-0.5, 6.5)
    ax1.set_ylim(-1.5, 5)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('Single 4×4 Block Structure', fontsize=14, fontweight='bold')
    
    # Right panel: Multiple blocks in array
    ax2 = axes[1]
    
    # Show 3×3 array of blocks (12×12 spots)
    array_blocks = 3
    array_size = array_blocks * block_size
    
    # Define colors for each block
    colors = plt.cm.Set3(np.linspace(0, 1, array_blocks * array_blocks))
    
    for block_row in range(array_blocks):
        for block_col in range(array_blocks):
            block_idx = block_row * array_blocks + block_col
            color = colors[block_idx]
            
            # Draw block rectangle
            x = block_col * block_size
            y = (array_blocks - block_row - 1) * block_size
            rect = Rectangle((x, y), block_size, block_size,
                           facecolor=color, alpha=0.6, 
                           edgecolor='black', linewidth=2)
            ax2.add_patch(rect)
            
            # Add block label
            center_x = x + block_size / 2
            center_y = y + block_size / 2
            ax2.text(center_x, center_y, f'Block {block_idx}', 
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Draw internal grid lines (lighter)
    for i in range(array_size + 1):
        ax2.axhline(i, color='gray', linewidth=0.5, alpha=0.5)
        ax2.axvline(i, color='gray', linewidth=0.5, alpha=0.5)
    
    # Draw block boundaries (heavier)
    for i in range(0, array_size + 1, block_size):
        ax2.axhline(i, color='black', linewidth=3)
        ax2.axvline(i, color='black', linewidth=3)
    
    ax2.set_xlim(-0.5, array_size + 0.5)
    ax2.set_ylim(-0.5, array_size + 0.5)
    ax2.set_aspect('equal')
    ax2.axis('off')
    ax2.set_title('Multiple Blocks in Array (12×12 spots)', 
                 fontsize=14, fontweight='bold')
    
    # Add legend
    legend_text = 'Each color = different bar_1\nWithin block: same bar_1, different bar_2'
    ax2.text(array_size / 2, -1.5, legend_text, 
            ha='center', fontsize=10, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig('figure_1a_block_structure.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_1a_block_structure.pdf', bbox_inches='tight')
    print("✅ Saved figure_1a_block_structure.png and .pdf")
    plt.close()

if __name__ == "__main__":
    create_block_structure_figure()
