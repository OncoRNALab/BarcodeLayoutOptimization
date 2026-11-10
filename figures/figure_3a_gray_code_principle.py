#!/usr/bin/env python3
"""
Figure 3A: Zipcodes Method - Spatial Barcode Encoding
Demonstrates how Gray codes ensure spatial similarity constraints
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def binary_to_gray(n):
    """Convert binary to Gray code"""
    return n ^ (n >> 1)

def gray_to_dna(gray_value, num_bits):
    """Convert Gray code to DNA using 2-bit encoding"""
    bases_map = {0b00: 'A', 0b01: 'C', 0b10: 'G', 0b11: 'T'}
    bases = []
    for bit_pos in range(0, num_bits, 2):
        two_bits = (gray_value >> bit_pos) & 0b11
        bases.append(bases_map[two_bits])
    return ''.join(bases)

def create_gray_code_figure():
    """Create visualization showing Zipcodes encoding method"""
    
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.4, top=0.95, bottom=0.05)
    
    # Panel 1: Example encoding for specific positions
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.axis('off')
    
    # Show examples
    examples = []
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (10, 10), (10, 11), (50, 100), (50, 101)]
    
    for row, col in positions:
        row_gray = binary_to_gray(row)
        col_gray = binary_to_gray(col)
        
        # Convert to DNA (using 4 bits for simplicity)
        row_dna = gray_to_dna(row_gray, 8)[:4]  # First 4 bases
        col_dna = gray_to_dna(col_gray, 8)[:4]  # First 4 bases
        
        barcode = row_dna + col_dna
        examples.append([f"({row},{col})", format(row_gray, '08b'), row_dna, 
                        format(col_gray, '08b'), col_dna, barcode])
    
    # Create table
    headers = ['Position', 'Row Gray', 'Row DNA', 'Col Gray', 'Col DNA', 'Barcode']
    
    table_text = "Position Encoding Examples:\n\n"
    table_text += f"{'Position':<12} {'Row Gray':<10} {'Row DNA':<10} {'Col Gray':<10} {'Col DNA':<10} {'Barcode':<12}\n"
    table_text += "-" * 80 + "\n"
    
    for ex in examples:
        table_text += f"{ex[0]:<12} {ex[1]:<10} {ex[2]:<10} {ex[3]:<10} {ex[4]:<10} {ex[5]:<12}\n"
        
        # Add Hamming distance calculation for adjacent positions
        if ex[0] == "(0,1)":
            prev_barcode = examples[0][5]
            curr_barcode = ex[5]
            ham_dist = sum(c1 != c2 for c1, c2 in zip(prev_barcode, curr_barcode))
            table_text += f"             ↑ Hamming distance from (0,0): {ham_dist} base\n"
        elif ex[0] == "(1,0)":
            prev_barcode = examples[0][5]
            curr_barcode = ex[5]
            ham_dist = sum(c1 != c2 for c1, c2 in zip(prev_barcode, curr_barcode))
            table_text += f"             ↑ Hamming distance from (0,0): {ham_dist} base\n"
        elif ex[0] == "(10,11)":
            prev_barcode = examples[4][5]
            curr_barcode = ex[5]
            ham_dist = sum(c1 != c2 for c1, c2 in zip(prev_barcode, curr_barcode))
            table_text += f"             ↑ Hamming distance from (10,10): {ham_dist} base\n"
    
    ax1.text(0.05, 0.95, table_text, ha='left', va='top', 
            fontsize=9, family='monospace', transform=ax1.transAxes)
    ax1.set_title('Encoding Examples', fontsize=13, fontweight='bold')
    
    # Panel 2: Visual spatial grid showing Hamming distances (EXAMPLE DATA)
    ax2 = fig.add_subplot(gs[0, 1])
    
    # Create small grid and calculate Hamming distances
    grid_size = 8
    ref_row, ref_col = 3, 3
    distances = np.zeros((grid_size, grid_size))
    
    ref_row_gray = binary_to_gray(ref_row)
    ref_col_gray = binary_to_gray(ref_col)
    ref_barcode = gray_to_dna(ref_row_gray, 8)[:4] + gray_to_dna(ref_col_gray, 8)[:4]
    
    for i in range(grid_size):
        for j in range(grid_size):
            row_gray = binary_to_gray(i)
            col_gray = binary_to_gray(j)
            barcode = gray_to_dna(row_gray, 8)[:4] + gray_to_dna(col_gray, 8)[:4]
            distances[i, j] = sum(c1 != c2 for c1, c2 in zip(ref_barcode, barcode))
    
    im = ax2.imshow(distances, cmap='RdYlGn_r', aspect='auto')
    ax2.plot(ref_col, ref_row, 'k*', markersize=20, markeredgewidth=2)
    
    plt.colorbar(im, ax=ax2, label='Hamming Distance', fraction=0.046, pad=0.04)
    ax2.set_xlabel('Column', fontsize=11)
    ax2.set_ylabel('Row', fontsize=11)
    ax2.set_title(f'Spatial Similarity Map (8×8 Example Grid)\nReference: position ({ref_row},{ref_col})', 
                 fontsize=13, fontweight='bold')
    ax2.text(0.5, -0.15, 'Note: Using example 8×8 grid for illustration', 
            ha='center', va='top', transform=ax2.transAxes, fontsize=9, style='italic',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    # Panel 3: Distance constraint verification
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.axis('off')
    
    constraint_text = """
Distance Constraint Verification:

Gray Code Property:
  Gray(n) XOR Gray(n+1) = single bit difference
  Example: Gray(5) = 0111, Gray(6) = 0101
           Difference: bit 1 changed (0111 → 0101)

DNA Encoding Preservation:
  1 bit change → affects 1 base (2-bit encoding)
  
  Gray(5) = 0111 → DNA: 01|11 → CG
  Gray(6) = 0101 → DNA: 01|01 → CC
                          ↑ Only 1 base changed (G→C)

Spatial Application:
  Position (r, c):
    - row_barcode encodes r using Gray(r)
    - col_barcode encodes c using Gray(c)
  
  Adjacent positions guarantee:
    (r, c) vs (r+1, c):  1 base difference in row_barcode
    (r, c) vs (r, c+1):  1 base difference in col_barcode
  
  Result: Hamming distance ∝ Manhattan distance
          (Smooth spatial gradient)
    """
    
    ax3.text(0.05, 0.95, constraint_text, ha='left', va='top', 
            fontsize=9, family='monospace', transform=ax3.transAxes,
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
    ax3.set_title('How Distance Constraints Are Ensured', fontsize=13, fontweight='bold')
    
    # Panel 4: Full barcode structure
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    structure_text = """
Complete Barcode Structure:

Configuration A: With Separator (40nt)
┌─────────────┬──────────┬─────────────┬─────────────┐
│ Row Barcode │Separator │ Col Barcode │   Library   │
│   (8 nt)    │ (4 nt)   │   (8 nt)    │   (20 nt)   │
└─────────────┴──────────┴─────────────┴─────────────┘
    ↓              ↓           ↓              ↓
  Spatial      Prevents    Spatial      High diversity
  encoding     run-on      encoding     for sequencing
  (Gray)       errors      (Gray)       (random)

Example for position (10, 20):
  Row(10) → Gray(10) = 00001111 → AACG (8nt)
  Sep              → GGGG (4nt)
  Col(20) → Gray(20) = 00011110 → AAGT (8nt)  
  Library          → ACGTACGT... (20nt random)
  
  Full barcode: AACGGGGGAAGTACGTACGT... (40nt)

Configuration B: No Separator (32nt)
  Same structure without GGGG separator

Trade-offs:
  + Separator: Reduces run-on errors, easier parsing
  - Separator: Longer barcode, potential homopolymer (GGGG)
    """
    
    ax4.text(0.05, 0.95, structure_text, ha='left', va='top', 
            fontsize=8.5, family='monospace', transform=ax4.transAxes,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    ax4.set_title('Full Barcode Architecture', fontsize=13, fontweight='bold')
    
    plt.savefig('figure_3a_gray_code_principle.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_3a_gray_code_principle.pdf', bbox_inches='tight')
    print("✅ Saved figure_3a_gray_code_principle.png and .pdf")
    plt.close()

if __name__ == "__main__":
    create_gray_code_figure()

if __name__ == "__main__":
    create_gray_code_figure()

