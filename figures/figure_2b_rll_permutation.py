#!/usr/bin/env python3
"""
Figure 2B: Quadtree RLL Permutation Selection
Shows how RLL algorithm selects base permutations to avoid homopolymers
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools

def create_rll_permutation_figure():
    """Visualize RLL permutation selection process using real quadtree data"""
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.35, top=0.95, bottom=0.05)
    
    # Load real quadtree data
    csv_path = '../Quadtree/outputs/no_homopylmers/layout_118x182_repeated_quadtree_24nt.csv'
    
    print("Loading real quadtree data...")
    data = []
    with open(csv_path, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split(',')
            row, col, barcode = int(parts[0]), int(parts[1]), parts[2]
            data.append((row, col, barcode))
    
    # Extract data for an early cycle (e.g., cycle 3)
    cycle_to_analyze = 3  # Early cycle where more violations expected
    max_run = 3  # Maximum allowed consecutive identical bases
    
    total_spots = len(data)
    quadrants = np.zeros(total_spots, dtype=int)
    prev_bases = np.zeros(total_spots, dtype=int)
    run_lengths = np.zeros(total_spots, dtype=int)
    
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    print(f"Analyzing cycle {cycle_to_analyze}...")
    for i, (row, col, barcode) in enumerate(data):
        if cycle_to_analyze > 0:
            prev_base = barcode[cycle_to_analyze - 1]
            prev_bases[i] = base_to_idx[prev_base]
            
            # Calculate run length
            run_len = 1
            for j in range(cycle_to_analyze - 1, 0, -1):
                if barcode[j] == prev_base:
                    run_len += 1
                else:
                    break
            run_lengths[i] = run_len
        
        # Calculate quadrant for this cycle from row/col coordinates
        # Scale coordinates to determine which quadrant at this level
        scale = 1 << cycle_to_analyze
        u = int(np.floor(row * scale / 118))
        v = int(np.floor(col * scale / 182))
        row_bit = (u >> (cycle_to_analyze - 1)) & 1 if cycle_to_analyze > 0 else 0
        col_bit = (v >> (cycle_to_analyze - 1)) & 1 if cycle_to_analyze > 0 else 0
        quadrants[i] = (row_bit << 1) | col_bit
    
    # Contingency table: [quadrant, prev_base] = count
    contingency_all = np.zeros((4, 4), dtype=int)
    for q, b in zip(quadrants, prev_bases):
        contingency_all[q, b] += 1
    
    # At-risk table: spots where run_length >= max_run (would violate if same base assigned)
    contingency_risk = np.zeros((4, 4), dtype=int)
    for q, b, r in zip(quadrants, prev_bases, run_lengths):
        if r >= max_run:  # Already at or beyond limit
            contingency_risk[q, b] += 1
    
    # Panel 1: Contingency table (all)
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(contingency_all, cmap='Blues', aspect='auto')
    
    ax1.set_xticks(range(4))
    ax1.set_yticks(range(4))
    ax1.set_xticklabels(['A', 'C', 'G', 'T'])
    ax1.set_yticklabels(['Quad 0', 'Quad 1', 'Quad 2', 'Quad 3'])
    ax1.set_xlabel('Previous Base', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Current Quadrant', fontsize=12, fontweight='bold')
    ax1.set_title(f'All Spots Distribution (Cycle {cycle_to_analyze})', fontsize=13, fontweight='bold')
    
    # Add values
    for i in range(4):
        for j in range(4):
            ax1.text(j, i, f'{contingency_all[i, j]}', 
                    ha='center', va='center', fontsize=10, fontweight='bold')
    
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    # Panel 2: At-risk spots (run_length >= max_run)
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(contingency_risk, cmap='Reds', aspect='auto')
    
    ax2.set_xticks(range(4))
    ax2.set_yticks(range(4))
    ax2.set_xticklabels(['A', 'C', 'G', 'T'])
    ax2.set_yticklabels(['Quad 0', 'Quad 1', 'Quad 2', 'Quad 3'])
    ax2.set_xlabel('Previous Base', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Current Quadrant', fontsize=12, fontweight='bold')
    ax2.set_title(f'At-Risk Spots (run_length ≥ {max_run}) - Cycle {cycle_to_analyze}', 
                 fontsize=13, fontweight='bold')
    
    # Add values
    for i in range(4):
        for j in range(4):
            color = 'white' if contingency_risk[i, j] > 500 else 'black'
            ax2.text(j, i, f'{contingency_risk[i, j]}', 
                    ha='center', va='center', fontsize=10, fontweight='bold', color=color)
    
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    
    # Add explanation
    total_at_risk = int(np.sum(contingency_risk))
    ax2.text(0.5, -0.15, 
            f'{total_at_risk:,} spots would violate max_run={max_run}\n' +
            'if assigned the same base as their previous base',
            ha='center', va='top', transform=ax2.transAxes, fontsize=9,
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    # Panel 3: Test all 24 permutations
    ax3 = fig.add_subplot(gs[1, :])
    
    bases = ['A', 'C', 'G', 'T']
    
    # Calculate violations for all permutations
    perms_results = []
    for perm in itertools.permutations(range(4)):
        # Violations = at-risk spots that would still get same base with this permutation
        violations = sum(contingency_risk[q, perm[q]] for q in range(4))
        repeats = sum(contingency_all[q, perm[q]] for q in range(4))
        perms_results.append((perm, violations, repeats))
    
    # Sort by violations, then repeats
    perms_results.sort(key=lambda x: (x[1], x[2]))
    
    # Show top 10 permutations
    top_n = 10
    perm_labels = []
    violation_counts = []
    repeat_counts = []
    
    for i, (perm, viol, reps) in enumerate(perms_results[:top_n]):
        label = ''.join([bases[perm[q]] for q in range(4)])
        perm_labels.append(f'{label}\n({i+1})')
        violation_counts.append(viol)
        repeat_counts.append(reps)
    
    x = np.arange(top_n)
    width = 0.35
    
    bars1 = ax3.bar(x - width/2, violation_counts, width, label='Violations', 
                    color='red', alpha=0.7)
    bars2 = ax3.bar(x + width/2, repeat_counts, width, label='Repeats (all)', 
                    color='blue', alpha=0.7)
    
    ax3.set_xlabel('Permutation (Quad0→Quad1→Quad2→Quad3)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax3.set_title('Top 10 Permutations (sorted by violations, then repeats)', 
                 fontsize=13, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(perm_labels, fontsize=9)
    ax3.legend(fontsize=11)
    ax3.grid(axis='y', alpha=0.3)
    
    # Highlight best permutation
    best_perm = perms_results[0][0]
    best_label = ''.join([bases[best_perm[q]] for q in range(4)])
    ax3.text(0, max(repeat_counts) * 1.1, 
            f'BEST: {best_label}\nViolations: {perms_results[0][1]}\nRepeats: {perms_results[0][2]}',
            fontsize=11, fontweight='bold', ha='center',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # Panel 4: Show how best permutation maps quadrants
    ax4 = fig.add_subplot(gs[2, 0])
    
    best_mapping = [bases[best_perm[q]] for q in range(4)]
    
    mapping_table = [
        ['Quadrant', 'Binary', 'Assigned Base'],
        ['0 (NW)', '00', best_mapping[0]],
        ['1 (NE)', '01', best_mapping[1]],
        ['2 (SW)', '10', best_mapping[2]],
        ['3 (SE)', '11', best_mapping[3]],
    ]
    
    # Draw table
    cell_height = 0.15
    cell_width = [0.3, 0.3, 0.4]
    
    for i, row in enumerate(mapping_table):
        y = 1 - (i + 1) * cell_height
        x_offset = 0
        for j, cell in enumerate(row):
            # Header row
            if i == 0:
                ax4.text(x_offset + cell_width[j]/2, y + cell_height/2, cell,
                        ha='center', va='center', fontsize=11, fontweight='bold',
                        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
            else:
                ax4.text(x_offset + cell_width[j]/2, y + cell_height/2, cell,
                        ha='center', va='center', fontsize=10)
            
            # Draw cell borders
            from matplotlib.patches import Rectangle
            rect = Rectangle((x_offset, y), cell_width[j], cell_height,
                           fill=False, edgecolor='black', linewidth=1)
            ax4.add_patch(rect)
            
            x_offset += cell_width[j]
    
    ax4.set_xlim(0, sum(cell_width))
    ax4.set_ylim(0, 1)
    ax4.axis('off')
    ax4.set_title('Selected Permutation Mapping', fontsize=13, fontweight='bold')
    
    # Panel 5: Explain forced flips
    ax5 = fig.add_subplot(gs[2, 1])
    
    explanation_text = """
RLL Algorithm Steps:

1. Build contingency tables:
   • C_all[q,b] = spots in quadrant q with prev_base b
   • C_risk[q,b] = spots at risk (run_length ≥ max_run)

2. Test all 24 permutations P:
   • For each P: violations = Σ C_risk[q, P[q]]
   • Select P* with minimum violations

3. Apply permutation:
   • Candidate base = P*[quadrant_index]

4. Forced flips (if still violates):
   • If (candidate == prev_base) AND (run_length ≥ max_run):
       offset = (row + col + cycle) % 3 + 1
       final_base = (prev_base + offset) % 4

5. Update state for next cycle:
   • If final_base == prev_base: run_length += 1
   • Else: run_length = 1
   • prev_base = final_base

Note: max_run = 3 means AAA is OK, AAAA is violation
    """
    
    ax5.text(0.05, 0.95, explanation_text, ha='left', va='top', 
            fontsize=9, family='monospace', transform=ax5.transAxes)
    ax5.axis('off')
    ax5.set_title('RLL Process Overview', fontsize=13, fontweight='bold')
    
    plt.savefig('figure_2b_rll_permutation.png', dpi=300, bbox_inches='tight')
    plt.savefig('figure_2b_rll_permutation.pdf', bbox_inches='tight')
    print("✅ Saved figure_2b_rll_permutation.png and .pdf")
    
    # Print best permutation
    print(f"\nBest permutation: {best_label}")
    print(f"  Violations (at-risk spots still violated): {perms_results[0][1]} / {int(np.sum(contingency_risk))}")
    print(f"  Repeats (total same as previous): {perms_results[0][2]} / {total_spots}")
    print(f"  max_run = {max_run}")
    
    plt.close()

if __name__ == "__main__":
    create_rll_permutation_figure()
