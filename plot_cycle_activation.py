#!/usr/bin/env python3
"""
Plot Cycle Activation Pattern

Visualizes which spots are activated at a specific synthesis cycle using
4-bit encoding (Cycle 0=A, 1=C, 2=G, 3=T, 4=A, ...).
Only spots that have the synthesized base at the current cycle position 
are shown as activated.

Usage:
    python plot_cycle_activation.py <layout.csv> <cycle_number> [output.png]
    
Example:
    python plot_cycle_activation.py layout_235x363_checkerboard_colmajor.csv 0
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches


# 4-bit one-hot encoding
BASE_TO_BIT = {
    'A': 0,  # 1000
    'C': 1,  # 0100
    'G': 2,  # 0010
    'T': 3   # 0001
}

BIT_TO_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

# Colors for each channel
CHANNEL_COLORS = {
    0: '#FF6B6B',  # Red for A
    1: '#4ECDC4',  # Cyan for C
    2: '#45B7D1',  # Blue for G
    3: '#FFA07A',  # Orange for T
    -1: '#CCCCCC'  # Gray for empty
}


def load_layout(filepath):
    """Load layout CSV file"""
    df = pd.read_csv(filepath)
    
    # Handle different column name possibilities
    barcode_col = None
    for col in ['barcode', 'zip_seq', 'sequence']:
        if col in df.columns:
            barcode_col = col
            break
    
    if barcode_col is None:
        raise ValueError(f"Could not find barcode column. Available: {df.columns.tolist()}")
    
    if 'row' not in df.columns or 'col' not in df.columns:
        raise ValueError(f"CSV must have 'row' and 'col' columns. Found: {df.columns.tolist()}")
    
    # Rename barcode column for consistency
    df = df.rename(columns={barcode_col: 'barcode'})
    
    return df[['row', 'col', 'barcode']]


def get_base_at_cycle(barcode, cycle):
    """Get the base at the specified cycle (0-indexed)"""
    if cycle < 0 or cycle >= len(barcode):
        raise ValueError(f"Cycle {cycle} out of range for barcode length {len(barcode)}")
    return barcode[cycle].upper()


def get_activated_channel(base):
    """Get the channel (bit) that is activated for this base"""
    return BASE_TO_BIT.get(base, -1)


def create_activation_grid(df, cycle):
    """
    Create a grid showing which spots are activated at 'cycle'.
    A spot is activated if its barcode at position 'cycle' matches 
    the base synthesized in that cycle (Cycle % 4: 0=A, 1=C, 2=G, 3=T).
    
    Returns:
        grid: 2D array with channel number for activated spots, -1 for others
        stats: Dictionary with activation statistics
    """
    max_row = df['row'].max() + 1
    max_col = df['col'].max() + 1
    
    # Initialize grid with -1 (inactive)
    grid = np.full((max_row, max_col), -1, dtype=int)
    
    # Determine target base for this synthesis cycle
    bases = ['A', 'C', 'G', 'T']
    target_base = bases[cycle % 4]
    target_channel = BASE_TO_BIT[target_base]
    
    # Count bases at this position
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    activated_count = 0
    
    # Fill grid only if base at this cycle position matches target_base
    for _, row in df.iterrows():
        r, c, barcode = int(row['row']), int(row['col']), row['barcode']
        
        # Get base at this cycle position
        base = get_base_at_cycle(barcode, cycle)
        base_counts[base] = base_counts.get(base, 0) + 1
        
        # If it matches the synthesized base for this cycle, activate it
        if base == target_base:
            grid[r, c] = target_channel
            activated_count += 1
        else:
            grid[r, c] = -1
    
    # Calculate statistics
    total_spots = len(df)
    stats = {
        'total_spots': total_spots,
        'cycle_position': cycle,
        'target_base': target_base,
        'target_channel': target_channel,
        'activated_count': activated_count,
        'activation_percentage': activated_count / total_spots * 100 if total_spots > 0 else 0,
        'base_counts': base_counts,
        'base_percentages': {base: count/total_spots*100 for base, count in base_counts.items()}
    }
    
    return grid, stats


def plot_activation_pattern(grid, stats, cycle, output_path=None, dpi=150):
    """
    Plot the activation pattern as a colored grid.
    
    Args:
        grid: 2D array with channel values (target_channel for active, -1 for inactive)
        stats: Statistics dictionary
        cycle: Cycle number being visualized
        output_path: Path to save figure (optional)
        dpi: Resolution for saved figure
    """
    fig, (ax_main, ax_stats) = plt.subplots(1, 2, figsize=(16, 8), 
                                             gridspec_kw={'width_ratios': [3, 1]})
    
    target_base = stats['target_base']
    target_channel = stats['target_channel']
    
    # Create colormap: gray for inactive, and the color for the target base
    # We use ListedColormap with 5 entries to match the vmin/vmax logic of the original
    # but only one channel will actually be used.
    colors = [CHANNEL_COLORS[-1]] + [CHANNEL_COLORS[i] for i in range(4)]
    cmap = ListedColormap(colors)
    
    # Plot grid
    im = ax_main.imshow(grid, cmap=cmap, vmin=-1, vmax=3, aspect='auto', interpolation='nearest')
    
    # Formatting
    ax_main.set_xlabel('Column', fontsize=12, fontweight='bold')
    ax_main.set_ylabel('Row', fontsize=12, fontweight='bold')
    ax_main.set_title(f'Synthesis Cycle {cycle} Activation Pattern\n' + 
                      f'Synthesized Base: {target_base} (4-bit encoding: A→C→G→T)', 
                      fontsize=14, fontweight='bold', pad=20)
    
    # Add grid lines for better readability
    rows, cols = grid.shape
    if rows <= 50 and cols <= 100:
        ax_main.set_xticks(np.arange(cols))
        ax_main.set_yticks(np.arange(rows))
        ax_main.grid(which='major', color='white', linestyle='-', linewidth=0.5, alpha=0.3)
    else:
        ax_main.set_xticks(np.arange(0, cols, max(1, cols//20)))
        ax_main.set_yticks(np.arange(0, rows, max(1, rows//20)))
    
    # Legend - Show only Active and Inactive
    legend_elements = [
        mpatches.Patch(facecolor=CHANNEL_COLORS[target_channel], 
                       label=f'Active: {target_base} (Channel {target_channel})'),
        mpatches.Patch(facecolor=CHANNEL_COLORS[-1], label='Inactive / Other Base')
    ]
    ax_main.legend(handles=legend_elements, loc='upper right', fontsize=10, 
                   framealpha=0.9, title=f'Synthesis: {target_base}')
    
    # Statistics panel
    ax_stats.axis('off')
    
    # Title
    ax_stats.text(0.5, 0.95, 'Activation Statistics', 
                  ha='center', va='top', fontsize=14, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    # Cycle info
    y_pos = 0.85
    ax_stats.text(0.1, y_pos, f'Cycle: {cycle}', 
                  ha='left', va='top', fontsize=12, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.08
    ax_stats.text(0.1, y_pos, f'Synthesized Base: {target_base}', 
                  ha='left', va='top', fontsize=11, fontweight='bold',
                  color=CHANNEL_COLORS[target_channel],
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.08
    ax_stats.text(0.1, y_pos, f'Total spots: {stats["total_spots"]:,}', 
                  ha='left', va='top', fontsize=11,
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.08
    ax_stats.text(0.1, y_pos, f'Activated spots: {stats["activated_count"]:,}', 
                  ha='left', va='top', fontsize=11,
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.06
    ax_stats.text(0.1, y_pos, f'Activation rate: {stats["activation_percentage"]:.2f}%', 
                  ha='left', va='top', fontsize=11,
                  transform=ax_stats.transAxes)
    
    # Base distribution at this position (context)
    y_pos -= 0.12
    ax_stats.text(0.1, y_pos, f'Full Distribution at Position {cycle}:', 
                  ha='left', va='top', fontsize=12, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    for base in ['A', 'C', 'G', 'T']:
        y_pos -= 0.08
        count = stats['base_counts'].get(base, 0)
        pct = stats['base_percentages'].get(base, 0)
        channel = BASE_TO_BIT[base]
        is_target = (base == target_base)
        
        # Color box
        rect = mpatches.Rectangle((0.1, y_pos - 0.03), 0.05, 0.05, 
                                   transform=ax_stats.transAxes,
                                   facecolor=CHANNEL_COLORS[channel] if is_target else CHANNEL_COLORS[-1],
                                   edgecolor='black', linewidth=1)
        ax_stats.add_patch(rect)
        
        # Text
        weight = 'bold' if is_target else 'normal'
        ax_stats.text(0.18, y_pos, f'{base}: {count:>6,} ({pct:5.1f}%)', 
                      ha='left', va='top', fontsize=11, family='monospace',
                      fontweight=weight,
                      transform=ax_stats.transAxes)
    
    # Balance check
    y_pos -= 0.12
    counts = [stats['base_counts'].get(base, 0) for base in ['A', 'C', 'G', 'T']]
    balance = (max(counts) - min(counts)) / stats["total_spots"] * 100 if stats["total_spots"] > 0 else 0
    
    ax_stats.text(0.1, y_pos, f'Balance:', 
                  ha='left', va='top', fontsize=12, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.08
    if balance < 2:
        balance_status = "Excellent"
        color = 'green'
    elif balance < 5:
        balance_status = "Good"
        color = 'orange'
    else:
        balance_status = "Poor"
        color = 'red'
    
    ax_stats.text(0.1, y_pos, f'Max-Min: {balance:.2f}%', 
                  ha='left', va='top', fontsize=11,
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.06
    ax_stats.text(0.1, y_pos, f'Status: {balance_status}', 
                  ha='left', va='top', fontsize=11, color=color, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    # Encoding info
    y_pos -= 0.12
    ax_stats.text(0.1, y_pos, '4-bit Encoding (A→C→G→T):', 
                  ha='left', va='top', fontsize=11, fontweight='bold',
                  transform=ax_stats.transAxes)
    
    y_pos -= 0.06
    for i, base in enumerate(['A', 'C', 'G', 'T']):
        ax_stats.text(0.1, y_pos, f'{base} = Channel {i}', 
                      ha='left', va='top', fontsize=10, family='monospace',
                      transform=ax_stats.transAxes)
        y_pos -= 0.05
    
    plt.tight_layout()
    
    # Save or show
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"  ✓ Saved: {output_path}")
    else:
        plt.show()
    
    plt.close()


def main():
    """Main execution"""
    print("=" * 80)
    print("CYCLE ACTIVATION PATTERN VISUALIZER")
    print("4-bit One-Hot Encoding (A→C→G→T)")
    print("=" * 80)
    print()
    
    # Parse arguments
    if len(sys.argv) < 3:
        print("Usage: python plot_cycle_activation.py <layout.csv> <cycle_number> [output.png]")
        print()
        print("Examples:")
        print("  python plot_cycle_activation.py layout.csv 0")
        print("  python plot_cycle_activation.py layout.csv 5 cycle5.png")
        sys.exit(1)
    
    layout_file = sys.argv[1]
    cycle = int(sys.argv[2])
    output_path = sys.argv[3] if len(sys.argv) > 3 else None
    
    print(f"Layout file: {layout_file}")
    print(f"Cycle: {cycle}")
    if output_path:
        print(f"Output: {output_path}")
    print()
    
    # Load layout
    print("[1/3] Loading layout...")
    try:
        df = load_layout(layout_file)
        print(f"  ✓ Loaded {len(df):,} spots")
        
        # Get barcode length
        barcode_len = len(df['barcode'].iloc[0])
        print(f"  ✓ Barcode length: {barcode_len} bases")
        
        if cycle >= barcode_len:
            print(f"  ✗ Error: Cycle {cycle} exceeds barcode length {barcode_len}")
            sys.exit(1)
            
    except Exception as e:
        print(f"  ✗ Error loading layout: {e}")
        sys.exit(1)
    
    # Create activation grid
    print(f"\n[2/3] Analyzing synthesis cycle {cycle}...")
    try:
        grid, stats = create_activation_grid(df, cycle)
        
        target_base = stats['target_base']
        target_channel = stats['target_channel']
        
        print(f"  ✓ Grid dimensions: {grid.shape[0]} × {grid.shape[1]}")
        print(f"  ✓ Synthesis cycle {cycle} uses base: {target_base}")
        print(f"  ✓ Activated spots: {stats['activated_count']:,} ({stats['activation_percentage']:.1f}%)")
        print(f"  ✓ Base distribution at position {cycle}:")
        for base in ['A', 'C', 'G', 'T']:
            count = stats['base_counts'].get(base, 0)
            pct = stats['base_percentages'].get(base, 0)
            marker = "<- ACTIVE" if base == target_base else ""
            print(f"      {base}: {count:>6,} spots ({pct:5.1f}%) {marker}")
        
    except Exception as e:
        print(f"  ✗ Error creating grid: {e}")
        sys.exit(1)
    
    # Plot
    print(f"\n[3/3] Generating visualization...")
    try:
        plot_activation_pattern(grid, stats, cycle, output_path)
        if not output_path:
            print("  ✓ Displaying plot...")
    except Exception as e:
        print(f"  ✗ Error plotting: {e}")
        sys.exit(1)
    
    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
