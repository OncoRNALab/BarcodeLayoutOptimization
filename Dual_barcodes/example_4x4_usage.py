"""
Example usage for the 4x4 Dual Barcode Layout Generator

Demonstrates how to use the modified dual barcode layout with 4x4 clusters
instead of the original 2x2 clusters.
"""

from dual_barcode_layout import DualBarcodeLayout, generate_random_barcodes, load_barcodes_from_file
import numpy as np


def example_1_minimum_grid():
    """Example 1: Minimum 4x4 grid with 16 barcodes"""
    print("\n" + "="*70)
    print("EXAMPLE 1: Minimum 4x4 Grid (1 block with 16 barcodes)")
    print("="*70)
    
    # Generate 16 barcodes (minimum for one 4x4 block)
    barcodes = generate_random_barcodes(n=16, length=6, seed=123)
    print(f"Using {len(barcodes)} barcodes")
    print(f"First 5: {barcodes[:5]}")
    
    layout = DualBarcodeLayout(rows=4, cols=4, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=123)
    
    print("\nbar_1 Grid (all same within 4x4 block):")
    print(grid_bar1)
    
    print("\nbar_2 Grid (all different within 4x4 block):")
    print(grid_bar2)
    
    stats = layout.get_statistics()
    print(f"\nTotal spots: {stats['total_spots']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    print(f"Number of 4x4 blocks: {stats['n_blocks']}")
    
    return layout


def example_2_medium_grid():
    """Example 2: 8x8 grid with 4 blocks"""
    print("\n" + "="*70)
    print("EXAMPLE 2: 8x8 Grid (4 blocks, 64 spots)")
    print("="*70)
    
    # Need at least 4 blocks worth = 4 barcodes minimum
    # But using 20 gives more combinations
    barcodes = generate_random_barcodes(n=20, length=8, seed=456)
    print(f"Using {len(barcodes)} barcodes of length 8")
    print(f"This allows: {len(barcodes)} blocks × 16 spots/block = {len(barcodes) * 16} maximum spots")
    
    layout = DualBarcodeLayout(rows=8, cols=8, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=456)
    
    stats = layout.get_statistics()
    print(f"\nGrid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']}")
    print(f"Number of 4x4 blocks: {stats['n_blocks']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    print(f"Theoretical max combinations: {stats['theoretical_max_combinations']}")
    
    # Export
    layout.export_layout('example_4x4_layout_8x8.csv')
    print("\n✓ Exported to example_4x4_layout_8x8.csv")
    
    return layout


def example_3_large_array():
    """Example 3: Large 16x16 array (16 blocks, 256 spots)"""
    print("\n" + "="*70)
    print("EXAMPLE 3: Large 16x16 Array (16 blocks, 256 spots)")
    print("="*70)
    
    # With 4x4 blocks, 16x16 grid = 16 blocks
    # Need at least 16 barcodes
    barcodes = generate_random_barcodes(n=30, length=10, seed=789)
    print(f"Using {len(barcodes)} barcodes of length 10")
    
    layout = DualBarcodeLayout(rows=16, cols=16, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=789)
    
    stats = layout.get_statistics()
    print(f"\nGrid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']}")
    print(f"Number of 4x4 blocks: {stats['n_blocks']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    
    # Check bar_1 distribution
    print("\nbar_1 distribution (16 spots per block):")
    bar1_usage = stats['bar1_usage']
    for bc, count in sorted(bar1_usage.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"  {bc}: {count} spots ({count/16:.0f} blocks)")
    
    # Export
    layout.export_layout('example_4x4_layout_16x16.tsv', format='tsv')
    print("\n✓ Exported to example_4x4_layout_16x16.tsv")
    
    # Visualize
    try:
        layout.visualize_layout('example_4x4_layout_16x16.png')
        print("✓ Visualization saved to example_4x4_layout_16x16.png")
    except:
        print("⚠ Matplotlib not available for visualization")
    
    return layout


def example_4_realistic_spatial_array():
    """Example 4: Realistic spatial transcriptomics array (20x20, 400 spots)"""
    print("\n" + "="*70)
    print("EXAMPLE 4: Realistic Spatial Array (20x20, 400 spots)")
    print("="*70)
    
    # 20x20 grid with 4x4 blocks = 25 blocks
    # Need at least 25 barcodes, using 40 for flexibility
    barcodes = generate_random_barcodes(n=40, length=12, seed=999)
    print(f"Using {len(barcodes)} barcodes of length 12")
    
    layout = DualBarcodeLayout(rows=20, cols=20, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=999)
    
    stats = layout.get_statistics()
    print(f"\nGrid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']}")
    print(f"Number of 4x4 blocks: {stats['n_blocks']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    print(f"Coverage: {stats['unique_combined_barcodes']}/{stats['total_spots']} "
          f"({100*stats['unique_combined_barcodes']/stats['total_spots']:.1f}%)")
    
    # Export
    layout.export_layout('example_4x4_layout_20x20.csv')
    print("\n✓ Exported to example_4x4_layout_20x20.csv")
    
    return layout


def example_5_with_barcode_file():
    """Example 5: Using barcodes from file"""
    print("\n" + "="*70)
    print("EXAMPLE 5: Using Barcodes from File")
    print("="*70)
    
    barcode_file = 'barcodes_1in4_half.txt'
    
    try:
        barcodes = load_barcodes_from_file(barcode_file)
        print(f"Loaded {len(barcodes)} barcodes from {barcode_file}")
        print(f"Barcode length: {len(barcodes[0])} bases")
        
        # Calculate what grid size we can support
        max_blocks = len(barcodes)
        max_grid_size = int(np.sqrt(max_blocks * 16))
        print(f"\nWith {len(barcodes)} barcodes and 4x4 clusters:")
        print(f"  Maximum blocks: {max_blocks}")
        print(f"  Maximum spots: {max_blocks * 16}")
        print(f"  Suggested grid size: up to {max_grid_size}x{max_grid_size}")
        
        # Create a 12x12 grid (9 blocks = 144 spots)
        rows, cols = 12, 12
        print(f"\nCreating {rows}x{cols} grid ({rows*cols} spots)...")
        
        layout = DualBarcodeLayout(rows=rows, cols=cols, barcodes=barcodes[:max_blocks])
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=123)
        
        stats = layout.get_statistics()
        print(f"\nGrid dimensions: {stats['grid_dimensions']}")
        print(f"Total spots: {stats['total_spots']}")
        print(f"Number of 4x4 blocks: {stats['n_blocks']}")
        print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
        
        layout.export_layout('example_4x4_from_file.csv')
        print("\n✓ Exported to example_4x4_from_file.csv")
        
        return layout
        
    except FileNotFoundError:
        print(f"⚠ Barcode file '{barcode_file}' not found in current directory")
        print("  This example requires a barcode file to run")
        return None


def comparison_2x2_vs_4x4():
    """Comparison: 2x2 vs 4x4 approach for same array size"""
    print("\n" + "="*70)
    print("COMPARISON: 2x2 vs 4x4 Cluster Approach")
    print("="*70)
    
    array_size = (16, 16)  # 256 spots
    
    print(f"\nFor a {array_size[0]}x{array_size[1]} array ({array_size[0]*array_size[1]} spots):\n")
    
    # 2x2 approach (theoretical)
    blocks_2x2 = ((array_size[0] + 1) // 2) * ((array_size[1] + 1) // 2)
    spots_per_block_2x2 = 4
    print("2x2 Cluster Approach:")
    print(f"  - Blocks: {blocks_2x2}")
    print(f"  - Spots per block: {spots_per_block_2x2}")
    print(f"  - Minimum barcodes needed: {spots_per_block_2x2}")
    print(f"  - Recommended barcodes: {blocks_2x2} (one unique bar_1 per block)")
    print(f"  - Maximum capacity with {blocks_2x2} barcodes: {blocks_2x2 * spots_per_block_2x2} spots")
    
    # 4x4 approach (current implementation)
    blocks_4x4 = ((array_size[0] + 3) // 4) * ((array_size[1] + 3) // 4)
    spots_per_block_4x4 = 16
    print("\n4x4 Cluster Approach:")
    print(f"  - Blocks: {blocks_4x4}")
    print(f"  - Spots per block: {spots_per_block_4x4}")
    print(f"  - Minimum barcodes needed: {spots_per_block_4x4}")
    print(f"  - Recommended barcodes: {blocks_4x4} (one unique bar_1 per block)")
    print(f"  - Maximum capacity with {blocks_4x4} barcodes: {blocks_4x4 * spots_per_block_4x4} spots")
    
    print("\n✅ Advantage of 4x4: Needs only {:.1f}% as many barcodes ({} vs {})".format(
        100 * blocks_4x4 / blocks_2x2, blocks_4x4, blocks_2x2))
    print("❌ Disadvantage of 4x4: Lower spatial resolution for error detection")


if __name__ == "__main__":
    print("="*70)
    print("4x4 DUAL BARCODE LAYOUT GENERATOR - EXAMPLES")
    print("="*70)
    
    # Run all examples
    try:
        example_1_minimum_grid()
        example_2_medium_grid()
        example_3_large_array()
        example_4_realistic_spatial_array()
        example_5_with_barcode_file()
        comparison_2x2_vs_4x4()
        
        print("\n" + "="*70)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("="*70)
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
