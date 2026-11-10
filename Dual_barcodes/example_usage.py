"""
Example usage scenarios for the Dual Barcode Layout Generator
"""

from dual_barcode_layout import DualBarcodeLayout, generate_random_barcodes
import numpy as np


def example_1_small_grid():
    """Example 1: Small 4x4 grid with 4 barcodes"""
    print("\n" + "="*70)
    print("EXAMPLE 1: Small 4x4 Grid with Minimum Barcodes")
    print("="*70)
    
    barcodes = ['AAAA', 'CCCC', 'GGGG', 'TTTT']
    layout = DualBarcodeLayout(rows=4, cols=4, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=123)
    
    print("\nbar_1 Grid:")
    print(grid_bar1)
    
    print("\nbar_2 Grid:")
    print(grid_bar2)
    
    print("\nCombined Barcode Grid:")
    print(grid_combined)
    
    stats = layout.get_statistics()
    print(f"\nUnique combined barcodes: {stats['unique_combined_barcodes']}")
    
    return layout


def example_2_rectangular_grid():
    """Example 2: Rectangular grid (6x10)"""
    print("\n" + "="*70)
    print("EXAMPLE 2: Rectangular 6x10 Grid")
    print("="*70)
    
    # Generate 8 random barcodes
    barcodes = generate_random_barcodes(n=8, length=6, seed=456)
    print(f"Using {len(barcodes)} barcodes: {barcodes}")
    
    layout = DualBarcodeLayout(rows=6, cols=10, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=456)
    
    stats = layout.get_statistics()
    print(f"\nGrid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    print(f"Coverage: {stats['unique_combined_barcodes']}/{stats['total_spots']} " 
          f"({100*stats['unique_combined_barcodes']/stats['total_spots']:.1f}%)")
    
    # Export
    layout.export_layout('example_2_layout.csv')
    print("\n✓ Exported to example_2_layout.csv")
    
    return layout


def example_3_large_array():
    """Example 3: Large array simulation (20x20)"""
    print("\n" + "="*70)
    print("EXAMPLE 3: Large 20x20 Array")
    print("="*70)
    
    # Use 10 barcodes for better diversity
    barcodes = generate_random_barcodes(n=10, length=8, seed=789)
    print(f"Using {len(barcodes)} barcodes of length 8")
    
    layout = DualBarcodeLayout(rows=20, cols=20, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=789)
    
    stats = layout.get_statistics()
    print(f"\nGrid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']}")
    print(f"Number of 2x2 blocks: {stats['n_blocks']}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    print(f"Theoretical max with {len(barcodes)} barcodes: {stats['theoretical_max_combinations']}")
    
    # Check bar_1 distribution
    print("\nbar_1 distribution across blocks:")
    bar1_usage = stats['bar1_usage']
    for bc, count in sorted(bar1_usage.items(), key=lambda x: x[1], reverse=True):
        print(f"  {bc}: {count} spots ({count/4:.0f} blocks)")
    
    # Export
    layout.export_layout('example_3_layout.tsv', format='tsv')
    print("\n✓ Exported to example_3_layout.tsv")
    
    return layout


def example_4_edge_cases():
    """Example 4: Edge cases with odd dimensions"""
    print("\n" + "="*70)
    print("EXAMPLE 4: Edge Cases - Odd Dimensions (5x7)")
    print("="*70)
    
    barcodes = ['ACGT', 'TGCA', 'GCTA', 'ATCG', 'CGAT', 'TACG']
    layout = DualBarcodeLayout(rows=5, cols=7, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=111)
    
    print("\nCombined Barcode Grid (5x7):")
    print(grid_combined)
    
    stats = layout.get_statistics()
    print(f"\nTotal spots: {stats['total_spots']}")
    print(f"Number of 2x2 blocks: {stats['n_blocks']} (some incomplete at edges)")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']}")
    
    # Verify that edge positions are properly handled
    print("\nEdge position examples:")
    print(f"  Position (0, 6): bar_1={grid_bar1[0, 6]}, bar_2={grid_bar2[0, 6]}, combined={grid_combined[0, 6]}")
    print(f"  Position (4, 0): bar_1={grid_bar1[4, 0]}, bar_2={grid_bar2[4, 0]}, combined={grid_combined[4, 0]}")
    print(f"  Position (4, 6): bar_1={grid_bar1[4, 6]}, bar_2={grid_bar2[4, 6]}, combined={grid_combined[4, 6]}")
    
    return layout


def example_5_validate_constraints():
    """Example 5: Detailed constraint validation"""
    print("\n" + "="*70)
    print("EXAMPLE 5: Constraint Validation")
    print("="*70)
    
    barcodes = ['AAA', 'CCC', 'GGG', 'TTT', 'ATC', 'GCA']
    layout = DualBarcodeLayout(rows=6, cols=6, barcodes=barcodes)
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=222)
    
    print("Checking constraints for each 2x2 block...\n")
    
    n_blocks_row = (layout.rows + 1) // 2
    n_blocks_col = (layout.cols + 1) // 2
    
    all_valid = True
    for block_row in range(n_blocks_row):
        for block_col in range(n_blocks_col):
            positions = []
            bar1_set = set()
            bar2_set = set()
            combined_set = set()
            
            for i in range(2):
                for j in range(2):
                    row = block_row * 2 + i
                    col = block_col * 2 + j
                    if row < layout.rows and col < layout.cols:
                        positions.append((row, col))
                        bar1_set.add(grid_bar1[row, col])
                        bar2_set.add(grid_bar2[row, col])
                        combined_set.add(grid_combined[row, col])
            
            # Validate constraints
            same_bar1 = len(bar1_set) == 1
            different_bar2 = len(bar2_set) == len(positions)
            unique_combined = len(combined_set) == len(positions)
            
            status = "✓" if (same_bar1 and different_bar2 and unique_combined) else "✗"
            
            if not (same_bar1 and different_bar2 and unique_combined):
                all_valid = False
                print(f"{status} Block ({block_row}, {block_col}): "
                      f"same_bar1={same_bar1}, different_bar2={different_bar2}, unique_combined={unique_combined}")
                print(f"    Positions: {positions}")
                print(f"    bar1: {bar1_set}")
                print(f"    bar2: {bar2_set}")
                print(f"    combined: {combined_set}")
    
    if all_valid:
        print("✓ All constraints satisfied!")
    else:
        print("\n✗ Some constraints violated!")
    
    return layout


def example_6_multiple_runs():
    """Example 6: Multiple runs with different seeds"""
    print("\n" + "="*70)
    print("EXAMPLE 6: Multiple Runs with Different Seeds")
    print("="*70)
    
    barcodes = generate_random_barcodes(n=5, length=5, seed=100)
    print(f"Using barcodes: {barcodes}\n")
    
    for seed in [100, 200, 300]:
        print(f"Run with seed={seed}:")
        layout = DualBarcodeLayout(rows=4, cols=4, barcodes=barcodes)
        grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=seed)
        
        stats = layout.get_statistics()
        print(f"  Unique combined barcodes: {stats['unique_combined_barcodes']}")
        
        # Show first row as example
        print(f"  First row combined: {grid_combined[0, :]}")
        print()


if __name__ == "__main__":
    print("="*70)
    print("DUAL BARCODE LAYOUT GENERATOR - EXAMPLE USAGE")
    print("="*70)
    
    # Run all examples
    example_1_small_grid()
    example_2_rectangular_grid()
    example_3_large_array()
    example_4_edge_cases()
    example_5_validate_constraints()
    example_6_multiple_runs()
    
    print("\n" + "="*70)
    print("ALL EXAMPLES COMPLETED")
    print("="*70)
