"""
Dual Barcode Layout Generator with Offset Block Pattern for Spatial Transcriptomics Arrays

This module implements an offset block pattern where both bar_1 and bar_2 use 2x2 blocks,
but bar_2 blocks are shifted by (1,1) to ensure unique combined barcodes at each spot.
"""

import numpy as np
import random
from typing import List, Tuple, Optional, Set
import itertools


def load_barcodes_from_file(filepath: str) -> List[str]:
    """
    Load barcodes from a text file (one barcode per line).
    
    Args:
        filepath: Path to the barcode file
    
    Returns:
        List of barcode sequences
    """
    barcodes = []
    with open(filepath, 'r') as f:
        for line in f:
            barcode = line.strip()
            if barcode:  # Skip empty lines
                barcodes.append(barcode)
    
    # Validate all barcodes have same length
    if barcodes:
        lengths = set(len(bc) for bc in barcodes)
        if len(lengths) > 1:
            raise ValueError(f"All barcodes must have the same length. Found lengths: {lengths}")
    
    return barcodes


class DualBarcodeOffsetLayout:
    """
    Generate dual barcode layouts using offset block pattern.
    
    bar_1: 2x2 blocks aligned at (0,0), (0,2), (2,0), etc.
    bar_2: 2x2 blocks offset at (1,1), (1,3), (3,1), etc.
    
    This ensures each spot has a unique bar_1+bar_2 combination.
    
    Attributes:
        rows (int): Number of rows in the grid
        cols (int): Number of columns in the grid
        barcodes (List[str]): List of available barcodes (used for both bar_1 and bar_2)
        grid_bar1 (np.ndarray): Grid storing bar_1 assignments
        grid_bar2 (np.ndarray): Grid storing bar_2 assignments
        grid_combined (np.ndarray): Grid storing combined bar_1+bar_2 barcodes
    """
    
    def __init__(self, rows: int, cols: int, barcodes: List[str]):
        """
        Initialize the dual barcode layout generator with offset pattern.
        
        Args:
            rows: Number of rows in the grid
            cols: Number of columns in the grid
            barcodes: List of barcode sequences
        
        Raises:
            ValueError: If insufficient barcodes provided
        """
        self.rows = rows
        self.cols = cols
        self.barcodes = barcodes
        self.n_barcodes = len(barcodes)
        
        # Initialize grids
        self.grid_bar1 = np.empty((rows, cols), dtype=object)
        self.grid_bar2 = np.empty((rows, cols), dtype=object)
        self.grid_combined = np.empty((rows, cols), dtype=object)
        
        # Calculate number of blocks
        self.n_blocks_bar1_row = (rows + 1) // 2
        self.n_blocks_bar1_col = (cols + 1) // 2
        self.total_blocks_bar1 = self.n_blocks_bar1_row * self.n_blocks_bar1_col
        
        # For bar_2 offset blocks starting at (1,1)
        self.n_blocks_bar2_row = rows // 2  # How many complete 2x2 blocks fit starting at row 1
        self.n_blocks_bar2_col = cols // 2  # How many complete 2x2 blocks fit starting at col 1
        self.total_blocks_bar2 = self.n_blocks_bar2_row * self.n_blocks_bar2_col
        
    def generate_layout(self, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate the complete dual barcode layout with offset pattern.
        
        Args:
            seed: Random seed for reproducibility
        
        Returns:
            Tuple of (grid_bar1, grid_bar2, grid_combined)
        """
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        
        # Step 1: Assign bar_1 to standard 2x2 blocks
        self._assign_bar1()
        
        # Step 2: Assign bar_2 to offset 2x2 blocks
        self._assign_bar2_offset()
        
        # Step 3: Create combined barcodes
        self._create_combined_barcodes()
        
        # Step 4: Validate the layout
        self._validate_layout()
        
        return self.grid_bar1, self.grid_bar2, self.grid_combined
    
    def _assign_bar1(self):
        """
        Assign bar_1 barcodes to standard 2x2 blocks starting at (0,0).
        """
        shuffled_barcodes = self.barcodes.copy()
        random.shuffle(shuffled_barcodes)
        barcode_cycle = itertools.cycle(shuffled_barcodes)
        
        for block_row in range(self.n_blocks_bar1_row):
            for block_col in range(self.n_blocks_bar1_col):
                bar1 = next(barcode_cycle)
                
                # Assign to all positions in this 2x2 block
                for i in range(2):
                    for j in range(2):
                        row = block_row * 2 + i
                        col = block_col * 2 + j
                        
                        if row < self.rows and col < self.cols:
                            self.grid_bar1[row, col] = bar1
    
    def _assign_bar2_offset(self):
        """
        Assign bar_2 barcodes to offset 2x2 blocks starting at (1,1).
        Uses backtracking to ensure no bar_1+bar_2 combination repeats.
        """
        # Track used combinations
        used_combinations = set()
        
        # Collect all blocks to fill
        blocks = []
        for block_row in range(self.n_blocks_bar2_row):
            for block_col in range(self.n_blocks_bar2_col):
                blocks.append((block_row, block_col))
        
        # Shuffle block order for randomization
        random.shuffle(blocks)
        
        # Try to assign bar_2 to each block
        for block_row, block_col in blocks:
            # Get positions in this offset block
            positions = []
            bar1_in_block = set()
            
            for i in range(2):
                for j in range(2):
                    row = 1 + block_row * 2 + i
                    col = 1 + block_col * 2 + j
                    
                    if row < self.rows and col < self.cols:
                        positions.append((row, col))
                        bar1_in_block.add(self.grid_bar1[row, col])
            
            # Find a bar_2 that doesn't create duplicate combinations with any bar_1 in this block
            available_bar2 = self.barcodes.copy()
            random.shuffle(available_bar2)
            
            bar2_assigned = None
            for bar2 in available_bar2:
                # Check if this bar_2 creates conflicts with any bar_1 in the block
                conflict = False
                for bar1 in bar1_in_block:
                    if (bar1, bar2) in used_combinations:
                        conflict = True
                        break
                
                if not conflict:
                    bar2_assigned = bar2
                    # Mark all bar_1+bar_2 combinations in this block as used
                    for bar1 in bar1_in_block:
                        used_combinations.add((bar1, bar2))
                    break
            
            if bar2_assigned is None:
                raise ValueError(
                    f"Cannot find valid bar_2 for offset block ({block_row}, {block_col}). "
                    f"bar_1 values in block: {bar1_in_block}. "
                    f"Try using more barcodes."
                )
            
            # Assign bar_2 to all positions in this block
            for row, col in positions:
                self.grid_bar2[row, col] = bar2_assigned
        
        # Fill edge positions (row 0, col 0, last row, last col) that aren't in offset blocks
        self._fill_edge_positions(used_combinations)
    
    def _fill_edge_positions(self, used_combinations: Set[Tuple[str, str]]):
        """
        Fill positions not covered by offset bar_2 blocks (edges).
        """
        # Find all positions that don't have bar_2 assigned yet
        edge_positions = []
        for i in range(self.rows):
            for j in range(self.cols):
                if self.grid_bar2[i, j] is None:
                    edge_positions.append((i, j))
        
        # Assign bar_2 to each edge position
        for row, col in edge_positions:
            bar1 = self.grid_bar1[row, col]
            
            # Find available bar_2 that doesn't create duplicate combination
            available_bar2 = [bc for bc in self.barcodes 
                             if (bar1, bc) not in used_combinations]
            
            if not available_bar2:
                raise ValueError(
                    f"Cannot find valid bar_2 for edge position ({row}, {col}) "
                    f"with bar_1='{bar1}'. Try using more barcodes."
                )
            
            bar2 = random.choice(available_bar2)
            self.grid_bar2[row, col] = bar2
            used_combinations.add((bar1, bar2))
    
    def _create_combined_barcodes(self):
        """
        Create combined barcodes by concatenating bar_1 and bar_2.
        """
        for i in range(self.rows):
            for j in range(self.cols):
                if self.grid_bar1[i, j] is not None and self.grid_bar2[i, j] is not None:
                    self.grid_combined[i, j] = self.grid_bar1[i, j] + self.grid_bar2[i, j]
    
    def _validate_layout(self) -> bool:
        """
        Validate that all combined barcodes are globally unique.
        
        Returns:
            True if valid, raises ValueError if invalid
        """
        all_combined = {}
        
        for i in range(self.rows):
            for j in range(self.cols):
                combined = self.grid_combined[i, j]
                
                if combined in all_combined:
                    raise ValueError(
                        f"Duplicate combined barcode '{combined}' found at "
                        f"position ({i}, {j}) and {all_combined[combined]}"
                    )
                all_combined[combined] = (i, j)
        
        return True
    
    def get_statistics(self) -> dict:
        """
        Get statistics about the generated layout.
        
        Returns:
            Dictionary with layout statistics
        """
        unique_combined = len(set(self.grid_combined.flatten()))
        total_spots = self.rows * self.cols
        
        # Count bar_1 usage
        bar1_counts = {}
        for barcode in self.grid_bar1.flatten():
            if barcode is not None:
                bar1_counts[barcode] = bar1_counts.get(barcode, 0) + 1
        
        # Count bar_2 usage
        bar2_counts = {}
        for barcode in self.grid_bar2.flatten():
            if barcode is not None:
                bar2_counts[barcode] = bar2_counts.get(barcode, 0) + 1
        
        # Count how many bar_1 blocks each bar_2 appears in
        bar2_block_diversity = {}
        for block_row in range(self.n_blocks_bar2_row):
            for block_col in range(self.n_blocks_bar2_col):
                row = 1 + block_row * 2
                col = 1 + block_col * 2
                if row < self.rows and col < self.cols:
                    bar2 = self.grid_bar2[row, col]
                    if bar2 is not None:
                        bar2_block_diversity[bar2] = bar2_block_diversity.get(bar2, 0) + 1
        
        return {
            'grid_dimensions': (self.rows, self.cols),
            'total_spots': total_spots,
            'n_barcodes_available': self.n_barcodes,
            'n_blocks_bar1': self.total_blocks_bar1,
            'n_blocks_bar2': self.total_blocks_bar2,
            'unique_combined_barcodes': unique_combined,
            'bar1_usage': bar1_counts,
            'bar2_usage': bar2_counts,
            'bar2_block_diversity': bar2_block_diversity,
            'theoretical_max_combinations': self.n_barcodes ** 2
        }
    
    def export_layout(self, filepath: str, format: str = 'csv'):
        """
        Export the layout to a file in column-major order.
        
        Args:
            filepath: Path to save the file
            format: Export format ('csv', 'tsv', or 'npy')
        """
        if format == 'csv':
            self._export_csv(filepath, delimiter=',')
        elif format == 'tsv':
            self._export_csv(filepath, delimiter='\t')
        elif format == 'npy':
            np.save(filepath, {
                'bar1': self.grid_bar1,
                'bar2': self.grid_bar2,
                'combined': self.grid_combined
            })
        else:
            raise ValueError(f"Unknown format: {format}")
    
    def _export_csv(self, filepath: str, delimiter: str = ','):
        """Export layout to CSV/TSV format in column-major order."""
        with open(filepath, 'w') as f:
            # Header
            f.write(f"row{delimiter}col{delimiter}bar_1{delimiter}bar_2{delimiter}combined\n")
            
            # Data rows in column-major order
            for j in range(self.cols):
                for i in range(self.rows):
                    f.write(f"{i}{delimiter}{j}{delimiter}")
                    f.write(f"{self.grid_bar1[i, j]}{delimiter}")
                    f.write(f"{self.grid_bar2[i, j]}{delimiter}")
                    f.write(f"{self.grid_combined[i, j]}\n")
    
    def visualize_layout(self, output_file: Optional[str] = None, include_combined: bool = True):
        """
        Create a visualization of the barcode layout showing block boundaries.
        
        Args:
            output_file: Optional path to save the figure
            include_combined: If True, show all 3 plots; if False, only bar_1 and bar_2
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
        except ImportError:
            print("Matplotlib not available. Install with: pip install matplotlib")
            return
        
        # Create color maps
        unique_bar1 = list(set(self.grid_bar1.flatten()))
        unique_bar2 = list(set(self.grid_bar2.flatten()))
        
        bar1_colors = {bc: i for i, bc in enumerate(unique_bar1)}
        bar2_colors = {bc: i for i, bc in enumerate(unique_bar2)}
        
        if include_combined:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            ax_bar1, ax_bar2, ax_combined = axes[0], axes[1], axes[2]
        else:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            ax_bar1, ax_bar2 = axes[0], axes[1]
        
        # Plot bar_1 with standard block boundaries
        bar1_numeric = np.array([[bar1_colors[self.grid_bar1[i, j]] 
                                  for j in range(self.cols)] 
                                 for i in range(self.rows)])
        ax_bar1.imshow(bar1_numeric, cmap='tab20', aspect='auto')
        ax_bar1.set_title('bar_1 Layout\n(2x2 blocks at (0,0))', fontsize=12, fontweight='bold')
        ax_bar1.set_xlabel('Column')
        ax_bar1.set_ylabel('Row')
        
        # Draw bar_1 block boundaries
        for i in range(0, self.rows, 2):
            ax_bar1.axhline(i - 0.5, color='black', linewidth=2)
        for j in range(0, self.cols, 2):
            ax_bar1.axvline(j - 0.5, color='black', linewidth=2)
        
        # Plot bar_2 with offset block boundaries
        bar2_numeric = np.array([[bar2_colors[self.grid_bar2[i, j]] 
                                  for j in range(self.cols)] 
                                 for i in range(self.rows)])
        ax_bar2.imshow(bar2_numeric, cmap='tab20', aspect='auto')
        ax_bar2.set_title('bar_2 Layout\n(2x2 blocks offset at (1,1))', fontsize=12, fontweight='bold')
        ax_bar2.set_xlabel('Column')
        ax_bar2.set_ylabel('Row')
        
        # Draw bar_2 offset block boundaries
        for i in range(1, self.rows, 2):
            ax_bar2.axhline(i - 0.5, color='red', linewidth=2)
        for j in range(1, self.cols, 2):
            ax_bar2.axvline(j - 0.5, color='red', linewidth=2)
        
        # Plot combined (if requested)
        if include_combined:
            unique_combined = list(set(self.grid_combined.flatten()))
            combined_colors = {bc: i for i, bc in enumerate(unique_combined)}
            
            combined_numeric = np.array([[combined_colors[self.grid_combined[i, j]] 
                                         for j in range(self.cols)] 
                                        for i in range(self.rows)])
            ax_combined.imshow(combined_numeric, cmap='tab20', aspect='auto')
            ax_combined.set_title('Combined Layout (bar_1 + bar_2)\n(All unique)', 
                                 fontsize=12, fontweight='bold')
            ax_combined.set_xlabel('Column')
            ax_combined.set_ylabel('Row')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Visualization saved to {output_file}")
        else:
            plt.show()


if __name__ == "__main__":
    print("Use run_offset_layout.py to generate layouts with custom parameters")
