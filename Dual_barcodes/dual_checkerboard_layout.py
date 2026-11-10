#!/usr/bin/env python3
"""
Dual Barcode Layout Generator for Checkerboard Arrays

This module implements a block-wise algorithm for checkerboard patterns where:
- The array consists of two interleaved sub-arrays (even/odd positions)
- Each 8x8 checkerboard block shares the same bar_1 barcode
- Each position within a block has a unique bar_2 barcode
- Final barcode = bar_1 + bar_2 (36nt = 18nt + 18nt)

Checkerboard Pattern:
- Row 0 starts at column 0 (even columns): (0,0), (0,2), (0,4), ...
- Row 1 starts at column 1 (odd columns): (1,1), (1,3), (1,5), ...
- This pattern repeats, creating two interleaved arrays
"""

import random
from typing import List, Tuple, Optional, Dict


class CheckerboardDualBarcodeLayout:
    """
    Generate dual barcode layouts for checkerboard spatial transcriptomics arrays.
    
    The checkerboard pattern consists of two interleaved arrays:
    - Array 1 (even rows): rows 0,2,4,... starting at column 0
    - Array 2 (odd rows): rows 1,3,5,... starting at column 1
    """
    
    def __init__(self, rows: int, cols: int, barcodes_set1: List[str], 
                 barcodes_set2: List[str], block_size: int = 8, split_barcodes: bool = True):
        """
        Initialize the checkerboard dual barcode layout generator.
        
        Args:
            rows: Number of rows in the grid (total, including both arrays)
            cols: Number of columns in the grid (total)
            barcodes_set1: List of barcodes for the first interleaved array (even rows)
            barcodes_set2: List of barcodes for the second interleaved array (odd rows)
            block_size: Size of blocks (default 8x8 in checkerboard pattern = 32 spots)
            split_barcodes: If True, split 36nt barcodes into 18nt halves (default True)
        
        Raises:
            ValueError: If insufficient barcodes provided
        """
        self.rows = rows
        self.cols = cols
        self.block_size = block_size
        
        # Split 36nt barcodes into 18nt halves if needed
        if split_barcodes:
            print(f"  Splitting {len(barcodes_set1)} x 36nt barcodes into 18nt halves...")
            self.barcodes_set1 = self._split_barcodes(barcodes_set1)
            self.barcodes_set2 = self._split_barcodes(barcodes_set2)
            print(f"  Result: {len(self.barcodes_set1)} unique 18nt barcodes from set 1")
            print(f"  Result: {len(self.barcodes_set2)} unique 18nt barcodes from set 2")
        else:
            self.barcodes_set1 = barcodes_set1
            self.barcodes_set2 = barcodes_set2
        
        # Calculate actual spot positions in checkerboard
        self.positions = self._generate_checkerboard_positions()
        self.total_spots = len(self.positions)
        
        # Separate positions into two interleaved arrays
        self.array1_positions = [(r, c) for r, c in self.positions if r % 2 == 0]
        self.array2_positions = [(r, c) for r, c in self.positions if r % 2 == 1]
        
        print(f"Checkerboard layout: {self.total_spots} total spots")
        print(f"  Array 1 (even rows): {len(self.array1_positions)} spots")
        print(f"  Array 2 (odd rows): {len(self.array2_positions)} spots")
        
        # Initialize dictionaries to store barcode assignments
        # Key: (row, col), Value: barcode
        self.grid_bar1 = {}
        self.grid_bar2 = {}
        self.grid_combined = {}
        
        # Validate barcode counts
        min_barcodes_needed = max(
            self._estimate_blocks_needed(self.array1_positions),
            self._estimate_blocks_needed(self.array2_positions)
        )
        
        if len(self.barcodes_set1) < min_barcodes_needed or len(self.barcodes_set2) < min_barcodes_needed:
            raise ValueError(
                f"Insufficient barcodes: need at least {min_barcodes_needed} "
                f"in each set, but got {len(self.barcodes_set1)} and {len(self.barcodes_set2)}"
            )
    
    def _split_barcodes(self, barcodes: List[str]) -> List[str]:
        """
        Split 36nt barcodes into 18nt halves and return unique 18nt sequences.
        
        Args:
            barcodes: List of 36nt barcodes
        
        Returns:
            List of unique 18nt barcode segments
        """
        segments = set()
        
        for bc in barcodes:
            if len(bc) != 36:
                # If not 36nt, assume it's already the right length
                segments.add(bc)
            else:
                # Split into two 18nt halves
                first_half = bc[:18]
                second_half = bc[18:]
                segments.add(first_half)
                segments.add(second_half)
        
        return list(segments)
    
    def _generate_checkerboard_positions(self) -> List[Tuple[int, int]]:
        """
        Generate all valid positions in the checkerboard pattern.
        
        Pattern:
        - Even rows (0,2,4,...): columns 0,2,4,... 
        - Odd rows (1,3,5,...): columns 1,3,5,...
        
        Returns:
            List of (row, col) tuples
        """
        positions = []
        
        for row in range(self.rows):
            if row % 2 == 0:
                # Even row: starts at column 0, uses even columns
                start_col = 0
                for col in range(start_col, self.cols, 2):
                    positions.append((row, col))
            else:
                # Odd row: starts at column 1, uses odd columns
                start_col = 1
                for col in range(start_col, self.cols, 2):
                    positions.append((row, col))
        
        return sorted(positions)
    
    def _estimate_blocks_needed(self, positions: List[Tuple[int, int]]) -> int:
        """
        Estimate number of blocks needed for a set of positions.
        
        For 8x8 checkerboard blocks, each block contains ~32 spots
        (half of an 8x8 grid due to checkerboard pattern).
        """
        # Conservative estimate: each block has at least 16 spots
        min_spots_per_block = 16
        return (len(positions) + min_spots_per_block - 1) // min_spots_per_block
    
    def _assign_to_checkerboard_blocks(self, positions: List[Tuple[int, int]]) -> Dict:
        """
        Assign positions to 8x8 checkerboard blocks.
        
        Args:
            positions: List of (row, col) positions
        
        Returns:
            Dictionary mapping block_id to list of positions
        """
        blocks = {}
        
        for row, col in positions:
            # Determine which 8x8 block this position belongs to
            # Block coordinates based on 8x8 grid
            block_row = row // self.block_size
            block_col = col // self.block_size
            block_id = (block_row, block_col)
            
            if block_id not in blocks:
                blocks[block_id] = []
            
            blocks[block_id].append((row, col))
        
        return blocks
    
    def generate_layout(self, seed: Optional[int] = None) -> Tuple[Dict, Dict, Dict]:
        """
        Generate the complete dual barcode layout for checkerboard pattern.
        
        Uses unified 8x8 spatial blocks where all spots (both arrays) within
        the same spatial region share the same bar_1.
        
        Args:
            seed: Random seed for reproducibility
        
        Returns:
            Tuple of (grid_bar1, grid_bar2, grid_combined) as dictionaries
        """
        if seed is not None:
            random.seed(seed)
        
        print("\n[1/3] Assigning barcodes using unified spatial blocks...")
        self._assign_unified_spatial_blocks()
        
        print("\n[2/3] Creating combined barcodes...")
        self._create_combined_barcodes()
        
        print("\n[3/3] Validating layout...")
        self._validate_layout()
        
        return self.grid_bar1, self.grid_bar2, self.grid_combined
    
    def _assign_unified_spatial_blocks(self):
        """
        Assign barcodes using unified 8x8 spatial blocks.
        
        All spots within the same 8x8 spatial region (both arrays) share the same bar_1,
        but each array uses its own barcode set.
        """
        # Group ALL positions into 8x8 spatial blocks
        all_blocks = {}
        for row, col in self.positions:
            block_row = row // self.block_size
            block_col = col // self.block_size
            block_id = (block_row, block_col)
            
            if block_id not in all_blocks:
                all_blocks[block_id] = {'array1': [], 'array2': []}
            
            # Separate by array
            if (row, col) in self.array1_positions:
                all_blocks[block_id]['array1'].append((row, col))
            else:
                all_blocks[block_id]['array2'].append((row, col))
        
        print(f"  Total unified spatial blocks: {len(all_blocks)}")
        print(f"  Total spots: {len(self.positions)}")
        
        # Combine both barcode sets (each block will use one bar_1 from either set)
        all_barcodes = list(set(self.barcodes_set1 + self.barcodes_set2))
        
        # Track barcode usage globally
        bar1_to_used_bar2 = {bc: set() for bc in all_barcodes}
        shuffled_barcodes = all_barcodes.copy()
        random.shuffle(shuffled_barcodes)
        
        barcode_usage = {bc: 0 for bc in all_barcodes}
        max_uses_per_barcode = len(all_barcodes)
        
        # Process each unified spatial block
        for block_id, arrays in sorted(all_blocks.items()):
            array1_positions = arrays['array1']
            array2_positions = arrays['array2']
            all_positions_in_block = array1_positions + array2_positions
            
            if not all_positions_in_block:
                continue
            
            # Assign ONE bar_1 for the entire spatial block
            available_bar1 = [bc for bc in shuffled_barcodes 
                             if barcode_usage[bc] < max_uses_per_barcode]
            
            if not available_bar1:
                raise ValueError(f"No available bar_1 barcodes for block {block_id}")
            
            bar1 = min(available_bar1, key=lambda bc: barcode_usage[bc])
            barcode_usage[bar1] += 1
            
            # Assign bar_2 for each position in this block
            n_positions = len(all_positions_in_block)
            used_bar2_for_this_bar1 = bar1_to_used_bar2[bar1]
            
            # Get all available bar_2 from both sets
            available_bar2 = [bc for bc in all_barcodes 
                             if bc not in used_bar2_for_this_bar1]
            
            if len(available_bar2) < n_positions:
                raise ValueError(
                    f"Cannot find enough unique bar_2 barcodes for block {block_id}. "
                    f"Need {n_positions}, but only {len(available_bar2)} available."
                )
            
            # Shuffle and select unique bar_2 for this block
            random.shuffle(available_bar2)
            selected_bar2 = available_bar2[:n_positions]
            
            # Mark bar_2 as used with this bar_1
            for bar2 in selected_bar2:
                used_bar2_for_this_bar1.add(bar2)
            
            # Assign to positions
            for pos, bar2 in zip(all_positions_in_block, selected_bar2):
                self.grid_bar1[pos] = bar1
                self.grid_bar2[pos] = bar2
    
    def _create_combined_barcodes(self):
        """Create combined barcodes by concatenating bar_1 and bar_2."""
        for pos in self.positions:
            self.grid_combined[pos] = self.grid_bar1[pos] + self.grid_bar2[pos]
        
        print(f"  Created {len(self.grid_combined)} combined barcodes")
    
    def _validate_layout(self) -> bool:
        """
        Validate that the layout meets all constraints.
        
        Returns:
            True if valid, raises ValueError if invalid
        """
        # Check all positions have barcodes
        missing = [pos for pos in self.positions 
                  if pos not in self.grid_combined]
        if missing:
            raise ValueError(f"Missing barcodes for {len(missing)} positions")
        
        # Check global uniqueness
        all_combined = list(self.grid_combined.values())
        unique_combined = set(all_combined)
        
        if len(unique_combined) != len(all_combined):
            duplicates = len(all_combined) - len(unique_combined)
            raise ValueError(f"Found {duplicates} duplicate combined barcodes")
        
        # Validate unified spatial blocks
        # Group ALL positions into 8x8 spatial blocks
        all_blocks = {}
        for row, col in self.positions:
            block_row = row // self.block_size
            block_col = col // self.block_size
            block_id = (block_row, block_col)
            
            if block_id not in all_blocks:
                all_blocks[block_id] = []
            all_blocks[block_id].append((row, col))
        
        for block_id, block_positions in all_blocks.items():
            block_bar1 = set(self.grid_bar1[pos] for pos in block_positions)
            block_bar2 = set(self.grid_bar2[pos] for pos in block_positions)
            
            # All positions in spatial block should share bar_1
            if len(block_bar1) != 1:
                raise ValueError(
                    f"Spatial block {block_id} has inconsistent bar_1 assignments: {block_bar1}"
                )
            
            # All positions in spatial block should have different bar_2
            if len(block_bar2) != len(block_positions):
                raise ValueError(
                    f"Spatial block {block_id} has duplicate bar_2 assignments"
                )
        
        print(f"  ✓ Validation passed")
        print(f"  ✓ All {len(unique_combined)} combined barcodes are unique")
        
        return True
    
    def get_statistics(self) -> dict:
        """
        Get statistics about the generated layout.
        
        Returns:
            Dictionary with layout statistics
        """
        unique_combined = len(set(self.grid_combined.values()))
        
        # Count bar_1 usage
        bar1_counts = {}
        for barcode in self.grid_bar1.values():
            bar1_counts[barcode] = bar1_counts.get(barcode, 0) + 1
        
        # Count bar_2 usage
        bar2_counts = {}
        for barcode in self.grid_bar2.values():
            bar2_counts[barcode] = bar2_counts.get(barcode, 0) + 1
        
        return {
            'grid_dimensions': (self.rows, self.cols),
            'total_spots': self.total_spots,
            'array1_spots': len(self.array1_positions),
            'array2_spots': len(self.array2_positions),
            'block_size': f'{self.block_size}x{self.block_size} checkerboard',
            'n_barcodes_set1': len(self.barcodes_set1),
            'n_barcodes_set2': len(self.barcodes_set2),
            'unique_combined_barcodes': unique_combined,
            'unique_bar1_used': len(bar1_counts),
            'unique_bar2_used': len(bar2_counts),
            'bar1_usage': bar1_counts,
            'bar2_usage': bar2_counts,
        }
    
    def export_layout(self, filepath: str, format: str = 'csv', column_major: bool = True):
        """
        Export the layout to a file.
        
        Args:
            filepath: Path to save the file
            format: Export format ('csv' or 'tsv')
            column_major: If True, export in column-major order (for photolithography)
        """
        if format == 'csv':
            self._export_csv(filepath, delimiter=',', column_major=column_major)
        elif format == 'tsv':
            self._export_csv(filepath, delimiter='\t', column_major=column_major)
        else:
            raise ValueError(f"Unknown format: {format}")
    
    def _export_csv(self, filepath: str, delimiter: str = ',', column_major: bool = True):
        """Export layout to CSV/TSV format."""
        with open(filepath, 'w') as f:
            # Header
            f.write(f"row{delimiter}col{delimiter}bar_1{delimiter}bar_2{delimiter}barcode\n")
            
            # Sort positions
            if column_major:
                # Column-major: iterate columns first, then rows
                sorted_positions = sorted(self.positions, key=lambda p: (p[1], p[0]))
            else:
                # Row-major: iterate rows first, then columns
                sorted_positions = sorted(self.positions, key=lambda p: (p[0], p[1]))
            
            # Data rows
            for (row, col) in sorted_positions:
                f.write(f"{row}{delimiter}{col}{delimiter}")
                f.write(f"{self.grid_bar1[(row, col)]}{delimiter}")
                f.write(f"{self.grid_bar2[(row, col)]}{delimiter}")
                f.write(f"{self.grid_combined[(row, col)]}\n")
        
        print(f"✓ Layout exported to {filepath}")
        print(f"  Format: {'column-major' if column_major else 'row-major'}")
        print(f"  Total rows: {len(self.positions) + 1} (including header)")


def load_barcodes_from_file(filepath: str) -> List[str]:
    """
    Load barcodes from a text file (one barcode per line).
    
    Args:
        filepath: Path to the text file containing barcodes
    
    Returns:
        List of barcode sequences
    """
    try:
        with open(filepath, 'r') as f:
            barcodes = [line.strip() for line in f if line.strip()]
        
        if not barcodes:
            raise ValueError(f"No barcodes found in file: {filepath}")
        
        print(f"Loaded {len(barcodes)} barcodes from {filepath}")
        return barcodes
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Barcode file not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error reading barcode file: {e}")


if __name__ == "__main__":
    import sys
    
    print("=" * 80)
    print("CHECKERBOARD DUAL BARCODE LAYOUT GENERATOR")
    print("=" * 80)
    print("\nThis script generates dual barcode layouts for checkerboard arrays.")
    print("Usage: python dual_checkerboard_layout.py <barcodes_set1.txt> <barcodes_set2.txt>")
    print("\nExample: python dual_checkerboard_layout.py barcodes_118x182_dual_set1.txt barcodes_118x182_dual_set2.txt")
    print("=" * 80)
