#!/usr/bin/env python3
"""
Get Layout Bounds

Simple script to extract min/max row and column coordinates from a layout CSV file.

Usage:
    python get_layout_bounds.py <layout.csv>
    
Example:
    python get_layout_bounds.py layout_235x363_checkerboard_colmajor.csv
"""

import sys
import csv


def get_layout_bounds(csv_file):
    """
    Read layout CSV and return min/max coordinates.
    
    Args:
        csv_file: Path to layout CSV file
    
    Returns:
        Dictionary with x_min, x_max, y_min, y_max (where x=row, y=col)
    """
    rows = []
    cols = []
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for line_num, row in enumerate(reader, start=2):  # Start at 2 (line 1 is header)
            try:
                row_val = int(row['row'])
                col_val = int(row['col'])
                rows.append(row_val)
                cols.append(col_val)
            except (KeyError, ValueError) as e:
                print(f"Warning: Skipping line {line_num}: {e}")
                continue
    
    if not rows or not cols:
        raise ValueError("No valid data found in CSV file")
    
    return {
        'x_min': min(rows),
        'x_max': max(rows),
        'y_min': min(cols),
        'y_max': max(cols),
        'total_spots': len(rows)
    }


def main():
    """Main execution"""
    if len(sys.argv) < 2:
        print("Usage: python get_layout_bounds.py <layout.csv>")
        print("\nExample:")
        print("  python get_layout_bounds.py layout_235x363_checkerboard_colmajor.csv")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    
    print("=" * 60)
    print("LAYOUT BOUNDS ANALYZER")
    print("=" * 60)
    print(f"\nReading: {csv_file}")
    
    try:
        bounds = get_layout_bounds(csv_file)
        
        print(f"\nResults:")
        print(f"  x_min (row min): {bounds['x_min']}")
        print(f"  x_max (row max): {bounds['x_max']}")
        print(f"  y_min (col min): {bounds['y_min']}")
        print(f"  y_max (col max): {bounds['y_max']}")
        print(f"  Total spots:     {bounds['total_spots']:,}")
        
        print(f"\nArray dimensions:")
        print(f"  Rows (x): {bounds['x_max'] - bounds['x_min'] + 1}")
        print(f"  Cols (y): {bounds['y_max'] - bounds['y_min'] + 1}")
        
        print("\n" + "=" * 60)
        
    except FileNotFoundError:
        print(f"\nError: File not found: {csv_file}")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()


# Example Python API usage:
"""
from get_layout_bounds import get_layout_bounds

# Get bounds
bounds = get_layout_bounds('layout_235x363_checkerboard_colmajor.csv')

print(f"Row range: {bounds['x_min']} to {bounds['x_max']}")
print(f"Col range: {bounds['y_min']} to {bounds['y_max']}")
print(f"Total spots: {bounds['total_spots']}")
"""
