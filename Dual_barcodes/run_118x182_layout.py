"""
Generate 4x4 dual barcode layout for 118x182 spot array
"""

from dual_barcode_layout import DualBarcodeLayout, load_barcodes_from_file
import sys

# Configuration
ROWS = 118
COLS = 182
TOTAL_SPOTS = ROWS * COLS
BARCODE_FILE = 'barcodes_1in4_half.txt'
OUTPUT_FILE = 'layout_118x182_4x4.csv'
SEED = 42

# Calculate requirements
N_BLOCKS_ROW = (ROWS + 3) // 4  # 30 blocks
N_BLOCKS_COL = (COLS + 3) // 4  # 46 blocks  
N_BLOCKS = N_BLOCKS_ROW * N_BLOCKS_COL  # 1,380 blocks
MIN_BARCODES = N_BLOCKS  # Need at least 1,380 barcodes

print("="*70)
print("DUAL BARCODE LAYOUT GENERATOR - 118x182 Array (4x4 Clusters)")
print("="*70)

print(f"\nArray Configuration:")
print(f"  Grid size: {ROWS} × {COLS} = {TOTAL_SPOTS:,} spots")
print(f"  Block size: 4 × 4 = 16 spots per block")
print(f"  Number of blocks: {N_BLOCKS_ROW} × {N_BLOCKS_COL} = {N_BLOCKS:,} blocks")
print(f"  Minimum barcodes needed: {MIN_BARCODES:,}")

# Load barcodes
print(f"\nLoading barcodes from: {BARCODE_FILE}")
try:
    barcodes = load_barcodes_from_file(BARCODE_FILE)
    print(f"✓ Loaded {len(barcodes):,} barcodes (length: {len(barcodes[0])} bases)")
    
    if len(barcodes) < MIN_BARCODES:
        print(f"\n✗ ERROR: Not enough barcodes!")
        print(f"  Need: {MIN_BARCODES:,}, Have: {len(barcodes):,}")
        print(f"  Missing: {MIN_BARCODES - len(barcodes):,} barcodes")
        sys.exit(1)
    
    print(f"✓ Sufficient barcodes available (using {MIN_BARCODES:,} of {len(barcodes):,})")
    
except FileNotFoundError:
    print(f"✗ ERROR: Barcode file '{BARCODE_FILE}' not found")
    print(f"  Please ensure the file is in the current directory")
    sys.exit(1)

# Create layout
print(f"\nGenerating {ROWS}×{COLS} layout...")
print(f"  This may take a moment...")

try:
    layout = DualBarcodeLayout(ROWS, COLS, barcodes[:MIN_BARCODES])
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=SEED)
    
    print("✓ Layout generated successfully!")
    
    # Get statistics
    stats = layout.get_statistics()
    print("\n" + "="*70)
    print("LAYOUT STATISTICS")
    print("="*70)
    print(f"Grid dimensions: {stats['grid_dimensions']}")
    print(f"Total spots: {stats['total_spots']:,}")
    print(f"Number of barcodes used: {stats['n_barcodes_available']:,}")
    print(f"Number of 4x4 blocks: {stats['n_blocks']:,}")
    print(f"Unique combined barcodes: {stats['unique_combined_barcodes']:,}")
    print(f"Coverage: {stats['unique_combined_barcodes']}/{stats['total_spots']} "
          f"({100*stats['unique_combined_barcodes']/stats['total_spots']:.1f}%)")
    
    # Export layout
    print("\n" + "="*70)
    print("EXPORTING LAYOUT")
    print("="*70)
    layout.export_layout(OUTPUT_FILE, format='csv')
    print(f"✓ Layout exported to: {OUTPUT_FILE}")
    print(f"  Format: CSV (column-major order)")
    print(f"  Columns: row, col, bar_1, bar_2, combined")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"✓ Successfully generated {ROWS}×{COLS} layout")
    print(f"✓ Total unique barcodes: {stats['unique_combined_barcodes']:,}")
    print(f"✓ Output file: {OUTPUT_FILE}")
    
    print("\nLayout Properties:")
    print(f"  • Each 4×4 block shares the same bar_1 barcode")
    print(f"  • Each spot in a block has a unique bar_2 barcode")
    print(f"  • Combined barcode = bar_1 + bar_2 (unique for each spot)")
    print(f"  • Block boundaries: every 4 rows and 4 columns")
    
    print("\n" + "="*70)
    print("DONE!")
    print("="*70)
    
except Exception as e:
    print(f"\n✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
