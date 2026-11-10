"""
Run dual barcode layout for large 118x182 grid
"""

from dual_barcode_layout import DualBarcodeLayout, load_barcodes_from_file

# Grid dimensions
ROWS = 118
COLS = 182

# Path to your barcode file (one barcode per line)
BARCODE_FILE = "barcodes_1in4_half.txt"

# Output files
OUTPUT_CSV = "layout_118x182.csv"
OUTPUT_VIZ = "layout_118x182.png"

print("=" * 70)
print(f"Dual Barcode Layout Generator: {ROWS}x{COLS} Grid")
print("=" * 70)

# Calculate requirements
n_blocks = ((ROWS + 1) // 2) * ((COLS + 1) // 2)
n_spots = ROWS * COLS
min_barcodes_needed = n_blocks

print(f"\nGrid requirements:")
print(f"  Total spots: {n_spots:,}")
print(f"  Number of 2x2 blocks: {n_blocks:,}")
print(f"  Minimum barcodes needed: {min_barcodes_needed:,}")
print(f"  (Each barcode can be used as bar_1 in exactly 1 block)")

# Load barcodes
print(f"\nLoading barcodes from: {BARCODE_FILE}")
try:
    barcodes = load_barcodes_from_file(BARCODE_FILE)
    print(f"✓ Loaded {len(barcodes):,} barcodes")
    
    if barcodes:
        print(f"  Barcode length: {len(barcodes[0])} nucleotides")
        print(f"  First barcode: {barcodes[0]}")
        print(f"  Last barcode: {barcodes[-1]}")
    
    # Check if we have enough barcodes
    if len(barcodes) < min_barcodes_needed:
        print(f"\n✗ ERROR: Not enough barcodes!")
        print(f"  Need: {min_barcodes_needed:,}")
        print(f"  Have: {len(barcodes):,}")
        print(f"  Missing: {min_barcodes_needed - len(barcodes):,}")
        print(f"\nPlease provide a file with at least {min_barcodes_needed:,} barcodes.")
        exit(1)
    
    print(f"✓ Sufficient barcodes available")
    
except FileNotFoundError:
    print(f"✗ ERROR: File not found: {BARCODE_FILE}")
    print(f"\nPlease create a text file with one barcode per line.")
    print(f"You need at least {min_barcodes_needed:,} barcodes for this grid size.")
    exit(1)

# Create layout
print(f"\nCreating {ROWS}x{COLS} layout...")
layout = DualBarcodeLayout(ROWS, COLS, barcodes)

# Generate layout
print("Generating barcode assignments...")
try:
    grid_bar1, grid_bar2, grid_combined = layout.generate_layout(seed=42)
    print("✓ Layout generated successfully!")
except ValueError as e:
    print(f"✗ ERROR: {e}")
    exit(1)

# Get statistics
print("\n" + "=" * 70)
print("LAYOUT STATISTICS")
print("=" * 70)
stats = layout.get_statistics()
print(f"Grid dimensions: {stats['grid_dimensions']}")
print(f"Total spots: {stats['total_spots']:,}")
print(f"Number of barcodes used: {stats['n_barcodes_available']:,}")
print(f"Number of blocks: {stats['n_blocks']:,}")
print(f"Unique combined barcodes: {stats['unique_combined_barcodes']:,}")
print(f"Theoretical max combinations: {stats['theoretical_max_combinations']:,}")

# Bar_1 usage distribution
print(f"\nbar_1 usage statistics:")
bar1_usage_counts = {}
for count in stats['bar1_usage'].values():
    bar1_usage_counts[count] = bar1_usage_counts.get(count, 0) + 1
print(f"  Barcodes used once: {bar1_usage_counts.get(4, 0):,}")
if any(c > 4 for c in stats['bar1_usage'].values()):
    print(f"  Warning: Some barcodes used in multiple blocks")

# Export
print("\n" + "=" * 70)
print("EXPORTING")
print("=" * 70)
print(f"Exporting layout to: {OUTPUT_CSV}")
layout.export_layout(OUTPUT_CSV, format='csv')
print("✓ CSV export complete")

# Generate visualization (only bar_1 and bar_2, skip combined)
print(f"\nGenerating visualization (bar_1 and bar_2 only)...")
print(f"Saving to: {OUTPUT_VIZ}")
layout.visualize_layout(OUTPUT_VIZ, include_combined=False)
print("✓ Visualization complete")

print("\n" + "=" * 70)
print("DONE!")
print("=" * 70)
print(f"\nOutput file: {OUTPUT_CSV}")
print(f"Total unique barcodes: {stats['unique_combined_barcodes']:,}")
