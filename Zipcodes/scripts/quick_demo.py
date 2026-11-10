#!/usr/bin/env python3
"""
Quick demo of spatial zipcode generation with a tiny example.
"""

import logging
from pathlib import Path

from spatial_zipcodes.axis_code_generator import generate_axis_codes
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.codebook import build_codebook

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

print("=" * 70)
print("Spatial Zipcode Generation - Quick Demo")
print("=" * 70)

# Configuration
rows, cols = 5, 8
target_distances = [2, 3, 4, 5]
far_min = 5

print(f"\nConfiguration:")
print(f"  Chip: {rows} rows x {cols} columns = {rows*cols} features")
print(f"  Target distances: {target_distances}")
print(f"  Far-features minimum: {far_min}")

# Generate X codes
print(f"\nGenerating X axis codes...")
x_codes, x_ver = generate_axis_codes(
    n_positions=cols,
    bit_length=10,
    target_distances=target_distances,
    far_features_min_distance=far_min,
    random_seed=42,
)
print(f"  ✓ {len(x_codes)} codes, {len(x_codes[0])} bits each")
print(f"  Hard constraint violations: {x_ver['hard_constraint_violations']}")

# Generate Y codes
print(f"\nGenerating Y axis codes...")
y_codes, y_ver = generate_axis_codes(
    n_positions=rows,
    bit_length=10,
    target_distances=target_distances,
    far_features_min_distance=far_min,
    random_seed=43,
)
print(f"  ✓ {len(y_codes)} codes, {len(y_codes[0])} bits each")
print(f"  Hard constraint violations: {y_ver['hard_constraint_violations']}")

# Convert to DNA
print(f"\nConverting to DNA sequences...")
encoder = BitToDNAEncoder(scheme="2-bits-per-base")
x_sequences = [encoder.encode(code) for code in x_codes]
y_sequences = [encoder.encode(code) for code in y_codes]
print(f"  X sequences: {len(x_sequences[0])} bp each")
print(f"  Y sequences: {len(y_sequences[0])} bp each")

# Build codebook
print(f"\nBuilding codebook...")
codebook = build_codebook(
    rows=rows,
    cols=cols,
    x_codes=x_codes,
    y_codes=y_codes,
    x_sequences=x_sequences,
    y_sequences=y_sequences,
    spacer="GGG",
)

stats = codebook.get_statistics()
print(f"  ✓ {stats['n_features']} features")
print(f"  Zip code length: {stats['zip_length']} bp")

# Show a few examples
print(f"\nExample zip codes:")
for i in range(min(5, len(codebook.features))):
    f = codebook.features[i]
    print(f"  {f.feature_id}: {f.zip_seq}")

# Export
output_file = Path("demo_codebook.csv")
codebook.export_to_csv(output_file)
print(f"\n✓ Codebook exported to {output_file}")

print("\n" + "=" * 70)
print("Demo complete!")
print("=" * 70)
