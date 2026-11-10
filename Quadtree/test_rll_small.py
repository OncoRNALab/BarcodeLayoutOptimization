#!/usr/bin/env python3
"""Test RLL implementation with small grid"""

import numpy as np
from quadtree_rll import generate_barcode_matrix_rll, verify_no_homopolymers, count_homopolymers

print("Testing RLL with small grid...")
print("="*60)

# Small test
ROWS = 10
COLS = 10
LENGTH = 24
MAX_RUN = 1

print(f"Grid: {ROWS}x{COLS} = {ROWS*COLS} spots")
print(f"Length: {LENGTH}nt")
print(f"Max run: {MAX_RUN}\n")

print("Generating barcodes...")
barcodes, stats = generate_barcode_matrix_rll(ROWS, COLS, LENGTH, MAX_RUN)

print(f"\n✅ Generation complete!")
print(f"  Shape: {barcodes.shape}")
print(f"  Sample barcode: {barcodes[0, 0]}")

print(f"\nRLL Statistics:")
print(f"  Total forced flips: {stats['total_forced_flips']}")
print(f"  Mean forced flip rate: {stats['mean_forced_flip_rate']:.4f}")
print(f"  Violations avoided: {stats['total_violations_avoided']}")

print(f"\nVerifying homopolymers...")
is_valid = verify_no_homopolymers(barcodes, MAX_RUN)
if is_valid:
    print(f"  ✅ No homopolymers found!")
else:
    print(f"  ⚠️  WARNING: Homopolymers detected!")
    homopolymers = count_homopolymers(barcodes, min_length=2)
    total = sum(sum(counts.values()) for counts in homopolymers.values())
    print(f"  Total homopolymer occurrences: {total}")
    for base in ['A', 'C', 'G', 'T']:
        for length in range(2, 5):
            count = homopolymers[base][length]
            if count > 0:
                print(f"    {base*length}: {count}")

print("\n" + "="*60)
print("Test complete!")
