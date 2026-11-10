#!/usr/bin/env python3
"""Greedy + 2-OPT decoys at the same decoy-to-real ratio as Dual and Zipcodes.

The Greedy decoys are pool barcodes that were not placed on the array. They
are unrelated to array position, so a random subset is as good as the full
set; subsampling only sets the decoy-to-real ratio (default 1:2).

Usage:
    python make_greedy_decoys.py \
        --real ../../data/1in1/layout_235x363_randombarcodes_opt.txt \
        --decoys ../../data/decoy_layouts/barcodes/1in1/greedy_decoybar_1in1_column_major.txt \
        --out-dir ../output/greedy
"""

import argparse
import json
import random
from pathlib import Path


def read_seqs(path):
    with open(path) as handle:
        return [line.strip().upper() for line in handle if line.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--real", required=True, help="Printed Greedy layout, one barcode per line")
    parser.add_argument("--decoys", required=True, help="Unplaced pool barcodes, one per line")
    parser.add_argument("--fraction", type=float, default=0.5, help="Decoys per real barcode")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    real = read_seqs(args.real)
    pool = read_seqs(args.decoys)
    if len(set(real)) != len(real):
        raise ValueError("Duplicate barcodes in the real layout")
    if len(set(pool)) != len(pool):
        raise ValueError("Duplicate barcodes in the decoy pool")
    if set(real) & set(pool):
        raise ValueError("Decoy pool overlaps the real layout")

    n_decoys = round(args.fraction * len(real))
    if n_decoys > len(pool):
        raise ValueError(f"Requested {n_decoys} decoys but the pool has {len(pool)}")
    keep = sorted(random.Random(args.seed).sample(range(len(pool)), n_decoys))
    decoys = [pool[i] for i in keep]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "greedy_decoys.txt", "w") as handle:
        handle.writelines(s + "\n" for s in decoys)
    with open(out_dir / "greedy_whitelist_real_plus_decoys.txt", "w") as handle:
        handle.writelines(s + "\n" for s in real)
        handle.writelines(s + "\n" for s in decoys)

    summary = {
        "real": str(args.real),
        "decoy_pool": str(args.decoys),
        "seed": args.seed,
        "n_real": len(real),
        "n_pool": len(pool),
        "n_decoys": len(decoys),
        "decoy_to_real": len(decoys) / len(real),
        "decoy_fraction_of_whitelist": len(decoys) / (len(decoys) + len(real)),
        "whitelist_order": "real barcodes in layout (column-major) order, then decoys in pool order",
    }
    with open(out_dir / "greedy_decoys_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
