#!/usr/bin/env python3
"""Structure-matched decoys for the printed Dual-barcode layout.

Every block keeps its real block identifier (bar_1) and receives extra
spot codes (bar_2) that were not printed with that bar_1. The decoy share is
therefore the same inside every block as in the whole whitelist, which is what
the normalised decoy-based precision estimator assumes.

Usage:
    python make_dual_decoys.py \
        --layout ../../data/1in1/layout_235x363_dual_4x4.csv \
        --library ../../inhouse_barcodes/bar_1in1_column_major_first18.txt \
        --out-dir ../output/dual
"""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path


def read_layout(path):
    layout = {}
    with open(path) as handle:
        for rec in csv.DictReader(handle):
            layout[(int(rec["row"]), int(rec["col"]))] = rec["combined"].strip().upper()
    return layout


def read_library(path):
    with open(path) as handle:
        return sorted({line.strip().upper() for line in handle if line.strip()})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--layout", required=True, help="Printed layout CSV with columns row,col,combined")
    parser.add_argument("--library", required=True, help="Segment library used to build the layout")
    parser.add_argument("--block-size", type=int, default=4)
    parser.add_argument("--segment-length", type=int, default=18)
    parser.add_argument("--fraction", type=float, default=0.5, help="Decoys per real barcode, per block")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    seg = args.segment_length
    layout = read_layout(args.layout)
    library = read_library(args.library)
    real_set = set(layout.values())

    blocks = defaultdict(list)
    for (row, col), seq in layout.items():
        blocks[(row // args.block_size, col // args.block_size)].append(seq)

    bar1_blocks = defaultdict(set)
    for key, seqs in blocks.items():
        ids = {s[:seg] for s in seqs}
        if len(ids) != 1:
            raise ValueError(f"Block {key} has {len(ids)} different bar_1 sequences")
        bar1_blocks[ids.pop()].add(key)
    shared = [b for b, keys in bar1_blocks.items() if len(keys) > 1]
    if shared:
        raise ValueError(f"{len(shared)} bar_1 sequences are used in more than one block")

    decoys = []
    seen = set()
    for key in sorted(blocks):
        seqs = blocks[key]
        bar1 = seqs[0][:seg]
        used_bar2 = {s[seg:] for s in seqs}
        n_decoys = round(len(seqs) * args.fraction)
        added = 0
        while added < n_decoys:
            bar2 = rng.choice(library)
            if bar2 in used_bar2:
                continue
            combined = bar1 + bar2
            if combined in real_set or combined in seen:
                continue
            used_bar2.add(bar2)
            seen.add(combined)
            decoys.append({
                "block_row": key[0],
                "block_col": key[1],
                "bar_1": bar1,
                "bar_2": bar2,
                "combined": combined,
            })
            added += 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "dual_decoys.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(decoys[0].keys()))
        writer.writeheader()
        writer.writerows(decoys)
    with open(out_dir / "dual_decoys.txt", "w") as handle:
        handle.writelines(d["combined"] + "\n" for d in decoys)
    with open(out_dir / "dual_whitelist_real_plus_decoys.txt", "w") as handle:
        handle.writelines(s + "\n" for s in sorted(real_set))
        handle.writelines(d["combined"] + "\n" for d in decoys)

    summary = {
        "layout": str(args.layout),
        "library": str(args.library),
        "seed": args.seed,
        "block_size": args.block_size,
        "n_blocks": len(blocks),
        "n_real": len(real_set),
        "n_decoys": len(decoys),
        "decoy_to_real": len(decoys) / len(real_set),
        "decoy_fraction_of_whitelist": len(decoys) / (len(decoys) + len(real_set)),
    }
    with open(out_dir / "dual_decoys_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
