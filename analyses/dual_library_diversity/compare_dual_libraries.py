#!/usr/bin/env python3
"""Does using one 18-nt library for both Dual-barcode parts reduce k-mer diversity?

Builds Dual layouts with the assignment rules of dual_barcode_layout.py:
  shared   : block identifiers and spot codes from the same library (as printed)
  separate : block identifiers from library 1, spot codes from library 2
and compares them, together with the printed layout, on metrics relevant to a
k-mer filtering decoder.

Usage:
    python compare_dual_libraries.py \
        --printed ../../data/1in1/layout_235x363_dual_4x4.csv \
        --lib1 ../../data/scripts/inhouse_barcodes/bar_1in1_column_major_first18.txt \
        --lib2 ../../data/scripts/inhouse_barcodes/bar_1in1_column_major_last18.txt \
        --out dual_library_diversity.json
"""

import argparse
import csv
import itertools
import json
import random
from collections import Counter, defaultdict

import numpy as np

ROWS, COLS, BLOCK, SEG = 235, 363, 4, 18
KMER_INDEX = {"".join(k): i for i, k in enumerate(itertools.product("ACGT", repeat=4))}


def read_lines(path):
    with open(path) as handle:
        return [line.strip().upper() for line in handle if line.strip()]


def build_layout(bar1_lib, bar2_lib, seed):
    """Block identifiers: least-used from a shuffled library; spot codes: random, unique per block."""
    rng = random.Random(seed)
    lib1 = bar1_lib[:]
    rng.shuffle(lib1)
    usage = Counter()
    layout = {}
    blocks = defaultdict(list)
    for r in range(ROWS):
        for c in range(COLS):
            blocks[(r // BLOCK, c // BLOCK)].append((r, c))
    order = iter(lib1)
    for key in sorted(blocks):
        bar1 = next(order)
        usage[bar1] += 1
        used = set()
        for pos in blocks[key]:
            while True:
                bar2 = rng.choice(bar2_lib)
                if bar2 not in used:
                    break
            used.add(bar2)
            layout[pos] = bar1 + bar2
    return layout


def read_printed(path):
    with open(path) as handle:
        return {(int(r["row"]), int(r["col"])): r["combined"].upper() for r in csv.DictReader(handle)}


def kmer_matrix(seqs):
    m = np.zeros((len(seqs), 256), dtype=np.float32)
    for i, s in enumerate(seqs):
        for j in range(len(s) - 3):
            m[i, KMER_INDEX[s[j:j + 4]]] += 1
    return m


def metrics(layout, sample, seed):
    positions = sorted(layout)
    seqs = [layout[p] for p in positions]
    block = np.array([(r // BLOCK) * 1000 + c // BLOCK for r, c in positions])
    bar1 = [s[:SEG] for s in seqs]
    bar2 = [s[SEG:] for s in seqs]

    both_roles = set(bar1) & set(bar2)
    cross_pairs = sum(Counter(bar1)[x] * Counter(bar2)[x] for x in both_roles)

    m = kmer_matrix(seqs)
    unit = m / np.linalg.norm(m, axis=1, keepdims=True)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(seqs), size=sample, replace=False)
    max_cos_out, shared_kmers_out = [], []
    for i in idx:
        cos = unit @ unit[i]
        cos[block == block[i]] = -1
        j = int(np.argmax(cos))
        max_cos_out.append(float(cos[j]))
        shared_kmers_out.append(float(np.minimum(m[i], m[j]).sum()))

    probs = m / m.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(np.where(probs > 0, probs * np.log2(probs), 0)).sum(1)
    library_kmers = m.sum(0) / m.sum()
    library_entropy = float(-(library_kmers[library_kmers > 0] * np.log2(library_kmers[library_kmers > 0])).sum())

    return {
        "n_barcodes": len(seqs),
        "distinct_bar1": len(set(bar1)),
        "distinct_bar2": len(set(bar2)),
        "sequences_in_both_roles": len(both_roles),
        "barcode_pairs_sharing_an_18mer_in_swapped_roles": cross_pairs,
        "nearest_outside_block_cosine_4mer_mean": float(np.mean(max_cos_out)),
        "nearest_outside_block_cosine_4mer_p95": float(np.percentile(max_cos_out, 95)),
        "nearest_outside_block_shared_4mers_mean": float(np.mean(shared_kmers_out)),
        "per_barcode_4mer_entropy_mean": float(ent.mean()),
        "library_4mer_entropy_bits": library_entropy,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--printed", required=True)
    parser.add_argument("--lib1", required=True)
    parser.add_argument("--lib2", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample", type=int, default=3000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    lib1 = sorted(set(read_lines(args.lib1)))
    lib2 = sorted(set(read_lines(args.lib2)) - set(lib1))
    designs = {
        "printed_shared": read_printed(args.printed),
        "simulated_shared": build_layout(lib1, lib1, args.seed),
        "simulated_separate": build_layout(lib1, lib2, args.seed),
    }
    report = {}
    for name, layout in designs.items():
        report[name] = metrics(layout, args.sample, args.seed)
        print(name, json.dumps(report[name], indent=2))
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=2)


if __name__ == "__main__":
    main()
