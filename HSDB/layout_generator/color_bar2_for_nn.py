#!/usr/bin/env python3
"""Raise HSDB combined NN by coloring blocks onto multiple bar_2 libraries.

Zipcode bar_1 can leave some non-adjacent block pairs at DNA Hamming < target
(e.g. 8). With a single reused bar_2 set, those pairs become combined NN
neighbours. This script:

  1. Builds a conflict graph on blocks (edge if hamm(bar_1)<target).
  2. Greedy-colors the graph.
  3. Generates one QC/min-Hamming bar_2 set of 16 per color.
  4. Writes a layout CSV with per-block bar_2 color.

Hierarchical decoding still works: bar_1 → block → bar_2 within that block's set.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from pathlib import Path
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Dual_barcodes"))
from sequence_qc import SequenceQC, hamming  # noqa: E402


def gen_bar2_set(n, length, min_h, qc, seed):
    rng = random.Random(seed)
    sel, seen = [], set()
    attempts = 0
    while len(sel) < n and attempts < 500_000:
        attempts += 1
        s = "".join(rng.choices("ACGT", k=length))
        if s in seen:
            continue
        seen.add(s)
        if not qc.passes_half(s):
            continue
        if all(hamming(s, t) >= min_h for t in sel):
            sel.append(s)
    if len(sel) < n:
        raise RuntimeError(f"only {len(sel)}/{n} bar_2 at minH={min_h}")
    return sel


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bar1", required=True)
    ap.add_argument("--array-rows", type=int, default=235)
    ap.add_argument("--array-cols", type=int, default=363)
    ap.add_argument("--block-size", type=int, default=4)
    ap.add_argument("--target-nn", type=int, default=8)
    ap.add_argument("--bar2-len", type=int, default=12)
    ap.add_argument("--bar2-min-hamming", type=int, default=8)
    ap.add_argument("--g-max", type=float, default=0.10)
    ap.add_argument("--max-homopolymer", type=int, default=5)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    bar1 = [l.strip().upper() for l in open(args.bar1) if l.strip()]
    B = args.block_size
    brs = (args.array_rows + B - 1) // B
    bcs = (args.array_cols + B - 1) // B
    if len(bar1) != brs * bcs:
        raise SystemExit(f"Need {brs*bcs} bar_1, got {len(bar1)}")

    print(f"Conflict graph (hamm(bar_1)<{args.target_nn})…")
    conflicts = [[] for _ in bar1]
    n_edges = 0
    for i, bi in enumerate(bar1):
        if i % 500 == 0:
            print(f"  {i}/{len(bar1)}")
        for j in range(i + 1, len(bar1)):
            if hamming(bi, bar1[j]) < args.target_nn:
                conflicts[i].append(j)
                conflicts[j].append(i)
                n_edges += 1
    print(f"  edges={n_edges} max_degree={max(map(len, conflicts))}")

    order = sorted(range(len(bar1)), key=lambda i: -len(conflicts[i]))
    color = [-1] * len(bar1)
    for i in order:
        used = {color[j] for j in conflicts[i] if color[j] >= 0}
        c = 0
        while c in used:
            c += 1
        color[i] = c
    n_colors = max(color) + 1
    counts = [color.count(k) for k in range(n_colors)]
    print(f"  n_colors={n_colors} sizes={counts}")

    qc = SequenceQC(
        g_content_max=args.g_max,
        max_homopolymer=args.max_homopolymer,
        max_leading=1,
    )
    bar2_sets = [
        gen_bar2_set(B * B, args.bar2_len, args.bar2_min_hamming, qc, args.seed + k)
        for k in range(n_colors)
    ]

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "block_bar2_colors.json", "w") as fh:
        json.dump(
            {
                "n_colors": n_colors,
                "min_bar1_hamming_threshold": args.target_nn,
                "n_conflict_edges": n_edges,
                "color_counts": counts,
                "bar2_sets": bar2_sets,
                "block_color": color,
            },
            fh,
            indent=2,
        )

    rng = random.Random(args.seed)
    spots = []
    for br in range(brs):
        for bc in range(bcs):
            bi = br * bcs + bc
            b1 = bar1[bi]
            lib = bar2_sets[color[bi]][:]
            rng.shuffle(lib)
            si = 0
            for di in range(B):
                for dj in range(B):
                    r, c = br * B + di, bc * B + dj
                    if r >= args.array_rows or c >= args.array_cols:
                        continue
                    b2 = lib[si]
                    si += 1
                    spots.append(
                        {
                            "row": r,
                            "col": c,
                            "block_row": br,
                            "block_col": bc,
                            "spot_in_block": si - 1,
                            "bar_1": b1,
                            "bar_2": b2,
                            "combined": b1 + b2,
                            "bar2_color": color[bi],
                        }
                    )

    n_unique = len({s["combined"] for s in spots})
    if n_unique != len(spots):
        raise RuntimeError(f"non-unique combined: {n_unique}/{len(spots)}")

    with open(out / "layout_235x363_hsdb.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(spots[0].keys()))
        w.writeheader()
        w.writerows(spots)
    with open(out / "barcodes_hsdb_235x363_36nt.txt", "w") as fh:
        fh.writelines(s["combined"] + "\n" for s in spots)

    # quick NN sample
    BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
    e = np.array([[BASE[b] for b in s["combined"]] for s in spots], dtype=np.int8)
    rng_np = np.random.default_rng(1)
    idx = rng_np.choice(len(e), min(1000, len(e)), replace=False)
    nn = []
    for q in e[idx]:
        d = (e != q).sum(1)
        d[d == 0] = 10**6
        nn.append(int(d.min()))
    print(f"Wrote {out}  NN sample mean={np.mean(nn):.3f} [{min(nn)}-{max(nn)}]")


if __name__ == "__main__":
    main()
