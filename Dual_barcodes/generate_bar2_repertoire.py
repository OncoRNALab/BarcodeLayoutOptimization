#!/usr/bin/env python3
"""
Generate a fixed Dual bar_2 repertoire sized by the block.

n_bar2 = block_size² (plus optional spares). Codes are selected for:
  - minimum pairwise Hamming distance (target, then fallback)
  - half-sequence QC (G, homopolymer, optional GC, leading run ≤ 1)

For lengths ≤ 12, candidates can be fully enumerated. For longer halves
(e.g. 18 nt), candidates are drawn at random.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from pathlib import Path
from typing import List, Optional

from sequence_qc import (
    SequenceQC,
    gc_content,
    g_fraction,
    hamming,
    leading_run,
    max_homopolymer,
    min_pairwise_hamming,
)

BASES = "ACGT"
ENUM_LENGTH_LIMIT = 12


def greedy_select(candidates: List[str], n: int, min_dist: int) -> List[str]:
    selected: List[str] = []
    for cand in candidates:
        if all(hamming(cand, s) >= min_dist for s in selected):
            selected.append(cand)
            if len(selected) == n:
                break
    return selected


def enumerate_candidates(length: int, qc: SequenceQC) -> List[str]:
    return [
        "".join(bases)
        for bases in itertools.product(BASES, repeat=length)
        if qc.passes_half("".join(bases))
    ]


def sample_candidates(
    length: int,
    qc: SequenceQC,
    n_candidates: int,
    seed: int,
) -> List[str]:
    rng = random.Random(seed)
    seen = set()
    out: List[str] = []
    attempts = 0
    max_attempts = max(n_candidates * 50, 100_000)
    while len(out) < n_candidates and attempts < max_attempts:
        attempts += 1
        seq = "".join(rng.choices(BASES, k=length))
        if seq in seen:
            continue
        seen.add(seq)
        if qc.passes_half(seq):
            out.append(seq)
    return out


def generate_bar2_repertoire(
    n: int,
    length: int,
    min_hamming: int = 8,
    g_max: float = 0.10,
    max_homopolymer_limit: int = 5,
    gc_min: Optional[float] = None,
    gc_max: Optional[float] = None,
    max_leading: int = 1,
    seed: int = 42,
    n_random_candidates: int = 200_000,
    min_hamming_floor: int = 4,
) -> tuple[List[str], int]:
    """
    Return (selected_barcodes, achieved_min_hamming).
    Tries min_hamming down to min_hamming_floor.
    """
    qc = SequenceQC(
        g_content_max=g_max,
        max_homopolymer=max_homopolymer_limit,
        gc_min=gc_min,
        gc_max=gc_max,
        max_leading=max_leading,
    )

    if length <= ENUM_LENGTH_LIMIT:
        print(f"Enumerating QC-passing {length}-mers...", flush=True)
        candidates = enumerate_candidates(length, qc)
    else:
        print(
            f"Sampling up to {n_random_candidates:,} QC-passing {length}-mers "
            f"(seed={seed})...",
            flush=True,
        )
        candidates = sample_candidates(length, qc, n_random_candidates, seed)

    rng = random.Random(seed)
    rng.shuffle(candidates)
    print(f"  {len(candidates):,} candidates pass half QC")

    selected: List[str] = []
    achieved = 0
    for dist in range(min_hamming, min_hamming_floor - 1, -1):
        selected = greedy_select(candidates, n, dist)
        if len(selected) == n:
            achieved = dist
            break
        print(f"  Only {len(selected)}/{n} at min_dist={dist}; retrying lower...")

    if len(selected) < n:
        raise RuntimeError(
            f"Could not select {n} bar_2 codes of length {length} with "
            f"min Hamming ≥ {min_hamming_floor} under the given QC"
        )

    return selected, min_pairwise_hamming(selected)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--block-size", type=int, default=4,
                    help="Block side length B; n_bar2 defaults to B²")
    ap.add_argument("--n", type=int, default=None,
                    help="Number of bar_2 codes (default block_size²)")
    ap.add_argument("--length", type=int, default=18)
    ap.add_argument("--min-hamming", type=int, default=8)
    ap.add_argument("--g-max", type=float, default=0.10)
    ap.add_argument("--max-homopolymer", type=int, default=5)
    ap.add_argument("--gc-min", type=float, default=None)
    ap.add_argument("--gc-max", type=float, default=None)
    ap.add_argument("--max-leading", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-candidates", type=int, default=200_000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats-out", default=None)
    args = ap.parse_args()

    n = args.n if args.n is not None else args.block_size ** 2
    print(f"Generating {n} bar_2 codes of length {args.length} "
          f"(block_size={args.block_size})", flush=True)

    try:
        selected, achieved = generate_bar2_repertoire(
            n=n,
            length=args.length,
            min_hamming=args.min_hamming,
            g_max=args.g_max,
            max_homopolymer_limit=args.max_homopolymer,
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            max_leading=args.max_leading,
            seed=args.seed,
            n_random_candidates=args.n_candidates,
        )
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nSelected {len(selected)} barcodes; min Hamming = {achieved}")
    print("  #   Sequence              GC    G     MaxHP  Lead")
    print("  " + "-" * 56)
    rows = []
    for i, seq in enumerate(selected):
        row = {
            "idx": i,
            "seq": seq,
            "gc_content": round(gc_content(seq), 4),
            "g_fraction": round(g_fraction(seq), 4),
            "max_homopolymer": max_homopolymer(seq),
            "leading_run": leading_run(seq),
        }
        rows.append(row)
        print(
            f"  {i+1:2d}  {seq}  {row['gc_content']:.2f}  "
            f"{row['g_fraction']:.2f}  {row['max_homopolymer']}      {row['leading_run']}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for seq in selected:
            fh.write(seq + "\n")
    print(f"\nWrote {len(selected)} barcodes → {out}")

    if args.stats_out:
        stats = {
            "n_barcodes": len(selected),
            "barcode_length": args.length,
            "block_size": args.block_size,
            "min_hamming_distance": achieved,
            "g_content_max": args.g_max,
            "max_homopolymer": args.max_homopolymer,
            "barcodes": rows,
        }
        with open(args.stats_out, "w") as fh:
            json.dump(stats, fh, indent=2)
        print(f"Wrote stats → {args.stats_out}")


if __name__ == "__main__":
    main()
