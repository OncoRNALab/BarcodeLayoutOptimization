#!/usr/bin/env python3
"""Structure-matched decoys for Dual / Dual-like (HSDB spot-level) layouts.

Every block keeps its real block identifier (bar_1) and receives extra
spot codes (bar_2) that were not printed with that bar_1. The decoy share is
therefore the same inside every block as in the whole whitelist, which is what
the normalised decoy-based precision estimator assumes.

Usage:
    python make_dual_decoys.py \\
        --layout .../layout_235x363_dual_shared_pool_g10.csv \\
        --library .../bar_1in1_column_major_first18.txt \\
        --block-size 4 --bar1-len 18 --bar2-len 18 \\
        --g-max 0.10 --max-homopolymer 5 \\
        --out-dir ../output/dual_shared_pool_g10
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    m = c = 1
    for i in range(1, len(seq)):
        c = c + 1 if seq[i] == seq[i - 1] else 1
        m = max(m, c)
    return m


def g_fraction(seq: str) -> float:
    return seq.count("G") / len(seq) if seq else 0.0


def passes_combined_qc(bar1: str, bar2: str, g_max, max_hp) -> bool:
    combined = bar1 + bar2
    if g_max is not None and g_fraction(combined) > g_max:
        return False
    if max_hp is not None and max_homopolymer(combined) > max_hp:
        return False
    return True


def mutate_bar2(parent: str, n_flips: int, rng: random.Random) -> str:
    """Flip n_flips positions in parent (structure-local decoy bar_2)."""
    bases = list(parent)
    positions = rng.sample(range(len(bases)), k=min(n_flips, len(bases)))
    alphabet = "ACGT"
    for i in positions:
        choices = [b for b in alphabet if b != bases[i]]
        bases[i] = rng.choice(choices)
    return "".join(bases)


def read_layout(path):
    layout = {}
    with open(path) as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        seq_key = "combined" if "combined" in fields else "barcode"
        for rec in reader:
            layout[(int(rec["row"]), int(rec["col"]))] = rec[seq_key].strip().upper()
    return layout


def read_library(path):
    with open(path) as handle:
        return sorted({line.strip().upper() for line in handle if line.strip()})


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--layout", required=True,
                        help="Layout CSV with row,col,combined (or barcode)")
    parser.add_argument("--library", default=None,
                        help="bar_2 candidate library (required for --bar2-mode library)")
    parser.add_argument("--block-size", type=int, default=4)
    parser.add_argument("--segment-length", type=int, default=None,
                        help="If set, use as both bar1-len and bar2-len")
    parser.add_argument("--bar1-len", type=int, default=None)
    parser.add_argument("--bar2-len", type=int, default=None)
    parser.add_argument("--fraction", type=float, default=0.5,
                        help="Decoys per real barcode, per block")
    parser.add_argument(
        "--bar2-mode",
        choices=("library", "mutate", "mutate_bar1", "mutate_mixed"),
        default="library",
        help="library: unused bar_2 from --library; "
             "mutate: flip in-block bar_2 (same bar_1); "
             "mutate_bar1: flip bar_1 (same in-block bar_2) — better when bar_2 is reused "
             "across blocks; mutate_mixed: alternate mutate / mutate_bar1",
    )
    parser.add_argument("--mutate-flips-min", type=int, default=2)
    parser.add_argument("--mutate-flips-max", type=int, default=6)
    parser.add_argument(
        "--target-nn",
        type=int,
        default=None,
        help="Match decoy→real NN to this Hamming distance: sets flip count to "
             "target-nn (± tol) and rejects candidates whose NN to relevant reals "
             "is outside [target-nn-tol, target-nn+tol]",
    )
    parser.add_argument("--nn-tol", type=int, default=1,
                        help="Allowed deviation from --target-nn (default 1)")
    parser.add_argument("--g-max", type=float, default=None,
                        help="Max G fraction on combined decoy (e.g. 0.10)")
    parser.add_argument("--max-homopolymer", type=int, default=None)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-attempts-per-decoy", type=int, default=50000)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    if args.segment_length is not None:
        bar1_len = bar2_len = args.segment_length
    else:
        bar1_len = args.bar1_len if args.bar1_len is not None else 18
        bar2_len = args.bar2_len if args.bar2_len is not None else 18

    flips_min, flips_max = args.mutate_flips_min, args.mutate_flips_max
    if args.target_nn is not None:
        if args.target_nn < 1:
            raise SystemExit("--target-nn must be ≥ 1")
        flips_min = max(1, args.target_nn - args.nn_tol)
        flips_max = args.target_nn + args.nn_tol

    rng = random.Random(args.seed)
    layout = read_layout(args.layout)
    library = []
    if args.bar2_mode == "library":
        if not args.library:
            raise SystemExit("--library is required with --bar2-mode library")
        library = [s for s in read_library(args.library) if len(s) == bar2_len]
        if not library:
            raise SystemExit(f"No bar_2 sequences of length {bar2_len} in {args.library}")
    real_set = set(layout.values())
    if flips_min < 1 or flips_max < flips_min:
        raise SystemExit("Need 1 <= mutate flips min <= max")

    blocks = defaultdict(list)
    for (row, col), seq in layout.items():
        if len(seq) != bar1_len + bar2_len:
            raise ValueError(
                f"Sequence length {len(seq)} at ({row},{col}) != "
                f"bar1+bar2={bar1_len}+{bar2_len}"
            )
        blocks[(row // args.block_size, col // args.block_size)].append(seq)

    bar1_blocks = defaultdict(set)
    for key, seqs in blocks.items():
        ids = {s[:bar1_len] for s in seqs}
        if len(ids) != 1:
            raise ValueError(f"Block {key} has {len(ids)} different bar_1 sequences")
        bar1_blocks[ids.pop()].add(key)
    shared = [b for b, keys in bar1_blocks.items() if len(keys) > 1]
    if shared:
        raise ValueError(f"{len(shared)} bar_1 sequences are used in more than one block")

    # Index for NN checks when targeting a distance shell
    bar1_by_bar2 = defaultdict(set)
    bar2_by_bar1 = defaultdict(set)
    for seq in real_set:
        b1, b2 = seq[:bar1_len], seq[bar1_len:]
        bar1_by_bar2[b2].add(b1)
        bar2_by_bar1[b1].add(b2)

    def hamming(a: str, b: str) -> int:
        return sum(x != y for x, y in zip(a, b))

    def nn_ok_mutate_bar1(decoy_bar1: str, bar2: str) -> bool:
        if args.target_nn is None:
            return True
        lo, hi = args.target_nn - args.nn_tol, args.target_nn + args.nn_tol
        # Same bar2: distance = hamm(decoy_bar1, real_bar1)
        dmin = min(
            (hamming(decoy_bar1, b1) for b1 in bar1_by_bar2.get(bar2, ())),
            default=10**9,
        )
        return lo <= dmin <= hi

    def nn_ok_mutate_bar2(bar1: str, decoy_bar2: str) -> bool:
        if args.target_nn is None:
            return True
        lo, hi = args.target_nn - args.nn_tol, args.target_nn + args.nn_tol
        # Same bar1: distance = hamm(decoy_bar2, real_bar2)
        dmin = min(
            (hamming(decoy_bar2, b2) for b2 in bar2_by_bar1.get(bar1, ())),
            default=10**9,
        )
        return lo <= dmin <= hi

    decoys = []
    seen = set()
    n_qc_rejects = 0
    n_nn_rejects = 0
    for key in sorted(blocks):
        seqs = blocks[key]
        bar1 = seqs[0][:bar1_len]
        used_bar2 = {s[bar1_len:] for s in seqs}
        n_decoys = round(len(seqs) * args.fraction)
        added = 0
        attempts = 0
        max_attempts = max(args.max_attempts_per_decoy * max(n_decoys, 1), 1000)
        real_bar2_list = list(used_bar2)
        real_pairs = [(s[:bar1_len], s[bar1_len:]) for s in seqs]
        while added < n_decoys:
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(
                    f"Block {key}: only added {added}/{n_decoys} decoys after "
                    f"{attempts} attempts (mode={args.bar2_mode}, library={len(library)}, "
                    f"QC rejects≈{n_qc_rejects}, NN rejects≈{n_nn_rejects}). "
                    f"Enlarge --library, change --bar2-mode/--target-nn, or relax QC."
                )
            mode = args.bar2_mode
            if mode == "mutate_mixed":
                mode = "mutate_bar1" if (added % 2 == 0) else "mutate"

            if mode == "mutate_bar1":
                # Same real bar_2, mutated bar_1 — matches cross-block same-bar_2 confusion.
                parent_b1, bar2 = rng.choice(real_pairs)
                n_flips = rng.randint(flips_min, flips_max)
                decoy_bar1 = mutate_bar2(parent_b1, n_flips, rng)
                if not passes_combined_qc(decoy_bar1, bar2, args.g_max, args.max_homopolymer):
                    n_qc_rejects += 1
                    continue
                if not nn_ok_mutate_bar1(decoy_bar1, bar2):
                    n_nn_rejects += 1
                    continue
                combined = decoy_bar1 + bar2
                if combined in real_set or combined in seen:
                    continue
                seen.add(combined)
                decoys.append({
                    "block_row": key[0],
                    "block_col": key[1],
                    "bar_1": decoy_bar1,
                    "bar_2": bar2,
                    "combined": combined,
                })
                added += 1
                continue

            if mode == "mutate":
                parent = rng.choice(real_bar2_list)
                n_flips = rng.randint(flips_min, flips_max)
                bar2 = mutate_bar2(parent, n_flips, rng)
            else:
                bar2 = rng.choice(library)
            if bar2 in used_bar2:
                continue
            if not passes_combined_qc(bar1, bar2, args.g_max, args.max_homopolymer):
                n_qc_rejects += 1
                continue
            if mode == "mutate" and not nn_ok_mutate_bar2(bar1, bar2):
                n_nn_rejects += 1
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
        "bar1_len": bar1_len,
        "bar2_len": bar2_len,
        "bar2_mode": args.bar2_mode,
        "mutate_flips_min": flips_min if args.bar2_mode != "library" else None,
        "mutate_flips_max": flips_max if args.bar2_mode != "library" else None,
        "target_nn": args.target_nn,
        "nn_tol": args.nn_tol if args.target_nn is not None else None,
        "g_max": args.g_max,
        "max_homopolymer": args.max_homopolymer,
        "n_blocks": len(blocks),
        "n_real": len(real_set),
        "n_decoys": len(decoys),
        "n_qc_rejects_during_sampling": n_qc_rejects,
        "n_nn_rejects_during_sampling": n_nn_rejects,
        "decoy_to_real": len(decoys) / len(real_set),
        "decoy_fraction_of_whitelist": len(decoys) / (len(decoys) + len(real_set)),
        "phi_design": len(decoys) / (len(decoys) + len(real_set)),
    }
    with open(out_dir / "dual_decoys_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
