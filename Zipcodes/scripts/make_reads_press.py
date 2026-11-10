#!/usr/bin/env python3
"""
Press-style simulation with a fixed TOTAL number of reads (Q), for structured barcodes:
  LEFT(16) + "GGGG" + RIGHT(16) => length 36

Edits allowed only outside the separator (no subs/ins/del within "GGGG").
Separator may shift due to indels in left/right.
Final output is forced to length 36 with pad/truncate, preserving the separator.

Sampling model (like Press notebook):
- Repeatedly pick a random barcode index
- Draw n ~ Poisson(nave) reads from it (capped so total reaches Q)
- For each read, apply makeerrors_structured_press

Usage:
  python press_structured_totalQ.py barcodes.txt reads.txt truth.txt \
      --total-reads 1000000 --nave 4 \
      --srate 0.0333 --irate 0.0333 --drate 0.0333 --seed 1

Notes:
- If you want EXACTLY one read per draw (uniform over barcodes), set --nave 1 and --poisson false
"""

import argparse
import random
import sys
from typing import List, Tuple
from collections import Counter

BASES = ["A", "C", "G", "T"]

L_PART = 16
SEP = "GGGG"
R_PART = 16
TOTAL_LEN = L_PART + len(SEP) + R_PART


def rand_base(rng: random.Random) -> str:
    return rng.choice(BASES)


def find_sep_span(seq: str, sep: str = SEP) -> Tuple[int, int]:
    pos = seq.find(sep)
    if pos == -1:
        raise RuntimeError("Separator not found in sequence (should never happen if protected).")
    return pos, pos + len(sep)


def binomial_py(rng: random.Random, n: int, p: float) -> int:
    if p <= 0.0:
        return 0
    if p >= 1.0:
        return n
    return sum(1 for _ in range(n) if rng.random() < p)


def poisson_py(rng: random.Random, lam: float) -> int:
    # Knuth algorithm; good for lam ~ small/moderate (Press uses small nave like 4)
    if lam <= 0:
        return 0
    import math
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return k - 1


def mutate_substitutions(chars: List[str], editable: List[bool], srate: float, rng: random.Random) -> None:
    n_edit = sum(editable)
    if n_edit == 0:
        return
    p = min(1.0, (4.0 / 3.0) * srate)
    ns = binomial_py(rng, n_edit, p)
    if ns <= 0:
        return
    edit_positions = [i for i, ok in enumerate(editable) if ok]
    for _ in range(ns):
        idx = edit_positions[rng.randrange(n_edit)]
        chars[idx] = rand_base(rng)  # Press assigns a random base (may equal original)


def apply_deletions(chars: List[str], editable: List[bool], drate: float, rng: random.Random) -> Tuple[List[str], List[bool]]:
    n_edit = sum(editable)
    if n_edit == 0:
        return chars, editable
    nd = binomial_py(rng, n_edit, drate)
    if nd <= 0:
        return chars, editable
    edit_positions = [i for i, ok in enumerate(editable) if ok]
    chosen = [edit_positions[rng.randrange(n_edit)] for _ in range(nd)]
    # collision behavior: duplicates effectively delete once
    for idx in sorted(set(chosen), reverse=True):
        del chars[idx]
        del editable[idx]
    return chars, editable


def apply_insertions(chars: List[str], editable: List[bool], irate: float, rng: random.Random) -> Tuple[List[str], List[bool]]:
    # Forbid insertions between SEP bases (inside the separator interior boundaries)
    seq_str = "".join(chars)
    s0, s1 = find_sep_span(seq_str, SEP)

    n = len(chars)
    valid_slots = []
    for slot in range(n + 1):
        # slots strictly inside SEP are s0+1 .. s1-1
        if (s0 + 1) <= slot <= (s1 - 1):
            continue
        valid_slots.append(slot)

    n_slots = len(valid_slots)
    if n_slots == 0:
        return chars, editable

    ni = binomial_py(rng, n_slots, irate)
    if ni <= 0:
        return chars, editable

    chosen_slots = [valid_slots[rng.randrange(n_slots)] for _ in range(ni)]
    counts = Counter(chosen_slots)

    offset = 0
    for slot in sorted(counts.keys()):
        c = counts[slot]
        insert_at = slot + offset
        for _ in range(c):
            chars.insert(insert_at, rand_base(rng))
            editable.insert(insert_at, True)
            insert_at += 1
            offset += 1
    return chars, editable


def enforce_length_preserve_sep(chars: List[str], target_len: int, rng: random.Random) -> List[str]:
    # pad
    while len(chars) < target_len:
        chars.append(rand_base(rng))
    # truncate while keeping SEP intact
    while len(chars) > target_len:
        seq = "".join(chars)
        s0 = seq.find(SEP)
        if s0 == -1:
            chars = chars[:target_len]
            break
        s1 = s0 + len(SEP)
        if s1 > target_len:
            chars.pop(0)  # drop from left until SEP fits
        else:
            chars.pop()   # drop from right
    return chars


def makeerrors_structured_press(barcode: str, srate: float, irate: float, drate: float, rng: random.Random) -> str:
    barcode = barcode.strip().upper()
    if len(barcode) != TOTAL_LEN:
        raise ValueError(f"Expected barcode length {TOTAL_LEN}, got {len(barcode)}")
    if barcode[L_PART:L_PART + len(SEP)] != SEP:
        raise ValueError(f"Expected separator '{SEP}' at position {L_PART}")

    chars = list(barcode)
    editable = [True] * L_PART + [False] * len(SEP) + [True] * R_PART

    mutate_substitutions(chars, editable, srate, rng)
    chars, editable = apply_deletions(chars, editable, drate, rng)
    chars, editable = apply_insertions(chars, editable, irate, rng)
    chars = enforce_length_preserve_sep(chars, TOTAL_LEN, rng)

    out = "".join(chars)
    if SEP not in out:
        raise RuntimeError("Separator lost; should not happen.")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("barcodes_in", help="Input: one structured barcode per line (16+GGGG+16)")
    ap.add_argument("reads_out", help="Output reads file")
    ap.add_argument("truth_out", help="Output truth file (barcode index per read)")
    ap.add_argument("--total-reads", type=int, required=True, help="Total number of reads to generate (Q)")
    ap.add_argument("--nave", type=float, default=4.0, help="Mean Poisson reads per selected barcode (Press-style)")
    ap.add_argument("--poisson", action="store_true",
                    help="Use Press-style Poisson burst sampling (default off unless set)")
    ap.add_argument("--srate", type=float, required=True, help="Substitution rate (effective; Press uses 4/3 correction internally)")
    ap.add_argument("--irate", type=float, required=True, help="Insertion rate")
    ap.add_argument("--drate", type=float, required=True, help="Deletion rate")
    ap.add_argument("--seed", type=int, default=None, help="Random seed")
    args = ap.parse_args()

    if args.total_reads <= 0:
        raise ValueError("--total-reads must be positive")

    rng = random.Random(args.seed)

    with open(args.barcodes_in, "r", encoding="utf-8") as f:
        barcodes = [line.strip().upper() for line in f if line.strip()]

    if not barcodes:
        print("No barcodes found.", file=sys.stderr)
        sys.exit(1)

    # Validate barcodes
    for i, b in enumerate(barcodes):
        if len(b) != TOTAL_LEN:
            raise ValueError(f"Line {i}: expected length {TOTAL_LEN}, got {len(b)}")
        if set(b) - set(BASES):
            raise ValueError(f"Line {i}: invalid characters {set(b) - set(BASES)}")
        if b[L_PART:L_PART + len(SEP)] != SEP:
            raise ValueError(f"Line {i}: expected separator '{SEP}' at position {L_PART}")

    Q = args.total_reads
    N = len(barcodes)

    with open(args.reads_out, "w", encoding="utf-8") as f_reads, \
         open(args.truth_out, "w", encoding="utf-8") as f_truth:

        q = 0
        while q < Q:
            idx = rng.randrange(N)
            n_to_make = 1
            if args.poisson:
                n_to_make = max(0, poisson_py(rng, args.nave))
                if n_to_make == 0:
                    continue
            # cap to reach exactly Q
            n_to_make = min(n_to_make, Q - q)

            bc = barcodes[idx]
            for _ in range(n_to_make):
                read = makeerrors_structured_press(bc, args.srate, args.irate, args.drate, rng)
                f_reads.write(read + "\n")
                f_truth.write(str(idx) + "\n")
            q += n_to_make


if __name__ == "__main__":
    main()
