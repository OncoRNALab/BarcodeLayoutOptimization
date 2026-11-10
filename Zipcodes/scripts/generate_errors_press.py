#!/usr/bin/env python3
"""
Press-style error simulation for structured barcodes:
  LEFT(16) + SEP("GGGG") + RIGHT(16)  => length 36

Edits are applied ONLY outside the separator (no subs/ins/del within SEP).
Separator may shift due to indels in LEFT/RIGHT.
Final output is forced to length 36 (truncate/pad), while preserving SEP.

Implements a close variant of Press's makeerrors():
- substitutions: Binomial(n_edit, (4/3)*srate), assign random bases at sampled indices
- deletions:    Binomial(n_edit, drate), delete sampled indices
- insertions:   Binomial(len(seq_after_del), irate), insert random bases at sampled positions
- pad/truncate to original length (36)

Usage:
  python press_structured_sim.py barcodes.txt reads.txt truth.txt \
      --srate 0.0333 --irate 0.0333 --drate 0.0333 --reads-per-barcode 5 --seed 1
"""

import argparse
import random
import sys
from typing import List, Tuple

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


def sample_indices_with_replacement(rng: random.Random, high: int, k: int) -> List[int]:
    # Like numpy randint(low=0, high=high, size=k)
    return [rng.randrange(high) for _ in range(k)]


def mutate_substitutions(chars: List[str], editable: List[bool], srate: float, rng: random.Random) -> None:
    """Apply Press-style substitutions on editable positions only."""
    n_edit = sum(editable)
    if n_edit == 0:
        return

    # Press uses Binomial(n, (4/3)*srate) to achieve effective srate
    p = min(1.0, (4.0 / 3.0) * srate)
    ns = rng.binomial(n_edit, p) if hasattr(rng, "binomial") else binomial_py(rng, n_edit, p)
    if ns <= 0:
        return

    # Map editable positions to their indices in chars
    edit_positions = [i for i, ok in enumerate(editable) if ok]
    # Sample indices WITH replacement among editable positions
    chosen = [edit_positions[rng.randrange(n_edit)] for _ in range(ns)]
    for idx in chosen:
        chars[idx] = rand_base(rng)  # Press assigns random base; can be same as original occasionally


def apply_deletions(chars: List[str], editable: List[bool], drate: float, rng: random.Random) -> Tuple[List[str], List[bool]]:
    """Apply Press-style deletions on editable positions only. Indices drawn with replacement."""
    n = len(chars)
    n_edit = sum(editable)
    if n_edit == 0:
        return chars, editable

    nd = rng.binomial(n_edit, drate) if hasattr(rng, "binomial") else binomial_py(rng, n_edit, drate)
    if nd <= 0:
        return chars, editable

    edit_positions = [i for i, ok in enumerate(editable) if ok]
    chosen = [edit_positions[rng.randrange(n_edit)] for _ in range(nd)]
    # Press/NumPy delete: duplicates effectively delete once; easiest is to delete unique indices
    to_delete = sorted(set(chosen), reverse=True)
    for idx in to_delete:
        del chars[idx]
        del editable[idx]
    return chars, editable


def apply_insertions(chars: List[str], editable: List[bool], irate: float, rng: random.Random) -> Tuple[List[str], List[bool]]:
    """
    Apply Press-style insertions on the post-deletion sequence.
    Insertions are allowed only at positions that are not inside SEP.
    We model insertion sites using the editable mask:
    - insertion at position i means "insert before i"
    - also allow insertion at end.
    """
    n = len(chars)
    # Define valid insertion slots: before any editable base, and before the first SEP base, and after SEP, etc.
    # Easiest rule: insertion is allowed at any boundary that is not *within* SEP.
    # We will forbid insertions between SEP bases by identifying SEP span in current chars.
    seq_str = "".join(chars)
    s0, s1 = find_sep_span(seq_str, SEP)

    # Slots are between characters: 0..n
    valid_slots = []
    for slot in range(n + 1):
        # slot between slot-1 and slot
        # It's "within SEP" if it falls strictly inside the SEP span boundaries.
        # SEP occupies indices [s0, s1-1]. Boundaries between bases are [0..n].
        # Slots inside SEP are s0+1 .. s1-1 (between SEP bases).
        if (s0 + 1) <= slot <= (s1 - 1):
            continue
        valid_slots.append(slot)

    n_slots = len(valid_slots)
    if n_slots == 0:
        return chars, editable

    ni = rng.binomial(n_slots, irate) if hasattr(rng, "binomial") else binomial_py(rng, n_slots, irate)
    if ni <= 0:
        return chars, editable

    chosen_slots = [valid_slots[rng.randrange(n_slots)] for _ in range(ni)]
    # Allow collisions => multiple insertions at same slot.
    # To apply correctly, count insertions per slot and insert in ascending slot order with offset tracking.
    from collections import Counter
    counts = Counter(chosen_slots)
    offset = 0
    for slot in sorted(counts.keys()):
        c = counts[slot]
        insert_at = slot + offset
        for _ in range(c):
            chars.insert(insert_at, rand_base(rng))
            editable.insert(insert_at, True)  # inserted bases are editable “error material”
            insert_at += 1
            offset += 1

    return chars, editable


def enforce_length_preserve_sep(chars: List[str], target_len: int, rng: random.Random) -> List[str]:
    """
    Press-style pad/truncate to target_len, but ensure SEP remains intact and not cut.
    If too long: trim from ends while keeping SEP fully present.
    If too short: pad with random bases (at end).
    """
    # Pad if short
    while len(chars) < target_len:
        chars.append(rand_base(rng))

    # Truncate if long, preserving SEP
    while len(chars) > target_len:
        seq = "".join(chars)
        s0 = seq.find(SEP)
        if s0 == -1:
            # Should never happen (SEP protected), fallback: truncate end
            chars = chars[:target_len]
            break
        s1 = s0 + len(SEP)

        # If SEP would not fit in first target_len, drop from left until it does
        if s1 > target_len:
            chars.pop(0)
            continue

        # Otherwise SEP is safe; trim from end
        chars.pop()

    return chars


def binomial_py(rng: random.Random, n: int, p: float) -> int:
    """Exact binomial via Bernoulli sum; OK for barcode lengths/slots ~ O(10^1-10^2)."""
    if p <= 0.0:
        return 0
    if p >= 1.0:
        return n
    return sum(1 for _ in range(n) if rng.random() < p)


def makeerrors_structured_press(barcode: str, srate: float, irate: float, drate: float, rng: random.Random) -> str:
    barcode = barcode.strip().upper()
    if len(barcode) != TOTAL_LEN:
        raise ValueError(f"Expected barcode length {TOTAL_LEN}, got {len(barcode)}")
    if barcode[L_PART:L_PART + len(SEP)] != SEP:
        raise ValueError(f"Barcode does not contain expected separator '{SEP}' at position {L_PART}")

    chars = list(barcode)

    # Editable mask: True for L and R positions, False for SEP positions
    editable = [True] * L_PART + [False] * len(SEP) + [True] * R_PART

    # Substitutions
    mutate_substitutions(chars, editable, srate, rng)

    # Deletions
    chars, editable = apply_deletions(chars, editable, drate, rng)

    # Insertions (rate applied to valid slots outside SEP interior)
    chars, editable = apply_insertions(chars, editable, irate, rng)

    # Pad/truncate to fixed length 36, preserving SEP
    chars = enforce_length_preserve_sep(chars, TOTAL_LEN, rng)

    out = "".join(chars)
    if SEP not in out:
        raise RuntimeError("Separator lost; should not happen with protected model.")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("barcodes_in", help="Input file: one structured barcode per line (16+GGGG+16)")
    ap.add_argument("reads_out", help="Output reads file")
    ap.add_argument("truth_out", help="Output truth file (0-based barcode index per read)")
    ap.add_argument("--srate", type=float, required=True, help="Effective substitution rate (Press-style corrected internally)")
    ap.add_argument("--irate", type=float, required=True, help="Insertion rate")
    ap.add_argument("--drate", type=float, required=True, help="Deletion rate")
    ap.add_argument("--reads-per-barcode", type=int, required=True, help="Reads to generate per barcode")
    ap.add_argument("--seed", type=int, default=None, help="Random seed")
    args = ap.parse_args()

    if args.reads_per_barcode <= 0:
        raise ValueError("--reads-per-barcode must be positive")

    rng = random.Random(args.seed)

    # Load barcodes
    with open(args.barcodes_in, "r", encoding="utf-8") as f:
        barcodes = [line.strip().upper() for line in f if line.strip()]

    if not barcodes:
        print("No barcodes found.", file=sys.stderr)
        sys.exit(1)

    # Validate
    for i, b in enumerate(barcodes):
        if len(b) != TOTAL_LEN:
            raise ValueError(f"Line {i}: expected length {TOTAL_LEN}, got {len(b)}")
        bad = set(b) - set(BASES)
        if bad:
            raise ValueError(f"Line {i}: invalid characters {bad}")
        if b[L_PART:L_PART + len(SEP)] != SEP:
            raise ValueError(f"Line {i}: expected separator '{SEP}' at position {L_PART}")

    with open(args.reads_out, "w", encoding="utf-8") as f_reads, \
         open(args.truth_out, "w", encoding="utf-8") as f_truth:
        for idx, bc in enumerate(barcodes):
            for _ in range(args.reads_per_barcode):
                read = makeerrors_structured_press(bc, args.srate, args.irate, args.drate, rng)
                f_reads.write(read + "\n")
                f_truth.write(str(idx) + "\n")


if __name__ == "__main__":
    main()
