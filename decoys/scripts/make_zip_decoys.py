#!/usr/bin/env python3
"""Structure-matched decoys for the printed Zipcodes layout.

A printed Zipcode is [column code (18 nt)][row code (18 nt)], where each axis
code is a 36-bit codeword encoded at 2 bits per base. Misassigned reads land
mostly on neighbouring rows or columns, so decoys are built from "virtual"
axis codes inserted between adjacent real codes:

  * virtual column code Dc between columns c and c+1 -> decoys Dc + row code
  * virtual row code    Dr between rows    r and r+1 -> decoys column code + Dr

A virtual code is chosen so that its bit distances to the real codes on both
sides match the distances observed between real codes at the same offsets,
and it is at least as far from every other code as real codes are from each
other. Virtual gaps are spread evenly along each axis, and the decoys are
split equally between the two axes.

Usage:
    python make_zip_decoys.py \
        --layout ../../data/1in1/layout_235x363_zipcodes_nosep.csv \
        --out-dir ../output/zipcodes
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

# Printed config: 2-bit pattern -> one-hot -> base (00->A, 01->C, 10->T, 11->G)
BASE_TO_BITS = {"A": (0, 0), "C": (0, 1), "T": (1, 0), "G": (1, 1)}
BITS_TO_BASE = {v: k for k, v in BASE_TO_BITS.items()}
NEAR_OFFSETS = (1, 2, 3, 4)


def dna_to_bits(seq):
    return np.array([b for base in seq for b in BASE_TO_BITS[base]], dtype=np.int8)


def bits_to_dna(bits):
    return "".join(BITS_TO_BASE[(int(bits[i]), int(bits[i + 1]))] for i in range(0, len(bits), 2))


def read_axis_codes(path, seg):
    cols, rows = {}, {}
    with open(path) as handle:
        for rec in csv.DictReader(handle):
            seq = rec["barcode"].strip().upper()
            cols.setdefault(int(rec["col"]), seq[:seg])
            rows.setdefault(int(rec["row"]), seq[seg:])
            if cols[int(rec["col"])] != seq[:seg] or rows[int(rec["row"])] != seq[seg:]:
                raise ValueError(f"Inconsistent axis code at row {rec['row']}, col {rec['col']}")
    col_codes = [cols[i] for i in range(len(cols))]
    row_codes = [rows[i] for i in range(len(rows))]
    return col_codes, row_codes


def axis_distance_stats(codes):
    """Bit Hamming distances between real axis codes, by position offset."""
    n = len(codes)
    dist = (codes[:, None, :] != codes[None, :, :]).sum(-1)
    stats = {}
    for off in NEAR_OFFSETS:
        d = np.array([dist[i, i + off] for i in range(n - off)])
        stats[off] = {"mean": float(d.mean()), "min": int(d.min()), "max": int(d.max())}
    far = dist[np.triu_indices(n, k=max(NEAR_OFFSETS) + 1)]
    stats["far"] = {"mean": float(far.mean()), "min": int(far.min())}
    return stats


def pick_gaps(n_codes, n_virtual, spacing="linspace"):
    """Gap indices g, where gap g lies between codes g and g+1.

    ``linspace``: unique rounded linspace on 0 .. n_codes-2 (manuscript decoys).
    ``bresenham``: integer Bresenham including both end gaps, so every window of
    real codes contains about the same number of virtual codes.
    """
    n_gaps = n_codes - 1
    if n_virtual < 1 or n_virtual > n_gaps:
        raise ValueError(f"Need 1..{n_gaps} virtual codes, got {n_virtual}")
    if spacing == "linspace":
        gaps = np.unique(np.round(np.linspace(0, n_codes - 2, n_virtual)).astype(int))
    elif spacing == "bresenham":
        if n_virtual == 1:
            gaps = np.array([0], dtype=int)
        else:
            gaps = np.array(
                [(i * (n_gaps - 1)) // (n_virtual - 1) for i in range(n_virtual)],
                dtype=int,
            )
    else:
        raise ValueError(f"Unknown spacing {spacing!r}")
    if len(np.unique(gaps)) != n_virtual:
        raise ValueError(f"Could not place {n_virtual} distinct gaps among {n_gaps}")
    return gaps


def make_virtual_code(codes, gap, stats, taken, rng, n_candidates):
    """Best random code for a virtual position between codes[gap] and codes[gap+1]."""
    n, n_bits = codes.shape
    near_min = stats[1]["min"]
    far_min = stats["far"]["min"]
    offsets = [stats[1]["min"] + k for k in range(stats[1]["max"] - stats[1]["min"] + 1)]

    # Real codes at virtual distance 1..4 on each side, with their target distances.
    near = []
    for off in NEAR_OFFSETS:
        for idx in (gap + 1 - off, gap + off):
            if 0 <= idx < n:
                near.append((idx, stats[off]["mean"]))
    near_idx = np.array([i for i, _ in near])
    near_target = np.array([t for _, t in near])
    far_idx = np.array([i for i in range(n) if i not in set(near_idx.tolist())])

    parents = codes[rng.choice([gap, gap + 1], size=n_candidates)]
    flips = rng.choice(offsets, size=n_candidates)
    cand = parents.copy()
    for k in range(n_candidates):
        pos = rng.choice(n_bits, size=flips[k], replace=False)
        cand[k, pos] ^= 1

    d_all = (cand[:, None, :] != codes[None, :, :]).sum(-1)
    ok = (d_all >= near_min).all(1)
    if len(far_idx):
        ok &= (d_all[:, far_idx] >= far_min).all(1)
    if taken:
        d_taken = (cand[:, None, :] != np.array(taken)[None, :, :]).sum(-1)
        ok &= (d_taken >= near_min).all(1)
    if not ok.any():
        raise RuntimeError(f"No valid virtual code found for gap {gap}; increase --candidates")

    score = ((d_all[:, near_idx] - near_target) ** 2).sum(1).astype(float)
    score[~ok] = np.inf
    best = int(np.argmin(score))
    return cand[best], d_all[best, near_idx].tolist(), near_target.tolist()


def build_axis(codes_dna, n_virtual, rng, n_candidates, label, spacing="linspace"):
    codes = np.array([dna_to_bits(s) for s in codes_dna], dtype=np.int8)
    for s in codes_dna:
        if bits_to_dna(dna_to_bits(s)) != s:
            raise ValueError(f"{label} code {s} does not round-trip through the bit mapping")
    stats = axis_distance_stats(codes)
    gaps = pick_gaps(len(codes), n_virtual, spacing=spacing)
    virtual, report, taken = [], [], []
    for i, gap in enumerate(gaps, start=1):
        print(f"  {label} virtual {i}/{n_virtual} at gap {int(gap)}", flush=True)
        code, dists, targets = make_virtual_code(codes, int(gap), stats, taken, rng, n_candidates)
        taken.append(code)
        virtual.append((int(gap), bits_to_dna(code)))
        report.append({"gap": int(gap), "distances_to_near_real": dists, "targets": targets})
    return virtual, stats, report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--layout", required=True, help="Printed layout CSV with columns row,col,barcode")
    parser.add_argument("--segment-length", type=int, default=18)
    parser.add_argument("--fraction", type=float, default=0.5, help="Decoys per real barcode (ignored if --n-virtual-* is set)")
    parser.add_argument("--n-virtual-columns", type=int, default=None, help="Override number of virtual column codes")
    parser.add_argument("--n-virtual-rows", type=int, default=None, help="Override number of virtual row codes")
    parser.add_argument("--spacing", choices=("linspace", "bresenham"), default="linspace")
    parser.add_argument("--candidates", type=int, default=4000, help="Random candidates per virtual code")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    col_codes, row_codes = read_axis_codes(args.layout, args.segment_length)
    n_cols, n_rows = len(col_codes), len(row_codes)
    n_real = n_cols * n_rows
    if args.n_virtual_columns is None or args.n_virtual_rows is None:
        if args.n_virtual_columns is not None or args.n_virtual_rows is not None:
            raise SystemExit("Set both --n-virtual-columns and --n-virtual-rows, or neither")
        half = args.fraction * n_real / 2
        n_virtual_cols = int(round(half / n_rows))
        n_virtual_rows = int(round(half / n_cols))
    else:
        n_virtual_cols = args.n_virtual_columns
        n_virtual_rows = args.n_virtual_rows

    virtual_cols, col_stats, col_report = build_axis(
        col_codes, n_virtual_cols, rng, args.candidates, "column", spacing=args.spacing)
    virtual_rows, row_stats, row_report = build_axis(
        row_codes, n_virtual_rows, rng, args.candidates, "row", spacing=args.spacing)

    real_set = {c + r for c in col_codes for r in row_codes}
    decoys = []
    for gap, code in virtual_cols:
        for row, rcode in enumerate(row_codes):
            decoys.append({"axis": "column", "gap": gap, "row": row, "col": gap + 0.5, "barcode": code + rcode})
    for gap, code in virtual_rows:
        for col, ccode in enumerate(col_codes):
            decoys.append({"axis": "row", "gap": gap, "row": gap + 0.5, "col": col, "barcode": ccode + code})

    seqs = [d["barcode"] for d in decoys]
    if len(set(seqs)) != len(seqs):
        raise RuntimeError("Duplicate decoy barcodes")
    if real_set & set(seqs):
        raise RuntimeError("Decoy barcodes overlap the real layout")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "zip_decoys.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(decoys[0].keys()))
        writer.writeheader()
        writer.writerows(decoys)
    with open(out_dir / "zip_decoys.txt", "w") as handle:
        handle.writelines(s + "\n" for s in seqs)
    with open(out_dir / "zip_whitelist_real_plus_decoys.txt", "w") as handle:
        handle.writelines(s + "\n" for s in sorted(real_set))
        handle.writelines(s + "\n" for s in seqs)
    with open(out_dir / "zip_virtual_axis_codes.json", "w") as handle:
        json.dump({"columns": col_report, "rows": row_report,
                   "virtual_columns": virtual_cols, "virtual_rows": virtual_rows}, handle, indent=2)

    summary = {
        "layout": str(args.layout),
        "seed": args.seed,
        "candidates_per_virtual_code": args.candidates,
        "n_real": len(real_set),
        "spacing": args.spacing,
        "n_virtual_columns": len(virtual_cols),
        "n_virtual_rows": len(virtual_rows),
        "n_decoys": len(decoys),
        "decoy_to_real": len(decoys) / len(real_set),
        "decoy_fraction_of_whitelist": len(decoys) / (len(decoys) + len(real_set)),
        "phi_design_same_column": len(virtual_rows) / (n_rows - 1 + len(virtual_rows)),
        "phi_design_same_row": len(virtual_cols) / (n_cols - 1 + len(virtual_cols)),
        "virtual_column_gaps": [g for g, _ in virtual_cols],
        "virtual_row_gaps": [g for g, _ in virtual_rows],
        "real_column_code_bit_distances": col_stats,
        "real_row_code_bit_distances": row_stats,
    }
    with open(out_dir / "zip_decoys_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if not k.startswith("real_")}, indent=2))
    print("column code distances:", json.dumps(col_stats))
    print("row code distances:   ", json.dumps(row_stats))


if __name__ == "__main__":
    main()
