#!/usr/bin/env python3
"""QC of a decoy set against the real barcodes of a printed layout.

Reports, using nucleotide Hamming distance as a fast proxy for edit distance:
  * size, uniqueness and overlap with the real barcodes
  * distance from each decoy to its nearest real barcode, compared with the
    distance from each real barcode to its nearest other real barcode
  * local decoy share: fraction of decoys among the K nearest whitelist
    entries of a real barcode, compared with the global decoy share. The
    decoy-based precision estimator is unbiased only if the two agree.
  * G content and longest homopolymer

Usage:
    python check_decoys.py --real REAL.txt --decoys NEW.txt --label new \
        --decoys OLD.txt --label old --match-size --out ../qc/dual_qc.json
"""

import argparse
import json
import re

import numpy as np

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_seqs(path):
    with open(path) as handle:
        return [line.strip().upper() for line in handle if line.strip()]


def encode(seqs):
    return np.array([[BASE_CODE[b] for b in s] for s in seqs], dtype=np.int8)


def nearest(queries, ref, exclude_self=False):
    out = np.empty(len(queries), dtype=int)
    for i, q in enumerate(queries):
        d = (ref != q).sum(1)
        if exclude_self:
            d[d == 0] = 10 ** 6
        out[i] = d.min()
    return out


def local_share(queries, whitelist, is_decoy, k):
    shares, tie_shares = [], []
    for q in queries:
        d = (whitelist != q).sum(1)
        d[d == 0] = 10 ** 6
        idx = np.argpartition(d, k)[:k]
        shares.append(is_decoy[idx].mean())
        tie_shares.append(is_decoy[d == d.min()].mean())
    return float(np.mean(shares)), float(np.mean(tie_shares))


def describe(values):
    values = np.asarray(values)
    return {
        "min": float(values.min()), "p05": float(np.percentile(values, 5)),
        "median": float(np.median(values)), "mean": float(values.mean()),
        "p95": float(np.percentile(values, 95)), "max": float(values.max()),
    }


def composition(seqs):
    g = [s.count("G") / len(s) for s in seqs]
    hp = [max(len(m.group(0)) for m in re.finditer(r"(.)\1*", s)) for s in seqs]
    return {"g_fraction": describe(g), "max_homopolymer": describe(hp)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--real", required=True)
    parser.add_argument("--decoys", action="append", required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--match-size", action="store_true",
                        help="Subsample every decoy set to the size of the first one")
    parser.add_argument("--sample", type=int, default=2000)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if len(args.decoys) != len(args.label):
        raise SystemExit("Give one --label per --decoys")

    rng = np.random.default_rng(args.seed)
    real_seqs = read_seqs(args.real)
    real = encode(real_seqs)
    q_idx = rng.choice(len(real), size=min(args.sample, len(real)), replace=False)
    real_nn = nearest(real[q_idx], real, exclude_self=True)
    report = {"real": {"n": len(real_seqs), "n_unique": len(set(real_seqs)),
                       "nearest_other_real": describe(real_nn),
                       "composition": composition(real_seqs)}}
    print(f"real: n={len(real_seqs)} nearest other real: {report['real']['nearest_other_real']}")

    target_n = None
    for path, label in zip(args.decoys, args.label):
        seqs = read_seqs(path)
        n_file = len(seqs)
        if args.match_size:
            if target_n is None:
                target_n = len(seqs)
            elif len(seqs) > target_n:
                seqs = [seqs[i] for i in sorted(rng.choice(len(seqs), size=target_n, replace=False))]
        decoy = encode(seqs)
        d_idx = rng.choice(len(decoy), size=min(args.sample, len(decoy)), replace=False)
        decoy_nn = nearest(decoy[d_idx], real)

        whitelist = np.vstack([real, decoy])
        is_decoy = np.r_[np.zeros(len(real), bool), np.ones(len(decoy), bool)]
        share_k, share_tie = local_share(real[q_idx], whitelist, is_decoy, args.k)
        global_share = len(decoy) / len(whitelist)

        entry = {
            "file": path, "n_in_file": n_file, "n_used": len(seqs),
            "n_unique": len(set(seqs)), "overlap_with_real": len(set(seqs) & set(real_seqs)),
            "global_decoy_share": global_share,
            f"local_decoy_share_k{args.k}": share_k,
            "local_decoy_share_nearest_ties": share_tie,
            "nearest_real": describe(decoy_nn),
            "composition": composition(seqs),
        }
        report[label] = entry
        print(f"\n{label}: used {len(seqs)} of {n_file}, unique {entry['n_unique']}, "
              f"overlap {entry['overlap_with_real']}")
        print(f"  decoy -> nearest real: {entry['nearest_real']}")
        print(f"  decoy share: global {global_share:.3f} | among {args.k} nearest of a real barcode "
              f"{share_k:.3f} | among nearest ties {share_tie:.3f}")
        print(f"  G fraction mean {entry['composition']['g_fraction']['mean']:.3f} "
              f"(real {report['real']['composition']['g_fraction']['mean']:.3f}); "
              f"max homopolymer max {entry['composition']['max_homopolymer']['max']:.0f} "
              f"(real {report['real']['composition']['max_homopolymer']['max']:.0f})")

    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=2)


if __name__ == "__main__":
    main()
