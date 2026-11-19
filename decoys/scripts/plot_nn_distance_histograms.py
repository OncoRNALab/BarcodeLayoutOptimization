#!/usr/bin/env python3
"""Nearest-neighbour distance histograms: real→real vs decoy→real.

For equal-length barcodes, nucleotide Hamming distance equals the
substitution-only edit distance. We sample ``--sample`` sequences and plot
both distributions so structure-matched decoys can be checked visually.

Usage:
    python plot_nn_distance_histograms.py --out-dir ../figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}

ARRAYS = [
    {
        "label": "Dual shared_pool minH6",
        "short": "dual_shared_pool_g10_minH6",
        "real": "real_dual_shared_pool_g10_minH6.txt",
        "decoys": "dual_shared_pool_g10_minH6/dual_decoys.txt",
    },
    {
        "label": "Dual fixed_bar2 8×8 NN≈8",
        "short": "dual_fixed_bar2_8x8_g10_nn8",
        "real": "real_dual_fixed_bar2_8x8_g10_nn8.txt",
        "decoys": "dual_fixed_bar2_8x8_g10_nn8/dual_decoys.txt",
    },
    {
        "label": "Zipcodes flip G10 NN≈8",
        "short": "zipcodes_flip_g10_nn8",
        "real": "real_zipcodes_flip_g10_nn8.txt",
        "decoys": "zipcodes_flip_g10_nn8/zip_decoys.txt",
    },
    {
        "label": "HSDB colored NN≈8",
        "short": "hsdb_g10_nn8_colored",
        "real": "real_hsdb_g10_nn8_colored.txt",
        "decoys": "hsdb_g10_nn8_colored/dual_decoys.txt",
    },
]


def read_seqs(path: Path):
    with open(path) as handle:
        return [line.strip().upper() for line in handle if line.strip()]


def encode(seqs):
    return np.array([[BASE_CODE[b] for b in s] for s in seqs], dtype=np.int8)


def nearest_hamming(queries, ref, exclude_self: bool = False):
    """NN Hamming distances for each query against ref (numpy, chunked)."""
    out = np.empty(len(queries), dtype=np.int32)
    chunk = 256
    for i in range(0, len(queries), chunk):
        q = queries[i : i + chunk]
        # (n_q, n_ref)
        d = (q[:, None, :] != ref[None, :, :]).sum(axis=2)
        if exclude_self:
            # identical rows: treat distance 0 as non-neighbour
            d = np.where(d == 0, 10**6, d)
        out[i : i + len(q)] = d.min(axis=1)
    return out


def describe(values: np.ndarray) -> dict:
    return {
        "n": int(len(values)),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "min": int(values.min()),
        "max": int(values.max()),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-dir", type=Path, default=None)
    parser.add_argument("--output-dir-decoys", type=Path, default=None,
                        help="Root decoys/output directory")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--xmin", type=int, default=0,
                        help="Shared x-axis lower bound (inclusive)")
    parser.add_argument("--xmax", type=int, default=12,
                        help="Shared x-axis upper bound (inclusive)")
    args = parser.parse_args()
    if args.xmax <= args.xmin:
        raise SystemExit("--xmax must be > --xmin")
    shared_bins = np.arange(args.xmin, args.xmax + 2) - 0.5

    root = Path(__file__).resolve().parents[1]
    qc_dir = args.qc_dir or (root / "qc")
    decoy_out = args.output_dir_decoys or (root / "output")
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    results = {}

    # Collect distances first
    for spec in ARRAYS:
        real_seqs = read_seqs(qc_dir / spec["real"])
        decoy_seqs = read_seqs(decoy_out / spec["decoys"])
        real = encode(real_seqs)
        decoy = encode(decoy_seqs)

        r_idx = rng.choice(len(real), size=min(args.sample, len(real)), replace=False)
        d_idx = rng.choice(len(decoy), size=min(args.sample, len(decoy)), replace=False)

        real_nn = nearest_hamming(real[r_idx], real, exclude_self=True)
        decoy_nn = nearest_hamming(decoy[d_idx], real, exclude_self=False)

        results[spec["short"]] = {
            "label": spec["label"],
            "real_nn": real_nn,
            "decoy_nn": decoy_nn,
            "real_stats": describe(real_nn),
            "decoy_stats": describe(decoy_nn),
        }
        print(
            f"{spec['label']}: real→real mean={real_nn.mean():.2f} "
            f"[{real_nn.min()}-{real_nn.max()}]; "
            f"decoy→real mean={decoy_nn.mean():.2f} "
            f"[{decoy_nn.min()}-{decoy_nn.max()}]"
        )

    n_arrays = len(ARRAYS)
    n_cols = 2
    n_rows = (n_arrays + n_cols - 1) // n_cols

    # --- Figure 1: overlay panels (one per array) ---
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(9.5, 3.6 * n_rows), sharex=True, sharey=False
    )
    axes_flat = np.atleast_1d(axes).ravel()
    for ax, spec in zip(axes_flat, ARRAYS):
        key = spec["short"]
        r = results[key]["real_nn"]
        d = results[key]["decoy_nn"]
        ax.hist(
            r, bins=shared_bins, density=True, alpha=0.55, color="#2c7bb6",
            edgecolor="white", linewidth=0.4, label="real → nearest real",
        )
        ax.hist(
            d, bins=shared_bins, density=True, alpha=0.55, color="#d7191c",
            edgecolor="white", linewidth=0.4, label="decoy → nearest real",
        )
        ax.set_xlim(args.xmin, args.xmax)
        ax.set_xticks(range(args.xmin, args.xmax + 1, 2))
        ax.set_title(spec["label"], fontsize=11)
        ax.set_xlabel("Nearest-neighbour Hamming distance (nt)")
        ax.set_ylabel("Density")
        ax.legend(fontsize=8, frameon=False, loc="upper right")
        rs, ds = results[key]["real_stats"], results[key]["decoy_stats"]
        ax.text(
            0.02, 0.98,
            f"real μ={rs['mean']:.1f}\ndecoy μ={ds['mean']:.1f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            family="monospace",
        )
    for ax in axes_flat[n_arrays:]:
        ax.set_visible(False)
    fig.suptitle(
        f"NN distance of {args.sample} sampled barcodes "
        "(Hamming = equal-length substitution edit distance)",
        fontsize=12,
    )
    fig.tight_layout()
    overlay_path = out_dir / "nn_distance_histograms_overlay.png"
    fig.savefig(overlay_path, dpi=160, bbox_inches="tight")
    fig.savefig(out_dir / "nn_distance_histograms_overlay.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", overlay_path)

    # --- Figure 2: n×2 separate panels (real | decoy) ---
    fig2, axes2 = plt.subplots(
        n_arrays, 2, figsize=(8.5, 2.6 * n_arrays), sharex=True, sharey=False
    )
    if n_arrays == 1:
        axes2 = np.array([axes2])
    for row, spec in enumerate(ARRAYS):
        key = spec["short"]
        r = results[key]["real_nn"]
        d = results[key]["decoy_nn"]
        axes2[row, 0].hist(
            r, bins=shared_bins, color="#2c7bb6", edgecolor="white", linewidth=0.4,
        )
        axes2[row, 1].hist(
            d, bins=shared_bins, color="#d7191c", edgecolor="white", linewidth=0.4,
        )
        for col in (0, 1):
            axes2[row, col].set_xlim(args.xmin, args.xmax)
            axes2[row, col].set_xticks(range(args.xmin, args.xmax + 1, 2))
        axes2[row, 0].set_ylabel(spec["label"].replace(" ", "\n"), fontsize=9)
        if row == 0:
            axes2[row, 0].set_title("Real → nearest other real")
            axes2[row, 1].set_title("Decoy → nearest real")
        if row == n_arrays - 1:
            axes2[row, 0].set_xlabel("Hamming distance (nt)")
            axes2[row, 1].set_xlabel("Hamming distance (nt)")
    fig2.suptitle(
        f"NN Hamming distances (sample={args.sample}, seed={args.seed})",
        fontsize=12,
    )
    fig2.tight_layout()
    split_path = out_dir / "nn_distance_histograms_split.png"
    fig2.savefig(split_path, dpi=160, bbox_inches="tight")
    fig2.savefig(out_dir / "nn_distance_histograms_split.pdf", bbox_inches="tight")
    plt.close(fig2)
    print("wrote", split_path)

    # Save numeric summary
    import json
    summary = {
        k: {"label": v["label"], "real": v["real_stats"], "decoy": v["decoy_stats"]}
        for k, v in results.items()
    }
    summary_path = out_dir / "nn_distance_histograms_summary.json"
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print("wrote", summary_path)


if __name__ == "__main__":
    main()
