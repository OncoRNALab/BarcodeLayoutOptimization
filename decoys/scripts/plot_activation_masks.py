#!/usr/bin/env python3
"""QC: photolithographic activation masks for the four NN≈8 layouts.

Uses the same leftmost-illuminated rule as the manuscript figures
(ACGT cycle order). Early cycles that act on shared sequence prefixes
should show the expected spatial structure:

  * Dual shared_pool minH6 (4×4): block tiles
  * Dual fixed_bar2 (8×8): larger block tiles
  * Zipcodes: full columns (shared column code = first 18 nt)
  * HSDB: block-columns (shared Zipcode X half of bar_1 = first 12 nt)

Usage:
    python plot_activation_masks.py --cycle 1 --out-dir ../figures
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]

ARRAYS = [
    {
        "label": "Dual shared_pool minH6 (4×4)",
        "path": REPO / "Dual_barcodes/layout_235x363_dual_shared_pool_g10_minH6.csv",
        "seq_key": "combined",
        "note": "expect 4×4 blocks",
    },
    {
        "label": "Dual fixed_bar2 (8×8)",
        "path": REPO / "Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10_nn8.csv",
        "seq_key": "combined",
        "note": "expect 8×8 blocks",
    },
    {
        "label": "Zipcodes flip",
        "path": REPO
        / "Zipcodes/output_235x363_flip_g10_h5_d12-16_far10/layout_235x363_zipcodes_nosep.csv",
        "seq_key": "barcode",
        "note": "expect full columns",
    },
    {
        "label": "HSDB colored",
        "path": REPO
        / "HSDB/outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored/layout_235x363_hsdb.csv",
        "seq_key": "combined",
        "note": "expect block-columns (width 4)",
    },
]


def leftmost_illuminated(seq: str, cycle: int, order: str = "ACGT") -> bool:
    """True if this spot is deprotected in synthesis cycle ``cycle`` (0-based)."""
    c = -1
    n = len(order)
    for base in seq:
        c += 1
        while order[c % n] != base:
            c += 1
        if c == cycle:
            return True
        if c > cycle:
            return False
    return False


def load_layout(path: Path, seq_key: str):
    rows, cols, seqs = [], [], []
    with open(path) as handle:
        reader = csv.DictReader(handle)
        for rec in reader:
            rows.append(int(rec["row"]))
            cols.append(int(rec["col"]))
            seqs.append(rec[seq_key].strip().upper())
    return np.asarray(rows), np.asarray(cols), seqs


def mask_grid(rows, cols, seqs, cycle: int) -> np.ndarray:
    n_rows = int(rows.max()) + 1
    n_cols = int(cols.max()) + 1
    g = np.zeros((n_rows, n_cols), dtype=np.uint8)
    for r, c, seq in zip(rows, cols, seqs):
        if leftmost_illuminated(seq, cycle):
            g[r, c] = 1
    return g


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cycle",
        type=int,
        default=1,
        help="1-based synthesis cycle to display (default 1 = first ACGT cycle)",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.cycle < 1:
        raise SystemExit("--cycle must be ≥ 1 (1-based)")
    cycle0 = args.cycle - 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, spec in zip(axes.ravel(), ARRAYS):
        print(f"mask {spec['label']} cycle {args.cycle}…", flush=True)
        rows, cols, seqs = load_layout(spec["path"], spec["seq_key"])
        g = mask_grid(rows, cols, seqs, cycle0)
        ax.imshow(g, cmap="gray", interpolation="nearest", aspect="equal", vmin=0, vmax=1)
        ax.set_title(f"{spec['label']}\n({spec['note']})", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.6)
        # Save individual mask too
        short = (
            spec["label"]
            .lower()
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("×", "x")
            .replace(",", "")
        )
        np.save(args.out_dir / f"activation_mask_{short}_cycle{args.cycle}.npy", g)

    fig.suptitle(
        f"Activation masks — synthesis cycle {args.cycle} (ACGT order; white = deprotected)",
        fontsize=13,
    )
    fig.tight_layout()
    out_png = args.out_dir / f"activation_masks_cycle{args.cycle}.png"
    out_pdf = args.out_dir / f"activation_masks_cycle{args.cycle}.pdf"
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out_png)

    # Zoomed Dual panels: 4×4 / 8×8 tiles are hard to see at full-array scale
    figz, axz = plt.subplots(1, 2, figsize=(10, 5))
    zoom = 64  # top-left crop
    for ax, spec in zip(axz, ARRAYS[:2]):
        rows, cols, seqs = load_layout(spec["path"], spec["seq_key"])
        g = mask_grid(rows, cols, seqs, cycle0)
        ax.imshow(
            g[:zoom, :zoom],
            cmap="gray",
            interpolation="nearest",
            aspect="equal",
            vmin=0,
            vmax=1,
        )
        ax.set_title(f"{spec['label']} — zoom {zoom}×{zoom}", fontsize=11)
        ax.set_xticks(range(0, zoom + 1, 8))
        ax.set_yticks(range(0, zoom + 1, 8))
        ax.grid(True, color="#666666", lw=0.4, alpha=0.5)
    figz.suptitle(
        f"Dual activation masks (zoom) — cycle {args.cycle}",
        fontsize=12,
    )
    figz.tight_layout()
    zoom_png = args.out_dir / f"activation_masks_dual_zoom_cycle{args.cycle}.png"
    figz.savefig(zoom_png, dpi=160, bbox_inches="tight")
    figz.savefig(
        args.out_dir / f"activation_masks_dual_zoom_cycle{args.cycle}.pdf",
        bbox_inches="tight",
    )
    plt.close(figz)
    print("wrote", zoom_png)


if __name__ == "__main__":
    main()
