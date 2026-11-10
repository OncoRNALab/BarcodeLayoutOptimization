#!/usr/bin/env python3
"""
Plot histograms of minimum Levenshtein distance per barcode for barcode sets.

Behavior:
- Reads barcode files for dual, quadtree, random, and zip sets.
- Samples a subset of barcodes from each set (default: 1000).
- For each barcode in the subset, computes its minimum Levenshtein distance
  to all other barcodes in the same subset.
- Saves a single 2x2 histogram figure (one panel per set).

This script is intentionally standalone and does not run automatically.
"""

from __future__ import annotations

import argparse
import json
import random
import re
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from statistics import mean, median
from typing import Dict, List

import matplotlib.pyplot as plt

try:
    from Levenshtein import distance as levenshtein_distance_fast  # type: ignore[import-not-found]
except ImportError:
    levenshtein_distance_fast = None


DEFAULT_FILES: Dict[str, str] = {
    "dual_barcodes": "/data/gent/vo/000/gvo00027/projects/Spatial/Spatial_experiments/BCSEQ_26018/Analysis/barcodes/layout_235x363_dual_4x4.txt",
    "quadtree_barcodes": "/data/gent/vo/000/gvo00027/projects/Spatial/Spatial_experiments/BCSEQ_26018/Analysis/barcodes/layout_235x363_quadtree_8x3_12nt.txt",
    "random_barcodes": "/data/gent/vo/000/gvo00027/projects/Spatial/Spatial_experiments/BCSEQ_26018/Analysis/barcodes/layout_235x363_randombarcodes_opt.txt",
    "zip_barcodes": "/data/gent/vo/000/gvo00027/projects/Spatial/Spatial_experiments/BCSEQ_26018/Analysis/barcodes/layout_235x363_zipcodes_nosep.txt",
}

DEFAULT_OUTPUT_DIR = Path(
    "/user/gent/446/vsc44685/ScratchVO_dir/layoutOPT_results/printed_layouts/plots"
)

DNA_TOKEN = re.compile(r"[ACGTNacgtn]+")


@dataclass
class MinDistanceSummary:
    n_barcodes: int
    min_of_mins: int
    max_of_mins: int
    mean_of_mins: float
    median_of_mins: float


def levenshtein_distance(a: str, b: str) -> int:
    """Compute Levenshtein distance using a 2-row dynamic programming table."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)

    if len(a) < len(b):
        a, b = b, a

    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        for j, cb in enumerate(b, start=1):
            insert_cost = current[j - 1] + 1
            delete_cost = previous[j] + 1
            sub_cost = previous[j - 1] + (0 if ca == cb else 1)
            current.append(min(insert_cost, delete_cost, sub_cost))
        previous = current
    return previous[-1]


def extract_barcode_from_line(line: str) -> str:
    """Extract barcode sequence from a text/CSV-like line."""
    line = line.strip()
    if not line:
        return ""

    dna_tokens = DNA_TOKEN.findall(line)
    if dna_tokens:
        return max(dna_tokens, key=len).upper()

    return ""


def load_barcodes(path: Path) -> List[str]:
    """Load barcodes from file, one extracted barcode per non-empty line."""
    barcodes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            seq = extract_barcode_from_line(raw_line)
            if seq:
                barcodes.append(seq)

    if not barcodes:
        raise ValueError(f"No barcode sequences could be parsed from: {path}")
    return barcodes


def sample_barcodes(barcodes: List[str], subset_size: int, rng: random.Random) -> List[str]:
    """Sample barcodes without replacement; use all if set is smaller than subset."""
    n = len(barcodes)
    if subset_size >= n:
        return list(barcodes)
    return rng.sample(barcodes, subset_size)


def min_distance_per_barcode(barcodes: List[str]) -> List[int]:
    """Compute nearest-neighbor Levenshtein distance for each barcode.

    Uses pairwise combinations and updates minima for both members of each pair.
    """
    n = len(barcodes)
    if n < 2:
        raise ValueError("Need at least 2 barcodes to compute nearest-neighbor distance.")

    distance_fn = levenshtein_distance_fast or levenshtein_distance
    mins = [10**9] * n

    for i, j in combinations(range(n), 2):
        d = distance_fn(barcodes[i], barcodes[j])
        if d < mins[i]:
            mins[i] = d
        if d < mins[j]:
            mins[j] = d

    return mins


def summarize_min_distances(min_distances: List[int]) -> MinDistanceSummary:
    """Summarize nearest-neighbor distance distribution statistics."""
    return MinDistanceSummary(
        n_barcodes=len(min_distances),
        min_of_mins=min(min_distances),
        max_of_mins=max(min_distances),
        mean_of_mins=float(mean(min_distances)),
        median_of_mins=float(median(min_distances)),
    )


def make_histograms(
    min_distance_map: Dict[str, List[int]],
    summary_map: Dict[str, MinDistanceSummary],
    output_png: Path,
) -> None:
    """Create a 2x2 histogram figure for the four barcode sets."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    axes = axes.flatten()

    for ax, (set_name, min_distances) in zip(axes, min_distance_map.items()):
        bins = range(min(min_distances), max(min_distances) + 2)
        ax.hist(min_distances, bins=bins, edgecolor="black", alpha=0.85)
        stats = summary_map[set_name]
        ax.set_title(
            f"{set_name}\n"
            f"n={stats.n_barcodes}, min={stats.min_of_mins}, max={stats.max_of_mins}, "
            f"mean={stats.mean_of_mins:.2f}",
            fontsize=11,
        )
        ax.set_xlabel("Minimum Levenshtein distance to any other barcode")
        ax.set_ylabel("Count")
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)

    fig.suptitle("Nearest-Neighbor Levenshtein Distance Distributions", fontsize=14)
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot histograms of per-barcode minimum Levenshtein distances.")
    parser.add_argument(
        "--subset-size",
        type=int,
        default=1000,
        help="Number of barcodes sampled per set (default: 1000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible subset sampling (default: 42).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where outputs are saved.",
    )
    parser.add_argument(
        "--dual-file",
        type=Path,
        default=Path(DEFAULT_FILES["dual_barcodes"]),
        help="Path to dual barcodes file.",
    )
    parser.add_argument(
        "--quadtree-file",
        type=Path,
        default=Path(DEFAULT_FILES["quadtree_barcodes"]),
        help="Path to quadtree barcodes file.",
    )
    parser.add_argument(
        "--random-file",
        type=Path,
        default=Path(DEFAULT_FILES["random_barcodes"]),
        help="Path to random barcodes file.",
    )
    parser.add_argument(
        "--zip-file",
        type=Path,
        default=Path(DEFAULT_FILES["zip_barcodes"]),
        help="Path to zip barcodes file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    backend = "python-Levenshtein" if levenshtein_distance_fast else "pure-python"
    print(f"Distance backend: {backend}")
    if not levenshtein_distance_fast:
        print("Tip: install python-Levenshtein for much faster pairwise distance computation.")

    files: Dict[str, Path] = {
        "dual_barcodes": args.dual_file,
        "quadtree_barcodes": args.quadtree_file,
        "random_barcodes": args.random_file,
        "zip_barcodes": args.zip_file,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)

    min_distance_map: Dict[str, List[int]] = {}
    summary_map: Dict[str, MinDistanceSummary] = {}

    for set_name, path in files.items():
        print(f"Loading {set_name}: {path}")
        barcodes = load_barcodes(path)
        subset = sample_barcodes(barcodes, args.subset_size, rng)
        print(
            f"  Parsed {len(barcodes):,} barcodes, using subset {len(subset):,}; "
            f"computing nearest-neighbor distances..."
        )
        min_distances = min_distance_per_barcode(subset)
        min_distance_map[set_name] = min_distances
        summary_map[set_name] = summarize_min_distances(min_distances)

    png_path = args.output_dir / "levenshtein_min_distance_histograms_subset.png"
    make_histograms(min_distance_map, summary_map, png_path)

    summary_path = args.output_dir / "levenshtein_min_distance_summary.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(
            {k: vars(v) for k, v in summary_map.items()},
            f,
            indent=2,
        )

    print("\nDone.")
    print(f"Histogram figure: {png_path}")
    print(f"Summary JSON:    {summary_path}")


if __name__ == "__main__":
    main()
