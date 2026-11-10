#!/usr/bin/env python3
"""
Plot Levenshtein distance histograms for barcode sets.

Behavior:
- Reads barcode files for dual, quadtree, random, and zip sets.
- Samples a subset of barcodes from each set (default: 1000).
- Computes all-to-all pairwise Levenshtein distances within each subset.
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
class DistanceSummary:
    n_barcodes: int
    n_pairs: int
    min_distance: int
    max_distance: int
    mean_distance: float
    median_distance: float


def levenshtein_distance(a: str, b: str) -> int:
    """Compute Levenshtein distance using a 2-row dynamic programming table."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)

    # Keep DP width small by making b the shorter string when possible.
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

    # Try direct DNA token extraction first (works for plain txt and mixed columns).
    dna_tokens = DNA_TOKEN.findall(line)
    if dna_tokens:
        # Pick the longest DNA-like token, which is typically the barcode.
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


def all_to_all_levenshtein(barcodes: List[str]) -> List[int]:
    """Compute all pairwise Levenshtein distances within one barcode set."""
    distances: List[int] = []
    distance_fn = levenshtein_distance_fast or levenshtein_distance
    for s1, s2 in combinations(barcodes, 2):
        distances.append(distance_fn(s1, s2))
    return distances


def summarize_distances(n_barcodes: int, distances: List[int]) -> DistanceSummary:
    """Summarize pairwise distance distribution statistics."""
    return DistanceSummary(
        n_barcodes=n_barcodes,
        n_pairs=len(distances),
        min_distance=min(distances),
        max_distance=max(distances),
        mean_distance=float(mean(distances)),
        median_distance=float(median(distances)),
    )


def make_histograms(
    distance_map: Dict[str, List[int]],
    summary_map: Dict[str, DistanceSummary],
    output_png: Path,
) -> None:
    """Create a 2x2 histogram figure for the four barcode sets."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    axes = axes.flatten()

    for ax, (set_name, distances) in zip(axes, distance_map.items()):
        bins = range(min(distances), max(distances) + 2)
        ax.hist(distances, bins=bins, edgecolor="black", alpha=0.85)
        stats = summary_map[set_name]
        ax.set_title(
            f"{set_name}\n"
            f"n={stats.n_barcodes}, pairs={stats.n_pairs:,}, "
            f"mean={stats.mean_distance:.2f}",
            fontsize=11,
        )
        ax.set_xlabel("Levenshtein distance")
        ax.set_ylabel("Count")
        ax.grid(alpha=0.25, linestyle="--", linewidth=0.5)

    fig.suptitle("Levenshtein Distance Distributions (Subset All-to-All)", fontsize=14)
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Levenshtein distance histograms from barcode sets.")
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

    distance_map: Dict[str, List[int]] = {}
    summary_map: Dict[str, DistanceSummary] = {}

    for set_name, path in files.items():
        print(f"Loading {set_name}: {path}")
        barcodes = load_barcodes(path)
        subset = sample_barcodes(barcodes, args.subset_size, rng)
        print(
            f"  Parsed {len(barcodes):,} barcodes, using subset {len(subset):,}; "
            f"computing all-to-all distances..."
        )
        distances = all_to_all_levenshtein(subset)
        distance_map[set_name] = distances
        summary_map[set_name] = summarize_distances(len(subset), distances)

    png_path = args.output_dir / "levenshtein_distance_histograms_subset.png"
    make_histograms(distance_map, summary_map, png_path)

    summary_path = args.output_dir / "levenshtein_distance_summary.json"
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
