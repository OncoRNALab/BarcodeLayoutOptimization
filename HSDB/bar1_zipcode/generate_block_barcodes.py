#!/usr/bin/env python3
"""
Phase 2 — Generate block barcodes (bar_1) for the HSDB design.

Runs the spatial_zipcodes axis-code generator on the block grid defined in
the config file, then combines X + Y into a bar_1 sequence per block. The
same script handles any array size by pointing it at a different config:

  118×182 array  →  block_zipcode_config.json           (30×46 blocks)
  235×363 array  →  block_zipcode_config_235x363.json   (59×91 blocks)
  253×363 array  →  block_zipcode_config_253x363.json   (64×91 blocks)

By default (and in the updated configs) axis codes are built in flip mode:
each candidate is the preceding code with a few bits flipped, and G content
and junction homopolymers are rejected inside the candidate loop. The old
post-hoc G→T / homopolymer rewrites are gone; they broke the X+Y structure
and the distance gradient.

Outputs (written to --out-dir):
  - bar1_barcodes_28nt.txt   : sequences in row-major block order
  - bar1_codebook.csv        : block_row, block_col, bar_1
  - bar1_qc_report.json      : QC and distance statistics
  - intermediate/axis_codes.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
ZIPCODES_DIR = REPO_ROOT / "Zipcodes"

if str(ZIPCODES_DIR) not in sys.path:
    sys.path.insert(0, str(ZIPCODES_DIR))

from spatial_zipcodes.axis_code_generator import generate_column_and_row_codes
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.config import load_config
from spatial_zipcodes.sequence_qc import (
    calculate_g_content,
    find_max_homopolymer,
    summarize_junction_homopolymers,
    terminal_run,
)

DEFAULT_BAR2 = REPO_ROOT / "HSDB" / "bar2_library" / "spot_barcodes_8nt.txt"


def dna_hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_spatial_distance_stats(
    barcodes: List[str],
    block_rows: int,
    block_cols: int,
) -> Dict[str, Any]:
    """DNA Hamming distances for adjacent blocks and blocks ≥5 steps apart."""
    adjacent_distances: List[int] = []
    far_distances: List[int] = []

    for br in range(block_rows):
        for bc in range(block_cols):
            idx = br * block_cols + bc
            seq = barcodes[idx]

            if bc + 1 < block_cols:
                adjacent_distances.append(
                    dna_hamming(seq, barcodes[br * block_cols + (bc + 1)])
                )
            if br + 1 < block_rows:
                adjacent_distances.append(
                    dna_hamming(seq, barcodes[(br + 1) * block_cols + bc])
                )

            for bc2 in range(bc + 5, min(bc + 15, block_cols)):
                far_distances.append(
                    dna_hamming(seq, barcodes[br * block_cols + bc2])
                )
            for br2 in range(br + 5, min(br + 15, block_rows)):
                far_distances.append(
                    dna_hamming(seq, barcodes[br2 * block_cols + bc])
                )

    stats: Dict[str, Any] = {}
    if adjacent_distances:
        stats["adjacent_min"] = min(adjacent_distances)
        stats["adjacent_mean"] = round(
            sum(adjacent_distances) / len(adjacent_distances), 2
        )
        stats["adjacent_max"] = max(adjacent_distances)
    if far_distances:
        stats["far_min"] = min(far_distances)
        stats["far_mean"] = round(sum(far_distances) / len(far_distances), 2)
        stats["far_max"] = max(far_distances)
    return stats


def soft_stats(verification: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for diff, stats in verification.get("soft_constraint_stats", {}).items():
        out[str(diff)] = {
            "target": stats["target"],
            "actual_mean": round(float(stats["actual_mean_distance"]), 3),
            "mean_abs_error": round(float(stats["mean_abs_error"]), 3),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--config",
        default=str(SCRIPT_DIR / "block_zipcode_config.json"),
        help="Path to Zipcode config JSON",
    )
    ap.add_argument(
        "--out-dir",
        default=str(SCRIPT_DIR),
        help="Output directory",
    )
    ap.add_argument(
        "--bar2",
        default=str(DEFAULT_BAR2),
        help="bar_2 library used after bar_1; trailing runs of row codes "
             "leave room for these leading runs",
    )
    args = ap.parse_args()

    config_path = Path(args.config)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    intermediate = out_dir / "intermediate"
    intermediate.mkdir(parents=True, exist_ok=True)

    print(f"Loading config from {config_path}", flush=True)
    config = load_config(config_path)

    block_rows: int = config.chip_geometry.rows
    block_cols: int = config.chip_geometry.cols
    n_blocks = block_rows * block_cols
    g_max = config.sequence_qc.g_content_max
    homo_max = config.sequence_qc.max_homopolymer_length
    spacer = config.codebook.spacer_sequence or ""
    enforce_junction = config.sequence_qc.reject_sequences_with_long_homopolymers
    candidate_mode = config.axis_code.candidate_mode
    flip_counts = tuple(config.axis_code.flip_counts)
    require_tolerance = config.axis_code.require_tolerance

    bar2_seqs = [line.strip().upper() for line in open(args.bar2) if line.strip()]
    downstream_start_runs = {base: 0 for base in "ACGT"}
    for seq in bar2_seqs:
        base, length = terminal_run(seq, at_end=False)
        downstream_start_runs[base] = max(downstream_start_runs[base], length)
    # Leave room for a unit start of every base, so a future bar_2 library
    # that uses all four starts still fits.
    for base in "ACGT":
        downstream_start_runs[base] = max(downstream_start_runs[base], 1)

    print(f"Block grid: {block_rows}×{block_cols} = {n_blocks:,} blocks")
    print(
        f"Generation: mode={candidate_mode}, flip_counts={list(flip_counts)}, "
        f"require_tolerance={require_tolerance}"
    )
    print(
        f"Composition: G≤{g_max:.0%} per axis code, "
        f"homopolymer≤{homo_max} nt"
        f"{' (junction enforced)' if enforce_junction else ''}"
    )
    print(f"bar_2 leading runs: {downstream_start_runs}")

    encoder = BitToDNAEncoder(
        scheme=config.bit_to_dna.scheme,
        mapping=config.bit_to_dna.mapping,
    )

    print("\n[Step 1] Generating column and row codes...")
    x_codes, y_codes, x_seqs, y_seqs, x_verif, y_verif = generate_column_and_row_codes(
        n_columns=block_cols,
        n_rows=block_rows,
        x_bit_length=config.x_bit_length,
        y_bit_length=config.y_bit_length,
        target_distances=config.error_correction.distances,
        far_features_min_distance=config.error_correction.far_features_min_distance,
        encoder=encoder,
        g_content_max=g_max,
        max_homopolymer_length=homo_max,
        spacer=spacer,
        enforce_junction=enforce_junction,
        max_attempts_per_position=config.axis_code.max_attempts_per_position,
        backtrack_depth=config.axis_code.backtrack_depth,
        tolerance=config.axis_code.tolerance,
        random_seed=config.axis_code.random_seed,
        candidate_mode=candidate_mode,
        flip_counts=flip_counts,
        require_tolerance=require_tolerance,
        downstream_start_runs=downstream_start_runs,
    )
    axis_nt = len(x_seqs[0])
    bar1_len = 2 * axis_nt + len(spacer)
    print(f"  X: {len(x_seqs)} × {axis_nt} nt")
    print(f"  Y: {len(y_seqs)} × {axis_nt} nt")
    print(
        f"  Hard violations: X={x_verif['hard_constraint_violations']}, "
        f"Y={y_verif['hard_constraint_violations']}"
    )
    print(
        f"  Outside tolerance: X={x_verif['positions_outside_tolerance']}, "
        f"Y={y_verif['positions_outside_tolerance']}"
    )

    with open(intermediate / "axis_codes.json", "w") as fh:
        json.dump(
            {
                "x_codes": [list(code) for code in x_codes],
                "y_codes": [list(code) for code in y_codes],
                "x_sequences": x_seqs,
                "y_sequences": y_seqs,
            },
            fh,
        )

    print(f"\n[Step 2] Combining into {bar1_len} nt bar_1 barcodes...")
    barcodes: List[str] = []
    for br in range(block_rows):
        for bc in range(block_cols):
            barcodes.append(x_seqs[bc] + spacer + y_seqs[br])

    assert len(barcodes) == n_blocks
    assert all(len(b) == bar1_len for b in barcodes)
    assert all(calculate_g_content(s) <= g_max + 1e-12 for s in x_seqs + y_seqs)

    # Full barcodes can exceed the per-axis G cap when both halves are near it.
    g_values = [calculate_g_content(b) for b in barcodes]
    n_g_fail = sum(1 for g in g_values if g > g_max)
    if n_g_fail:
        print(
            f"  NOTE: {n_g_fail}/{n_blocks} full bar_1 barcodes exceed "
            f"G {g_max:.0%} (axis halves are within the cap)."
        )

    junction = summarize_junction_homopolymers(
        x_seqs, y_seqs, homo_max, spacer=spacer
    )
    print(
        f"  X–Y junction homopolymers: max={junction['max_junction_run']}, "
        f"violations={junction['n_violations']}"
    )
    if junction["n_violations"]:
        raise RuntimeError(
            f"{junction['n_violations']} barcodes violate the junction "
            f"homopolymer limit of {homo_max}"
        )

    n_bar2_bad = 0
    for bar1 in barcodes:
        end_base, end_len = terminal_run(bar1, at_end=True)
        for bar2 in bar2_seqs:
            start_base, start_len = terminal_run(bar2, at_end=False)
            if end_base == start_base and end_len + start_len > homo_max:
                n_bar2_bad += 1
    if n_bar2_bad:
        raise RuntimeError(
            f"{n_bar2_bad} bar_1×bar_2 pairs would exceed the homopolymer "
            f"limit of {homo_max} at the bar_1/bar_2 junction"
        )
    print(f"  bar_1×bar_2 junctions: all {len(barcodes)*len(bar2_seqs):,} pairs OK")

    print("\n[Step 3] Validating uniqueness and X+Y structure...")
    n_unique = len(set(barcodes))
    if n_unique != n_blocks:
        raise RuntimeError(f"{n_blocks - n_unique} duplicate bar_1 barcodes")
    structure_ok = all(
        barcodes[br * block_cols + bc] == x_seqs[bc] + spacer + y_seqs[br]
        for br in range(block_rows)
        for bc in range(block_cols)
    )
    if not structure_ok:
        raise RuntimeError("X+Y structure broken")
    print(f"  OK — {n_blocks:,} unique barcodes, X+Y structure intact")

    print("\n[Step 4] Spatial DNA distances...")
    dist_stats = compute_spatial_distance_stats(barcodes, block_rows, block_cols)
    print(
        f"  Adjacent: min={dist_stats.get('adjacent_min')}  "
        f"mean={dist_stats.get('adjacent_mean')}  "
        f"max={dist_stats.get('adjacent_max')}"
    )
    print(
        f"  Far (≥5): min={dist_stats.get('far_min')}  "
        f"mean={dist_stats.get('far_mean')}  "
        f"max={dist_stats.get('far_max')}"
    )

    gc_values = [gc_fraction(b) for b in barcodes]
    hp_values = [find_max_homopolymer(b)[0] for b in barcodes]
    print(
        f"\n  G content:  mean={sum(g_values)/len(g_values):.3f}  "
        f"max={max(g_values):.3f}  >{g_max:.0%}: {n_g_fail}"
    )
    print(
        f"  GC content: mean={sum(gc_values)/len(gc_values):.3f}  "
        f"range [{min(gc_values):.3f}, {max(gc_values):.3f}]"
    )
    print(f"  Max homopolymer: {max(hp_values)}")

    txt_path = out_dir / "bar1_barcodes_28nt.txt"
    with open(txt_path, "w") as fh:
        for seq in barcodes:
            fh.write(seq + "\n")
    print(f"\nWrote barcodes       → {txt_path}")

    csv_path = out_dir / "bar1_codebook.csv"
    with open(csv_path, "w") as fh:
        fh.write("block_row,block_col,bar_1\n")
        for br in range(block_rows):
            for bc in range(block_cols):
                fh.write(f"{br},{bc},{barcodes[br * block_cols + bc]}\n")
    print(f"Wrote codebook CSV   → {csv_path}")

    qc_path = out_dir / "bar1_qc_report.json"
    qc_report = {
        "n_blocks": n_blocks,
        "block_grid": f"{block_rows}x{block_cols}",
        "bar1_length_nt": bar1_len,
        "axis_length_nt": axis_nt,
        "n_globally_unique": n_unique,
        "x_y_structure_intact": True,
        "candidate_mode": candidate_mode,
        "flip_counts": list(flip_counts),
        "require_tolerance": require_tolerance,
        "g_content_max_per_axis": g_max,
        "g_content": {
            "mean": round(sum(g_values) / len(g_values), 4),
            "max": round(max(g_values), 4),
            "n_above_cap": n_g_fail,
        },
        "gc_content": {
            "mean": round(sum(gc_values) / len(gc_values), 4),
            "min": round(min(gc_values), 4),
            "max": round(max(gc_values), 4),
        },
        "max_homopolymer": int(max(hp_values)),
        "junction_homopolymers": junction,
        "downstream_start_runs": downstream_start_runs,
        "bar2_file": str(Path(args.bar2)),
        "n_bar2": len(bar2_seqs),
        "spatial_hamming_distance": dist_stats,
        "axis_code_generation": {
            "x_hard_violations": x_verif["hard_constraint_violations"],
            "y_hard_violations": y_verif["hard_constraint_violations"],
            "x_outside_tolerance": x_verif["positions_outside_tolerance"],
            "y_outside_tolerance": y_verif["positions_outside_tolerance"],
            "x_soft_stats": soft_stats(x_verif),
            "y_soft_stats": soft_stats(y_verif),
        },
    }
    with open(qc_path, "w") as fh:
        json.dump(qc_report, fh, indent=2)
    print(f"Wrote QC report      → {qc_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
