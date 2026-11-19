#!/usr/bin/env python3
"""
Generate a Dual barcode layout.

Modes (barcode_configuration.mode):
  - shared_pool (default / legacy): one predefined half library for bar_1 and bar_2
  - fixed_bar2: distance-designed bar_2 repertoire (size ≥ B²) + bar_1 pool
      bar_1_source: "file" | "generate"

Optional sequence_qc enforces G / homopolymer (and optional GC) on combined barcodes.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

from dual_barcode_layout import (
    DualBarcodeLayout,
    generate_compatible_bar1_pool,
    load_barcodes_from_file,
)
from generate_bar2_repertoire import generate_bar2_repertoire
from sequence_qc import SequenceQC, min_pairwise_hamming


def load_config(config_path: str) -> dict:
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(config_path, "r") as f:
        return json.load(f)


def resolve_path(path: str, config_dir: Path) -> str:
    p = Path(path)
    if p.is_file():
        return str(p)
    alt = config_dir / path
    if alt.is_file():
        return str(alt)
    return path


def build_qc(cfg: dict | None) -> SequenceQC | None:
    if not cfg:
        return None
    return SequenceQC(
        g_content_max=cfg.get("g_content_max"),
        max_homopolymer=cfg.get("max_homopolymer_length"),
        gc_min=cfg.get("gc_content_min"),
        gc_max=cfg.get("gc_content_max"),
        max_leading=cfg.get("max_leading"),
    )


def save_statistics(stats: dict, config: dict, filepath: str) -> None:
    output = {
        "generation_timestamp": datetime.now().isoformat(),
        "configuration": config,
        "layout_statistics": {
            k: stats[k]
            for k in (
                "grid_dimensions",
                "total_spots",
                "block_size",
                "mode",
                "n_barcodes_available",
                "n_bar2_repertoire",
                "n_blocks",
                "unique_combined_barcodes",
                "theoretical_max_combinations",
                "g_content_combined",
                "max_homopolymer",
                "sequence_qc",
            )
            if k in stats
        },
        "barcode_usage": {
            "bar_1_unique_count": len(stats["bar1_usage"]),
            "bar_2_unique_count": len(stats["bar2_usage"]),
            "bar_1_usage_distribution": {
                "min": min(stats["bar1_usage"].values()) if stats["bar1_usage"] else 0,
                "max": max(stats["bar1_usage"].values()) if stats["bar1_usage"] else 0,
                "mean": float(np.mean(list(stats["bar1_usage"].values())))
                if stats["bar1_usage"]
                else 0,
            },
            "bar_2_usage_distribution": {
                "min": min(stats["bar2_usage"].values()) if stats["bar2_usage"] else 0,
                "max": max(stats["bar2_usage"].values()) if stats["bar2_usage"] else 0,
                "mean": float(np.mean(list(stats["bar2_usage"].values())))
                if stats["bar2_usage"]
                else 0,
            },
        },
    }
    with open(filepath, "w") as f:
        json.dump(output, f, indent=2)
    print(f"✓ Statistics saved to {filepath}")


def main() -> None:
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config_235x363.json"
    config_dir = Path(config_file).resolve().parent

    print("=" * 80)
    print("DUAL BARCODE LAYOUT GENERATOR")
    print("=" * 80)
    print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\n[1/6] Loading configuration from {config_file}...")
    config = load_config(config_file)
    print("✓ Configuration loaded")

    rows = config["array_dimensions"]["rows"]
    cols = config["array_dimensions"]["columns"]
    total_spots = config["array_dimensions"].get("total_spots", rows * cols)
    block_size = int(config.get("grid_configuration", {}).get("block_size", 4))
    n_blocks = ((rows + block_size - 1) // block_size) * (
        (cols + block_size - 1) // block_size
    )
    bc_cfg = config["barcode_configuration"]
    mode = bc_cfg.get("mode", "shared_pool")
    half_len = int(bc_cfg.get("barcode_length_per_segment", 18))
    seed = config.get("random_seed", 42)
    qc = build_qc(config.get("sequence_qc"))
    out_cfg = config["output"]
    output_layout = out_cfg["layout_file"]
    output_stats = out_cfg["statistics_file"]
    output_viz = out_cfg.get("visualization_file")

    print("\nArray configuration:")
    print(f"  Dimensions: {rows} × {cols} = {total_spots} spots")
    print(f"  Block size: {block_size}×{block_size} → {n_blocks} blocks")
    print(f"  Mode: {mode}")
    print(f"  Half length: {half_len} nt → combined {2 * half_len} nt")
    if qc is not None:
        print(
            f"  QC: G≤{qc.g_content_max}, HP≤{qc.max_homopolymer}, "
            f"leading≤{qc.max_leading}"
        )
    else:
        print("  QC: none (legacy behaviour)")

    bar2_set = None
    bar1_pool: list[str]

    print(f"\n[2/6] Preparing barcode libraries ({mode})...")
    if mode == "shared_pool":
        barcode_file = resolve_path(bc_cfg["barcode_file"], config_dir)
        # Also try repo-relative inhouse path
        if not os.path.exists(barcode_file):
            repo_alt = (
                Path(__file__).resolve().parent.parent
                / "inhouse_barcodes"
                / Path(bc_cfg["barcode_file"]).name
            )
            if repo_alt.is_file():
                barcode_file = str(repo_alt)
        bar1_pool = load_barcodes_from_file(barcode_file)
        print(f"✓ Shared pool: {len(bar1_pool)} halves from {barcode_file}")

    elif mode == "fixed_bar2":
        # --- bar_2 repertoire ---
        bar2_cfg = bc_cfg.get("bar2", {})
        n_bar2 = int(bar2_cfg.get("n", block_size * block_size))
        bar2_path = bar2_cfg.get("file")
        if bar2_path:
            bar2_path = resolve_path(bar2_path, config_dir)
            bar2_set = load_barcodes_from_file(bar2_path)
            print(
                f"✓ bar_2 repertoire: {len(bar2_set)} from {bar2_path} "
                f"(min Hamming={min_pairwise_hamming(bar2_set)})"
            )
        else:
            print(f"  Generating bar_2 repertoire (n={n_bar2}, len={half_len})...")
            bar2_set, achieved = generate_bar2_repertoire(
                n=n_bar2,
                length=half_len,
                min_hamming=int(bar2_cfg.get("min_hamming", 8)),
                g_max=float(
                    (config.get("sequence_qc") or {}).get("g_content_max", 0.10)
                ),
                max_homopolymer_limit=int(
                    (config.get("sequence_qc") or {}).get("max_homopolymer_length", 5)
                ),
                gc_min=(config.get("sequence_qc") or {}).get("gc_content_min"),
                gc_max=(config.get("sequence_qc") or {}).get("gc_content_max"),
                max_leading=int(
                    (config.get("sequence_qc") or {}).get("max_leading", 1)
                ),
                seed=seed,
                n_random_candidates=int(bar2_cfg.get("n_candidates", 200_000)),
            )
            print(f"✓ Generated {len(bar2_set)} bar_2 codes (min Hamming={achieved})")
            bar2_out = bar2_cfg.get("out", f"bar2_repertoire_{block_size}x{block_size}_{half_len}nt.txt")
            with open(bar2_out, "w") as fh:
                for s in bar2_set:
                    fh.write(s + "\n")
            print(f"  Wrote {bar2_out}")

        # --- bar_1 pool ---
        bar1_cfg = bc_cfg.get("bar1", {})
        bar1_source = bar1_cfg.get("source", "generate")
        if bar1_source == "file":
            bar1_file = resolve_path(bar1_cfg["file"], config_dir)
            if not os.path.exists(bar1_file):
                repo_alt = (
                    Path(__file__).resolve().parent.parent
                    / "inhouse_barcodes"
                    / Path(bar1_cfg["file"]).name
                )
                if repo_alt.is_file():
                    bar1_file = str(repo_alt)
            raw = load_barcodes_from_file(bar1_file)
            if qc is None:
                bar1_pool = raw
            else:
                # Filter to compatible bar_1s
                strict = bc_cfg.get("require_full_bar2_compatibility", True)
                need = len(bar2_set) if strict else block_size * block_size
                bar1_pool = []
                for b1 in raw:
                    if not qc.passes_half(b1):
                        continue
                    n_safe = sum(1 for b2 in bar2_set if qc.passes_combined(b1, b2))
                    if n_safe >= need:
                        bar1_pool.append(b1)
                print(
                    f"✓ bar_1 from file: {len(bar1_pool)}/{len(raw)} pass "
                    f"compatibility with bar_2 set"
                )
        else:
            if qc is None:
                raise ValueError("fixed_bar2 generate mode requires sequence_qc")
            n_bar1 = int(bar1_cfg.get("n", n_blocks + 100))
            bar1_min_h = int(bar1_cfg.get("min_hamming", 0))
            print(
                f"  Generating {n_bar1} compatible bar_1 codes"
                + (f" (min Hamming {bar1_min_h})" if bar1_min_h else "")
                + "..."
            )
            bar1_pool = generate_compatible_bar1_pool(
                n=n_bar1,
                length=half_len,
                bar2_set=bar2_set,
                qc=qc,
                seed=seed,
                require_full_compatibility=bc_cfg.get(
                    "require_full_bar2_compatibility", True
                ),
                min_safe_bar2=block_size * block_size,
                min_hamming=bar1_min_h,
                max_attempts_factor=int(bar1_cfg.get("max_attempts_factor", 2000)),
                candidate_pool_factor=int(bar1_cfg.get("candidate_pool_factor", 40)),
            )
            print(f"✓ Generated {len(bar1_pool)} compatible bar_1 codes")
            if bar1_min_h:
                print(
                    f"  bar_1 min pairwise Hamming="
                    f"{min_pairwise_hamming(bar1_pool)}"
                )
            bar1_out = bar1_cfg.get("out", f"bar1_pool_{half_len}nt_compatible.txt")
            with open(bar1_out, "w") as fh:
                for s in bar1_pool:
                    fh.write(s + "\n")
            print(f"  Wrote {bar1_out}")
    else:
        raise ValueError(f"Unknown mode: {mode}")

    print("\n[3/6] Validating libraries...")
    if len(bar1_pool) < n_blocks:
        raise ValueError(
            f"Need ≥ {n_blocks} bar_1 codes for this array, got {len(bar1_pool)}"
        )
    print(f"✓ bar_1 pool size OK ({len(bar1_pool)} ≥ {n_blocks} blocks)")

    print("\n[4/6] Generating dual barcode layout...")
    layout = DualBarcodeLayout(
        rows=rows,
        cols=cols,
        barcodes=bar1_pool,
        block_size=block_size,
        bar2_barcodes=bar2_set,
        sequence_qc=qc,
        require_full_bar2_compatibility=bc_cfg.get(
            "require_full_bar2_compatibility", True
        ),
    )
    layout.generate_layout(seed=seed)
    n_assigned = int(np.sum(layout.grid_combined != None))
    print(f"✓ Layout generated ({n_assigned}/{total_spots} spots)")

    print("\n[5/6] Computing statistics...")
    stats = layout.get_statistics()
    print(f"  Unique combined: {stats['unique_combined_barcodes']:,}")
    print(f"  Unique bar_1: {len(stats['bar1_usage']):,}")
    print(f"  Unique bar_2: {len(stats['bar2_usage']):,}")
    g = stats["g_content_combined"]
    print(f"  Combined G: mean={g['mean']:.4f} max={g['max']:.4f}")
    print(f"  Max homopolymer: {stats['max_homopolymer']}")
    save_statistics(stats, config, output_stats)

    print("\n[6/6] Exporting results...")
    layout.export_layout(output_layout, format="csv")
    print(f"✓ Layout exported to {output_layout}")
    if output_viz:
        try:
            layout.visualize_layout(output_viz, include_combined=False)
        except Exception as e:
            print(f"⚠ Visualization skipped: {e}")

    print("\n" + "=" * 80)
    print("GENERATION COMPLETE!")
    print("=" * 80)
    print(f"  Mode: {stats['mode']}")
    print(f"  Blocks: {block_size}×{block_size} ({stats['n_blocks']})")
    print(f"  Unique barcodes: {stats['unique_combined_barcodes']:,}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n✗ Interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
