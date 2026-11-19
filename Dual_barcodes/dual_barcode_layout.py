"""
Dual Barcode Layout Generator for Spatial Transcriptomics Arrays

Two library modes:

1. shared_pool (legacy / printed): one predefined half-barcode library is used
   for both bar_1 and bar_2. bar_2 may be reused across blocks; uniqueness is
   the (bar_1, bar_2) pair. Block size defaults to 4×4.

2. fixed_bar2: a small distance-designed bar_2 repertoire of size ≥ B² is reused
   in every block; bar_1 comes from a separate pool (file or generated). Optional
   combined QC (G, homopolymer) filters which pairs are allowed.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from sequence_qc import SequenceQC, g_fraction, max_homopolymer, summarize_g


class DualBarcodeLayout:
    """
    Generate dual barcode layouts.

    Args:
        rows, cols: array dimensions
        barcodes: bar_1 pool (and, in shared_pool mode, also the bar_2 pool)
        block_size: side length B of each block (default 4)
        bar2_barcodes: if set, fixed bar_2 repertoire (fixed_bar2 mode)
        sequence_qc: optional combined/half QC; None disables filtering
        require_full_bar2_compatibility: in fixed_bar2 mode, each bar_1 must
            form a QC-passing combination with every bar_2 (strict). If False,
            a bar_1 is usable when at least as many safe bar_2s exist as spots
            in its block.
    """

    def __init__(
        self,
        rows: int,
        cols: int,
        barcodes: List[str],
        block_size: int = 4,
        bar2_barcodes: Optional[List[str]] = None,
        sequence_qc: Optional[SequenceQC] = None,
        require_full_bar2_compatibility: bool = True,
    ):
        if block_size < 1:
            raise ValueError(f"block_size must be ≥ 1, got {block_size}")
        spots_per_full_block = block_size * block_size
        if bar2_barcodes is None:
            if len(barcodes) < spots_per_full_block:
                raise ValueError(
                    f"shared_pool mode needs ≥ {spots_per_full_block} barcodes "
                    f"for {block_size}×{block_size} blocks, got {len(barcodes)}"
                )
        else:
            if len(bar2_barcodes) < spots_per_full_block:
                raise ValueError(
                    f"fixed_bar2 mode needs ≥ {spots_per_full_block} bar_2 codes, "
                    f"got {len(bar2_barcodes)}"
                )
            if len(barcodes) < 1:
                raise ValueError("bar_1 pool is empty")

        self.rows = rows
        self.cols = cols
        self.block_size = block_size
        self.spots_per_full_block = spots_per_full_block
        self.barcodes = list(barcodes)
        self.n_barcodes = len(self.barcodes)
        self.bar2_barcodes = list(bar2_barcodes) if bar2_barcodes is not None else None
        self.mode = "fixed_bar2" if bar2_barcodes is not None else "shared_pool"
        self.sequence_qc = sequence_qc
        self.require_full_bar2_compatibility = require_full_bar2_compatibility

        self.grid_bar1 = np.empty((rows, cols), dtype=object)
        self.grid_bar2 = np.empty((rows, cols), dtype=object)
        self.grid_combined = np.empty((rows, cols), dtype=object)

        self.n_blocks_row = (rows + block_size - 1) // block_size
        self.n_blocks_col = (cols + block_size - 1) // block_size
        self.total_blocks = self.n_blocks_row * self.n_blocks_col

        self._qc_rejects = 0
        self._compatible_bar1: List[str] = []

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate_layout(
        self, seed: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        if self.mode == "fixed_bar2":
            self._prepare_compatible_bar1()
            self._assign_bar1_fixed()
            self._assign_bar2_fixed()
        else:
            self._assign_bar1_shared()
            self._assign_bar2_shared()

        self._create_combined_barcodes()
        self._validate_layout()
        return self.grid_bar1, self.grid_bar2, self.grid_combined

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    def _block_positions(self, block_row: int, block_col: int) -> List[Tuple[int, int]]:
        B = self.block_size
        positions = []
        for i in range(B):
            for j in range(B):
                row = block_row * B + i
                col = block_col * B + j
                if row < self.rows and col < self.cols:
                    positions.append((row, col))
        return positions

    def _safe_bar2_for(self, bar1: str, pool: Sequence[str]) -> List[str]:
        if self.sequence_qc is None:
            return list(pool)
        return [b2 for b2 in pool if self.sequence_qc.passes_combined(bar1, b2)]

    # ------------------------------------------------------------------
    # Legacy shared-pool assignment
    # ------------------------------------------------------------------

    def _assign_bar1_shared(self) -> None:
        max_uses_per_barcode = self.n_barcodes
        if self.total_blocks > self.n_barcodes * max_uses_per_barcode:
            raise ValueError(
                f"Grid too large: {self.total_blocks} blocks but only "
                f"{self.n_barcodes} barcodes in the shared pool"
            )

        # With QC, keep bar_1 candidates that can possibly form enough safe pairs.
        # Combined G≤g_max on length 2L implies ≤ floor(g_max*2L) G's total, so a
        # bar_1 with more G's than that can never pass. Full pair checks happen in
        # _assign_bar2_shared (avoid O(n²) here).
        if self.sequence_qc is None:
            candidates = list(self.barcodes)
        else:
            half_len = len(self.barcodes[0])
            max_g_total = None
            if self.sequence_qc.g_content_max is not None:
                max_g_total = int(self.sequence_qc.g_content_max * (2 * half_len) + 1e-9)
            candidates = []
            for b1 in self.barcodes:
                if max_g_total is not None and b1.count("G") > max_g_total:
                    continue
                if (
                    self.sequence_qc.max_homopolymer is not None
                    and max_homopolymer(b1) > self.sequence_qc.max_homopolymer
                ):
                    continue
                candidates.append(b1)
            if len(candidates) < self.total_blocks:
                raise ValueError(
                    f"shared_pool + QC: only {len(candidates)} plausible bar_1 "
                    f"codes; need {self.total_blocks} blocks. Relax QC or enlarge "
                    f"the pool."
                )

        shuffled = candidates.copy()
        random.shuffle(shuffled)
        usage = {bc: 0 for bc in candidates}
        # Prefer distinct bar_1 per block when the pool is large enough.
        prefer_unique = len(candidates) >= self.total_blocks

        for br in range(self.n_blocks_row):
            for bc in range(self.n_blocks_col):
                if prefer_unique:
                    available = [b for b in shuffled if usage[b] == 0]
                else:
                    available = [
                        b for b in shuffled if usage[b] < max_uses_per_barcode
                    ]
                if not available:
                    available = [
                        b for b in shuffled if usage[b] < max_uses_per_barcode
                    ]
                if not available:
                    raise ValueError("No available bar_1 barcodes")
                bar1 = min(available, key=lambda b: usage[b])
                usage[bar1] += 1
                for row, col in self._block_positions(br, bc):
                    self.grid_bar1[row, col] = bar1

    def _assign_bar2_shared(self) -> None:
        used: Dict[str, Set[str]] = {bc: set() for bc in self.barcodes}
        # Bucket by G count to prune combined-G checks quickly.
        by_g: Dict[int, List[str]] = {}
        for b2 in self.barcodes:
            by_g.setdefault(b2.count("G"), []).append(b2)
        max_g_total = None
        if self.sequence_qc is not None and self.sequence_qc.g_content_max is not None:
            half_len = len(self.barcodes[0])
            max_g_total = int(self.sequence_qc.g_content_max * (2 * half_len) + 1e-9)

        for br in range(self.n_blocks_row):
            for bc in range(self.n_blocks_col):
                positions = self._block_positions(br, bc)
                if not positions:
                    continue
                bar1 = self.grid_bar1[positions[0]]
                n_pos = len(positions)
                used_here = used[bar1]
                if self.sequence_qc is None:
                    pool = self.barcodes
                else:
                    g1 = bar1.count("G")
                    pool = []
                    g2_max = (
                        max(by_g)
                        if max_g_total is None
                        else max_g_total - g1
                    )
                    for g2, bucket in by_g.items():
                        if g2 > g2_max:
                            continue
                        pool.extend(bucket)
                available = [
                    b2
                    for b2 in pool
                    if b2 not in used_here
                    and (
                        self.sequence_qc is None
                        or self.sequence_qc.passes_combined(bar1, b2)
                    )
                ]
                if len(available) < n_pos:
                    raise ValueError(
                        f"Cannot fill block ({br},{bc}) for bar_1={bar1}: "
                        f"need {n_pos} unused QC-safe bar_2, have {len(available)}"
                    )
                random.shuffle(available)
                selected = available[:n_pos]
                for b2 in selected:
                    used_here.add(b2)
                for (row, col), b2 in zip(positions, selected):
                    self.grid_bar2[row, col] = b2

    # ------------------------------------------------------------------
    # Fixed bar_2 repertoire assignment
    # ------------------------------------------------------------------

    def _prepare_compatible_bar1(self) -> None:
        assert self.bar2_barcodes is not None
        compatible = []
        for bar1 in self.barcodes:
            safe = self._safe_bar2_for(bar1, self.bar2_barcodes)
            if self.require_full_bar2_compatibility:
                if len(safe) == len(self.bar2_barcodes):
                    compatible.append(bar1)
            else:
                if len(safe) >= self.spots_per_full_block:
                    compatible.append(bar1)

        if len(compatible) < self.total_blocks:
            raise ValueError(
                f"Need {self.total_blocks} QC-compatible bar_1 codes for "
                f"{self.n_blocks_row}×{self.n_blocks_col} blocks, found "
                f"{len(compatible)} in a pool of {len(self.barcodes)}. "
                f"Generate/filter a larger bar_1 pool or relax QC."
            )
        self._compatible_bar1 = compatible

    def _assign_bar1_fixed(self) -> None:
        # Prefer distinct bar_1 per block (one use each).
        pool = self._compatible_bar1.copy()
        random.shuffle(pool)
        if len(pool) < self.total_blocks:
            raise ValueError("Internal error: insufficient compatible bar_1")

        idx = 0
        for br in range(self.n_blocks_row):
            for bc in range(self.n_blocks_col):
                bar1 = pool[idx]
                idx += 1
                for row, col in self._block_positions(br, bc):
                    self.grid_bar1[row, col] = bar1

    def _assign_bar2_fixed(self) -> None:
        assert self.bar2_barcodes is not None
        for br in range(self.n_blocks_row):
            for bc in range(self.n_blocks_col):
                positions = self._block_positions(br, bc)
                if not positions:
                    continue
                bar1 = self.grid_bar1[positions[0]]
                safe = self._safe_bar2_for(bar1, self.bar2_barcodes)
                if len(safe) < len(positions):
                    raise ValueError(
                        f"Block ({br},{bc}) bar_1={bar1} has only "
                        f"{len(safe)} safe bar_2 for {len(positions)} spots"
                    )
                random.shuffle(safe)
                for (row, col), b2 in zip(positions, safe[: len(positions)]):
                    self.grid_bar2[row, col] = b2

    # ------------------------------------------------------------------
    # Combine / validate / stats
    # ------------------------------------------------------------------

    def _create_combined_barcodes(self) -> None:
        for i in range(self.rows):
            for j in range(self.cols):
                if self.grid_bar1[i, j] is not None and self.grid_bar2[i, j] is not None:
                    self.grid_combined[i, j] = self.grid_bar1[i, j] + self.grid_bar2[i, j]

    def _validate_layout(self) -> bool:
        all_combined: Dict[str, Tuple[int, int]] = {}
        B = self.block_size

        for br in range(self.n_blocks_row):
            for bc in range(self.n_blocks_col):
                block_bar1 = set()
                block_bar2 = set()
                block_combined = set()
                n_positions = 0

                for i in range(B):
                    for j in range(B):
                        row = br * B + i
                        col = bc * B + j
                        if row >= self.rows or col >= self.cols:
                            continue
                        n_positions += 1
                        block_bar1.add(self.grid_bar1[row, col])
                        block_bar2.add(self.grid_bar2[row, col])
                        combined = self.grid_combined[row, col]
                        block_combined.add(combined)
                        if combined in all_combined:
                            raise ValueError(
                                f"Duplicate combined barcode '{combined}' at "
                                f"({row},{col}) and {all_combined[combined]}"
                            )
                        all_combined[combined] = (row, col)

                        if self.sequence_qc is not None:
                            if not self.sequence_qc.passes_combined(
                                self.grid_bar1[row, col], self.grid_bar2[row, col]
                            ):
                                raise ValueError(
                                    f"Combined QC failed at ({row},{col}): {combined}"
                                )

                if n_positions == 0:
                    continue
                if len(block_bar1) != 1:
                    raise ValueError(f"Block ({br},{bc}) has inconsistent bar_1")
                if len(block_bar2) != n_positions:
                    raise ValueError(f"Block ({br},{bc}) has duplicate bar_2")
                if len(block_combined) != n_positions:
                    raise ValueError(f"Block ({br},{bc}) has duplicate combined barcodes")
        return True

    def get_statistics(self) -> dict:
        combined = [
            self.grid_combined[i, j]
            for i in range(self.rows)
            for j in range(self.cols)
            if self.grid_combined[i, j] is not None
        ]
        bar1_counts: Dict[str, int] = {}
        bar2_counts: Dict[str, int] = {}
        for i in range(self.rows):
            for j in range(self.cols):
                b1 = self.grid_bar1[i, j]
                b2 = self.grid_bar2[i, j]
                if b1 is not None:
                    bar1_counts[b1] = bar1_counts.get(b1, 0) + 1
                if b2 is not None:
                    bar2_counts[b2] = bar2_counts.get(b2, 0) + 1

        hp_vals = [max_homopolymer(c) for c in combined]
        return {
            "grid_dimensions": (self.rows, self.cols),
            "total_spots": self.rows * self.cols,
            "block_size": self.block_size,
            "mode": self.mode,
            "n_barcodes_available": self.n_barcodes,
            "n_bar2_repertoire": (
                len(self.bar2_barcodes) if self.bar2_barcodes is not None else None
            ),
            "n_blocks": self.total_blocks,
            "unique_combined_barcodes": len(set(combined)),
            "bar1_usage": bar1_counts,
            "bar2_usage": bar2_counts,
            "theoretical_max_combinations": (
                self.n_barcodes * len(self.bar2_barcodes)
                if self.bar2_barcodes is not None
                else self.n_barcodes ** 2
            ),
            "g_content_combined": summarize_g(combined),
            "max_homopolymer": max(hp_vals) if hp_vals else 0,
            "sequence_qc": None
            if self.sequence_qc is None
            else {
                "g_content_max": self.sequence_qc.g_content_max,
                "max_homopolymer": self.sequence_qc.max_homopolymer,
                "gc_min": self.sequence_qc.gc_min,
                "gc_max": self.sequence_qc.gc_max,
            },
        }

    def export_layout(self, filepath: str, format: str = "csv"):
        if format == "csv":
            self._export_csv(filepath, delimiter=",")
        elif format == "tsv":
            self._export_csv(filepath, delimiter="\t")
        elif format == "npy":
            np.save(
                filepath,
                {
                    "bar1": self.grid_bar1,
                    "bar2": self.grid_bar2,
                    "combined": self.grid_combined,
                },
            )
        else:
            raise ValueError(f"Unknown format: {format}")

    def _export_csv(self, filepath: str, delimiter: str = ","):
        with open(filepath, "w") as f:
            f.write(
                f"row{delimiter}col{delimiter}bar_1{delimiter}bar_2{delimiter}combined\n"
            )
            for j in range(self.cols):
                for i in range(self.rows):
                    f.write(f"{i}{delimiter}{j}{delimiter}")
                    f.write(f"{self.grid_bar1[i, j]}{delimiter}")
                    f.write(f"{self.grid_bar2[i, j]}{delimiter}")
                    f.write(f"{self.grid_combined[i, j]}\n")

    def visualize_layout(
        self, output_file: Optional[str] = None, include_combined: bool = True
    ):
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
        except ImportError:
            print("Matplotlib not available. Install with: pip install matplotlib")
            return

        unique_bar1 = list({b for b in self.grid_bar1.flatten() if b is not None})
        unique_bar2 = list({b for b in self.grid_bar2.flatten() if b is not None})
        bar1_colors = {bc: i for i, bc in enumerate(unique_bar1)}
        bar2_colors = {bc: i for i, bc in enumerate(unique_bar2)}
        B = self.block_size

        if include_combined:
            fig, axes = plt.subplots(1, 3, figsize=(18, 6))
            ax_bar1, ax_bar2, ax_combined = axes
        else:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6))
            ax_bar1, ax_bar2 = axes

        bar1_numeric = np.array(
            [
                [bar1_colors[self.grid_bar1[i, j]] for j in range(self.cols)]
                for i in range(self.rows)
            ]
        )
        ax_bar1.imshow(bar1_numeric, cmap="tab20", aspect="auto")
        ax_bar1.set_title(
            f"bar_1 Layout\n(Same within {B}×{B} blocks)",
            fontsize=12,
            fontweight="bold",
        )
        for i in range(0, self.rows, B):
            ax_bar1.axhline(i - 0.5, color="black", linewidth=2)
        for j in range(0, self.cols, B):
            ax_bar1.axvline(j - 0.5, color="black", linewidth=2)

        bar2_numeric = np.array(
            [
                [bar2_colors[self.grid_bar2[i, j]] for j in range(self.cols)]
                for i in range(self.rows)
            ]
        )
        ax_bar2.imshow(bar2_numeric, cmap="tab20", aspect="auto")
        ax_bar2.set_title(
            f"bar_2 Layout\n(Different within {B}×{B} blocks)",
            fontsize=12,
            fontweight="bold",
        )
        for i in range(0, self.rows, B):
            ax_bar2.axhline(i - 0.5, color="black", linewidth=2)
        for j in range(0, self.cols, B):
            ax_bar2.axvline(j - 0.5, color="black", linewidth=2)

        if include_combined:
            ax_combined.set_xlim(-0.5, self.cols - 0.5)
            ax_combined.set_ylim(self.rows - 0.5, -0.5)
            ax_combined.set_aspect("equal")
            for i in range(self.rows):
                for j in range(self.cols):
                    color_idx = bar1_colors[self.grid_bar1[i, j]]
                    color = plt.cm.tab20(color_idx / max(len(unique_bar1), 1))
                    rect = mpatches.Rectangle(
                        (j - 0.4, i - 0.4),
                        0.8,
                        0.8,
                        facecolor=color,
                        edgecolor="black",
                        linewidth=0.5,
                    )
                    ax_combined.add_patch(rect)
                    fontsize = max(6, min(10, 100 // max(self.rows, self.cols)))
                    ax_combined.text(
                        j,
                        i,
                        self.grid_combined[i, j],
                        ha="center",
                        va="center",
                        fontsize=fontsize,
                        fontweight="bold",
                    )
            for i in range(0, self.rows, B):
                ax_combined.axhline(i - 0.5, color="red", linewidth=2)
            for j in range(0, self.cols, B):
                ax_combined.axvline(j - 0.5, color="red", linewidth=2)
            ax_combined.set_title("Combined Layout (bar_1 + bar_2)")

        plt.tight_layout()
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            print(f"Visualization saved to {output_file}")
        else:
            plt.show()


def generate_random_barcodes(
    n: int,
    length: int = 8,
    seed: Optional[int] = None,
    qc: Optional[SequenceQC] = None,
    max_attempts_factor: int = 100,
) -> List[str]:
    """Generate random DNA barcodes, optionally filtered by half QC."""
    if seed is not None:
        random.seed(seed)
    bases = ["A", "C", "G", "T"]
    barcodes: List[str] = []
    seen = set()
    attempts = 0
    max_attempts = max(n * max_attempts_factor, n)
    while len(barcodes) < n and attempts < max_attempts:
        attempts += 1
        barcode = "".join(random.choices(bases, k=length))
        if barcode in seen:
            continue
        if qc is not None and not qc.passes_half(barcode):
            continue
        seen.add(barcode)
        barcodes.append(barcode)
    if len(barcodes) < n:
        raise RuntimeError(
            f"Only generated {len(barcodes)}/{n} barcodes of length {length} "
            f"under QC after {attempts} attempts"
        )
    return barcodes


def generate_compatible_bar1_pool(
    n: int,
    length: int,
    bar2_set: Sequence[str],
    qc: SequenceQC,
    seed: int = 42,
    require_full_compatibility: bool = True,
    min_safe_bar2: Optional[int] = None,
    max_attempts_factor: int = 500,
    min_hamming: int = 0,
    candidate_pool_factor: int = 40,
) -> List[str]:
    """
    Generate bar_1 halves that are QC-compatible with a fixed bar_2 set.

    If require_full_compatibility, every bar_2 must form a valid combined
    barcode with bar_1. Otherwise at least min_safe_bar2 combinations must pass
    (default: len(bar2_set)).

    If min_hamming > 0, first build a large QC-compatible pool, then greedily
    pack n codes at that minimum Hamming (falls back by 1 if needed down to
    max(min_hamming-2, 4)).
    """
    from sequence_qc import hamming

    rng = random.Random(seed)
    bases = ["A", "C", "G", "T"]
    need = (
        len(bar2_set)
        if require_full_compatibility
        else (min_safe_bar2 if min_safe_bar2 is not None else len(bar2_set))
    )

    def is_compatible(bar1: str) -> bool:
        if not qc.passes_half(bar1):
            return False
        n_safe = sum(1 for b2 in bar2_set if qc.passes_combined(bar1, b2))
        return n_safe >= need

    if min_hamming <= 0:
        out: List[str] = []
        seen = set()
        attempts = 0
        max_attempts = max(n * max_attempts_factor, n)
        while len(out) < n and attempts < max_attempts:
            attempts += 1
            bar1 = "".join(rng.choices(bases, k=length))
            if bar1 in seen:
                continue
            seen.add(bar1)
            if is_compatible(bar1):
                out.append(bar1)
        if len(out) < n:
            raise RuntimeError(
                f"Only generated {len(out)}/{n} bar_1 codes compatible with the "
                f"bar_2 set under QC after {attempts} attempts"
            )
        return out

    # Phase 1: large compatible pool (no mutual distance yet)
    pool_target = max(n * candidate_pool_factor, n + 500)
    pool: List[str] = []
    seen = set()
    attempts = 0
    max_attempts = max(pool_target * max_attempts_factor // 10, pool_target * 50)
    while len(pool) < pool_target and attempts < max_attempts:
        attempts += 1
        bar1 = "".join(rng.choices(bases, k=length))
        if bar1 in seen:
            continue
        seen.add(bar1)
        if is_compatible(bar1):
            pool.append(bar1)
    if len(pool) < n:
        raise RuntimeError(
            f"Only generated {len(pool)}/{n} QC-compatible bar_1 candidates "
            f"after {attempts} attempts (need a larger search)"
        )
    rng.shuffle(pool)

    # Phase 2: greedy pack at min_hamming, falling back if needed
    floor = max(min_hamming - 2, 4)
    for mh in range(min_hamming, floor - 1, -1):
        selected: List[str] = []
        for cand in pool:
            if all(hamming(cand, s) >= mh for s in selected):
                selected.append(cand)
                if len(selected) == n:
                    if mh < min_hamming:
                        print(
                            f"  bar_1 pack: achieved min Hamming {mh} "
                            f"(requested {min_hamming})"
                        )
                    return selected
        print(
            f"  bar_1 pack: only {len(selected)}/{n} at min Hamming {mh} "
            f"from pool of {len(pool)}; retrying lower…"
        )
    raise RuntimeError(
        f"Could not pack {n} bar_1 codes at min Hamming ≥ {floor} "
        f"from a pool of {len(pool)}"
    )


def load_barcodes_from_file(filepath: str) -> List[str]:
    """Load barcodes from a text file (one barcode per line)."""
    try:
        with open(filepath, "r") as f:
            barcodes = [line.strip() for line in f if line.strip()]
        if not barcodes:
            raise ValueError(f"No barcodes found in file: {filepath}")
        for i, bc in enumerate(barcodes):
            if not bc:
                raise ValueError(f"Empty barcode at line {i + 1}")
        print(f"Loaded {len(barcodes)} barcodes from {filepath}")
        return barcodes
    except FileNotFoundError:
        raise FileNotFoundError(f"Barcode file not found: {filepath}")
    except Exception as e:
        raise ValueError(f"Error reading barcode file: {e}")


if __name__ == "__main__":
    import os
    import sys

    print("=" * 70)
    print("Dual Barcode Layout Generator")
    print("=" * 70)

    ROWS = 8
    COLS = 8
    SEED = 42
    BARCODE_FILE = None

    if len(sys.argv) > 1:
        BARCODE_FILE = sys.argv[1]
        if len(sys.argv) > 2:
            ROWS = int(sys.argv[2])
        if len(sys.argv) > 3:
            COLS = int(sys.argv[3])

    N_BLOCKS = ((ROWS + 3) // 4) * ((COLS + 3) // 4)
    print(f"\nGrid size: {ROWS}x{COLS} ({ROWS * COLS} spots)")
    print(f"Number of 4x4 blocks: {N_BLOCKS}")

    if BARCODE_FILE and os.path.exists(BARCODE_FILE):
        barcodes = load_barcodes_from_file(BARCODE_FILE)
    else:
        if BARCODE_FILE:
            print(f"⚠ File '{BARCODE_FILE}' not found; generating random barcodes.")
        barcodes = generate_random_barcodes(max(N_BLOCKS, 16), 6, seed=SEED)

    if len(barcodes) < N_BLOCKS:
        print(f"✗ Need ≥ {N_BLOCKS} barcodes, got {len(barcodes)}")
        sys.exit(1)

    layout = DualBarcodeLayout(ROWS, COLS, barcodes)
    try:
        layout.generate_layout(seed=SEED)
        stats = layout.get_statistics()
        print("✓ Layout generated")
        print(f"  Unique combined: {stats['unique_combined_barcodes']}")
        layout.export_layout("barcode_layout.csv")
        print("✓ Exported barcode_layout.csv")
    except ValueError as e:
        print(f"✗ Error: {e}")
        sys.exit(1)
