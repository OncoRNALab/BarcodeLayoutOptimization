"""Sequence QC helpers for Dual barcode halves and combined barcodes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def g_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return seq.count("G") / len(seq)


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    max_run = cur = 1
    for i in range(1, len(seq)):
        cur = cur + 1 if seq[i] == seq[i - 1] else 1
        max_run = max(max_run, cur)
    return max_run


def leading_run(seq: str) -> int:
    if not seq:
        return 0
    i = 1
    while i < len(seq) and seq[i] == seq[0]:
        i += 1
    return i


def trailing_run(seq: str) -> int:
    if not seq:
        return 0
    i = 1
    while i < len(seq) and seq[-1 - i] == seq[-1]:
        i += 1
    return i


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def min_pairwise_hamming(seqs: list[str]) -> int:
    if len(seqs) < 2:
        return 0
    return min(hamming(a, b) for i, a in enumerate(seqs) for b in seqs[i + 1 :])


@dataclass
class SequenceQC:
    """QC thresholds. None means 'do not enforce' for that field."""

    g_content_max: Optional[float] = None
    max_homopolymer: Optional[int] = None
    gc_min: Optional[float] = None
    gc_max: Optional[float] = None
    max_leading: Optional[int] = None  # for bar_2 halves only

    def passes_half(self, seq: str) -> bool:
        if self.g_content_max is not None and g_fraction(seq) > self.g_content_max:
            return False
        if self.max_homopolymer is not None and max_homopolymer(seq) > self.max_homopolymer:
            return False
        if self.gc_min is not None and gc_content(seq) < self.gc_min:
            return False
        if self.gc_max is not None and gc_content(seq) > self.gc_max:
            return False
        if self.max_leading is not None and leading_run(seq) > self.max_leading:
            return False
        return True

    def passes_combined(self, bar1: str, bar2: str) -> bool:
        combined = bar1 + bar2
        if self.g_content_max is not None and g_fraction(combined) > self.g_content_max:
            return False
        if self.max_homopolymer is not None and max_homopolymer(combined) > self.max_homopolymer:
            return False
        if self.gc_min is not None and gc_content(combined) < self.gc_min:
            return False
        if self.gc_max is not None and gc_content(combined) > self.gc_max:
            return False
        return True


def summarize_g(seqs: list[str]) -> dict:
    if not seqs:
        return {"mean": 0.0, "max": 0.0, "n": 0}
    vals = [g_fraction(s) for s in seqs]
    return {
        "mean": round(sum(vals) / len(vals), 4),
        "max": round(max(vals), 4),
        "n": len(seqs),
    }
