#!/usr/bin/env python3
"""Per-barcode recall in 1M simulations, and a subsample rescue of unassigned real reads.

Writes simulation_recall.csv, decoding_rescue_subsample.csv and
half_rescue.csv.
"""

from __future__ import annotations

import gzip
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein

from io_utils import (
    ARRAYS,
    BARCODE_LEN,
    OUT_DIR,
    SEG,
    SPATIAL_R,
    UMI_LEN,
    load_layout,
    raw_fastq,
    run_prefix,
    sim_prefix,
)
from read_accounting import insert_class, insert_end

N_FEATURES = 235 * 363
SIM_READS = 1_000_000
SUBSAMPLE = 15_000
HAMMING_PREFILTER = 12
RNG = np.random.default_rng(2026)


def pack_barcodes(seqs: list[str]) -> np.ndarray:
    enc = np.empty((len(seqs), BARCODE_LEN), dtype=np.uint8)
    tbl = {b: i for i, b in enumerate("ACGT")}
    for i, s in enumerate(seqs):
        enc[i] = [tbl.get(c, 4) for c in s]
    return enc


def hamming_to_library(query: str, lib: np.ndarray) -> np.ndarray:
    tbl = {b: i for i, b in enumerate("ACGT")}
    q = np.array([tbl.get(c, 4) for c in query[:BARCODE_LEN]], dtype=np.uint8)
    if len(q) < BARCODE_LEN:
        q = np.pad(q, (0, BARCODE_LEN - len(q)), constant_values=4)
    return (lib != q).sum(axis=1)


def simulation_recall(name: str, layout: pd.DataFrame) -> pd.DataFrame:
    real = np.array(layout.barcode)
    n_correct = np.zeros(len(layout), dtype=np.int32)
    n_from = np.zeros(len(layout), dtype=np.int32)
    n_assigned_as = np.zeros(len(layout), dtype=np.int32)
    lookup = {b: i for i, b in enumerate(real)}
    path = f"{sim_prefix(name, SPATIAL_R)}_R1_filtered.fastq"
    with open(path) as handle:
        for i, line in enumerate(handle):
            if i % 4:
                continue
            parts = line.strip().split("_")
            # @read_n_barcode_<true>_calledidx_<idx>_<SEQ>
            true_idx = int(parts[3])
            called = parts[-1].split()[0]
            n_from[true_idx] += 1
            if called == real[true_idx]:
                n_correct[true_idx] += 1
            idx = lookup.get(called)
            if idx is not None:
                n_assigned_as[idx] += 1
    depth = SIM_READS / len(layout)
    return pd.DataFrame({
        "barcode": real,
        "sim_n_from_assigned": n_from,
        "sim_n_correct": n_correct,
        "sim_n_assigned_as": n_assigned_as,
        "sim_recall": n_correct / depth,
        "sim_assignment_from": n_from / depth,
        "sim_precision_from": np.divide(n_correct, n_from, out=np.full(len(layout), np.nan), where=n_from > 0),
    })


def rescue_subsample(name: str, layout: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    real = set(layout.barcode)
    assigned_ids = set()
    with open(f"{run_prefix(name, SPATIAL_R)}_R1_filtered.fastq") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 0:
                assigned_ids.add(line[1:].split()[0].partition("_calledidx_")[0])

    candidates = []
    with gzip.open(raw_fastq(name, 1), "rt") as handle:
        rid = None
        for i, line in enumerate(handle):
            if i % 4 == 0:
                rid = line.split()[0][1:]
            elif i % 4 == 1:
                if rid in assigned_ids:
                    continue
                seq = line.strip()
                if insert_class(insert_end(seq)) != "full":
                    continue
                candidates.append(seq[UMI_LEN : UMI_LEN + BARCODE_LEN])
    n_unassigned_full = len(candidates)
    if n_unassigned_full > SUBSAMPLE:
        pick = RNG.choice(n_unassigned_full, SUBSAMPLE, replace=False)
        queries = [candidates[i] for i in pick]
    else:
        queries = candidates

    queries = [q for q in queries if len(q) == BARCODE_LEN and "N" not in q]
    lib = pack_barcodes(lib_seqs)
    qenc = pack_barcodes(queries)
    half1 = defaultdict(list)
    half2 = defaultdict(list)
    for i, s in enumerate(lib_seqs):
        half1[s[:SEG]].append(i)
        half2[s[SEG:]].append(i)

    rows = []
    half_counts = Counter()
    chunk = 32
    for start in range(0, len(queries), chunk):
        qs = qenc[start : start + chunk]
        d_h = (qs[:, None, :] != lib[None, :, :]).sum(axis=2)
        for k, q in enumerate(queries[start : start + chunk]):
            dh = d_h[k]
            keep = np.flatnonzero(dh <= HAMMING_PREFILTER)
            if len(keep) == 0:
                keep = np.array([int(dh.argmin())])
            elif len(keep) > 200:
                keep = keep[np.argpartition(dh[keep], 200)[:200]]
            best_i, best_d = None, 99
            for i in keep:
                d = Levenshtein.distance(q, lib_seqs[int(i)])
                if d < best_d:
                    best_d, best_i = d, int(i)
            rows.append({
                "nearest_barcode": lib_seqs[best_i],
                "levenshtein": best_d,
                "hamming": int(dh[best_i]),
                "rescued_r7": best_d <= 7,
                "rescued_r10": best_d <= 10,
            })
            h1, h2 = q[:SEG], q[SEG:]
            m1, m2 = h1 in half1, h2 in half2
            if m1 and m2:
                half_counts["both_halves"] += 1
            elif m1:
                half_counts["half1_only"] += 1
            elif m2:
                half_counts["half2_only"] += 1
            else:
                half_counts["neither_half"] += 1

    rescue = pd.DataFrame(rows)
    rescue["layout"] = name
    summary = {
        "layout": name,
        "unassigned_full_length": n_unassigned_full,
        "subsample": len(rescue),
        "frac_lev_le7": float((rescue.levenshtein <= 7).mean()) if len(rescue) else np.nan,
        "frac_lev_le10": float((rescue.levenshtein <= 10).mean()) if len(rescue) else np.nan,
        "median_levenshtein": float(rescue.levenshtein.median()) if len(rescue) else np.nan,
        **{k: half_counts[k] / max(len(rescue), 1) for k in ("both_halves", "half1_only", "half2_only", "neither_half")},
    }
    return rescue, summary


def main():
    sim_frames, rescue_frames, summaries = [], [], []
    for name in ARRAYS:
        print("simulation recall", name, flush=True)
        layout = load_layout(name)
        sim = simulation_recall(name, layout)
        sim["layout"] = name
        sim_frames.append(sim)
        print("rescue subsample", name, flush=True)
        rescue, summary = rescue_subsample(name, layout)
        rescue_frames.append(rescue)
        summaries.append(summary)
        print(summary, flush=True)
    pd.concat(sim_frames, ignore_index=True).to_csv(OUT_DIR / "simulation_recall.csv", index=False)
    pd.concat(rescue_frames, ignore_index=True).to_csv(OUT_DIR / "decoding_rescue_subsample.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUT_DIR / "decoding_rescue_summary.csv", index=False)


if __name__ == "__main__":
    main()
