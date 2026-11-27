#!/usr/bin/env python3
"""Per-barcode sequence, synthesis and decoder features for the three printed layouts.

Writes barcode_features.csv (one row per feature, all layouts) and greedy_bubbles.csv.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import json

import numpy as np
import pandas as pd

from io_utils import (
    ARRAYS,
    BARCODE_LEN,
    CYCLE_ORDER,
    NEXTERA_R1,
    NEXTERA_R2,
    OUT_DIR,
    ROWS,
    COLS,
    SEG,
    SPATIAL_R,
    THRESHOLDS,
    UMI_LEN,
    dead_band_start,
    eight_neighbours,
    grid_from,
    load_layout,
    run_prefix,
)

ORDER = {b: i for i, b in enumerate(CYCLE_ORDER)}
BASE_TO_INT = np.full(256, 255, dtype=np.uint8)
BASE_TO_INT[ord("A")] = 0
BASE_TO_INT[ord("C")] = 1
BASE_TO_INT[ord("G")] = 2
BASE_TO_INT[ord("T")] = 3


def encode_seqs(seqs: pd.Series) -> np.ndarray:
    arr = np.frombuffer("".join(seqs.tolist()).encode("ascii"), dtype=np.uint8)
    return BASE_TO_INT[arr].reshape(len(seqs), BARCODE_LEN)


def nearest_hamming(layout: pd.DataFrame, enc: np.ndarray, sample: int = 5000, seed: int = 2026) -> tuple[np.ndarray, np.ndarray]:
    """Exact Hamming to spatial neighbours and same-half partners, plus a library sample."""
    n = len(enc)
    mins = np.full(n, 99, dtype=np.int16)
    n_le8 = np.zeros(n, dtype=np.int32)
    index = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(layout.row, layout.col))}
    for i, (r, c) in enumerate(zip(layout.row.to_numpy(), layout.col.to_numpy())):
        js = [index[(rr, cc)] for rr, cc in eight_neighbours(int(r), int(c))]
        if js:
            d = (enc[i] != enc[js]).sum(axis=1)
            mins[i] = min(mins[i], int(d.min()))
            n_le8[i] += int((d <= 8).sum())

    def update_groups(labels):
        for _, idx in pd.Series(np.arange(n)).groupby(labels, sort=False):
            ix = idx.to_numpy()
            if len(ix) < 2 or len(ix) > 400:
                continue
            sub = enc[ix]
            d = (sub[:, None, :] != sub[None, :, :]).sum(axis=2, dtype=np.int16)
            np.fill_diagonal(d, 99)
            mins[ix] = np.minimum(mins[ix], d.min(axis=1))
            n_le8[ix] += (d <= 8).sum(axis=1)

    update_groups(layout.bar_1.to_numpy())
    update_groups(layout.bar_2.to_numpy())

    rng = np.random.default_rng(seed)
    ref = rng.choice(n, size=min(sample, n), replace=False)
    ref_enc = enc[ref]
    chunk = 200
    for i in range(0, n, chunk):
        a = enc[i : i + chunk]
        d = (a[:, None, :] != ref_enc[None, :, :]).sum(axis=2, dtype=np.int16)
        for k, gi in enumerate(range(i, i + len(a))):
            d[k, ref == gi] = 99
        mins[i : i + len(a)] = np.minimum(mins[i : i + len(a)], d.min(axis=1))
        n_le8[i : i + len(a)] += (d <= 8).sum(axis=1)
    return mins, n_le8


def composition_features(seq: str) -> dict:
    n = len(seq)
    counts = Counter(seq)
    runs, run = [], 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            run += 1
        else:
            runs.append(run)
            run = 1
    runs.append(run)
    g_run = 0
    g_score = 0.0
    i = 0
    while i < n:
        if seq[i] == "G":
            j = i
            while j < n and seq[j] == "G":
                j += 1
            L = j - i
            g_run = max(g_run, L)
            g_score += L * L
            i = j
        else:
            i += 1
    return {
        "A_frac": counts["A"] / n,
        "C_frac": counts["C"] / n,
        "G_frac": counts["G"] / n,
        "T_frac": counts["T"] / n,
        "GC_frac": (counts["G"] + counts["C"]) / n,
        "purine_frac": (counts["A"] + counts["G"]) / n,
        "G_frac_half1": seq[:SEG].count("G") / SEG,
        "G_frac_half2": seq[SEG:].count("G") / SEG,
        "longest_homopolymer": max(runs),
        "GG_count": seq.count("GG"),
        "GGG_count": seq.count("GGG"),
        "G4_score": g_score / n,
        "longest_G_run": g_run,
    }


def schedule_cycles(seq: str) -> list[int]:
    """Left-most embedding on the repeating A, C, G, T deposition order.

    Returns the 0-based cycle indices at which the barcode is extended.
    """
    t = 0
    hits = []
    for base in seq:
        need = ORDER[base]
        wait = (need - (t % 4)) % 4
        t = t + wait
        hits.append(t)
        t = t + 1
    return hits


def schedule_bitset(seq: str, n_cycles: int) -> np.ndarray:
    bits = np.zeros(n_cycles, dtype=np.uint8)
    bits[schedule_cycles(seq)] = 1
    return bits


def adapter_complementarity(seq: str) -> dict:
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    adapters = NEXTERA_R1 + "N" + NEXTERA_R2
    best = 0
    for k in range(6, 13):
        hits = any(rc[i : i + k] in adapters for i in range(0, len(rc) - k + 1))
        if hits:
            best = k
    self_best = 0
    for k in range(4, 9):
        if any(seq[i : i + k] in rc for i in range(0, len(seq) - k + 1)):
            self_best = k
    return {"adapter_rc_k": best, "self_rc_k": self_best}


def kmer_sharing(seqs: pd.Series, k: int = 4) -> np.ndarray:
    freq = Counter()
    windows = []
    for s in seqs:
        kmers = [s[i : i + k] for i in range(len(s) - k + 1)]
        windows.append(kmers)
        freq.update(set(kmers))
    return np.array([np.mean([freq[m] for m in kmers]) for kmers in windows])


def local_median(grid: np.ndarray, k: int = 5) -> np.ndarray:
    from scipy.ndimage import median_filter

    loc = median_filter(grid.astype(float), size=k, mode="nearest")
    # remove the contribution of the centre cell when it is an outlier of the window
    return loc


def connected_components(mask: np.ndarray) -> list[np.ndarray]:
    seen = np.zeros_like(mask, dtype=bool)
    comps = []
    h, w = mask.shape
    for i in range(h):
        for j in range(w):
            if not mask[i, j] or seen[i, j]:
                continue
            stack = [(i, j)]
            seen[i, j] = True
            cells = []
            while stack:
                r, c = stack.pop()
                cells.append((r, c))
                for rr, cc in eight_neighbours(r, c):
                    if mask[rr, cc] and not seen[rr, cc]:
                        seen[rr, cc] = True
                        stack.append((rr, cc))
            comps.append(np.array(cells))
    return comps


def greedy_bubbles(grid: np.ndarray, printed_cols: int, min_size: int = 80) -> np.ndarray:
    """Dark connected regions in the printed area of the Greedy array."""
    area = grid[:, :printed_cols].astype(float)
    med = np.median(area)
    from scipy.ndimage import uniform_filter

    local = uniform_filter(area, size=11, mode="nearest")
    mask = np.zeros_like(grid, dtype=bool)
    mask[:, :printed_cols] = local < 0.15 * max(med, 1e-6)
    keep = np.zeros_like(mask)
    for comp in connected_components(mask):
        if len(comp) >= min_size:
            keep[comp[:, 0], comp[:, 1]] = True
    return keep


def count_assigned(name: str, layout: pd.DataFrame) -> pd.DataFrame:
    real = set(layout.barcode)
    lookup = {b: i for i, b in enumerate(layout.barcode)}
    out = layout[["row", "col", "barcode"]].copy()
    for r in THRESHOLDS:
        reads = np.zeros(len(layout), dtype=np.int32)
        umis = [set() for _ in range(len(layout))]
        short = np.zeros(len(layout), dtype=np.int32)
        path = f"{run_prefix(name, r)}_R1_filtered.fastq"
        with open(path) as handle:
            for i, line in enumerate(handle):
                if i % 4 == 0:
                    seq = line.rsplit("_", 1)[-1].split()[0]
                elif i % 4 == 1:
                    if seq not in real:
                        continue
                    idx = lookup[seq]
                    reads[idx] += 1
                    umis[idx].add(line.strip()[:UMI_LEN])
                    if ADAPTER_IN_INSERT(line):
                        short[idx] += 1
        out[f"reads_r{r}"] = reads
        out[f"umis_r{r}"] = [len(s) for s in umis]
        out[f"short_assigned_r{r}"] = short
    return out


def ADAPTER_IN_INSERT(seq_line: str) -> bool:
    return "CTGTCTCTTATAC" in seq_line


def neighbour_schedule_distance(layout: pd.DataFrame, bits: np.ndarray) -> np.ndarray:
    index = {(r, c): i for i, (r, c) in enumerate(zip(layout.row, layout.col))}
    dist = np.zeros(len(layout), dtype=float)
    for i, (r, c) in enumerate(zip(layout.row.to_numpy(), layout.col.to_numpy())):
        neigh = [index[(rr, cc)] for rr, cc in eight_neighbours(r, c)]
        dist[i] = np.mean([np.abs(bits[i] - bits[j]).sum() for j in neigh])
    return dist


def build_layout_features(name: str) -> pd.DataFrame:
    layout = load_layout(name)
    print(f"{name}: counting assigned reads", flush=True)
    counts = count_assigned(name, layout)
    df = layout.merge(counts.drop(columns=["row", "col"]), on="barcode")
    grid7 = grid_from(df, df.reads_r7)
    start = dead_band_start(grid7)
    if name == "Dual barcodes":
        start = 4 * (start // 4)
    df["dead_start"] = start
    df["printed"] = df.col < start
    df["border"] = (df.row < 2) | (df.row >= ROWS - 2) | (df.col < 2) | (df.col >= start - 2)
    bubble = np.zeros(len(df), dtype=bool)
    if name == "Greedy + 2-OPT":
        mask = greedy_bubbles(grid7, start)
        bubble = mask[df.row.to_numpy(), df.col.to_numpy()]
        pd.DataFrame({"row": np.where(mask)[0], "col": np.where(mask)[1]}).to_csv(
            OUT_DIR / "greedy_bubbles.csv", index=False
        )
    df["bubble"] = bubble
    df["analysis_set"] = df.printed & ~df.border & ~df.bubble

    loc = local_median(grid7)
    df["local_median_r7"] = loc[df.row.to_numpy(), df.col.to_numpy()]
    df["local_ratio_r7"] = df.reads_r7 / df.local_median_r7.replace(0, np.nan)

    print(f"{name}: nearest Hamming", flush=True)
    enc = encode_seqs(df.barcode)
    nn, n8 = nearest_hamming(df, enc)
    df["nn_hamming"] = nn
    df["n_hamming_le8"] = n8
    df["mean_4mer_share"] = kmer_sharing(df.barcode)

    cycles = df.barcode.map(lambda s: schedule_cycles(s))
    df["n_cycles"] = cycles.map(len)  # always 36
    df["n_synthesis_cycles"] = cycles.map(lambda xs: xs[-1] + 1)
    df["cycles_per_nt"] = df.n_synthesis_cycles / BARCODE_LEN
    max_t = int(df.n_synthesis_cycles.max())
    bits = np.stack([schedule_bitset(s, max_t) for s in df.barcode])
    df["neighbour_schedule_hamming"] = neighbour_schedule_distance(df, bits)

    comp = pd.DataFrame(df.barcode.map(composition_features).tolist())
    adp = pd.DataFrame(df.barcode.map(adapter_complementarity).tolist())
    df = pd.concat([df.reset_index(drop=True), comp, adp], axis=1)
    df["layout"] = name
    return df


def main():
    frames = [build_layout_features(name) for name in ARRAYS]
    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_csv(OUT_DIR / "barcode_features.csv", index=False)
    meta = {
        "n_rows": int(len(all_df)),
        "dead_start": {n: int(all_df.loc[all_df.layout == n, "dead_start"].iloc[0]) for n in ARRAYS},
        "analysis_set": {n: int(all_df.loc[all_df.layout == n, "analysis_set"].sum()) for n in ARRAYS},
        "bubbles": int(all_df.bubble.sum()),
    }
    (OUT_DIR / "barcode_features_summary.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
