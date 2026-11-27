#!/usr/bin/env python3
"""Error profiles of assigned reads and R1/R2 concordance.

Writes error_profiles_by_group.csv, positional_errors.csv and
r1_r2_concordance.csv. PDFs are produced by error_profiles.ipynb.
"""

from __future__ import annotations

import gzip
from collections import defaultdict

import numpy as np
import pandas as pd

from io_utils import (
    ADAPTER_SEED,
    ARRAYS,
    BARCODE_LEN,
    OUT_DIR,
    SPATIAL_R,
    UMI_LEN,
    load_layout,
    raw_fastq,
    run_prefix,
)
from read_accounting import insert_end

MAX_ALIGN = 40_000
MAX_PAIRS = 20_000
RNG = np.random.default_rng(2026)


def nw_counts(ref: str, query: str) -> dict:
    """Return substitution, insertion and deletion counts plus a per-ref-position substitution mask."""
    n, m = len(ref), len(query)
    dp = np.zeros((n + 1, m + 1), dtype=np.int16)
    dp[:, 0] = np.arange(n + 1)
    dp[0, :] = np.arange(m + 1)
    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            dp[i, j] = min(
                dp[i - 1, j - 1] + (ri != query[j - 1]),
                dp[i - 1, j] + 1,
                dp[i, j - 1] + 1,
            )
    i, j = n, m
    subst = ins = dele = 0
    subst_pos = np.zeros(n, dtype=np.int8)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i, j] == dp[i - 1, j - 1] + (ref[i - 1] != query[j - 1]):
            if ref[i - 1] != query[j - 1]:
                subst += 1
                subst_pos[i - 1] = 1
            i -= 1
            j -= 1
        elif i > 0 and dp[i, j] == dp[i - 1, j] + 1:
            dele += 1
            i -= 1
        else:
            ins += 1
            j -= 1
    return {"substitutions": subst, "insertions": ins, "deletions": dele, "subst_pos": subst_pos}


def rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def analyse_layout(name: str, features: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    feat = features.loc[features.layout == name].set_index("barcode")
    groups = feat["count_group"]
    real = set(feat.index)

    n_by_group = defaultdict(int)
    short_by_group = defaultdict(int)
    subst = defaultdict(int)
    ins = defaultdict(int)
    dele = defaultdict(int)
    pos_subst = {g: np.zeros(BARCODE_LEN, dtype=int) for g in ("zero", "low", "medium", "high")}
    pos_n = {g: np.zeros(BARCODE_LEN, dtype=int) for g in pos_subst}
    aligned = 0

    records = []
    with open(f"{run_prefix(name, SPATIAL_R)}_R1_filtered.fastq") as handle:
        called = None
        for i, line in enumerate(handle):
            if i % 4 == 0:
                called = line.rsplit("_", 1)[-1].split()[0]
            elif i % 4 == 1:
                if called not in real:
                    continue
                g = groups.get(called, "medium")
                if g not in pos_subst:
                    g = "medium"
                n_by_group[g] += 1
                seq = line.strip()
                p = insert_end(seq)
                if p is not None and p < UMI_LEN + BARCODE_LEN:
                    short_by_group[g] += 1
                query = seq[UMI_LEN : UMI_LEN + BARCODE_LEN]
                if aligned < MAX_ALIGN and RNG.random() < 0.25:
                    c = nw_counts(called, query)
                    subst[g] += c["substitutions"]
                    ins[g] += c["insertions"]
                    dele[g] += c["deletions"]
                    pos_subst[g] += c["subst_pos"]
                    pos_n[g] += 1
                    aligned += 1
                    records.append({"group": g, **{k: c[k] for k in ("substitutions", "insertions", "deletions")}})

    rows = []
    for g in ("low", "medium", "high"):
        n = n_by_group[g]
        a = max(sum(1 for rec in records if rec["group"] == g), 1)
        rows.append({
            "layout": name,
            "group": g,
            "assigned_reads": n,
            "short_insert_frac": short_by_group[g] / n if n else np.nan,
            "mean_substitutions": subst[g] / a,
            "mean_insertions": ins[g] / a,
            "mean_deletions": dele[g] / a,
        })

    pos_rows = []
    for g, arr in pos_subst.items():
        denom = np.maximum(pos_n[g], 1)
        for i, (s, d) in enumerate(zip(arr, denom)):
            pos_rows.append({"layout": name, "group": g, "position": i + 1, "substitution_rate": s / d, "n": int(d)})

    # R1 vs R2 concordance on a subsample of full-length assigned reads
    assigned_call = {}
    with open(f"{run_prefix(name, SPATIAL_R)}_R1_filtered.fastq") as handle:
        for i, line in enumerate(handle):
            if i % 4:
                continue
            body = line[1:].split()[0]
            rid, _, rest = body.partition("_calledidx_")
            seq = rest.rsplit("_", 1)[-1]
            if seq in real:
                assigned_call[rid] = seq

    dists = []
    with gzip.open(raw_fastq(name, 1), "rt") as f1, gzip.open(raw_fastq(name, 2), "rt") as f2:
        n_pairs = 0
        for i, (a, b) in enumerate(zip(f1, f2)):
            if i % 4 == 0:
                rid = a.split()[0][1:]
            elif i % 4 == 1:
                if rid not in assigned_call:
                    continue
                r1, r2 = a.strip(), b.strip()
                if insert_end(r1) is not None:
                    continue
                q1 = r1[UMI_LEN : UMI_LEN + BARCODE_LEN]
                r2rc = rc(r2)
                best = min(
                    sum(x != y for x, y in zip(q1, r2rc[k : k + BARCODE_LEN]))
                    for k in range(0, 8)
                    if k + BARCODE_LEN <= len(r2rc)
                )
                dists.append(best)
                n_pairs += 1
                if n_pairs >= MAX_PAIRS:
                    break
    conc = {
        "layout": name,
        "n_pairs": len(dists),
        "median_r1_r2_hamming": float(np.median(dists)) if dists else np.nan,
        "frac_hamming_le3": float(np.mean(np.array(dists) <= 3)) if dists else np.nan,
        "frac_hamming_le6": float(np.mean(np.array(dists) <= 6)) if dists else np.nan,
    }
    return pd.DataFrame(rows), pd.DataFrame(pos_rows), conc


def add_count_groups(features: pd.DataFrame) -> pd.DataFrame:
    def group(n):
        if n == 0:
            return "zero"
        if n == 1:
            return "low"
        if n <= 4:
            return "medium"
        return "high"

    features = features.copy()
    features["count_group"] = features.reads_r7.map(group)
    return features


def main():
    features = add_count_groups(pd.read_csv(OUT_DIR / "barcode_features.csv"))
    profile_frames, pos_frames, concs = [], [], []
    for name in ARRAYS:
        print("error profiles", name, flush=True)
        prof, pos, conc = analyse_layout(name, features)
        profile_frames.append(prof)
        pos_frames.append(pos)
        concs.append(conc)
        print(conc, flush=True)
    pd.concat(profile_frames, ignore_index=True).to_csv(OUT_DIR / "error_profiles_by_group.csv", index=False)
    pd.concat(pos_frames, ignore_index=True).to_csv(OUT_DIR / "positional_errors.csv", index=False)
    pd.DataFrame(concs).to_csv(OUT_DIR / "r1_r2_concordance.csv", index=False)


if __name__ == "__main__":
    main()
