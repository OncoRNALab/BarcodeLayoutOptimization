#!/usr/bin/env python3
"""Classify every raw read by insert length and match QUIK assignments.

Writes read_accounting.csv (per layout and threshold) and
short_insert_attribution.csv (prefix/suffix matches of short inserts).
"""

from __future__ import annotations

import gzip
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from io_utils import (
    ADAPTER_SEED,
    ARRAYS,
    NEXTERA_R1,
    NEXTERA_R2,
    OUT_DIR,
    SEG,
    THRESHOLDS,
    UMI_LEN,
    load_layout,
    raw_fastq,
    run_prefix,
)

RC = str.maketrans("ACGTN", "TGCAN")


def insert_end(seq: str) -> int | None:
    """0-based start of the Read 2 handle in Read 1, or None if not found."""
    p = seq.find(ADAPTER_SEED)
    return None if p < 0 else p


def insert_class(p: int | None) -> str:
    if p is None:
        return "full"
    if p <= 10:
        return "le10"
    if p <= 25:
        return "11_25"
    return "26_34"


def load_assignments(name: str, r: int, real: set[str]) -> dict[str, tuple[str, bool]]:
    """read_id (Illumina) -> (called_seq, is_real)."""
    out = {}
    with open(f"{run_prefix(name, r)}_R1_filtered.fastq") as handle:
        for i, line in enumerate(handle):
            if i % 4:
                continue
            body = line[1:].split()[0]
            rid, _, rest = body.partition("_calledidx_")
            seq = rest.rsplit("_", 1)[-1]
            out[rid] = (seq, seq in real)
    return out


def is_adapter_dimer(seq: str) -> bool:
    return seq.startswith(NEXTERA_R2[:12]) or seq.startswith(NEXTERA_R1[:12])


def unique_prefix_map(barcodes: list[str], length: int) -> dict[str, str]:
    buckets: dict[str, list[str]] = defaultdict(list)
    for b in barcodes:
        buckets[b[:length]].append(b)
    return {p: v[0] for p, v in buckets.items() if len(v) == 1}


def unique_suffix_map(barcodes: list[str], length: int) -> dict[str, str]:
    buckets: dict[str, list[str]] = defaultdict(list)
    for b in barcodes:
        buckets[b[-length:]].append(b)
    return {p: v[0] for p, v in buckets.items() if len(v) == 1}


def analyse_layout(name: str) -> tuple[list[dict], dict]:
    layout = load_layout(name)
    real = set(layout.barcode)
    barcodes = layout.barcode.tolist()
    c1 = layout.bar_1.value_counts()
    unique_h1 = set(c1[c1 == 1].index)

    assigned = {r: load_assignments(name, r, real) for r in THRESHOLDS}

    class_counts = Counter()
    dimer = 0
    n = 0
    short_assigned = {r: Counter() for r in THRESHOLDS}  # real/decoy/unassigned
    full_assigned = {r: Counter() for r in THRESHOLDS}
    prefix_hit = Counter()
    suffix_hit = Counter()
    half_hit = Counter()
    short_n = 0
    umi_entropy_bins = Counter()

    prefix_maps = {L: unique_prefix_map(barcodes, L) for L in (8, 12, 18)}
    suffix_maps = {L: unique_suffix_map(barcodes, L) for L in (8, 12, 18)}

    with gzip.open(raw_fastq(name, 1), "rt") as handle:
        rid = None
        for i, line in enumerate(handle):
            if i % 4 == 0:
                rid = line.split()[0][1:]
            elif i % 4 == 1:
                seq = line.strip()
                n += 1
                p = insert_end(seq)
                cls = insert_class(p)
                class_counts[cls] += 1
                if is_adapter_dimer(seq):
                    dimer += 1
                for r in THRESHOLDS:
                    hit = assigned[r].get(rid)
                    bucket = full_assigned if cls == "full" else short_assigned
                    if hit is None:
                        bucket[r]["unassigned"] += 1
                    elif hit[1]:
                        bucket[r]["real"] += 1
                    else:
                        bucket[r]["decoy"] += 1
                if cls == "full":
                    continue
                short_n += 1
                frag = seq[UMI_LEN:p] if p is not None and p > UMI_LEN else ""
                umi = seq[: min(UMI_LEN, p if p is not None else UMI_LEN)]
                umi_entropy_bins["umi_has_N"] += ("N" in umi)
                if frag:
                    L = len(frag)
                    for key, mp in (("prefix8", prefix_maps[8]), ("prefix12", prefix_maps[12]), ("prefix18", prefix_maps[18])):
                        if L >= int(key.replace("prefix", "")) and frag[: int(key.replace("prefix", ""))] in mp:
                            prefix_hit[key] += 1
                    for key, mp in (("suffix8", suffix_maps[8]), ("suffix12", suffix_maps[12]), ("suffix18", suffix_maps[18])):
                        k = int(key.replace("suffix", ""))
                        if L >= k and frag[-k:] in mp:
                            suffix_hit[key] += 1
                    if L >= SEG:
                        if frag[:SEG] in unique_h1:
                            half_hit["unique_half1"] += 1
                        if frag[:SEG] in set(c1.index):
                            half_hit["any_half1"] += 1

    rows = []
    full_n = class_counts["full"]
    for r in THRESHOLDS:
        real_all = full_assigned[r]["real"] + short_assigned[r]["real"]
        real_full = full_assigned[r]["real"]
        decoy_all = full_assigned[r]["decoy"] + short_assigned[r]["decoy"]
        decoy_full = full_assigned[r]["decoy"]
        rows.append({
            "layout": name,
            "threshold": r,
            "reads": n,
            "full_length": full_n,
            "short_le10": class_counts["le10"],
            "short_11_25": class_counts["11_25"],
            "short_26_34": class_counts["26_34"],
            "adapter_dimer_like": dimer,
            "assigned_real_all": real_all,
            "assigned_decoy_all": decoy_all,
            "assigned_real_full": real_full,
            "assigned_decoy_full": decoy_full,
            "assigned_real_short": short_assigned[r]["real"],
            "assigned_decoy_short": short_assigned[r]["decoy"],
            "assignment_rate_all": real_all / n,
            "assignment_rate_full": real_full / full_n,
            "decoy_fraction_all": decoy_all / (real_all + decoy_all) if real_all + decoy_all else np.nan,
            "decoy_fraction_full": decoy_full / (real_full + decoy_full) if real_full + decoy_full else np.nan,
            "decoy_fraction_short": (
                short_assigned[r]["decoy"] / (short_assigned[r]["real"] + short_assigned[r]["decoy"])
                if short_assigned[r]["real"] + short_assigned[r]["decoy"] else np.nan
            ),
        })
    extra = {
        "layout": name,
        "short_inserts": short_n,
        "prefix8_unique": prefix_hit["prefix8"] / short_n if short_n else 0,
        "prefix12_unique": prefix_hit["prefix12"] / short_n if short_n else 0,
        "prefix18_unique": prefix_hit["prefix18"] / short_n if short_n else 0,
        "suffix8_unique": suffix_hit["suffix8"] / short_n if short_n else 0,
        "suffix12_unique": suffix_hit["suffix12"] / short_n if short_n else 0,
        "half1_any": half_hit["any_half1"] / short_n if short_n else 0,
        "half1_unique": half_hit["unique_half1"] / short_n if short_n else 0,
        "adapter_dimer_like": dimer / n,
    }
    return rows, extra


def main():
    rows, extras = [], []
    for name in ARRAYS:
        print("accounting", name, flush=True)
        r, e = analyse_layout(name)
        rows.extend(r)
        extras.append(e)
    pd.DataFrame(rows).to_csv(OUT_DIR / "read_accounting.csv", index=False)
    pd.DataFrame(extras).to_csv(OUT_DIR / "short_insert_attribution.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(pd.DataFrame(extras).to_string(index=False))


if __name__ == "__main__":
    main()
