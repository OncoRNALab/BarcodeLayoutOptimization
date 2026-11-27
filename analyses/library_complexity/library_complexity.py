#!/usr/bin/env python3
"""Library complexity of the three printed layouts (Methods §2.8).

For each layout this script:

1. Samples ``--subset-size`` barcodes (default 10,000) and, for each query,
   computes the Levenshtein distance to its nearest *other* barcode in the
   **full** printed library of 85,305 sequences.
2. Embeds every barcode as a 768-dimensional positional 4-mer count vector
   (three consecutive bins along the sequence × 256 tetramers), L2-normalises
   it, and records the maximum cosine similarity to any other barcode of the
   same library.
3. Records the Shannon entropy of each barcode's unbinned 4-mer counts.

Copy this file and the three layout ``.txt`` files to the cluster, then e.g.:

    python library_complexity.py \\
        --data-dir /path/to/1in1 \\
        --out-dir  $VSC_SCRATCH/library_complexity \\
        --workers 32

Dependencies: numpy, pandas, matplotlib, rapidfuzz (see analyses/requirements.txt).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist

LAYOUTS = {
    "Dual barcodes": "layout_235x363_dual_4x4.txt",
    "Greedy + 2-OPT": "layout_235x363_randombarcodes_opt.txt",
    "Zipcodes": "layout_235x363_zipcodes_nosep.txt",
}
COLORS = {
    "Dual barcodes": "tab:blue",
    "Greedy + 2-OPT": "tab:orange",
    "Zipcodes": "tab:green",
}
BARCODE_LEN = 36
K = 4
N_KMER = 4 ** K
N_BINS = 3
N_KMER_POS = BARCODE_LEN - K + 1  # 33
BIN_SIZE = N_KMER_POS // N_BINS  # 11
BASE_TO_INT = np.full(256, 255, dtype=np.uint8)
BASE_TO_INT[ord("A")] = 0
BASE_TO_INT[ord("C")] = 1
BASE_TO_INT[ord("G")] = 2
BASE_TO_INT[ord("T")] = 3


def default_data_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data" / "1in1"


def read_seqs(path: Path) -> list[str]:
    seqs = [line.strip().upper() for line in path.open() if line.strip()]
    if not seqs:
        raise ValueError(f"No sequences in {path}")
    bad = [s for s in seqs if len(s) != BARCODE_LEN or any(c not in "ACGT" for c in s)]
    if bad:
        raise ValueError(f"{path}: {len(bad)} sequences are not {BARCODE_LEN}-nt ACGT")
    if len(set(seqs)) != len(seqs):
        raise ValueError(f"{path}: duplicate barcodes")
    return seqs


def encode(seqs: list[str]) -> np.ndarray:
    raw = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)
    return BASE_TO_INT[raw].reshape(len(seqs), BARCODE_LEN)


def kmer_codes(enc: np.ndarray) -> np.ndarray:
    """Integer 4-mer codes, shape (n_barcodes, 33), values in 0..255."""
    return (
        enc[:, 0:N_KMER_POS] * 64
        + enc[:, 1 : N_KMER_POS + 1] * 16
        + enc[:, 2 : N_KMER_POS + 2] * 4
        + enc[:, 3 : N_KMER_POS + 3]
    ).astype(np.int32)


def positional_4mer_counts(codes: np.ndarray) -> np.ndarray:
    """Raw positional counts, shape (n, 3, 256)."""
    n = codes.shape[0]
    bins = np.minimum(np.arange(N_KMER_POS) // BIN_SIZE, N_BINS - 1)
    counts = np.zeros((n, N_BINS, N_KMER), dtype=np.float32)
    rows = np.arange(n)
    for j in range(N_KMER_POS):
        np.add.at(counts, (rows, bins[j], codes[:, j]), 1.0)
    return counts


def unbinned_4mer_counts(codes: np.ndarray) -> np.ndarray:
    n = codes.shape[0]
    counts = np.zeros((n, N_KMER), dtype=np.float32)
    rows = np.arange(n)
    for j in range(N_KMER_POS):
        np.add.at(counts, (rows, codes[:, j]), 1.0)
    return counts


def l2_normalize(x: np.ndarray, axis: int = 1) -> np.ndarray:
    nrm = np.linalg.norm(x, axis=axis, keepdims=True)
    nrm = np.maximum(nrm, 1e-12)
    return x / nrm


def shannon_entropy(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    p = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.log2(p, where=p > 0, out=np.zeros_like(p))
    return -(p * logp).sum(axis=1)


def nearest_levenshtein(
    queries: list[str],
    query_idx: np.ndarray,
    refs: list[str],
    *,
    batch: int,
    workers: int,
) -> np.ndarray:
    """Min Levenshtein from each query to the full library, excluding self."""
    mins = np.empty(len(queries), dtype=np.int32)
    for start in range(0, len(queries), batch):
        stop = min(start + batch, len(queries))
        d = cdist(
            queries[start:stop],
            refs,
            scorer=Levenshtein.distance,
            dtype=np.int32,
            workers=workers,
        )
        d[np.arange(stop - start), query_idx[start:stop]] = np.iinfo(np.int32).max
        mins[start:stop] = d.min(axis=1)
    return mins


def max_cosine_similarity(unit: np.ndarray, *, batch: int) -> np.ndarray:
    """Max cosine of each row of ``unit`` to any other row."""
    n = unit.shape[0]
    nn = np.empty(n, dtype=np.float32)
    for start in range(0, n, batch):
        stop = min(start + batch, n)
        sim = unit[start:stop] @ unit.T
        sim[np.arange(stop - start), np.arange(start, stop)] = -np.inf
        nn[start:stop] = sim.max(axis=1)
    return nn


def summarise(values: np.ndarray) -> dict:
    return {
        "n": int(values.size),
        "min": float(np.min(values)),
        "p25": float(np.percentile(values, 25)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "p75": float(np.percentile(values, 75)),
        "max": float(np.max(values)),
    }


def plot_figure(
    lev: dict[str, np.ndarray],
    cosine: dict[str, np.ndarray],
    entropy: dict[str, np.ndarray],
    out_pdf: Path,
) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12,
        "axes.titlesize": 13, "axes.labelsize": 12,
        "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
    })
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.4))

    ax = axes[0]
    max_d = max(int(v.max()) for v in lev.values())
    bins = np.arange(0, max_d + 2) - 0.5
    for name, values in lev.items():
        ax.hist(
            values,
            bins=bins,
            histtype="step",
            density=True,
            color=COLORS[name],
            label=name,
            linewidth=1.4,
        )
    ax.set_xlabel("Nearest-neighbour Levenshtein distance")
    ax.set_ylabel("Density")
    ax.set_title("(a)", loc="left", fontweight="bold")
    ax.legend(frameon=False)

    ax = axes[1]
    for name, values in cosine.items():
        ax.hist(
            values,
            bins=40,
            histtype="step",
            density=True,
            color=COLORS[name],
            label=name,
            linewidth=1.4,
        )
    ax.set_xlabel("Max positional 4-mer cosine similarity")
    ax.set_ylabel("Density")
    ax.set_title("(b)", loc="left", fontweight="bold")

    ax = axes[2]
    for name, values in entropy.items():
        ax.hist(
            values,
            bins=40,
            histtype="step",
            density=True,
            color=COLORS[name],
            label=name,
            linewidth=1.4,
        )
    ax.set_xlabel("Per-barcode 4-mer Shannon entropy (bits)")
    ax.set_ylabel("Density")
    ax.set_title("(c)", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=Path, default=default_data_dir(), help="Directory of printed layout .txt files")
    p.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent, help="Output directory")
    p.add_argument("--subset-size", type=int, default=10_000, help="Queries for Levenshtein NN (default 10000)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--batch", type=int, default=500, help="Query batch size for rapidfuzz cdist and cosine")
    p.add_argument("--workers", type=int, default=-1, help="rapidfuzz threads (-1 = all cores)")
    p.add_argument("--skip-levenshtein", action="store_true")
    p.add_argument("--skip-kmer", action="store_true")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument(
        "--plot-only",
        action="store_true",
        help="Rebuild the PDF from nn_levenshtein_subset.csv and kmer_complexity_full.csv in --out-dir",
    )
    return p.parse_args()


def load_results(out_dir: Path) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    lev_df = pd.read_csv(out_dir / "nn_levenshtein_subset.csv")
    kmer_df = pd.read_csv(out_dir / "kmer_complexity_full.csv")
    lev = {name: g.nn_levenshtein.to_numpy() for name, g in lev_df.groupby("layout", sort=False)}
    cosine = {name: g.max_positional_4mer_cosine.to_numpy() for name, g in kmer_df.groupby("layout", sort=False)}
    entropy = {name: g.shannon_entropy_4mer_bits.to_numpy() for name, g in kmer_df.groupby("layout", sort=False)}
    return lev, cosine, entropy


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        lev, cosine, entropy = load_results(args.out_dir)
        pdf = args.out_dir / "library_complexity.pdf"
        plot_figure(lev, cosine, entropy, pdf)
        print(f"Figure: {pdf}", flush=True)
        return
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    libraries = {}
    for name, fname in LAYOUTS.items():
        path = args.data_dir / fname
        print(f"Loading {name}: {path}", flush=True)
        libraries[name] = read_seqs(path)
        print(f"  {len(libraries[name]):,} barcodes", flush=True)

    summary: dict = {
        "subset_size": args.subset_size,
        "seed": args.seed,
        "data_dir": str(args.data_dir.resolve()),
        "layouts": {},
    }
    lev_dist: dict[str, np.ndarray] = {}
    cosine_dist: dict[str, np.ndarray] = {}
    entropy_dist: dict[str, np.ndarray] = {}
    lev_rows = []
    kmer_rows = []

    for name, seqs in libraries.items():
        n = len(seqs)
        layout_summary: dict = {"n_library": n}

        if not args.skip_levenshtein:
            n_query = min(args.subset_size, n)
            qidx = rng.choice(n, size=n_query, replace=False)
            qidx.sort()
            queries = [seqs[i] for i in qidx]
            print(
                f"{name}: Levenshtein NN, {n_query:,} queries vs {n:,} references, "
                f"batch={args.batch}, workers={args.workers}",
                flush=True,
            )
            t = time.time()
            nn = nearest_levenshtein(
                queries, qidx, seqs, batch=args.batch, workers=args.workers
            )
            print(f"  done in {time.time() - t:.1f}s; mean={nn.mean():.2f}", flush=True)
            lev_dist[name] = nn
            layout_summary["levenshtein_nn"] = summarise(nn)
            for i, d in zip(qidx, nn):
                lev_rows.append({"layout": name, "index": int(i), "barcode": seqs[i], "nn_levenshtein": int(d)})

        if not args.skip_kmer:
            print(f"{name}: positional 4-mer embedding and entropy", flush=True)
            t = time.time()
            codes = kmer_codes(encode(seqs))
            pos = positional_4mer_counts(codes).reshape(n, N_BINS * N_KMER)
            unit = l2_normalize(pos)
            cos = max_cosine_similarity(unit, batch=args.batch)
            ent = shannon_entropy(unbinned_4mer_counts(codes))
            print(
                f"  done in {time.time() - t:.1f}s; mean max cosine={cos.mean():.3f}, "
                f"mean entropy={ent.mean():.3f} bits",
                flush=True,
            )
            cosine_dist[name] = cos
            entropy_dist[name] = ent
            layout_summary["max_positional_4mer_cosine"] = summarise(cos)
            layout_summary["shannon_entropy_4mer_bits"] = summarise(ent)
            for i in range(n):
                kmer_rows.append(
                    {
                        "layout": name,
                        "index": i,
                        "barcode": seqs[i],
                        "max_positional_4mer_cosine": float(cos[i]),
                        "shannon_entropy_4mer_bits": float(ent[i]),
                    }
                )

        summary["layouts"][name] = layout_summary

    if lev_rows:
        pd.DataFrame(lev_rows).to_csv(args.out_dir / "nn_levenshtein_subset.csv", index=False)
    if kmer_rows:
        pd.DataFrame(kmer_rows).to_csv(args.out_dir / "kmer_complexity_full.csv", index=False)

    summary_path = args.out_dir / "library_complexity_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)

    if not args.no_plot and lev_dist and cosine_dist and entropy_dist:
        pdf = args.out_dir / "library_complexity.pdf"
        plot_figure(lev_dist, cosine_dist, entropy_dist, pdf)
        print(f"Figure: {pdf}", flush=True)

    print(f"Summary: {summary_path}", flush=True)
    print(f"Total wall time {time.time() - t0:.1f}s", flush=True)
    print(json.dumps({k: {m: v.get("mean") for m, v in s.items() if isinstance(v, dict)}
                      for k, s in summary["layouts"].items()}, indent=2))


if __name__ == "__main__":
    main()
