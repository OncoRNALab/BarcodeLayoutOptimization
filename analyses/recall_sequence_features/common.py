"""Shared paths, loaders and helpers for the recall / sequence-feature analyses."""
from functools import lru_cache
from pathlib import Path
import gzip
import re

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "data"
LAYOUT_DIR = DATA / "1in1"
QUIK_DIR = DATA / "quik_out_manuscript"
RAW_DIR = DATA / "raw_fastq_files"
DECOY_DIR = DATA / "decoys" / "output"
OUT = HERE / "output"
OUT.mkdir(exist_ok=True)

ROWS, COLS, SEG, BC_LEN, UMI_LEN = 235, 363, 18, 36, 10
THRESHOLDS = [5, 6, 7]
CYCLE_ORDER = "ACGT"
READ2_HANDLE = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
READ1_HANDLE = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ADAPTER_PROBE = READ2_HANDLE[:13]
FULL_LENGTH_MIN_INSERT = 35

ARRAYS = {
    "Dual barcodes": {"key": "dual", "sample": "S1", "quik": "dual_combined_decoy_thr{r}",
                      "sim": "dual", "whitelist": "layout_235x363_dual_4x4.txt"},
    "Greedy + 2-OPT": {"key": "greedy", "sample": "S3", "quik": "greedy_combined_decoy_thr{r}",
                       "sim": "greedy", "whitelist": "layout_235x363_randombarcodes_opt.txt"},
    "Zipcodes": {"key": "zipcodes", "sample": "S4", "quik": "zipcodesV2_combined_decoy_thr{r}",
                 "sim": "zipV2", "whitelist": "layout_235x363_zipcodes_nosep.txt"},
}
COLORS = {"Dual barcodes": "tab:blue", "Greedy + 2-OPT": "tab:orange", "Zipcodes": "tab:green"}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def setup_matplotlib():
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    return plt


@lru_cache(maxsize=None)
def load_layout(name):
    """Printed layout as (row, col, barcode), barcodes upper case, in layout-file order."""
    if name == "Dual barcodes":
        df = pd.read_csv(LAYOUT_DIR / "layout_235x363_dual_4x4.csv").rename(columns={"combined": "barcode"})
    elif name == "Zipcodes":
        df = pd.read_csv(LAYOUT_DIR / "layout_235x363_zipcodes_nosep.csv")
    else:
        df = pd.read_csv(LAYOUT_DIR / "layout_235x363_randombarcodes_opt_coords.csv", header=None,
                         names=["barcode", "col", "row"])
    df["barcode"] = df.barcode.str.upper()
    assert df.barcode.is_unique and len(df) == ROWS * COLS
    return df[["row", "col", "barcode"]].reset_index(drop=True)


@lru_cache(maxsize=None)
def load_whitelist(name, r=7):
    """Combined whitelist (real + decoys) used by QUIK, read from the path in its log."""
    text = (quik_prefix(name, r).parent / f"26_018_{ARRAYS[name]['sample']}_quik_barcode_calling_stats.txt").read_text()
    size = int(re.search(r"Loaded (\d+) barcodes", text)[1])
    real = set(load_layout(name).barcode)
    decoy_file = {"Dual barcodes": DECOY_DIR / "dual" / "dual_decoys.txt",
                  "Zipcodes": DECOY_DIR / "zipcodes_balanced_khalf" / "zip_decoysv2.txt",
                  "Greedy + 2-OPT": DECOY_DIR / "greedy" / "greedy_decoys.txt"}[name]
    decoys = [l.strip().upper() for l in open(decoy_file) if l.strip()]
    assert len(real) + len(decoys) == size, (name, len(real), len(decoys), size)
    return sorted(real), decoys


def quik_prefix(name, r):
    cfg = ARRAYS[name]
    return QUIK_DIR / cfg["quik"].format(r=r) / f"26_018_{cfg['sample']}_quik"


@lru_cache(maxsize=None)
def load_raw(name):
    """All raw read pairs: read_id, r1, q1, r2, q2, umi, insert_len (-1 if the adapter is not in Read 1)."""
    s = ARRAYS[name]["sample"]
    rec = {"read_id": [], "r1": [], "q1": [], "r2": [], "q2": []}
    with gzip.open(RAW_DIR / f"26_018_{s}_{s}_L001_R1_001.fastq.gz", "rt") as f1, \
            gzip.open(RAW_DIR / f"26_018_{s}_{s}_L001_R2_001.fastq.gz", "rt") as f2:
        for i, (a, b) in enumerate(zip(f1, f2)):
            k = i % 4
            if k == 0:
                rid = a.split()[0][1:]
                assert rid == b.split()[0][1:]
                rec["read_id"].append(rid)
            elif k == 1:
                rec["r1"].append(a.strip())
                rec["r2"].append(b.strip())
            elif k == 3:
                rec["q1"].append(a.strip())
                rec["q2"].append(b.strip())
    df = pd.DataFrame(rec)
    df["umi"] = df.r1.str[:UMI_LEN]
    df["insert_len"] = df.r1.str.find(ADAPTER_PROBE)
    df["full_length"] = (df.insert_len < 0) | (df.insert_len >= FULL_LENGTH_MIN_INSERT)
    return df


@lru_cache(maxsize=None)
def load_assignments(name, r):
    """Reads assigned by QUIK: read_id, called_seq, is_real."""
    ids, seqs = [], []
    with open(f"{quik_prefix(name, r)}_R1_filtered.fastq") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                head = line.split()[0][1:]
                rid, _, rest = head.partition("_calledidx_")
                ids.append(rid)
                seqs.append(rest.split("_", 1)[1])
    df = pd.DataFrame({"read_id": ids, "called_seq": seqs})
    df["is_real"] = df.called_seq.isin(set(load_layout(name).barcode))
    return df


def dead_band_start(grid, ref_cols=300, frac_below=0.95, rel_threshold=0.1):
    col_mean = grid.mean(axis=0)
    below = col_mean < rel_threshold * np.median(col_mean[:ref_cols])
    for c in range(grid.shape[1]):
        if below[c:].mean() >= frac_below:
            return c
    return grid.shape[1]


@lru_cache(maxsize=None)
def feature_counts(name, r=7):
    """Per-feature read counts and UMI-deduplicated counts (exact UMI), with printed-area flags."""
    lay = load_layout(name).copy()
    asg = load_assignments(name, r)
    asg = asg[asg.is_real].merge(load_raw(name)[["read_id", "umi", "full_length"]], on="read_id", how="left")
    lay["reads"] = asg.called_seq.value_counts().reindex(lay.barcode).fillna(0).astype(int).values
    lay["reads_full_length"] = (asg[asg.full_length].called_seq.value_counts()
                                .reindex(lay.barcode).fillna(0).astype(int).values)
    lay["umis"] = (asg.drop_duplicates(["called_seq", "umi"]).called_seq.value_counts()
                   .reindex(lay.barcode).fillna(0).astype(int).values)
    grid = to_grid(lay, "reads")
    start = dead_band_start(grid)
    if name == "Dual barcodes":
        start = 4 * (start // 4)
    lay["printed"] = lay.col < start
    return lay, start


def to_grid(df, column, fill=0):
    grid = np.full((ROWS, COLS), fill, dtype=float)
    grid[df.row.values, df.col.values] = df[column].values
    return grid


def local_median(grid, mask, half=2):
    """Median of the (2*half+1)^2 window around each feature, excluding the feature itself and masked-out cells."""
    out = np.full(grid.shape, np.nan)
    padded = np.pad(np.where(mask, grid, np.nan), half, constant_values=np.nan)
    windows = np.lib.stride_tricks.sliding_window_view(padded, (2 * half + 1, 2 * half + 1)).copy()
    windows[:, :, half, half] = np.nan
    out[:] = np.nanmedian(windows.reshape(*grid.shape, -1), axis=2)
    return out


def local_mean(grid, mask, half=2):
    padded = np.pad(np.where(mask, grid, np.nan), half, constant_values=np.nan)
    windows = np.lib.stride_tricks.sliding_window_view(padded, (2 * half + 1, 2 * half + 1)).copy()
    windows[:, :, half, half] = np.nan
    return np.nanmean(windows.reshape(*grid.shape, -1), axis=2)


def leftmost_schedule(seq, order=CYCLE_ORDER):
    """Cycles (0-based) in which each base of seq is coupled under a repeating deposition order."""
    cycles, c = [], -1
    n = len(order)
    for base in seq:
        c += 1
        while order[c % n] != base:
            c += 1
        cycles.append(c)
    return cycles


def fold_batch(seqs):
    """MFE (kcal/mol) of each barcode alone and inside the probe (Read 1 handle + N10 + barcode + Read 2 handle)."""
    import RNA
    RNA.params_load_DNA_Mathews2004()
    return [(RNA.fold(s)[1], RNA.fold(READ1_HANDLE + "N" * UMI_LEN + s + READ2_HANDLE)[1]) for s in seqs]


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out
