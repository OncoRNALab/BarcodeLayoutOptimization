"""Shared paths, layout loading and dead-band detection for the recall analyses."""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
LAYOUT_DIR = DATA / "1in1"
QUIK_DIR = DATA / "quik_out_manuscript"
RAW_DIR = DATA / "raw_fastq_files"
SIM_DIR = QUIK_DIR
OUT_DIR = Path(__file__).resolve().parent

ROWS, COLS, SEG, BARCODE_LEN = 235, 363, 18, 36
THRESHOLDS = [5, 6, 7]
SPATIAL_R = 7
CYCLE_ORDER = "ACGT"
NEXTERA_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
NEXTERA_R2 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
ADAPTER_SEED = "CTGTCTCTTATAC"
UMI_LEN = 10

ARRAYS = {
    "Dual barcodes": {
        "key": "dual",
        "sample": "26_018_S1",
        "raw": "26_018_S1_S1",
        "dir": "dual_combined_decoy_thr{r}",
        "sim": "dual",
        "whitelist": "layout_235x363_dual_4x4.txt",
    },
    "Greedy + 2-OPT": {
        "key": "greedy",
        "sample": "26_018_S3",
        "raw": "26_018_S3_S3",
        "dir": "greedy_combined_decoy_thr{r}",
        "sim": "greedy",
        "whitelist": "layout_235x363_randombarcodes_opt.txt",
    },
    "Zipcodes": {
        "key": "zipcodes",
        "sample": "26_018_S4",
        "raw": "26_018_S4_S4",
        "dir": "zipcodesV2_combined_decoy_thr{r}",
        "sim": "zipV2",
        "whitelist": "layout_235x363_zipcodes_nosep.txt",
    },
}

CALLED_HEADER = re.compile(r"_calledidx_(\d+)_([ACGTN]+)\s")
SIM_HEADER = re.compile(r"^@read_(\d+)_barcode_(\d+)_calledidx_(\d+)_([ACGTN]+)\s*$")


def run_prefix(name: str, r: int) -> Path:
    cfg = ARRAYS[name]
    return QUIK_DIR / cfg["dir"].format(r=r) / f"{cfg['sample']}_quik"


def sim_prefix(name: str, r: int) -> Path:
    cfg = ARRAYS[name]
    return QUIK_DIR / f"simreads_thr{r}_1M" / f"quik_{cfg['sim']}_combined_1M"


def raw_fastq(name: str, read: int) -> Path:
    stem = ARRAYS[name]["raw"]
    return RAW_DIR / f"{stem}_L001_R{read}_001.fastq.gz"


def load_layout(name: str) -> pd.DataFrame:
    if name == "Dual barcodes":
        df = pd.read_csv(LAYOUT_DIR / "layout_235x363_dual_4x4.csv").rename(columns={"combined": "barcode"})
    elif name == "Zipcodes":
        df = pd.read_csv(LAYOUT_DIR / "layout_235x363_zipcodes_nosep.csv")
    else:
        df = pd.read_csv(
            LAYOUT_DIR / "layout_235x363_randombarcodes_opt_coords.csv",
            header=None,
            names=["barcode", "col", "row"],
        )
    df["barcode"] = df.barcode.str.upper()
    df["bar_1"] = df.barcode.str[:SEG]
    df["bar_2"] = df.barcode.str[SEG:]
    if name == "Dual barcodes":
        df["unit_1"] = df.row // 4
        df["unit_2"] = df.col // 4
    elif name == "Zipcodes":
        df["unit_1"] = df.col
        df["unit_2"] = df.row
    else:
        df["unit_1"] = np.arange(len(df))
        df["unit_2"] = np.arange(len(df))
    assert df.barcode.is_unique and len(df) == ROWS * COLS
    return df[["row", "col", "barcode", "bar_1", "bar_2", "unit_1", "unit_2"]]


def grid_from(df: pd.DataFrame, values: pd.Series) -> np.ndarray:
    grid = np.zeros((ROWS, COLS), dtype=values.dtype if values.dtype != object else float)
    grid[df.row.to_numpy(), df.col.to_numpy()] = values.to_numpy()
    return grid


def dead_band_start(grid: np.ndarray, ref_cols: int = 300, frac_below: float = 0.95, rel_threshold: float = 0.1) -> int:
    col_mean = grid.mean(axis=0)
    below = col_mean < rel_threshold * np.median(col_mean[:ref_cols])
    for c in range(grid.shape[1]):
        if below[c:].mean() >= frac_below:
            return c
    return grid.shape[1]


def eight_neighbours(row: int, col: int):
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            rr, cc = row + dr, col + dc
            if 0 <= rr < ROWS and 0 <= cc < COLS:
                yield rr, cc
