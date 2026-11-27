#!/usr/bin/env python3
"""Recompute Zipcodes metrics with the balanced 1:1 decoy whitelist (V2)."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, spearmanr

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "recall_sequence_features"))
from io_utils import dead_band_start  # noqa: E402

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent
QUIK = DATA / "quik_out_manuscript"
LAYOUT = DATA / "1in1" / "layout_235x363_zipcodes_nosep.csv"
REAL_TXT = DATA / "1in1" / "layout_235x363_zipcodes_nosep.txt"
DECOY_TXT = DATA / "decoys" / "output" / "zipcodes_balanced_khalf" / "zip_decoysv2.txt"
DECOY_CSV = DATA / "decoys" / "output" / "zipcodes_balanced_khalf" / "zip_decoysv2.csv"
RAW_R1 = DATA / "raw_fastq_files" / "26_018_S4_S4_L001_R1_001.fastq.gz"
ADAPTER = "CTGTCTCTTATAC"
SEG, ROWS, COLS, N, UMI_LEN = 18, 235, 363, 395_770, 10
PHI_DESIGN = 1 / 3
HEADER_REAL = re.compile(r"_calledidx_(\d+)_([ACGTN]+)\s")
HEADER_SIM = re.compile(r"^@read_(\d+)_barcode_(\d+)_calledidx_(\d+)_([ACGTN]+)\s*$")


def read_lines(path):
    return [line.strip().upper() for line in open(path) if line.strip()]


def precision(D, A, phi):
    p = 1 - (D / A) / phi
    D_lo = chi2.ppf(0.025, 2 * D) / 2 if D > 0 else 0.0
    D_hi = chi2.ppf(0.975, 2 * (D + 1)) / 2
    return p, 1 - (D_hi / A) / phi, 1 - (D_lo / A) / phi


def read_called_real(path):
    seqs, rids = [], []
    with open(path) as handle:
        for i, line in enumerate(handle):
            if i % 4:
                continue
            m = HEADER_REAL.search(line)
            if m is None:
                raise ValueError(line[:160])
            seqs.append(m[2])
            rids.append(line[1:].split()[0].partition("_calledidx_")[0])
    return pd.DataFrame({"rid": rids, "called_seq": seqs})


def read_called_sim(path):
    true_idx, called = [], []
    with open(path) as handle:
        for i, line in enumerate(handle):
            if i % 4:
                continue
            m = HEADER_SIM.match(line)
            if m is None:
                raise ValueError(line[:160])
            true_idx.append(int(m[2]))
            called.append(m[4])
    return pd.DataFrame({"true_idx": true_idx, "called_seq": called})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    layout = pd.read_csv(LAYOUT)
    layout["barcode"] = layout.barcode.str.upper()
    real = np.array(read_lines(REAL_TXT))
    real_set = set(real)
    decoy_set = set(read_lines(DECOY_TXT))
    col_codes = {s[:SEG] for s in real}
    row_codes = {s[SEG:] for s in real}
    n_real, n_decoy = len(real), len(decoy_set)
    phi_global = n_decoy / (n_real + n_decoy)
    bc_to_pos = dict(zip(layout.barcode, zip(layout.row.astype(int), layout.col.astype(int))))

    sim_rows, real_rows, axis_rows, land_rows = [], [], [], []
    counts_r7 = layout.copy()
    counts_r7["reads"] = 0
    assigned_by_r = {}

    for r in (5, 6, 7):
        sim = read_called_sim(QUIK / f"simreads_thr{r}_1M" / "quik_zipV2_combined_1M_R1_filtered.fastq")
        sim["true_seq"] = real[sim.true_idx.values]
        sim["is_decoy"] = sim.called_seq.isin(decoy_set)
        sim["in_real"] = sim.called_seq.isin(real_set)
        sim["correct"] = sim.called_seq == sim.true_seq
        A, D, C = len(sim), int(sim.is_decoy.sum()), int(sim.correct.sum())
        F = int(sim.in_real.sum()) - C
        true_p = C / A
        hat = 1 - (D / A) / PHI_DESIGN
        phi_eff = D / (D + F)
        sim_rows.append({
            "threshold": r, "assigned_any": A, "D": D, "F": F, "C": C,
            "true_precision": true_p, "hat_P_phi_design": hat, "phi_eff": phi_eff,
            "bias_pp": (hat - true_p) * 100,
        })
        wrong = sim[~sim.correct]
        first = wrong.called_seq.str[:SEG] == wrong.true_seq.str[:SEG]
        last = wrong.called_seq.str[SEG:] == wrong.true_seq.str[SEG:]
        keep_col = int((first & ~last).sum())
        keep_row = int((last & ~first).sum())
        land_rows.append({"threshold": r, "kept_column": keep_col, "kept_row": keep_row,
                          "neither_or_both": len(wrong) - keep_col - keep_row})

        rec = read_called_real(QUIK / f"zipcodesV2_combined_decoy_thr{r}" / "26_018_S4_quik_R1_filtered.fastq")
        rec["is_decoy"] = rec.called_seq.isin(decoy_set)
        rec["in_real"] = rec.called_seq.isin(real_set)
        A, D, Rreal = len(rec), int(rec.is_decoy.sum()), int(rec.in_real.sum())
        p, lo, hi = precision(D, A, PHI_DESIGN)
        p_g, _, _ = precision(D, A, phi_global)
        real_rows.append({
            "threshold": r, "reads": N, "assigned_any": A, "assigned_real": Rreal, "assigned_decoy": D,
            "assignment_rate_real": Rreal / N, "precision_phi_design": p,
            "precision_lo": lo, "precision_hi": hi, "precision_phi_global": p_g,
            "recall": A * p / N,
        })
        assigned_by_r[r] = rec
        dec = rec.loc[rec.is_decoy, "called_seq"]
        vr = int((dec.str[:SEG].isin(col_codes) & ~dec.str[SEG:].isin(row_codes)).sum())
        vc = int((~dec.str[:SEG].isin(col_codes) & dec.str[SEG:].isin(row_codes)).sum())
        sim_dec = sim.loc[sim.is_decoy, "called_seq"]
        svr = int((sim_dec.str[:SEG].isin(col_codes) & ~sim_dec.str[SEG:].isin(row_codes)).sum())
        svc = int((~sim_dec.str[:SEG].isin(col_codes) & sim_dec.str[SEG:].isin(row_codes)).sum())
        axis_rows.append({
            "threshold": r, "real_virtual_row": vr / D, "real_virtual_column": vc / D,
            "sim_virtual_row": svr / max(len(sim_dec), 1), "sim_virtual_column": svc / max(len(sim_dec), 1),
        })
        if r == 7:
            vc_counts = rec.loc[rec.in_real, "called_seq"].value_counts()
            counts_r7["reads"] = counts_r7.barcode.map(vc_counts).fillna(0).astype(int)

    grid = np.zeros((ROWS, COLS))
    for row, col, n in zip(counts_r7.row, counts_r7.col, counts_r7.reads):
        grid[int(row), int(col)] = n
    start = int(dead_band_start(grid))
    printed = counts_r7.col < start
    det_all = float((counts_r7.reads > 0).mean())
    det_print = float((counts_r7.loc[printed, "reads"] > 0).mean())

    decoy_table = pd.read_csv(DECOY_CSV)
    keep = np.where(decoy_table.axis == "column", decoy_table.gap + 1 < start, decoy_table.col < start)
    kept_decoys = set(decoy_table.barcode[keep].str.upper())
    printed_real = set(layout.loc[printed, "barcode"])

    printed_rows = []
    for r, rec in assigned_by_r.items():
        recp = rec[rec.called_seq.isin(printed_real) | rec.called_seq.isin(kept_decoys)]
        Ap, Dp, Rp = len(recp), int(recp.is_decoy.sum()), int(recp.in_real.sum())
        pp, lo, hi = precision(Dp, Ap, PHI_DESIGN)
        printed_rows.append({
            "threshold": r, "assigned_any": Ap, "assigned_real": Rp, "assigned_decoy": Dp,
            "assignment_rate_real": Rp / N, "precision_phi_design": pp,
            "precision_lo": lo, "precision_hi": hi, "recall": Ap * pp / N,
        })
    rec7p = printed_rows[2]

    col_g = layout.groupby("col").barcode.first().str[:SEG].map(lambda s: s.count("G") / SEG)
    row_g = layout.groupby("row").barcode.first().str[SEG:].map(lambda s: s.count("G") / SEG)
    col_mean_reads = counts_r7.loc[printed].groupby("col").reads.mean()
    row_mean_reads = counts_r7.loc[printed].groupby("row").reads.mean()
    rho_c = spearmanr(col_g.loc[col_mean_reads.index], col_mean_reads).correlation
    rho_r = spearmanr(row_g.loc[row_mean_reads.index], row_mean_reads).correlation

    edge = np.minimum.reduce([
        counts_r7.row, ROWS - 1 - counts_r7.row, counts_r7.col, COLS - 1 - counts_r7.col,
    ])
    area = printed & (edge >= 2)
    sub = counts_r7.loc[area]
    mean_reads = float(sub.reads.mean())
    frac_zero = float((sub.reads == 0).mean())
    pois = float(np.exp(-mean_reads))
    gfrac = sub.barcode.str.count("G") / 36.0
    zeros = gfrac[sub.reads == 0].to_numpy()
    high = gfrac[sub.reads > 4].to_numpy()
    # Cliff's delta: P(zero > high) - P(zero < high)
    n_cmp = len(zeros) * len(high)
    n_gt = int((zeros[:, None] > high).sum()) if len(zeros) and len(high) else 0
    n_lt = int((zeros[:, None] < high).sum()) if len(zeros) and len(high) else 0
    cliffs = (n_gt - n_lt) / n_cmp if n_cmp else np.nan
    has_ggg = sub.barcode.str.contains("GGG")
    zero_mask = sub.reads == 0
    high_mask = sub.reads > 4
    table = np.array([
        [int((has_ggg & zero_mask).sum()), int((~has_ggg & zero_mask).sum())],
        [int((has_ggg & high_mask).sum()), int((~has_ggg & high_mask).sum())],
    ])
    ggg_or, ggg_p = fisher_exact(table)
    ggg_zero = table[0, 0] / table[0].sum() if table[0].sum() else np.nan
    ggg_high = table[1, 0] / table[1].sum() if table[1].sum() else np.nan

    # Full-length vs short inserts
    assigned = {r: dict(zip(df.rid, zip(df.called_seq, df.in_real))) for r, df in assigned_by_r.items()}
    class_n = {"full": 0, "short": 0}
    full = {r: {"real": 0, "decoy": 0, "unassigned": 0} for r in (5, 6, 7)}
    short = {r: {"real": 0, "decoy": 0, "unassigned": 0} for r in (5, 6, 7)}
    n_raw = 0
    with gzip.open(RAW_R1, "rt") as handle:
        rid = None
        for i, line in enumerate(handle):
            if i % 4 == 0:
                rid = line.split()[0][1:]
            elif i % 4 == 1:
                n_raw += 1
                p = line.find(ADAPTER)
                bucket = full if p < 0 or p >= 35 else short
                class_n["full" if bucket is full else "short"] += 1
                for r in (5, 6, 7):
                    hit = assigned[r].get(rid)
                    if hit is None:
                        bucket[r]["unassigned"] += 1
                    elif hit[1]:
                        bucket[r]["real"] += 1
                    else:
                        bucket[r]["decoy"] += 1
    fl_rows = []
    for r in (5, 6, 7):
        nf = class_n["full"]
        ar = full[r]["real"]
        ashort = short[r]["real"]
        dshort = short[r]["decoy"]
        dfull = full[r]["decoy"]
        any_short = ashort + dshort
        any_all = real_rows[r - 5]["assigned_any"]
        fl_rows.append({
            "threshold": r, "full_length_reads": nf,
            "assigned_real_full": ar, "assigned_real_short": ashort,
            "assignment_rate_full_length": ar / nf,
            "decoy_fraction_short": dshort / any_short if any_short else np.nan,
            "decoy_fraction_all": real_rows[r - 5]["assigned_decoy"] / any_all,
        })

    report = {
        "n_real": n_real, "n_decoy": n_decoy, "phi_global": phi_global, "phi_design": PHI_DESIGN,
        "dead_band_start": int(start), "n_kept_decoys": len(kept_decoys),
        "detection_whole": det_all, "detection_printed": det_print,
        "precision_r7_printed_area": rec7p["precision_phi_design"],
        "assignment_r7_printed_real": rec7p["assignment_rate_real"],
        "recall_r7_printed": rec7p["recall"],
        "spearman_G_column": float(rho_c), "spearman_G_row": float(rho_r),
        "printed_mean_reads_r7": mean_reads, "printed_frac_zero": frac_zero,
        "poisson_zero": pois, "cliffs_delta_G_zero_vs_high": float(cliffs),
        "ggg_frac_zero": float(ggg_zero), "ggg_frac_high": float(ggg_high),
        "ggg_odds_ratio": float(ggg_or), "ggg_p": float(ggg_p),
        "raw_reads": n_raw, "full_length_reads": class_n["full"],
        "simulation": sim_rows, "real": real_rows, "printed_area": printed_rows,
        "axis": axis_rows, "landing": land_rows, "full_length": fl_rows,
    }
    (OUT / "zipcodes_v2_report.json").write_text(json.dumps(report, indent=2))
    pd.DataFrame(sim_rows).to_csv(OUT / "sim_precision.csv", index=False)
    pd.DataFrame(real_rows).to_csv(OUT / "real_precision.csv", index=False)
    pd.DataFrame(printed_rows).to_csv(OUT / "printed_area_precision.csv", index=False)
    pd.DataFrame(axis_rows).to_csv(OUT / "axis_split.csv", index=False)
    pd.DataFrame(fl_rows).to_csv(OUT / "full_length.csv", index=False)
    counts_r7.to_csv(OUT / "reads_per_feature_r7.csv", index=False)
    skip = {"simulation", "real", "axis", "landing", "full_length", "printed_area"}
    print(json.dumps({k: report[k] for k in report if k not in skip}, indent=2))
    print("sim", pd.DataFrame(sim_rows).to_string(index=False))
    print("real", pd.DataFrame(real_rows).to_string(index=False))
    print("printed", pd.DataFrame(printed_rows).to_string(index=False))
    print("axis", pd.DataFrame(axis_rows).to_string(index=False))
    print("full", pd.DataFrame(fl_rows).to_string(index=False))
    print("landing", pd.DataFrame(land_rows).to_string(index=False))


if __name__ == "__main__":
    main()
