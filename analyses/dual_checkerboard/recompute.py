#!/usr/bin/env python3
"""Metrics for the half-density Dual checkerboard array (printed + medium-error sim)."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
QUIK = DATA / "quik_out_manuscript"
OUT = Path(__file__).resolve().parent
WL_PATH = DATA / "decoys" / "output" / "dual_1in2" / "dual_whitelist_real_plus_decoys.txt"
LAY_PATH = DATA / "1in2" / "layout_1in2_dual_4x4.csv"
TXT_PATH = DATA / "1in2" / "bar_1in2_dual_4x4.txt"
RAW_R1 = DATA / "raw_fastq_files" / "26_022_S1_S1_L001_R1_001.fastq.gz"
ADAPTER = "CTGTCTCTTATAC"
FULL_MIN = 35
N_SIM = 1_000_000
HEADER_REAL = re.compile(r"_calledidx_(\d+)_([ACGTN]+)")
HEADER_SIM = re.compile(r"^@read_(\d+)_barcode_(\d+)_calledidx_(\d+)_([ACGTN]+)")


def hat_p(D, A, phi):
    p = 1 - (D / A) / phi
    lo = 0.0 if D == 0 else chi2.ppf(0.025, 2 * D) / 2
    hi = chi2.ppf(0.975, 2 * (D + 1)) / 2
    return p, 1 - (hi / A) / phi, 1 - (lo / A) / phi


def main():
    wl = [line.strip().upper() for line in open(WL_PATH) if line.strip()]
    layout = pd.read_csv(LAY_PATH)
    layout["combined"] = layout.combined.str.upper()
    n_real = len(layout)
    n_decoy = len(wl) - n_real
    phi = n_decoy / len(wl)
    real_set = set(wl[:n_real])
    decoy_set = set(wl[n_real:])
    assert set(wl[:n_real]) == set(layout.combined)

    def from_pb(path):
        counts = pd.read_csv(path).read_count.to_numpy()
        assert len(counts) == len(wl)
        R = int(counts[:n_real].sum())
        D = int(counts[n_real:].sum())
        det = int((counts[:n_real] > 0).sum())
        return R, D, det, counts

    # raw read count and insert classes
    n_raw = 0
    insert_of = {}
    with gzip.open(RAW_R1, "rt") as handle:
        rid, seq = None, None
        for i, line in enumerate(handle):
            k = i % 4
            if k == 0:
                rid = line.split()[0][1:]
            elif k == 1:
                seq = line.strip()
                n_raw += 1
                p = seq.find(ADAPTER)
                insert_of[rid] = "full" if p < 0 or p >= FULL_MIN else "short"

    real_rows, fl_rows = [], []
    assigned_by_r = {}
    for r in (6, 7, 8):
        pb = QUIK / f"dual_checkerboard_decoy_thr{r}" / "26_022_S1_Dual_quik_per_barcode.csv"
        R, D, det, counts = from_pb(pb)
        A = R + D
        p, lo, hi = hat_p(D, A, phi)
        rec = (A - D / phi) / n_raw
        real_rows.append({
            "threshold": r, "N": n_raw, "assigned_any": A, "assigned_real": R, "assigned_decoy": D,
            "assignment_rate_real": R / n_raw, "assignment_rate_any": A / n_raw,
            "precision": p, "precision_lo": lo, "precision_hi": hi, "recall": rec,
            "detected_real": det, "detection": det / n_real, "decoy_fraction": D / A,
        })
        fq = QUIK / f"dual_checkerboard_decoy_thr{r}" / "26_022_S1_Dual_quik_R1_filtered.fastq"
        if fq.stat().st_size == 0:
            assigned_by_r[r] = None
            continue
        recs = {}
        with fq.open() as handle:
            for i, line in enumerate(handle):
                if i % 4:
                    continue
                m = HEADER_REAL.search(line)
                rid = line[1:].split()[0].partition("_calledidx_")[0]
                recs[rid] = m[2] in real_set
        assigned_by_r[r] = recs
        assert len(recs) == A

    for r in (6, 7, 8):
        recs = assigned_by_r[r]
        if recs is None:
            fl_rows.append({"threshold": r, "note": "filtered FASTQ empty"})
            continue
        full = {"real": 0, "decoy": 0, "unassigned": 0}
        short = {"real": 0, "decoy": 0, "unassigned": 0}
        n_full = n_short = 0
        for rid, kind in insert_of.items():
            bucket = full if kind == "full" else short
            if kind == "full":
                n_full += 1
            else:
                n_short += 1
            hit = recs.get(rid)
            if hit is None:
                bucket["unassigned"] += 1
            elif hit:
                bucket["real"] += 1
            else:
                bucket["decoy"] += 1
        any_short = short["real"] + short["decoy"]
        any_full = full["real"] + full["decoy"]
        fl_rows.append({
            "threshold": r, "full_length_reads": n_full, "short_reads": n_short,
            "assigned_real_full": full["real"], "assigned_decoy_full": full["decoy"],
            "assigned_real_short": short["real"], "assigned_decoy_short": short["decoy"],
            "assignment_rate_full": full["real"] / n_full,
            "assignment_rate_short": short["real"] / n_short,
            "decoy_fraction_full": full["decoy"] / any_full if any_full else None,
            "decoy_fraction_short": short["decoy"] / any_short if any_short else None,
        })

    real_txt = np.array([line.strip().upper() for line in open(TXT_PATH) if line.strip()])
    decoy_set = set(wl[n_real:])
    sim_rows = []
    for r in (6, 7, 8):
        pb = QUIK / f"simreads_thr{r}_1M_checkerboard" / "quik_dual_combined_1M_per_barcode.csv"
        R, D, det, _ = from_pb(pb)
        A = R + D
        p, lo, hi = hat_p(D, A, phi)
        fq = QUIK / f"simreads_thr{r}_1M_checkerboard" / "quik_dual_combined_1M_R1_filtered.fastq"
        true_idx, called = [], []
        with fq.open() as handle:
            for i, line in enumerate(handle):
                if i % 4:
                    continue
                m = HEADER_SIM.match(line.strip())
                if m is None:
                    raise ValueError(line[:160])
                true_idx.append(int(m[2]))
                called.append(m[4])
        called = np.array(called)
        true = real_txt[np.array(true_idx)]
        assert len(called) == A, (r, len(called), A)
        correct = called == true
        C = int(correct.sum())
        is_decoy = np.fromiter((s in decoy_set for s in called), dtype=bool, count=A)
        in_real = ~is_decoy
        D_fq = int(is_decoy.sum())
        F = int(in_real.sum()) - C
        phi_eff = D_fq / (D_fq + F) if D_fq + F else float("nan")
        true_p = C / A
        hat = 1 - (D_fq / A) / phi
        wrong_real = in_real & ~correct
        same_bar1 = int(sum(called[i][:18] == true[i][:18] for i in np.flatnonzero(wrong_real)))
        sim_rows.append({
            "threshold": r, "N": N_SIM, "assigned_any": A, "assigned_real": R, "assigned_decoy": D,
            "assignment_rate_real": R / N_SIM, "hat_P": p, "hat_P_lo": lo, "hat_P_hi": hi,
            "recall_cal": (A - D / phi) / N_SIM, "decoy_fraction": D / A, "detected_real": det,
            "C": C, "F": F, "D_fastq": D_fq, "true_precision": true_p, "hat_P_fastq": hat,
            "phi_eff": phi_eff, "bias_pp": (hat - true_p) * 100, "recall_true": C / N_SIM,
            "wrong_real": int(wrong_real.sum()), "wrong_real_same_bar1": same_bar1,
        })

    pd.DataFrame(real_rows).to_csv(OUT / "real_precision.csv", index=False)
    pd.DataFrame(sim_rows).to_csv(OUT / "sim_precision.csv", index=False)
    pd.DataFrame(fl_rows).to_csv(OUT / "full_length.csv", index=False)
    report = {
        "n_real": n_real, "n_decoy": n_decoy, "phi": phi, "raw_reads": n_raw,
        "error_model_sim": {"sub": 0.04, "ins": 0.05, "del": 0.14},
        "real": real_rows, "simulation": sim_rows, "full_length": fl_rows,
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, default=str))
    print("raw_reads", n_raw, "phi", round(phi, 6))
    print("REAL")
    print(pd.DataFrame(real_rows).to_string(index=False))
    print("SIM")
    cols = ["threshold", "assigned_any", "assignment_rate_real", "true_precision", "hat_P",
            "bias_pp", "phi_eff", "recall_true", "wrong_real_same_bar1", "wrong_real"]
    print(pd.DataFrame(sim_rows)[cols].to_string(index=False))
    print("FULL-LENGTH")
    print(pd.DataFrame(fl_rows).to_string(index=False))


if __name__ == "__main__":
    main()
