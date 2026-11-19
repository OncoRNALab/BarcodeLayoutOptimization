#!/usr/bin/env python3
"""Generate labeled garbled reads (txt + paired FASTQ) from a layout CSV.

Channel: unequal walk (sub/ins/del). R2 = reverse_complement(R1).
Headers carry true_barcode for decode-precision scoring.

Usage:
    python make_reads.py --layout .../layout.csv --out-dir output/dual_shared \\
        --n-reads 100000 --seed 1
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path

from error_model import (
    DEFAULT_DRATE,
    DEFAULT_IRATE,
    DEFAULT_SRATE,
    reverse_complement,
    simulate_phred,
    simulate_read,
)

REPO = Path(__file__).resolve().parents[1]

PRESETS = {
    "dual_shared_pool_g10_minH6": REPO
    / "Dual_barcodes/layout_235x363_dual_shared_pool_g10_minH6.csv",
    "dual_fixed_bar2_8x8_g10_nn8": REPO
    / "Dual_barcodes/layout_235x363_dual_fixed_bar2_8x8_g10_nn8.csv",
    "zipcodes_flip_g10_nn8": REPO
    / "Zipcodes/output_235x363_flip_g10_h5_d12-16_far10/layout_235x363_zipcodes_nosep.csv",
    "hsdb_g10_nn8_colored": REPO
    / "HSDB/outputs/235x363_g10_h5_d12-16_far8_24+12_nn8_colored/layout_235x363_hsdb.csv",
}


def load_layout(path: Path):
    rows, cols, barcodes = [], [], []
    with open(path) as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        seq_key = "combined" if "combined" in fields else "barcode"
        for rec in reader:
            rows.append(int(rec["row"]))
            cols.append(int(rec["col"]))
            barcodes.append(rec[seq_key].strip().upper())
    return rows, cols, barcodes


def write_fastq_record(handle, header: str, seq: str, qual: str) -> None:
    if len(seq) != len(qual):
        raise ValueError(f"seq/qual length mismatch: {len(seq)} vs {len(qual)}")
    handle.write(f"@{header}\n{seq}\n+\n{qual}\n")


def run_qc(truth_rows: list[dict], barcode_len: int = 36) -> dict:
    n = len(truth_rows)
    mean_sub = sum(r["n_sub"] for r in truth_rows) / n
    mean_ins = sum(r["n_ins"] for r in truth_rows) / n
    mean_del = sum(r["n_del"] for r in truth_rows) / n
    mean_len = sum(r["len_R1"] for r in truth_rows) / n
    # spot-check R2 == RC(R1)
    for r in truth_rows[: min(100, n)]:
        if r["read_R2"] != reverse_complement(r["read_R1"]):
            raise RuntimeError(f"R2 != RC(R1) for read_id={r['read_id']}")
    return {
        "n_reads": n,
        "mean_n_sub": mean_sub,
        "mean_n_ins": mean_ins,
        "mean_n_del": mean_del,
        "mean_sub_per_ref_base": mean_sub / barcode_len,
        "mean_del_per_ref_base": mean_del / barcode_len,
        "mean_len_R1": mean_len,
        "r2_equals_rc_r1_checked": min(100, n),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", type=Path, default=None)
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default=None,
        help="Named NN≈8 layout (alternative to --layout)",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-reads", type=int, default=100_000)
    parser.add_argument("--srate", type=float, default=DEFAULT_SRATE)
    parser.add_argument("--irate", type=float, default=DEFAULT_IRATE)
    parser.add_argument("--drate", type=float, default=DEFAULT_DRATE)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mean-q", type=float, default=30.0)
    parser.add_argument("--sd-q", type=float, default=3.0)
    parser.add_argument("--qmin", type=int, default=10)
    parser.add_argument("--qmax", type=int, default=40)
    args = parser.parse_args()

    if (args.layout is None) == (args.preset is None):
        raise SystemExit("Provide exactly one of --layout or --preset")
    layout_path = args.layout if args.layout is not None else PRESETS[args.preset]
    if not layout_path.is_file():
        raise SystemExit(f"Layout not found: {layout_path}")

    rng = random.Random(args.seed)
    rows, cols, barcodes = load_layout(layout_path)
    n_spots = len(barcodes)
    print(f"Loaded {n_spots} spots from {layout_path}", flush=True)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    path_r1_txt = out / "reads_R1.txt"
    path_r2_txt = out / "reads_R2.txt"
    path_r1_fq = out / "reads_R1.fastq"
    path_r2_fq = out / "reads_R2.fastq"
    path_truth = out / "truth.csv"

    truth_rows: list[dict] = []
    with (
        open(path_r1_txt, "w") as f1t,
        open(path_r2_txt, "w") as f2t,
        open(path_r1_fq, "w") as f1q,
        open(path_r2_fq, "w") as f2q,
        open(path_truth, "w", newline="") as ft,
    ):
        writer = csv.DictWriter(
            ft,
            fieldnames=[
                "read_id",
                "row",
                "col",
                "true_barcode",
                "read_R1",
                "read_R2",
                "n_sub",
                "n_ins",
                "n_del",
                "len_R1",
            ],
        )
        writer.writeheader()

        for read_id in range(args.n_reads):
            idx = rng.randrange(n_spots)
            true = barcodes[idx]
            row, col = rows[idx], cols[idx]
            r1, n_sub, n_ins, n_del = simulate_read(
                true, args.srate, args.irate, args.drate, rng
            )
            r2 = reverse_complement(r1)
            q1 = simulate_phred(
                len(r1), rng, args.mean_q, args.sd_q, args.qmin, args.qmax
            )
            q2 = simulate_phred(
                len(r2), rng, args.mean_q, args.sd_q, args.qmin, args.qmax
            )
            header_base = (
                f"read{read_id} true_barcode={true} row={row} col={col}"
            )
            write_fastq_record(f1q, f"{header_base} mate=1", r1, q1)
            write_fastq_record(f2q, f"{header_base} mate=2", r2, q2)
            f1t.write(r1 + "\n")
            f2t.write(r2 + "\n")
            rec = {
                "read_id": read_id,
                "row": row,
                "col": col,
                "true_barcode": true,
                "read_R1": r1,
                "read_R2": r2,
                "n_sub": n_sub,
                "n_ins": n_ins,
                "n_del": n_del,
                "len_R1": len(r1),
            }
            writer.writerow(rec)
            truth_rows.append(rec)
            if (read_id + 1) % 20000 == 0:
                print(f"  {read_id + 1}/{args.n_reads}", flush=True)

    qc = run_qc(truth_rows, barcode_len=len(barcodes[0]))
    params = {
        "layout": str(layout_path),
        "preset": args.preset,
        "n_reads": args.n_reads,
        "seed": args.seed,
        "srate": args.srate,
        "irate": args.irate,
        "drate": args.drate,
        "r2_mode": "reverse_complement_of_R1",
        "phred": {
            "mean_q": args.mean_q,
            "sd_q": args.sd_q,
            "qmin": args.qmin,
            "qmax": args.qmax,
            "encoding": "phred+33",
        },
        "header_contract": (
            "@read{id} true_barcode={TRUE} row={r} col={c} mate={1|2}"
        ),
        "qc": qc,
    }
    with open(out / "params.json", "w") as handle:
        json.dump(params, handle, indent=2)
        handle.write("\n")

    print(json.dumps(qc, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
