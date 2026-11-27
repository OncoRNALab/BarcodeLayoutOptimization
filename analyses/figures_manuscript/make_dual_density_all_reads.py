#!/usr/bin/env python3
"""Dual barcodes, full vs half density, counting all reads: assignment vs precision and recall."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
FULL = ROOT / "analyses" / "real_arrays" / "real_arrays_summary.csv"
HALF = ROOT / "analyses" / "dual_checkerboard" / "real_precision.csv"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
})
C_FULL, C_HALF = "tab:blue", "#08306b"


def main() -> None:
    full = pd.read_csv(FULL)
    full = full[full.layout == "Dual barcodes"]
    half = pd.read_csv(HALF)
    series = [
        ("Full density", "o-", C_FULL, full.threshold.to_numpy(), full["Assignment rate (%)"].to_numpy(),
         full["Precision, calibrated phi (%)"].to_numpy(), full["Recall, calibrated (%)"].to_numpy(),
         {5: (6, 6), 6: (6, 6), 7: (6, -14)}),
        ("Half density", "s-", C_HALF, half.threshold.to_numpy(), 100 * half.assignment_rate_real.to_numpy(),
         100 * half.precision.to_numpy(), 100 * half.recall.to_numpy(),
         {6: (-30, 8), 7: (6, 8), 8: (6, -14)}),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0))
    for label, fmt, color, r, x, y, rec, off in series:
        axes[0].plot(x, y, fmt, color=color, label=label)
        for ri, xi, yi in zip(r, x, y):
            axes[0].annotate(f"$r={int(ri)}$", (xi, yi), textcoords="offset points",
                             xytext=off[int(ri)], fontsize=10, color=color)
        axes[1].plot(r, rec, fmt, color=color, label=label)
        print(label, [f"r={int(a)}: {b:.1f}% at {c:.1f}%, recall {d:.1f}%" for a, b, c, d in zip(r, x, y, rec)])

    axes[0].set_title("(a) all reads")
    axes[0].set_xlabel("assignment rate (% of all reads)")
    axes[0].set_ylabel("calibrated precision (%)")
    axes[0].set_xlim(25, 85)
    axes[0].set_ylim(84, 101)
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].set_title("(b) all reads")
    axes[1].set_xlabel("rejection threshold $r$")
    axes[1].set_ylabel("calibrated recall (% of all reads)")
    axes[1].set_xticks([5, 6, 7, 8])
    axes[1].set_xlim(4.6, 8.4)
    axes[1].set_ylim(25, 80)
    axes[1].legend(frameon=False, loc="center right")
    fig.tight_layout()
    out = FIG / "dual_full_vs_half_all_reads_assignment_vs_precision"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
