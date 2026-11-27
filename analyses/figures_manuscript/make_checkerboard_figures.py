#!/usr/bin/env python3
"""Half-density Dual checkerboard: reads-per-spot heatmaps and assignment vs precision."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
FIG = ROOT / "figures"
ANA = ROOT / "analyses" / "dual_checkerboard"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
})

LAYOUT = DATA / "1in2" / "layout_1in2_dual_4x4.csv"
WL = DATA / "decoys" / "output" / "dual_1in2" / "dual_whitelist_real_plus_decoys.txt"
PB7 = DATA / "quik_out_manuscript" / "dual_checkerboard_decoy_thr7" / "26_022_S1_Dual_quik_per_barcode.csv"
FULL = ROOT / "analyses" / "recall_sequence_features" / "full_length_summary.csv"
HALF = ANA / "full_length.csv"


COLORS = {
    "Dual barcodes": "tab:blue",
    "Greedy + 2-OPT": "tab:orange",
    "Zipcodes": "tab:green",
    "Dual barcodes, half density": "#08306b",
}
# Zoom on printer coordinates (inclusive start, exclusive end), away from the array edge.
ZOOM_X = (48, 80)
ZOOM_Y = (80, 128)


def reads_at_r7() -> pd.DataFrame:
    layout = pd.read_csv(LAYOUT)
    layout["combined"] = layout.combined.str.upper()
    wl = [line.strip().upper() for line in WL.open() if line.strip()]
    counts = pd.read_csv(PB7).read_count.to_numpy()
    n_real = len(layout)
    real = pd.DataFrame({"combined": wl[:n_real], "reads": counts[:n_real]})
    df = layout.merge(real, on="combined", how="left")
    if df.reads.isna().any():
        raise SystemExit(f"{int(df.reads.isna().sum())} barcodes missing from the whitelist")
    df = df.rename(columns={"combined": "barcode"})
    df["reads"] = df.reads.astype(int)
    df.to_csv(ANA / "reads_per_feature_r7.csv", index=False)
    return df


def printer_grid(df: pd.DataFrame) -> tuple[np.ndarray, int, int]:
    x0, y0 = int(df.x.min()), int(df.y.min())
    g = np.zeros((int(df.x.max()) - x0 + 1, int(df.y.max()) - y0 + 1))
    g[df.x.to_numpy() - x0, df.y.to_numpy() - y0] = df.reads.to_numpy()
    return g, x0, y0


def spatial(df: pd.DataFrame) -> None:
    g, x0, y0 = printer_grid(df)
    vmax = np.percentile(df.loc[df.reads > 0, "reads"], 99)
    fig, axes = plt.subplots(
        1, 2, figsize=(12.6, 4.6), gridspec_kw={"width_ratios": [2.35, 1.0]},
    )
    im = axes[0].imshow(
        np.ma.masked_equal(g, 0), cmap="viridis", norm=LogNorm(vmin=1, vmax=vmax),
        interpolation="nearest", aspect="equal",
        extent=(y0 - 0.5, y0 + g.shape[1] - 0.5, x0 + g.shape[0] - 0.5, x0 - 0.5),
    )
    axes[0].add_patch(plt.Rectangle(
        (ZOOM_Y[0] - 0.5, ZOOM_X[0] - 0.5),
        ZOOM_Y[1] - ZOOM_Y[0], ZOOM_X[1] - ZOOM_X[0],
        fill=False, edgecolor="white", lw=0.9, ls="--",
    ))
    axes[0].set_facecolor("black")
    axes[0].set_title("(a) Dual barcodes")
    axes[0].set_xlabel("column")
    axes[0].set_ylabel("row")

    x1, x2 = ZOOM_X
    y1, y2 = ZOOM_Y
    z = g[x1 - x0:x2 - x0, y1 - y0:y2 - y0]
    axes[1].imshow(
        np.ma.masked_equal(z, 0), cmap="viridis", norm=LogNorm(vmin=1, vmax=vmax),
        interpolation="nearest", aspect="equal",
        extent=(y1 - 0.5, y2 - 0.5, x2 - 0.5, x1 - 0.5),
    )
    axes[1].set_facecolor("black")
    axes[1].set_title("(b) checkerboard detail")
    axes[1].set_xlabel("column")
    axes[1].set_ylabel("row")
    fig.colorbar(im, ax=axes, shrink=0.82, label="reads per feature (black = empty or 0)")
    fig.savefig(FIG / "checkerboard_spatial.pdf", bbox_inches="tight")
    plt.close(fig)


def assignment() -> None:
    full = pd.read_csv(FULL)
    half = pd.read_csv(HALF)
    d = full[full.layout == "Dual barcodes"]
    phi = 1 / 3
    hx = 100 * half.assignment_rate_full
    hy = 100 * (1 - half.decoy_fraction_full / phi)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2))
    axes[0].plot(
        d["assignment_rate_full (%)"], d["precision_full (%)"], "o-",
        color=COLORS["Dual barcodes"], label="full density",
    )
    for r, x, y in zip(d.threshold, d["assignment_rate_full (%)"], d["precision_full (%)"]):
        dx, dy = (5, -14) if int(r) == 7 else (6, 7)
        axes[0].annotate(f"$r={r}$", (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=10)
    axes[0].plot(
        hx, hy, "s-", color=COLORS["Dual barcodes, half density"], label="half density",
    )
    half_off = {6: (-22, 8), 7: (6, 8), 8: (6, -12)}
    for r, x, y in zip(half.threshold, hx, hy):
        dx, dy = half_off[int(r)]
        axes[0].annotate(f"$r={r}$", (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=10)
    axes[1].plot(d.threshold, d["assignment_rate_full (%)"], "o-",
                 color=COLORS["Dual barcodes"], label="full density")
    axes[1].plot(half.threshold, hx, "s-",
                 color=COLORS["Dual barcodes, half density"], label="half density")
    axes[0].set_title("(a)")
    axes[0].set_xlabel("assignment rate (% of full-length inserts)")
    axes[0].set_ylabel("precision (%)")
    axes[0].legend(frameon=False, loc="lower right")
    xmin, xmax = axes[0].get_xlim()
    ymin, ymax = axes[0].get_ylim()
    dx, dy = xmax - xmin, ymax - ymin
    axes[0].set_xlim(xmin - 0.12 * dx, xmax + 0.14 * dx)
    axes[0].set_ylim(ymin - 0.12 * dy, ymax + 0.16 * dy)
    axes[1].set_title("(b)")
    axes[1].set_xlabel("rejection threshold $r$")
    axes[1].set_ylabel("assignment rate (% of full-length inserts)")
    axes[1].set_xticks([5, 6, 7, 8])
    axes[1].legend(frameon=False, loc="center right")
    xmin, xmax = axes[1].get_xlim()
    ymin, ymax = axes[1].get_ylim()
    dx, dy = xmax - xmin, ymax - ymin
    axes[1].set_xlim(xmin - 0.08 * dx, xmax + 0.08 * dx)
    axes[1].set_ylim(ymin - 0.12 * dy, ymax + 0.12 * dy)
    fig.tight_layout()
    fig.savefig(FIG / "checkerboard_assignment.pdf")
    plt.close(fig)


def main() -> None:
    df = reads_at_r7()
    print(
        f"features {len(df)} zeros {int((df.reads == 0).sum())} "
        f"mean {df.reads.mean():.2f} detection {(df.reads > 0).mean():.4f}",
        flush=True,
    )
    spatial(df)
    assignment()
    print("wrote", FIG / "checkerboard_spatial.pdf", FIG / "checkerboard_assignment.pdf")


if __name__ == "__main__":
    main()
