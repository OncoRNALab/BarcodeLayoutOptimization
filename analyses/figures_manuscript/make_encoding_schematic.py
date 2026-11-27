#!/usr/bin/env python3
"""Schematic of Dual and Zipcodes barcode encoding (sequence halves on the grid)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12})

DUAL_BLOCKS = ["#4C78A8", "#F58518", "#54A24B", "#E45756"]
BAR2 = "#D9D9D9"
COL_C, ROW_C = "#F58518", "#54A24B"


def draw_molecule(ax, x0, y0, w, h, c_left, c_right, lab_left=None, lab_right=None,
                  ink_right="0.15", end_5=True, end_3=True, nt=True):
    ax.add_patch(Rectangle((x0, y0), w / 2, h, facecolor=c_left, edgecolor="0.2", lw=0.8))
    ax.add_patch(Rectangle((x0 + w / 2, y0), w / 2, h, facecolor=c_right, edgecolor="0.2", lw=0.8))
    if lab_left:
        ax.text(x0 + w / 4, y0 + h / 2, lab_left, ha="center", va="center", fontsize=11,
                color="white", fontweight="bold")
    if lab_right:
        ax.text(x0 + 3 * w / 4, y0 + h / 2, lab_right, ha="center", va="center", fontsize=11,
                color=ink_right, fontweight="bold")
    if end_5:
        ax.text(x0 - 0.08, y0 + h / 2, "5′", ha="right", va="center", fontsize=11)
    if end_3:
        ax.text(x0 + w + 0.08, y0 + h / 2, "3′", ha="left", va="center", fontsize=11)
    if nt:
        ax.text(x0 + w / 4, y0 - 0.12, "18 nt", ha="center", va="top", fontsize=10, color="0.3")
        ax.text(x0 + 3 * w / 4, y0 - 0.12, "18 nt", ha="center", va="top", fontsize=10, color="0.3")


def dual_panel(ax):
    ax.set_title("(a) Dual barcodes", loc="left", fontsize=13, pad=8)
    ax.set_xlim(-0.35, 11.0)
    ax.set_ylim(-0.55, 10.4)
    ax.set_aspect("equal")
    ax.axis("off")

    n, bs = 8, 4
    x0, y0, s = 0.15, 2.15, 0.60
    for r in range(n):
        for c in range(n):
            b = (r // bs) * 2 + (c // bs)
            ax.add_patch(Rectangle((x0 + c * s, y0 + (n - 1 - r) * s), s, s,
                                   facecolor=DUAL_BLOCKS[b], edgecolor="white", lw=0.6))
    for br in range(2):
        for bc in range(2):
            ax.add_patch(Rectangle(
                (x0 + bc * bs * s, y0 + (n - (br + 1) * bs) * s),
                bs * s, bs * s, fill=False, edgecolor="0.15", lw=1.3,
            ))

    mw, mh = 2.35, 0.62
    gap = 0.55
    mx = [x0 + n * s + 0.55, x0 + n * s + 0.55 + mw + gap]
    grid_mid = y0 + 4 * s
    my = [grid_mid + 0.55, grid_mid - 0.55 - mh]
    coords = [(0, 0), (1, 0), (0, 1), (1, 1)]
    for i, ((cx, cy), color) in enumerate(zip(coords, DUAL_BLOCKS)):
        draw_molecule(
            ax, mx[cx], my[cy], mw, mh, color, BAR2,
            lab_left=r"$\mathrm{bar}_1$" if i == 0 else None,
            lab_right=r"$\mathrm{bar}_2$" if i == 0 else None,
            end_5=cx == 0, end_3=cx == 1, nt=False,
        )
    ax.text((mx[0] + mx[1] + mw) / 2, my[0] + mh + 0.35, "one barcode per block",
            ha="center", va="bottom", fontsize=11)
    ax.text(x0 + 4 * s, y0 - 0.45,
            r"colour = $\mathrm{bar}_1$ (shared); $\mathrm{bar}_2$ is random",
            ha="center", va="top", fontsize=11, color="0.25")


def zip_panel(ax):
    ax.set_title("(b) Zipcodes", loc="left", fontsize=13, pad=8)
    ax.set_xlim(-1.35, 8.8)
    ax.set_ylim(-1.35, 10.4)
    ax.set_aspect("equal")
    ax.axis("off")

    draw_molecule(ax, 1.15, 9.15, 6.2, 0.85, COL_C, ROW_C,
                  "column", "row", ink_right="white")
    ax.text(4.25, 10.22, "column code          row code", ha="center", va="bottom", fontsize=11)

    nr, nc, hc, hr = 5, 6, 1, 2
    x0, y0, s = 0.55, 0.55, 1.05
    for r in range(nr):
        for c in range(nc):
            x, y = x0 + c * s, y0 + (nr - 1 - r) * s
            ax.add_patch(Rectangle((x, y), s, s, facecolor="white", edgecolor="0.75", lw=0.7))
            if c == hc and r != hr:
                ax.add_patch(Rectangle((x, y), s, s, facecolor=COL_C, edgecolor="0.75", lw=0.7))
            elif r == hr and c != hc:
                ax.add_patch(Rectangle((x, y), s, s, facecolor=ROW_C, edgecolor="0.75", lw=0.7))
            elif c == hc and r == hr:
                ax.add_patch(Rectangle((x, y), s / 2, s, facecolor=COL_C, edgecolor="white", lw=0.4))
                ax.add_patch(Rectangle((x + s / 2, y), s / 2, s, facecolor=ROW_C, edgecolor="white", lw=0.4))
                ax.add_patch(Rectangle((x, y), s, s, fill=False, edgecolor="0.1", lw=1.6))
    ax.annotate(
        "column code",
        xy=(x0 + (hc + 0.5) * s, y0 + nr * s),
        xytext=(x0 + (hc + 0.5) * s, y0 + nr * s + 0.45),
        fontsize=11, ha="center", va="bottom",
        arrowprops=dict(arrowstyle="-|>", color="0.2", lw=0.8),
    )
    ax.annotate(
        "row code",
        xy=(x0 + nc * s, y0 + (nr - 1 - hr + 0.5) * s),
        xytext=(x0 + nc * s + 0.3, y0 + (nr - 1 - hr + 0.5) * s),
        fontsize=11, ha="left", va="center",
        arrowprops=dict(arrowstyle="-|>", color="0.2", lw=0.8),
    )
    xi, yi = x0 + hc * s + s / 2, y0 + (nr - 1 - hr) * s
    ax.annotate(
        "intersection:\nfull barcode",
        xy=(xi, yi), xytext=(x0 + nc * s + 0.3, yi - 1.55),
        fontsize=11, ha="left", va="center",
        arrowprops=dict(arrowstyle="-|>", color="0.2", lw=0.8),
    )


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.6))
    dual_panel(axes[0])
    zip_panel(axes[1])
    fig.tight_layout(w_pad=1.2)
    out = FIG / "encoding_schematic.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
