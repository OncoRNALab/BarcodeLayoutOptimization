#!/usr/bin/env python3
"""Schematic of one greedy placement step (BarLay).

At the target feature the unused barcode with the smallest summed Hamming
distance of synthesis schedules to already placed 8-neighbours is chosen.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12,
    "axes.titlesize": 13,
})

TARGET = "#F58518"
NEIGH = "#4C78A8"
PLACED = "#D9DEE7"
EMPTY = "#FFFFFF"
MISMATCH = "#E45756"
INK = "#1A1A1A"
MUTED = "#4A4A4A"
ON = "#243B53"
OFF = "#EEF1F5"


def cell(ax, c, r, s, x0, y0, n, **kw):
    x = x0 + c * s
    y = y0 + (n - 1 - r) * s
    ax.add_patch(Rectangle((x, y), s, s, **kw))
    return x + s / 2, y + s / 2


def legend_swatch(ax, x, y, fc, label, ec="0.25", lw=0.8, hatch=None):
    ax.add_patch(Rectangle((x, y), 0.38, 0.38, facecolor=fc, edgecolor=ec,
                           lw=lw, hatch=hatch))
    ax.text(x + 0.48, y + 0.19, label, ha="left", va="center", fontsize=11, color=INK)


def greedy_panel(ax):
    ax.set_xlim(-0.15, 10.6)
    ax.set_ylim(-0.15, 8.55)
    ax.set_aspect("equal")
    ax.axis("off")

    n, s, x0, y0 = 7, 0.88, 0.05, 1.85
    ax.text(x0, y0 + n * s + 0.28, "(a) One greedy step", ha="left", va="bottom",
            fontsize=13, color=INK)
    tr, tc = 3, 3
    placed = {(r, c) for r in range(n) for c in range(n) if (r, c) < (tr, tc)}
    neigh = {(tr + dr, tc + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)
             if not (dr == 0 and dc == 0) and 0 <= tr + dr < n and 0 <= tc + dc < n}
    placed_n = sorted(neigh & placed)
    empty_n = neigh - placed

    # 8-neighbourhood frame
    ax.add_patch(Rectangle(
        (x0 + (tc - 1) * s - 0.07, y0 + (n - 1 - (tr + 1)) * s - 0.07),
        3 * s + 0.14, 3 * s + 0.14, fill=False, edgecolor="0.15", lw=1.6, zorder=3,
    ))

    for r in range(n):
        for c in range(n):
            if (r, c) == (tr, tc):
                fc, ec, lw = TARGET, "#7A3E00", 1.6
            elif (r, c) in placed_n:
                fc, ec, lw = NEIGH, "#1F4E79", 1.2
            elif (r, c) in empty_n:
                fc, ec, lw = EMPTY, NEIGH, 1.2
            elif (r, c) in placed:
                fc, ec, lw = PLACED, "0.55", 0.55
            else:
                fc, ec, lw = EMPTY, "0.75", 0.55
            cell(ax, c, r, s, x0, y0, n, facecolor=fc, edgecolor=ec, lw=lw, zorder=2)

    ax.text(x0 + tc * s + s / 2, y0 + (n - 1 - tr) * s + s / 2, "T",
            ha="center", va="center", fontsize=15, fontweight="bold", color="white", zorder=4)

    for i, (r, c) in enumerate(placed_n, start=1):
        ax.text(x0 + c * s + s / 2, y0 + (n - 1 - r) * s + s / 2, fr"$d_{i}$",
                ha="center", va="center", fontsize=13, fontweight="bold",
                color="white", zorder=4)

    for r, c in empty_n:
        x, y = x0 + c * s, y0 + (n - 1 - r) * s
        ax.plot([x + 0.20, x + s - 0.20], [y + 0.20, y + s - 0.20],
                color=NEIGH, lw=1.0, alpha=0.5, zorder=3)
        ax.plot([x + s - 0.20, x + 0.20], [y + 0.20, y + s - 0.20],
                color=NEIGH, lw=1.0, alpha=0.5, zorder=3)

    ax.text(
        x0 + n * s / 2, 1.28,
        r"choose unused barcode minimizing  $d_1+d_2+d_3+d_4$",
        ha="center", va="center", fontsize=12, color=INK,
    )
    ax.text(
        x0 + n * s / 2, 0.78,
        "empty neighbours (outlined, not yet placed) are ignored",
        ha="center", va="center", fontsize=11, color=MUTED,
    )

    items = [
        (TARGET, "target", "#7A3E00", 1.2, None),
        (NEIGH, "placed neighbour", "#1F4E79", 1.0, None),
        (EMPTY, "empty neighbour", NEIGH, 1.2, None),
        (PLACED, "already placed", "0.55", 0.8, None),
        (EMPTY, "unplaced", "0.75", 0.8, None),
    ]
    ly = 7.55
    for fc, lab, ec, lw, _ in items:
        legend_swatch(ax, 6.75, ly, fc, lab, ec=ec, lw=lw)
        ly -= 0.58


def bits_row(ax, x, y, bits, w=0.42, h=0.58, mismatches=None):
    for i, b in enumerate(bits):
        fc = ON if b else OFF
        hit = mismatches and i in mismatches
        ax.add_patch(Rectangle(
            (x + i * w, y), w, h, facecolor=fc,
            edgecolor=MISMATCH if hit else "0.35", lw=1.6 if hit else 0.55,
        ))
        ax.text(x + i * w + w / 2, y + h / 2, str(b), ha="center", va="center",
                fontsize=10, color="white" if b else MUTED, fontweight="bold")


def schedule_panel(ax):
    ax.set_xlim(-0.15, 11.4)
    ax.set_ylim(-0.15, 8.55)
    ax.axis("off")
    ax.text(0.1, 8.25, "(b) Distance is schedule Hamming distance",
            ha="left", va="top", fontsize=13, color=INK)

    cand = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1]
    neigh = [1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1]
    mism = {i for i, (a, b) in enumerate(zip(cand, neigh)) if a != b}

    ax.text(0.1, 7.55, "Each barcode is a binary synthesis schedule.",
            ha="left", va="top", fontsize=12, color=INK)
    ax.text(0.1, 7.05, "1 = illuminated in that cycle;  0 = remains protected.",
            ha="left", va="top", fontsize=12, color=MUTED)

    bits_x = 0.1
    ax.text(0.1, 6.55, "candidate at T", ha="left", va="bottom", fontsize=11, color=INK)
    bits_row(ax, bits_x, 5.85, cand, mismatches=mism)
    ax.text(0.1, 5.52, "placed neighbour", ha="left", va="bottom", fontsize=11, color=NEIGH)
    bits_row(ax, bits_x, 4.82, neigh, mismatches=mism)

    ax.text(
        0.1, 4.38,
        "red outline: cycles in which one feature is on and the other is off",
        ha="left", va="top", fontsize=11, color=MISMATCH,
    )
    ax.text(
        0.1, 3.88,
        fr"$d_i$ = Hamming distance = {len(mism)}   (here, {len(mism)} of {len(cand)} cycles)",
        ha="left", va="top", fontsize=12, color=INK,
    )

    ax.text(0.1, 3.05, "Placement cost at T", ha="left", va="top",
            fontsize=12, fontweight="bold", color=INK)
    ax.text(
        0.1, 2.52,
        r"cost $= d_1 + d_2 + \cdots + d_k$  over the $k$ already placed neighbours.",
        ha="left", va="top", fontsize=12, color=INK,
    )
    ax.text(
        0.1, 1.90,
        "Every unused barcode in the pool is scored; the one with the\n"
        "smallest cost is written at T. Filling then continues to the\n"
        "next empty feature. 2-OPT later swaps pairs if that lowers the\n"
        "sum of neighbour distances over the whole array.",
        ha="left", va="top", fontsize=12, color=MUTED, linespacing=1.38,
    )


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 6.2))
    greedy_panel(axes[0])
    schedule_panel(axes[1])
    fig.tight_layout(w_pad=1.4)
    out = FIG / "greedy_principle.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
