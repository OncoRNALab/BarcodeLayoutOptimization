#!/usr/bin/env python3
"""Build manuscript Results figures from current (Zipcodes V2) local data."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
FIG = ROOT / "figures"
ANA = ROOT / "analyses"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
})
COLORS = {"Dual barcodes": "tab:blue", "Greedy + 2-OPT": "tab:orange", "Zipcodes": "tab:green"}
ORDER = ["Dual barcodes", "Greedy + 2-OPT", "Zipcodes"]
THRESHOLDS = [5, 6, 7]
ROWS, COLS, SEG = 235, 363, 18
DEAD = {"Dual barcodes": 308, "Greedy + 2-OPT": 317, "Zipcodes": 327}
PHI_DESIGN = 1 / 3
SAMPLE, SEED, WORKERS = 3000, 1, -1


def read_seqs(path: Path) -> list[str]:
    return [line.strip().upper() for line in path.open() if line.strip()]


def nn_lev(queries: list[str], targets: list[str], self_match: bool) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(queries), size=min(SAMPLE, len(queries)), replace=False)
    q = [queries[i] for i in idx]
    d = cdist(q, targets, scorer=Levenshtein.distance, workers=WORKERS)
    if self_match:
        d = np.where(d == 0, np.inf, d)
    return d.min(axis=1)


def grid_from(df: pd.DataFrame) -> np.ndarray:
    g = np.zeros((ROWS, COLS))
    g[df.row.to_numpy(), df.col.to_numpy()] = df.reads.to_numpy()
    return g


def leftmost_illuminated(seq: str, cycle: int, order: str = "ACGT") -> bool:
    c = -1
    n = len(order)
    for base in seq:
        c += 1
        while order[c % n] != base:
            c += 1
        if c == cycle:
            return True
        if c > cycle:
            return False
    return False


def mask_grid(df: pd.DataFrame, cycle: int) -> np.ndarray:
    g = np.zeros((ROWS, COLS), dtype=np.uint8)
    for row, col, bc in zip(df.row.to_numpy(), df.col.to_numpy(), df.barcode.to_numpy()):
        if leftmost_illuminated(bc, int(cycle)):
            g[int(row), int(col)] = 1
    return g


def masks():
    print("activation masks", flush=True)
    dual = pd.read_csv(DATA / "1in1" / "layout_235x363_dual_4x4.csv").rename(columns={"combined": "barcode"})
    zipc = pd.read_csv(DATA / "1in1" / "layout_235x363_zipcodes_nosep.csv")
    greedy = pd.read_csv(DATA / "1in1" / "layout_235x363_randombarcodes_opt_coords.csv",
                         header=None, names=["barcode", "col", "row"])
    for df in (dual, zipc, greedy):
        df["barcode"] = df.barcode.str.upper()
    random = greedy.copy()
    random["barcode"] = np.random.default_rng(1).permutation(random.barcode.to_numpy())
    cycle = 0
    panels = [
        ("Random assignment", mask_grid(random, cycle)),
        ("Greedy + 2-OPT", mask_grid(greedy, cycle)),
        ("Dual barcodes", mask_grid(dual, cycle)),
        ("Zipcodes", mask_grid(zipc, cycle)),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.4))
    for ax, (title, g) in zip(axes.ravel(), panels):
        ax.imshow(g, cmap="gray", interpolation="nearest", aspect="equal", vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.6)
    fig.tight_layout()
    fig.savefig(FIG / "masks.pdf")
    plt.close(fig)


def decoy_distances():
    print("decoy distances", flush=True)
    sets = {
        "Dual barcodes": [
            ("real barcodes", DATA / "1in1" / "layout_235x363_dual_4x4.txt", True,
             dict(color="black", lw=2)),
            ("layout-matched decoys", DATA / "decoys" / "output" / "dual" / "dual_decoys.txt", False,
             dict(color="tab:blue", lw=2)),
            ("independent decoys", DATA / "decoy_layouts" / "barcodes" / "1in1" / "layout_235x363_dual_4x4.txt", False,
             dict(color="tab:red", lw=1.5, ls="--")),
        ],
        "Zipcodes": [
            ("real barcodes", DATA / "1in1" / "layout_235x363_zipcodes_nosep.txt", True,
             dict(color="black", lw=2)),
            ("layout-matched decoys", DATA / "decoys" / "output" / "zipcodes_balanced_khalf" / "zip_decoysv2.txt", False,
             dict(color="tab:blue", lw=2)),
            ("independent decoys", DATA / "decoy_layouts" / "barcodes" / "1in1" / "layout_235x363_zipcodes_nosep.txt", False,
             dict(color="tab:red", lw=1.5, ls="--")),
        ],
        "Greedy + 2-OPT": [
            ("real barcodes", DATA / "1in1" / "layout_235x363_randombarcodes_opt.txt", True,
             dict(color="black", lw=2)),
            ("unplaced pool decoys", DATA / "decoys" / "output" / "greedy" / "greedy_decoys.txt", False,
             dict(color="tab:blue", lw=2)),
        ],
    }
    distances, rows = {}, []
    for layout, series in sets.items():
        real = read_seqs(series[0][1])
        for label, path, self, style in series:
            seqs = read_seqs(path)
            d = nn_lev(seqs, real, self_match=self)
            distances[layout, label] = (d, style)
            rows.append({"layout": layout, "distribution": label, "n": len(d),
                         "min": float(d.min()), "median": float(np.median(d)),
                         "mean": float(d.mean()), "max": float(d.max())})
            print(f"  {layout}: {label} mean {d.mean():.2f}", flush=True)
    pd.DataFrame(rows).to_csv(ANA / "decoy_quality" / "decoy_nearest_real_levenshtein.csv", index=False)

    fig, axes = plt.subplots(3, 1, figsize=(7.4, 9.6), sharex=True)
    bins = np.arange(0, 25) - 0.5
    for i, (ax, layout) in enumerate(zip(axes, ORDER)):
        for (lay, label), (d, style) in distances.items():
            if lay != layout:
                continue
            ax.hist(d, bins=bins, density=True, histtype="step",
                    label=f"{label} ({d.mean():.1f})", **style)
        ax.set_title(f"({chr(97 + i)}) {layout}", loc="left")
        ax.set_ylabel("fraction of barcodes")
        ax.legend(frameon=False, loc="upper right")
    axes[-1].set_xlabel("Levenshtein distance to nearest real")
    fig.tight_layout()
    fig.savefig(FIG / "decoy_distances.pdf")
    plt.close(fig)


def calibration():
    print("calibration", flush=True)
    dual_g = pd.read_csv(ANA / "precision_calibration" / "precision_calibration_1M.csv")
    zipv2 = pd.read_csv(ANA / "zipcodes_v2" / "sim_precision.csv")
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 7.4), sharex=True)
    kw = dict(markersize=7, lw=1.8)
    for j, name in enumerate(ORDER):
        ax = axes[0, j]
        ax2 = axes[1, j]
        if name == "Zipcodes":
            r = zipv2.threshold
            true = zipv2.true_precision * 100
            phi_g = 85006 / (85305 + 85006)
            hat_g = (1 - (zipv2.D / zipv2.assigned_any) / phi_g) * 100
            hat_d = zipv2.hat_P_phi_design * 100
            ax.plot(r, true, "o-", color="black", label="true $C/A$", **kw)
            ax.plot(r, hat_g, "s--", color="tab:red",
                    label=r"$\hat P$, $\phi_{\mathrm{global}}\approx 1/2$", **kw)
            ax.plot(r, hat_d, "D-", color="tab:blue", label=r"$\hat P$, $\phi=1/3$", **kw)
            ax2.plot(r, zipv2.phi_eff, "o-", color="black", label=r"$\phi_{\mathrm{eff}}$", **kw)
            ax2.axhline(0.499, color="grey", ls=":", label=r"$\phi_{\mathrm{global}}$")
            ax2.axhline(PHI_DESIGN, color="tab:blue", ls="--", label=r"$\phi=1/3$")
        else:
            df = dual_g[dual_g.layout == name]
            ax.plot(df.threshold, df.true_precision * 100, "o-", color="black",
                    label="true $C/A$", **kw)
            ax.plot(df.threshold, df.decoy_precision * 100, "s--", color="tab:blue",
                    label=r"$\hat P$, $\phi_{\mathrm{global}}$", **kw)
            ax2.plot(df.threshold, df.phi_eff, "o-", color="black",
                     label=r"$\phi_{\mathrm{eff}}$", **kw)
            ax2.axhline(df.phi_global.iloc[0], color="grey", ls=":", label=r"$\phi_{\mathrm{global}}$")
        ax.set_title(f"({chr(97 + j)}) {name}")
        ax.set_ylabel("precision (%)" if j == 0 else "")
        ax.set_ylim(80, 100.5)
        ax.legend(frameon=False, loc="center left" if name == "Zipcodes" else "lower left")
        ax2.set_title(f"({chr(100 + j)})")
        ax2.set_ylim(0, 0.62)
        ax2.set_xlabel("rejection threshold $r$")
        ax2.set_ylabel("decoy share of misassignments" if j == 0 else "")
        ax2.set_xticks(THRESHOLDS)
        ax2.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "precision_calibration.pdf")
    plt.close(fig)


def expand_limits(ax, x=0.12, y=0.12):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.set_xlim(xmin - x * (xmax - xmin), xmax + x * (xmax - xmin))
    ax.set_ylim(ymin - y * (ymax - ymin), ymax + y * (ymax - ymin))


def assignment():
    print("assignment vs precision (full-length)", flush=True)
    df = pd.read_csv(ANA / "recall_sequence_features" / "full_length_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2))
    for name in ORDER:
        d = df[df.layout == name]
        axes[0].plot(d["assignment_rate_full (%)"], d["precision_full (%)"], "o-",
                     color=COLORS[name], label=name)
        for r, x, y in zip(d.threshold, d["assignment_rate_full (%)"], d["precision_full (%)"]):
            dx, dy = (4, -14) if name == "Greedy + 2-OPT" else (6, 7)
            axes[0].annotate(f"$r={r}$", (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=10)
        axes[1].plot(d.threshold, d["assignment_rate_full (%)"], "o-", color=COLORS[name], label=name)
    axes[0].set_title("(a)")
    axes[0].set_xlabel("assignment rate (% of full-length inserts)")
    axes[0].set_ylabel("precision (%)")
    axes[0].legend(frameon=False, loc="lower left")
    expand_limits(axes[0], x=0.14, y=0.10)
    axes[0].set_ylim(axes[0].get_ylim()[0], 102)
    axes[1].set_title("(b)")
    axes[1].set_xlabel("rejection threshold $r$")
    axes[1].set_ylabel("assignment rate (% of full-length inserts)")
    axes[1].set_xticks(THRESHOLDS)
    axes[1].legend(frameon=False)
    expand_limits(axes[1], x=0.08, y=0.12)
    fig.tight_layout()
    fig.savefig(FIG / "assignment_vs_precision.pdf")
    plt.close(fig)


def spatial():
    print("heatmaps and Zipcodes G", flush=True)
    counts = {
        "Dual barcodes": pd.read_csv(ANA / "real_arrays" / "reads_per_feature_dual_thr7.csv"),
        "Greedy + 2-OPT": pd.read_csv(ANA / "real_arrays" / "reads_per_feature_greedy_thr7.csv"),
        "Zipcodes": pd.read_csv(ANA / "zipcodes_v2" / "reads_per_feature_r7.csv"),
    }
    grids = {name: grid_from(df) for name, df in counts.items()}
    vmax = max(np.percentile(g[g > 0], 99) for g in grids.values())

    fig = plt.figure(figsize=(13.8, 7.8))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.25, 1.0], hspace=0.42, wspace=0.35)
    axes_hm = [fig.add_subplot(gs[0, 2 * i:2 * i + 2]) for i in range(3)]
    for i, (ax, name) in enumerate(zip(axes_hm, ORDER)):
        g = grids[name]
        im = ax.imshow(np.ma.masked_equal(g, 0), cmap="viridis", norm=LogNorm(vmin=1, vmax=vmax),
                       interpolation="nearest", aspect="equal")
        ax.axvline(DEAD[name] - 0.5, color="white", lw=0.8, ls="--")
        ax.set_facecolor("black")
        ax.set_title(f"({chr(97 + i)}) {name}")
        ax.set_xlabel("column")
        ax.set_ylabel("row")
        ax.set_anchor("C")
    cbar = fig.colorbar(im, ax=axes_hm, orientation="horizontal", fraction=0.055, pad=0.22,
                        label="reads per feature (black = 0)")
    cbar.ax.tick_params(labelsize=11)

    zip_df = counts["Zipcodes"]
    printed = zip_df.col < DEAD["Zipcodes"]
    col_g = zip_df.groupby("col").barcode.first().str[:SEG].map(lambda s: 100 * s.count("G") / SEG)
    row_g = zip_df.groupby("row").barcode.first().str[SEG:].map(lambda s: 100 * s.count("G") / SEG)
    col_m = zip_df.loc[printed].groupby("col").reads.mean()
    row_m = zip_df.loc[printed].groupby("row").reads.mean()
    axc = fig.add_subplot(gs[1, 1:3])
    axr = fig.add_subplot(gs[1, 3:5])
    axc.scatter(col_g.loc[col_m.index], col_m, s=8, color=COLORS["Zipcodes"], alpha=0.7)
    axr.scatter(row_g.loc[row_m.index], row_m, s=8, color=COLORS["Zipcodes"], alpha=0.7)
    rc, _ = spearmanr(col_g.loc[col_m.index], col_m)
    rr, _ = spearmanr(row_g.loc[row_m.index], row_m)
    axc.set_xlabel("G in column code (%)")
    axc.set_ylabel("mean reads per column")
    axc.set_title(fr"(d) Zipcodes columns, $\rho={rc:.2f}$")
    axr.set_xlabel("G in row code (%)")
    axr.set_ylabel("mean reads per row")
    axr.set_title(fr"(e) Zipcodes rows, $\rho={rr:.2f}$")
    fig.savefig(FIG / "spatial_reads.pdf", bbox_inches="tight")
    plt.close(fig)
    return zip_df


def g_content(zip_df: pd.DataFrame):
    print("G content", flush=True)
    feats = pd.read_csv(ANA / "recall_sequence_features" / "barcode_features.csv")
    dual_g = feats[feats.layout.isin(["Dual barcodes", "Greedy + 2-OPT"])].copy()
    dual_g = dual_g[dual_g.analysis_set]
    dual_g["reads"] = dual_g.reads_r7
    dual_g["G_frac"] = dual_g.G_frac

    zip_df = zip_df.copy()
    zip_df["edge"] = np.minimum.reduce([
        zip_df.row, ROWS - 1 - zip_df.row, zip_df.col, COLS - 1 - zip_df.col,
    ])
    z = zip_df[(zip_df.col < DEAD["Zipcodes"]) & (zip_df.edge >= 2)].copy()
    z["layout"] = "Zipcodes"
    z["G_frac"] = z.barcode.str.count("G") / 36

    frames = [dual_g[["layout", "barcode", "reads", "G_frac"]], z[["layout", "barcode", "reads", "G_frac"]]]
    allg = pd.concat(frames, ignore_index=True)
    allg["group"] = np.select([allg.reads == 0, allg.reads == 1, allg.reads <= 4],
                              ["zero", "low", "medium"], default="high")

    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.8))
    ax = axes[0]
    positions, data, facecolors = [], [], []
    pos = 1
    xticks, xlabels = [], []
    for name in ORDER:
        for g, color in [("zero", "#bdbdbd"), ("high", "0.35")]:
            data.append(100 * allg.loc[(allg.layout == name) & (allg.group == g), "G_frac"])
            positions.append(pos)
            facecolors.append(color)
            pos += 1
        xticks.append(pos - 1.5)
        xlabels.append(name.replace(" + 2-OPT", "\n+ 2-OPT"))
        pos += 0.6
    bp = ax.boxplot(data, positions=positions, widths=0.7, showfliers=False, patch_artist=True)
    for patch, color in zip(bp["boxes"], facecolors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel("G content (%)")
    ax.set_title("(a) zero vs high-count barcodes")
    ax.plot([], [], color="#bdbdbd", lw=8, label="zero")
    ax.plot([], [], color="0.35", lw=8, label="high")
    ax.legend(frameon=False)

    ax = axes[1]
    x = np.arange(3)
    zeros, highs = [], []
    for name in ORDER:
        sub = allg[allg.layout == name]
        zeros.append(100 * sub.loc[sub.group == "zero", "barcode"].str.contains("GGG").mean())
        highs.append(100 * sub.loc[sub.group == "high", "barcode"].str.contains("GGG").mean())
        print(f"  {name} GGG zero {zeros[-1]:.1f}% high {highs[-1]:.1f}%", flush=True)
    ax.bar(x - 0.18, zeros, 0.36, color="#bdbdbd", label="zero")
    ax.bar(x + 0.18, highs, 0.36, color="0.35", label="high")
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace(" + 2-OPT", "\n+ 2-OPT") for n in ORDER])
    ax.set_ylabel("barcodes containing GGG (%)")
    ax.set_title("(b) GGG 3-mer")
    ax.legend(frameon=False)

    ax = axes[2]
    for name in ORDER:
        sub = allg[allg.layout == name]
        rhos = [spearmanr(sub.barcode.str[pos].eq("G").astype(int), sub.reads)[0] for pos in range(36)]
        ax.plot(np.arange(1, 37), rhos, color=COLORS[name], lw=1.3, label=name)
        print(f"  {name} mean rho {np.mean(rhos):.3f} pos18 {rhos[17]:.3f}", flush=True)
    ax.axhline(0, color="grey", lw=0.6)
    ax.axvline(18.5, color="grey", ls=":", lw=0.8)
    ax.set_xlim(1, 36)
    ax.set_xlabel("position (5′ → 3′)")
    ax.set_ylabel(r"Spearman $\rho$ of G vs reads")
    ax.set_title("(c) G along the barcode")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "g_content.pdf")
    plt.close(fig)


def main():
    masks()
    decoy_distances()
    calibration()
    assignment()
    zip_df = spatial()
    g_content(zip_df)
    print("wrote", sorted(p.name for p in FIG.glob("*.pdf")))


if __name__ == "__main__":
    main()
