#!/usr/bin/env python3
"""G context (isolated vs runs) and positional G / 3-mer analysis.

Answers:
  4.   Does 3-mer / G position along the barcode (5' = synthesis start) correlate with yield?
  4.1  Are the extra G's in zero-count barcodes clustered as runs (kmers) or scattered?
  4.2  Per-position base content of zero/low vs high-count barcodes.

Requires barcode_features.csv. Writes CSVs and editable PDFs.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from analyse import add_groups
from io_utils import ARRAYS, BARCODE_LEN, OUT_DIR

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
G_KMERS = ["GGG", "GGC", "GCG", "CGG", "TGG", "GGT", "GGA", "AGG"]
ORDER = ["zero", "low", "medium", "high"]
COLORS = {"zero": "#4c4c4c", "low": "#6baed6", "medium": "#74c476", "high": "#238b45"}


def g_context(seq: str) -> dict:
    isolated = in_gg = in_run3 = 0
    n_gg_runs = n_ggg_runs = 0
    i, n = 0, len(seq)
    while i < n:
        if seq[i] != "G":
            i += 1
            continue
        j = i
        while j < n and seq[j] == "G":
            j += 1
        length = j - i
        if length == 1:
            isolated += 1
        elif length == 2:
            in_gg += 2
            n_gg_runs += 1
        else:
            in_run3 += length
            n_ggg_runs += 1
        i = j
    n_g = isolated + in_gg + in_run3
    return {
        "n_G": n_g,
        "n_G_isolated": isolated,
        "n_G_in_GG": in_gg,
        "n_G_in_run3plus": in_run3,
        "n_GG_runs": n_gg_runs,
        "n_GGG_runs": n_ggg_runs,
        "frac_G_isolated": isolated / n_g if n_g else np.nan,
        "frac_G_clustered": (in_gg + in_run3) / n_g if n_g else np.nan,
        "has_GGG": int(n_ggg_runs > 0),
        "has_GG": int(n_gg_runs + n_ggg_runs > 0),
    }


def positional_g(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        enc = np.array([list(s) for s in a.barcode])
        is_g = enc == "G"
        reads = a.reads_r7.to_numpy()
        groups = a.count_group.to_numpy()
        for pos in range(BARCODE_LEN):
            x = is_g[:, pos]
            rho, p = stats.spearmanr(x, reads)
            rows.append({
                "layout": name, "feature": "G", "position": pos + 1,
                "spearman": rho, "p_raw": p,
                "frac_zero": x[groups == "zero"].mean() if (groups == "zero").any() else np.nan,
                "frac_low": x[groups == "low"].mean() if (groups == "low").any() else np.nan,
                "frac_medium": x[groups == "medium"].mean() if (groups == "medium").any() else np.nan,
                "frac_high": x[groups == "high"].mean() if (groups == "high").any() else np.nan,
                "mean_reads_present": reads[x].mean() if x.any() else np.nan,
                "mean_reads_absent": reads[~x].mean() if (~x).any() else np.nan,
            })
        for mer in G_KMERS + ["GG"]:
            k = len(mer)
            for pos in range(BARCODE_LEN - k + 1):
                x = np.fromiter((s[pos : pos + k] == mer for s in a.barcode), dtype=bool, count=len(a))
                if x.sum() < 20:
                    continue
                rho, p = stats.spearmanr(x, reads)
                rows.append({
                    "layout": name, "feature": mer, "position": pos + 1,
                    "spearman": rho, "p_raw": p,
                    "frac_zero": x[groups == "zero"].mean() if (groups == "zero").any() else np.nan,
                    "frac_high": x[groups == "high"].mean() if (groups == "high").any() else np.nan,
                    "n_present": int(x.sum()),
                    "mean_reads_present": reads[x].mean(),
                    "mean_reads_absent": reads[~x].mean(),
                })
    out = pd.DataFrame(rows)
    out["p_adj"] = multipletests(out.p_raw.fillna(1), method="fdr_bh")[1]
    return out


def matched_g_table(df: pd.DataFrame) -> pd.DataFrame:
    """At fixed G count, do clustered G's (GGG) predict lower yield?"""
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        for n_g, sub in a.groupby("n_G"):
            if len(sub) < 80 or n_g < 2:
                continue
            with_run = sub[sub.has_GGG == 1]
            without = sub[sub.has_GGG == 0]
            if len(with_run) < 20 or len(without) < 20:
                continue
            u, p = stats.mannwhitneyu(with_run.reads_r7, without.reads_r7, alternative="two-sided")
            rows.append({
                "layout": name, "n_G": int(n_g),
                "n_with_GGG": len(with_run), "n_without_GGG": len(without),
                "mean_reads_with_GGG": with_run.reads_r7.mean(),
                "mean_reads_without_GGG": without.reads_r7.mean(),
                "p_raw": p,
            })
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj"] = multipletests(out.p_raw, method="fdr_bh")[1]
    return out


def context_by_group(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        for grp, sub in a.groupby("count_group"):
            rows.append({
                "layout": name, "group": grp, "n": len(sub),
                "mean_n_G": sub.n_G.mean(),
                "mean_isolated": sub.n_G_isolated.mean(),
                "mean_in_GG": sub.n_G_in_GG.mean(),
                "mean_in_run3plus": sub.n_G_in_run3plus.mean(),
                "frac_G_isolated": sub.frac_G_isolated.mean(),
                "frac_has_GGG": sub.has_GGG.mean(),
            })
    return pd.DataFrame(rows)


def plot_g_by_position(pos: pd.DataFrame, path):
    g = pos[pos.feature == "G"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for ax, name in zip(axes, ARRAYS):
        d = g[g.layout == name]
        for grp, col in [("frac_zero", COLORS["zero"]), ("frac_low", COLORS["low"]),
                          ("frac_medium", COLORS["medium"]), ("frac_high", COLORS["high"])]:
            label = grp.replace("frac_", "")
            ax.plot(d.position, d[grp] * 100, color=col, lw=1.2, label=label)
        ax.axvline(18.5, color="grey", ls=":", lw=0.6)
        ax.set_title(name)
        ax.set_xlabel("Position (5′ → 3′)")
        ax.set_xlim(1, 36)
    axes[0].set_ylabel("% of barcodes with G")
    axes[0].legend(frameon=False, fontsize=8, title="count group")
    fig.suptitle("G content along the barcode, by read-count group (r = 7)", y=1.02)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_bases_zero_vs_high(comp: pd.DataFrame, path):
    fig, axes = plt.subplots(3, 2, figsize=(10, 7.2), sharex=True, sharey=True)
    for i, name in enumerate(ARRAYS):
        for j, grp in enumerate(("zero", "high")):
            ax = axes[i, j]
            d = comp[(comp.layout == name) & (comp.group == grp)]
            wide = d.pivot(index="position", columns="base", values="frequency")
            bottom = np.zeros(len(wide))
            for b, c in zip("ACGT", ["#74c476", "#6baed6", "#e6550d", "#d62728"]):
                ax.bar(wide.index, wide[b] * 100, bottom=bottom, color=c, width=1.0, linewidth=0, label=b)
                bottom = bottom + wide[b].to_numpy() * 100
            ax.axvline(18.5, color="white", ls=":", lw=0.6)
            if i == 0:
                ax.set_title(f"{grp}-count barcodes")
            if j == 0:
                ax.set_ylabel(f"{name}\n% of barcodes")
    axes[0, 1].legend(frameon=False, ncol=4, fontsize=8, loc="upper right")
    axes[-1, 0].set_xlabel("Position (5′ → 3′)")
    axes[-1, 1].set_xlabel("Position (5′ → 3′)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_kmer_position(pos: pd.DataFrame, path):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for ax, name in zip(axes, ARRAYS):
        d = pos[(pos.layout == name) & (pos.feature == "GGG")]
        ax.plot(d.position, d.frac_zero * 100, color=COLORS["zero"], lw=1.3, label="zero")
        ax.plot(d.position, d.frac_high * 100, color=COLORS["high"], lw=1.3, label="high")
        ax.axvline(18.5, color="grey", ls=":", lw=0.6)
        ax.set_title(name)
        ax.set_xlabel("GGG start position")
        ax.set_xlim(1, 34)
    axes[0].set_ylabel("% of barcodes containing GGG at this start")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("GGG position: zero- vs high-count barcodes", y=1.02)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_position_spearman(pos: pd.DataFrame, path):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    for ax, name in zip(axes, ARRAYS):
        g = pos[(pos.layout == name) & (pos.feature == "G")]
        ggg = pos[(pos.layout == name) & (pos.feature == "GGG")]
        ax.plot(g.position, g.spearman, color="#e6550d", lw=1.2, label="G at this position")
        ax.plot(ggg.position, ggg.spearman, color="black", lw=1.2, label="GGG starting here")
        ax.axhline(0, color="grey", lw=0.5)
        ax.axvline(18.5, color="grey", ls=":", lw=0.6)
        ax.set_title(name)
        ax.set_xlabel("Position (5′ → 3′)")
        ax.set_xlim(1, 36)
    axes[0].set_ylabel("Spearman ρ vs reads (r = 7)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_g_context(ctx: pd.DataFrame, path):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, name in zip(axes, ARRAYS):
        d = ctx[ctx.layout == name].set_index("group").reindex(ORDER)
        x = np.arange(len(ORDER))
        ax.bar(x, d.mean_isolated, color="#74c476", label="isolated G")
        ax.bar(x, d.mean_in_GG, bottom=d.mean_isolated, color="#fdae6b", label="G in GG")
        ax.bar(x, d.mean_in_run3plus, bottom=d.mean_isolated + d.mean_in_GG, color="#e6550d", label="G in GGG+")
        ax.set_xticks(x)
        ax.set_xticklabels(["0", "1", "2–4", ">4"])
        ax.set_title(name)
        ax.set_xlabel("Reads at r = 7")
    axes[0].set_ylabel("Mean G nucleotides per barcode")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_matched_g(matched: pd.DataFrame, path):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, name in zip(axes, ARRAYS):
        d = matched[matched.layout == name]
        if d.empty:
            ax.set_title(name)
            continue
        ax.plot(d.n_G, d.mean_reads_without_GGG, "o-", color="#238b45", label="no GGG")
        ax.plot(d.n_G, d.mean_reads_with_GGG, "s--", color="#e6550d", label="has GGG")
        ax.set_title(name)
        ax.set_xlabel("Number of G nucleotides")
    axes[0].set_ylabel("Mean reads (r = 7)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Yield at matched G content: clustered (GGG) vs scattered G", y=1.02)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    df = add_groups(pd.read_csv(OUT_DIR / "barcode_features.csv"))
    ctx = pd.DataFrame(df.barcode.map(g_context).tolist())
    df = pd.concat([df.reset_index(drop=True), ctx], axis=1)
    pos = positional_g(df)
    matched = matched_g_table(df)
    grp = context_by_group(df)
    comp = pd.read_csv(OUT_DIR / "position_composition.csv")

    pos.to_csv(OUT_DIR / "g_position.csv", index=False)
    matched.to_csv(OUT_DIR / "g_matched_yield.csv", index=False)
    grp.to_csv(OUT_DIR / "g_context_by_group.csv", index=False)

    plot_g_by_position(pos, OUT_DIR / "g_content_by_position.pdf")
    plot_bases_zero_vs_high(comp, OUT_DIR / "base_content_zero_vs_high.pdf")
    plot_kmer_position(pos, OUT_DIR / "ggg_by_position.pdf")
    plot_position_spearman(pos, OUT_DIR / "g_position_vs_reads.pdf")
    plot_g_context(grp, OUT_DIR / "g_isolated_vs_runs.pdf")
    plot_matched_g(matched, OUT_DIR / "g_matched_clustered_vs_scattered.pdf")

    print("G context by group")
    print(grp.round(3).to_string(index=False))
    print("\nMatched G (has GGG vs not), FDR < 0.05")
    sig = matched[matched.p_adj < 0.05] if len(matched) else matched
    print(sig.round(3).to_string(index=False) if len(sig) else "none")
    print("\nSpearman of G-at-position vs reads: mean over positions")
    g = pos[pos.feature == "G"]
    print(g.groupby("layout").spearman.agg(["mean", "min", "max"]).round(3).to_string())
    print("\nSpearman of GGG-at-position vs reads: mean over positions")
    ggg = pos[pos.feature == "GGG"]
    print(ggg.groupby("layout").spearman.agg(["mean", "min", "max"]).round(3).to_string())


if __name__ == "__main__":
    main()
