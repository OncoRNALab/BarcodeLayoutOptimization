# %% [markdown]
# # Low- and high-count features: sampling baseline and sequence differences
#
# Steps 0 and 2 of `results_and_analysis_plan.md`.
#
# * **Analysis area:** printed area (dead band excluded), at least 2 features from the edge, and outside dark regions (11 x 11 local mean below 25% of the median local mean; this captures the bubble-like spots in Greedy + 2-OPT).
# * **Counts:** UMI-deduplicated counts at r = 7 (exact UMI match), because PCR duplicates differ between arrays (22 cycles for Dual barcodes, 18 for the others). Read counts are used as a sensitivity check.
# * **Groups:** zero (0), low (1), medium (2-4) and high (> 4).
# * **Local expectation:** mean count of the 5 x 5 neighbourhood within the analysis area, excluding the feature itself.

# %%
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu, nbinom, poisson, spearmanr

from common import (ARRAYS, COLORS, OUT, THRESHOLDS, BC_LEN, bh, feature_counts, local_mean, setup_matplotlib, to_grid)

plt = setup_matplotlib()
COUNT = "umis"
GROUPS = ["zero", "low", "medium", "high"]
FEATURES = ["frac_A", "frac_C", "frac_G", "frac_T", "GC", "G_first18", "G_last18", "max_homopolymer", "n_GGG",
            "g4hunter", "n_cycles", "hamming_neighbours", "stray_events", "stray_cycles", "mfe_barcode", "mfe_probe",
            "nn_dist", "n_within_7", "kmer4_mean_log_freq", "kmer4_min_freq"]


def group_of(x):
    return pd.Series(np.select([x == 0, x == 1, x <= 4], GROUPS[:3], "high"), index=x.index)


data = {}
for name, cfg in ARRAYS.items():
    feats = pd.read_parquet(OUT / f"features_{cfg['key']}.parquet")
    lay, start = feature_counts(name, 7)
    df = feats.merge(lay[["row", "col", "reads", "reads_full_length", "umis", "printed"]], on=["row", "col"])
    for r in [5, 6]:
        other, _ = feature_counts(name, r)
        df = df.merge(other[["row", "col", "umis"]].rename(columns={"umis": f"umis_r{r}"}), on=["row", "col"])
    grid = to_grid(df, COUNT)
    printed = to_grid(df, "printed").astype(bool)
    smooth = local_mean(grid, printed, half=5)
    dark = smooth < 0.25 * np.nanmedian(smooth[printed])
    df["dark_region"] = dark[df.row, df.col]
    df["analysis"] = df.printed & (df.edge_dist >= 2) & ~df.dark_region
    area = to_grid(df, "analysis").astype(bool)
    df["local_mean"] = local_mean(grid, area, half=2)[df.row, df.col]
    df["local_log2_ratio"] = np.log2((df[COUNT] + 0.5) / (df.local_mean + 0.5))
    df["group"] = group_of(df[COUNT])
    data[name] = df

area_table = pd.DataFrame({name: {
    "features": len(df), "printed": int(df.printed.sum()), "dark regions (printed)": int((df.dark_region & df.printed).sum()),
    "analysis area": int(df.analysis.sum()), "mean UMIs per feature": round(df[COUNT][df.analysis].mean(), 3),
    "mean reads per feature": round(df.reads[df.analysis].mean(), 3)} for name, df in data.items()}).T
area_table

# %% [markdown]
# ## Step 0: how many zeros does sampling explain?
#
# * **Poisson, global mean:** every feature has the same expected count (mean of the analysis area).
# * **Poisson, local mean:** each feature's expected count is its neighbourhood mean; spatial variation is absorbed.
# * **Negative binomial, local mean:** adds over-dispersion with a single dispersion parameter, estimated by the method of moments from the residual variance around the local mean; this includes PCR jackpotting and any feature-to-feature variation, so it is an upper bound for "noise".

# %%
rows = []
for name, df in data.items():
    a = df[df.analysis]
    x, mu = a[COUNT].values, np.clip(a.local_mean.values, 1e-6, None)
    lam = x.mean()
    alpha = max(((x - mu) ** 2 - mu).sum() / (mu ** 2).sum(), 1e-6)
    n_nb = 1 / alpha
    rows.append({
        "layout": name, "mean count": lam,
        "zero observed (%)": 100 * (x == 0).mean(),
        "zero, Poisson global (%)": 100 * poisson.pmf(0, lam),
        "zero, Poisson local (%)": 100 * poisson.pmf(0, mu).mean(),
        "zero, NB local (%)": 100 * nbinom.pmf(0, n_nb, n_nb / (n_nb + mu)).mean(),
        "NB dispersion alpha": alpha,
        "one observed (%)": 100 * (x == 1).mean(),
        "one, Poisson local (%)": 100 * poisson.pmf(1, mu).mean(),
    })
sampling = pd.DataFrame(rows).set_index("layout")
sampling["excess zeros over local Poisson (pp)"] = sampling["zero observed (%)"] - sampling["zero, Poisson local (%)"]
sampling.round(2).to_csv(OUT / "sampling_baseline.csv")
sampling.round(2)

# %% [markdown]
# Group sizes in the analysis area, and stability of the groups across rejection thresholds (Spearman correlation of UMI counts between thresholds).

# %%
sizes = pd.DataFrame({name: df[df.analysis].group.value_counts().reindex(GROUPS) for name, df in data.items()}).T
stability = pd.DataFrame({name: {
    "rho r5 vs r7": spearmanr(df[df.analysis].umis_r5, df[df.analysis][COUNT])[0],
    "rho r6 vs r7": spearmanr(df[df.analysis].umis_r6, df[df.analysis][COUNT])[0],
    "rho reads vs UMIs (r7)": spearmanr(df[df.analysis].reads, df[df.analysis][COUNT])[0]} for name, df in data.items()}).T
pd.concat([sizes, stability.round(3)], axis=1)

# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
for ax, (name, df) in zip(axes, data.items()):
    a = df[df.analysis]
    x, mu = a[COUNT].values, np.clip(a.local_mean.values, 1e-6, None)
    k = np.arange(0, 13)
    obs = np.array([(x == i).mean() for i in k]) * 100
    pois = np.array([poisson.pmf(i, mu).mean() for i in k]) * 100
    pois_g = poisson.pmf(k, x.mean()) * 100
    ax.bar(k, obs, color=COLORS[name], alpha=0.6, label="observed")
    ax.plot(k, pois, "ko-", ms=3, label="Poisson, local mean")
    ax.plot(k, pois_g, "k:", label="Poisson, global mean")
    ax.set_title(name)
    ax.set_xlabel(f"UMIs per feature (r = 7)")
    ax.set_ylabel("Features (%)")
    ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "count_distribution_vs_poisson.pdf")
plt.show()

# %% [markdown]
# ## Step 2a: feature distributions per group
#
# For each feature: zero vs high and low vs high (Mann-Whitney U), with Cliff's delta as effect size (positive: larger values in the first group). Benjamini-Hochberg across all tests. In addition, the Spearman correlation of each feature with the local log2 ratio (count relative to the neighbourhood), which removes smooth spatial effects.
#
# Caveat: in Dual barcodes and Zipcodes, barcodes share 18-nt halves (block identifiers and column codes are shared by 16 and 235 features), so features are not independent and p-values are anti-conservative. Effect sizes are the quantity to interpret; the unit-level analysis is in `04_models`.

# %%
rows = []
for name, df in data.items():
    a = df[df.analysis]
    for f in FEATURES:
        for g1 in ["zero", "low"]:
            x, y = a[f][a.group == g1], a[f][a.group == "high"]
            if len(x) == 0 or len(y) == 0:
                continue
            u, p = mannwhitneyu(x, y)
            rows.append({"layout": name, "feature": f, "comparison": f"{g1} vs high", "n1": len(x), "n2": len(y),
                         "median_1": x.median(), "median_2": y.median(), "cliffs_delta": 2 * u / (len(x) * len(y)) - 1, "p": p})
        rho, p = spearmanr(a[f], a.local_log2_ratio)
        rows.append({"layout": name, "feature": f, "comparison": "rho with local log2 ratio", "n1": len(a),
                     "cliffs_delta": rho, "p": p})
tests = pd.DataFrame(rows)
tests["q"] = bh(tests.p)
tests.to_csv(OUT / "group_tests.csv", index=False)
effect = tests.pivot_table(index="feature", columns=["comparison", "layout"], values="cliffs_delta").reindex(FEATURES)
effect.round(3)

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
for ax, comp in zip(axes, ["zero vs high", "rho with local log2 ratio"]):
    m = effect[comp][list(ARRAYS)]
    lim = np.nanmax(np.abs(m.values))
    im = ax.imshow(m.values, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(ARRAYS)), list(ARRAYS), rotation=20)
    ax.set_yticks(range(len(FEATURES)), FEATURES)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title("Cliff's delta, zero vs high" if comp == "zero vs high" else "Spearman rho with local log2 ratio")
    fig.colorbar(im, ax=ax, shrink=0.6)
fig.tight_layout()
fig.savefig(OUT / "feature_effects.pdf")
plt.show()

# %% [markdown]
# Distributions of the features with the largest effects, per group.

# %%
top = (effect["zero vs high"].abs().max(axis=1).sort_values(ascending=False).index[:6]).tolist()
fig, axes = plt.subplots(len(top), 3, figsize=(12, 2.2 * len(top)), sharey="row")
for i, f in enumerate(top):
    for j, (name, df) in enumerate(data.items()):
        a = df[df.analysis]
        vals = [a[f][a.group == g].values for g in GROUPS]
        ax = axes[i, j]
        ax.violinplot([v for v in vals if len(v)], showmedians=True)
        ax.set_xticks(range(1, len(GROUPS) + 1), GROUPS)
        if i == 0:
            ax.set_title(name)
        if j == 0:
            ax.set_ylabel(f)
fig.tight_layout()
fig.savefig(OUT / "feature_violins.pdf")
plt.show()

# %% [markdown]
# ## Step 2b: k-mer enrichment, zero vs high
#
# Presence of each 2-, 3- and 4-mer per barcode; Fisher exact test of presence in zero- vs high-count barcodes; log2 odds ratio with 0.5 pseudocount; Benjamini-Hochberg within layout.

# %%
from itertools import product

KMERS = ["".join(p) for k in (2, 3, 4) for p in product("ACGT", repeat=k)]
rows = []
for name, df in data.items():
    a = df[df.analysis & df.group.isin(["zero", "high"])]
    zero = (a.group == "zero").values
    for k in KMERS:
        has = a.barcode.str.contains(k, regex=False).values
        t = np.array([[np.sum(has & zero), np.sum(~has & zero)], [np.sum(has & ~zero), np.sum(~has & ~zero)]])
        _, p = fisher_exact(t)
        lor = np.log2(((t[0, 0] + .5) * (t[1, 1] + .5)) / ((t[0, 1] + .5) * (t[1, 0] + .5)))
        rows.append({"layout": name, "kmer": k, "k": len(k), "log2_odds_ratio": lor, "p": p,
                     "frac_zero": t[0, 0] / zero.sum(), "frac_high": t[1, 0] / (~zero).sum()})
kmers = pd.DataFrame(rows)
kmers["q"] = kmers.groupby("layout").p.transform(bh)
kmers.to_csv(OUT / "kmer_enrichment.csv", index=False)
kmers.sort_values("p").groupby("layout").head(8).round(4)

# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, name in zip(axes, ARRAYS):
    d = kmers[kmers.layout == name]
    y = -np.log10(d.q.clip(lower=1e-300))
    ax.scatter(d.log2_odds_ratio, y, s=6, c=d.k, cmap="viridis", alpha=0.7)
    for _, row in d.sort_values("p").head(8).iterrows():
        ax.annotate(row.kmer, (row.log2_odds_ratio, -np.log10(max(row.q, 1e-300))), fontsize=6)
    ax.axhline(-np.log10(0.05), color="grey", ls=":", lw=0.8)
    ax.set_title(name)
    ax.set_xlabel("log2 odds ratio (zero vs high)")
    ax.set_ylabel("-log10 q")
fig.tight_layout()
fig.savefig(OUT / "kmer_enrichment_volcano.pdf")
plt.show()

# %% [markdown]
# ## Step 2c: positional base composition, zero vs high
#
# Difference in base frequency (zero minus high) at each of the 36 barcode positions, with a chi-square test per position (BH within layout). Positions are in synthesis order: position 1 is coupled first (5' to 3' synthesis).

# %%
fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
BASE_COLORS = {"A": "tab:green", "C": "tab:blue", "G": "orange", "T": "tab:red"}
rows = []
for ax, (name, df) in zip(axes, data.items()):
    a = df[df.analysis]
    mats = {}
    for g in ["zero", "high"]:
        seqs = np.array([list(s) for s in a.barcode[a.group == g]])
        mats[g] = pd.DataFrame({b: (seqs == b).sum(axis=0) for b in "ACGT"})
    diff = mats["zero"].div(mats["zero"].sum(axis=1), axis=0) - mats["high"].div(mats["high"].sum(axis=1), axis=0)
    pvals = [chi2_contingency(np.vstack([mats["zero"].iloc[i], mats["high"].iloc[i]]))[1] for i in range(BC_LEN)]
    q = bh(pvals)
    pos_bottom, neg_bottom = np.zeros(BC_LEN), np.zeros(BC_LEN)
    for b in "ACGT":
        v = diff[b].values * 100
        ax.bar(np.arange(1, BC_LEN + 1), np.where(v > 0, v, 0), bottom=pos_bottom, color=BASE_COLORS[b], label=b)
        ax.bar(np.arange(1, BC_LEN + 1), np.where(v < 0, v, 0), bottom=neg_bottom, color=BASE_COLORS[b])
        pos_bottom += np.where(v > 0, v, 0)
        neg_bottom += np.where(v < 0, v, 0)
    for i in np.where(q < 0.05)[0]:
        ax.text(i + 1, pos_bottom[i] + 0.2, "*", ha="center", fontsize=8)
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(18.5, color="grey", ls=":", lw=0.8)
    ax.set_ylabel("Zero minus high (pp)")
    ax.set_title(name, fontsize=10)
    for i in range(BC_LEN):
        rows.append({"layout": name, "position": i + 1, **{f"diff_{b}": diff[b].iloc[i] for b in "ACGT"}, "q": q[i]})
axes[0].legend(frameon=False, ncol=4, fontsize=8)
axes[-1].set_xlabel("Barcode position (synthesis order)")
fig.tight_layout()
fig.savefig(OUT / "positional_composition_zero_vs_high.pdf")
plt.show()
positional = pd.DataFrame(rows)
positional.to_csv(OUT / "positional_composition_zero_vs_high.csv", index=False)

# %%
for name, df in data.items():
    df.to_parquet(OUT / f"analysis_table_{ARRAYS[name]['key']}.parquet")
