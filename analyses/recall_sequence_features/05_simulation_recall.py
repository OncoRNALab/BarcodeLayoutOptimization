# %% [markdown]
# # Intrinsic decodability: per-barcode recall in simulation vs real counts
#
# Step 4, parts 1-2. In the 1M-read simulations every barcode received reads with the same error model, so differences in per-barcode recall reflect only how hard the barcode is for QUIK (close whitelist neighbours, shared 4-mers). If real counts follow simulated recall, low counts are a decoding problem; if not, they come from the molecules.
#
# The number of reads simulated per barcode is not stored in the local files (the local simulated FASTQs contain 100,000 reads), so recall is estimated as correct assignments divided by the expected 1,000,000 / 85,305 = 11.72 reads per barcode. This adds Poisson noise (about 29% per barcode) but no bias.

# %%
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import ARRAYS, COLORS, LAYOUT_DIR, OUT, QUIK_DIR, THRESHOLDS, setup_matplotlib

plt = setup_matplotlib()
COUNT = "umis"
N_SIM, N_BC = 1_000_000, 85_305
HEADER = re.compile(r"^@read_(\d+)_barcode_(\d+)_calledidx_(\d+)_([ACGTN]+)")
data = {name: pd.read_parquet(OUT / f"analysis_table_{cfg['key']}.parquet") for name, cfg in ARRAYS.items()}
FEATURES = ["frac_G", "max_homopolymer", "n_GGG", "mfe_probe", "n_cycles", "stray_events", "nn_dist", "n_within_7",
            "kmer4_mean_log_freq", "kmer4_min_freq"]


def sim_recall(name, r):
    cfg = ARRAYS[name]
    order = [l.strip().upper() for l in open(LAYOUT_DIR / cfg["whitelist"]) if l.strip()]
    path = QUIK_DIR / f"simreads_thr{r}_1M" / f"quik_{cfg['sim']}_combined_1M_R1_filtered.fastq"
    true, called = [], []
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                m = HEADER.match(line)
                true.append(int(m[2]))
                called.append(m[4])
    df = pd.DataFrame({"true_seq": np.array(order)[true], "called_seq": called})
    correct = df[df.true_seq == df.called_seq].true_seq.value_counts()
    assigned = df.true_seq.value_counts()
    out = pd.DataFrame({"barcode": order})
    out[f"sim_correct_r{r}"] = correct.reindex(order).fillna(0).values
    out[f"sim_recall_r{r}"] = out[f"sim_correct_r{r}"] / (N_SIM / N_BC)
    out[f"sim_assigned_r{r}"] = assigned.reindex(order).fillna(0).values
    return out


for name in ARRAYS:
    for r in THRESHOLDS:
        data[name] = data[name].merge(sim_recall(name, r), on="barcode", how="left")
    data[name].to_parquet(OUT / f"analysis_table_sim_{ARRAYS[name]['key']}.parquet")

pd.DataFrame({name: df[[f"sim_recall_r{r}" for r in THRESHOLDS]].mean() for name, df in data.items()}).T.round(3)

# %% [markdown]
# ## Does simulated recall predict real counts?
#
# Spearman correlations in the analysis area at r = 7. For reference, the split-half reliability of simulated recall itself is limited by Poisson noise; with about 4-5 correct reads per barcode at r = 7, the maximum attainable correlation is well below 1.

# %%
rows = []
for name, df in data.items():
    a = df[df.analysis]
    for target in [COUNT, "local_log2_ratio"]:
        rho, p = spearmanr(a.sim_recall_r7, a[target])
        rows.append({"layout": name, "real quantity": target, "spearman_rho": rho, "p": p, "n": len(a)})
    rho, p = spearmanr(a.sim_recall_r5, a.sim_recall_r7)
    rows.append({"layout": name, "real quantity": "(sim recall r5 vs r7)", "spearman_rho": rho, "p": p, "n": len(a)})
sim_vs_real = pd.DataFrame(rows)
sim_vs_real.to_csv(OUT / "sim_recall_vs_real.csv", index=False)
sim_vs_real.round(4)

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for name, df in data.items():
    a = df[df.analysis].copy()
    a["bin"] = pd.qcut(a.sim_recall_r7.rank(method="first"), 10, labels=False)
    g = a.groupby("bin").agg(sim=("sim_recall_r7", "mean"), real=(COUNT, "mean"), zero=(COUNT, lambda x: (x == 0).mean()))
    axes[0].plot(g.sim, g.real / a[COUNT].mean(), "o-", color=COLORS[name], label=name)
    axes[1].plot(g.sim, g.zero * 100, "o-", color=COLORS[name], label=name)
axes[0].axhline(1, color="grey", lw=0.8)
axes[0].set_xlabel("Simulated recall at r = 7 (decile mean)")
axes[0].set_ylabel("Real UMI count relative to layout mean")
axes[1].set_xlabel("Simulated recall at r = 7 (decile mean)")
axes[1].set_ylabel("Real features with zero UMIs (%)")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "sim_recall_vs_real.pdf")
plt.show()

# %% [markdown]
# ## Which features determine simulated recall?
#
# Negative control: in simulation the error model is the same for every barcode, so only decoding features (nearest-neighbour distance, crowding, 4-mer frequency) should matter. A synthesis feature such as `stray_events` can only correlate with simulated recall through its correlation with sequence.

# %%
rows = []
for name, df in data.items():
    a = df[df.analysis]
    for f in FEATURES:
        rows.append({"layout": name, "feature": f,
                     "rho with sim recall": spearmanr(a[f], a.sim_recall_r7)[0],
                     "rho with real local log2 ratio": spearmanr(a[f], a.local_log2_ratio)[0]})
feat_sim = pd.DataFrame(rows)
feat_sim.to_csv(OUT / "features_vs_sim_and_real.csv", index=False)
feat_sim.pivot(index="feature", columns="layout").round(3)
