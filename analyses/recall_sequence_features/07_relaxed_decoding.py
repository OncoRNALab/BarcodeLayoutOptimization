# %% [markdown]
# # Relaxed decoding of unassigned full-length reads
#
# Step 4, parts 3-5. Reads with a full-length insert that QUIK did not assign at r = 7 are matched exhaustively (Levenshtein distance of Read 1 positions 11-46 to every entry of the combined whitelist). If zero-count features gain reads at relaxed distances at a rate similar to other features, their molecules exist but carry too many errors; if they gain few, the molecules are absent.
#
# The precision of relaxed assignments is estimated from the decoys in the whitelist ($1 - d/\phi_{\text{eff}}$, with $\phi_{\text{eff}}$ from the r = 7 simulation; an approximation, because the calibration was measured at r ≤ 7).
#
# Results are cached in `output/relaxed_<layout>.parquet` (about 20 minutes on 8 cores for the three arrays).

# %%
import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist

from common import ARRAYS, COLORS, OUT, SEG, UMI_LEN, BC_LEN, load_assignments, load_raw, load_whitelist, setup_matplotlib

plt = setup_matplotlib()
COUNT = "umis"
CUTOFF = 16
data = {name: pd.read_parquet(OUT / f"analysis_table_{cfg['key']}.parquet") for name, cfg in ARRAYS.items()}
calib = pd.read_csv(OUT.parents[1] / "precision_calibration" / "precision_calibration_1M.csv", index_col=["layout", "threshold"])


def best_matches(queries, wl):
    best, second, idx = (np.empty(len(queries), dtype=int) for _ in range(3))
    for s in range(0, len(queries), 2000):
        d = cdist(queries[s:s + 2000], wl, scorer=Levenshtein.distance, dtype=np.uint8, workers=-1, score_cutoff=CUTOFF)
        part = np.argpartition(d, 1, axis=1)[:, :2]
        two = np.take_along_axis(d, part, axis=1)
        o = np.argsort(two, axis=1)
        best[s:s + 2000] = np.take_along_axis(two, o, axis=1)[:, 0]
        second[s:s + 2000] = np.take_along_axis(two, o, axis=1)[:, 1]
        idx[s:s + 2000] = np.take_along_axis(part, o, axis=1)[:, 0]
    return best, second, idx


relaxed = {}
for name, cfg in ARRAYS.items():
    path = OUT / f"relaxed_{cfg['key']}.parquet"
    if not path.exists():
        raw = load_raw(name)
        assigned = set(load_assignments(name, 7).read_id)
        un = raw[raw.full_length & ~raw.read_id.isin(assigned)].copy()
        real, decoys = load_whitelist(name)
        wl = real + decoys
        best, second, idx = best_matches(list(un.r1.str[UMI_LEN:UMI_LEN + BC_LEN]), wl)
        un = un[["read_id", "umi"]].assign(best=best, second=second, match=np.array(wl)[idx],
                                           is_decoy=idx >= len(real))
        un.to_parquet(path)
    relaxed[name] = pd.read_parquet(path)
pd.DataFrame({name: {"unassigned full-length reads": len(df), "median best distance": df.best.median()}
              for name, df in relaxed.items()}).T

# %% [markdown]
# ## Distance of unassigned reads to the whitelist
#
# A unique best match is required (best < second best).

# %%
rows = []
for name, df in relaxed.items():
    phi = calib.loc[(name, 7), "phi_eff"]
    for d in range(6, 14):
        sel = df[(df.best <= d) & (df.best < df.second)]
        dec = sel.is_decoy.mean() if len(sel) else np.nan
        rows.append({"layout": name, "max distance": d, "unique matches": len(sel),
                     "% of unassigned": 100 * len(sel) / len(df), "decoy fraction (%)": 100 * dec,
                     "estimated precision (%)": 100 * max(0, 1 - dec / phi),
                     "estimated correct reads": len(sel) * max(0, 1 - dec / phi)})
rescue = pd.DataFrame(rows).set_index(["layout", "max distance"])
rescue.round(2).to_csv(OUT / "relaxed_decoding_rates.csv")
rescue.round(2)

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
for name, df in relaxed.items():
    h = df.best.clip(upper=CUTOFF).value_counts(normalize=True).sort_index() * 100
    axes[0].plot(h.index, h.values, "o-", color=COLORS[name], label=name)
    r = rescue.loc[name]
    axes[1].plot(r.index, r["estimated precision (%)"], "o-", color=COLORS[name], label=name)
axes[0].set_xlabel("Levenshtein distance to best whitelist entry")
axes[0].set_ylabel("Unassigned full-length reads (%)")
axes[0].legend(frameon=False, fontsize=8)
axes[1].set_xlabel("Maximum distance for a relaxed assignment")
axes[1].set_ylabel("Estimated precision (%)")
fig.tight_layout()
fig.savefig(OUT / "relaxed_decoding.pdf")
plt.show()

# %% [markdown]
# ## Do zero-count features have molecules with too many errors?
#
# Relaxed assignments to real barcodes (unique best match, distance 8-10) are counted per feature (UMI-deduplicated) and compared between count groups. Under the decoding hypothesis, zero-count features receive rescued reads at a similar rate per feature as low- and medium-count features. The expected number of false rescues per feature, from the decoy rate, is shown for reference.

# %%
rows = []
for name, df in relaxed.items():
    phi = calib.loc[(name, 7), "phi_eff"]
    sel = df[(df.best >= 8) & (df.best <= 10) & (df.best < df.second)]
    real_hits = sel[~sel.is_decoy].drop_duplicates(["match", "umi"]).match.value_counts()
    n_dec = len(load_whitelist(name)[1])
    false_per_real = sel.is_decoy.sum() / phi * (1 - phi) / len(data[name])
    a = data[name][data[name].analysis].copy()
    a["rescued"] = real_hits.reindex(a.barcode).fillna(0).values
    g = a.groupby("group").rescued.agg(["mean", "size"]).rename(columns={"mean": "rescued UMIs per feature", "size": "features"})
    g["expected false rescues per feature"] = false_per_real
    g["layout"] = name
    rows.append(g)
by_group = pd.concat(rows).reset_index().set_index(["layout", "group"]).reindex(
    pd.MultiIndex.from_product([list(ARRAYS), ["zero", "low", "medium", "high"]]))
by_group.round(4).to_csv(OUT / "relaxed_rescue_by_group.csv")
by_group.round(4)

# %% [markdown]
# ## Half matches (Dual barcodes and Zipcodes)
#
# For unassigned full-length reads, the first 18 nt of the barcode window are compared with all first halves (Dual: block identifiers; Zipcodes: column codes) and the last 18 nt with all second halves (Dual: spot codes; Zipcodes: row codes). A read with one good half and one bad half comes from a molecule that is corrupted in one half only. Assigned reads are shown as reference.

# %%
def half_distances(seqs, halves):
    best = np.empty(len(seqs), dtype=int)
    for s in range(0, len(seqs), 5000):
        d = cdist(seqs[s:s + 5000], halves, scorer=Levenshtein.distance, dtype=np.uint8, workers=-1)
        best[s:s + 5000] = d.min(axis=1)
    return best


rng = np.random.default_rng(0)
rows = []
for name in ["Dual barcodes", "Zipcodes"]:
    lay = data[name]
    first, second = sorted(set(lay.barcode.str[:SEG])), sorted(set(lay.barcode.str[SEG:]))
    raw = load_raw(name).set_index("read_id")
    asg = set(load_assignments(name, 7).read_id)
    full = raw[raw.full_length]
    for label, ids in [("unassigned", [i for i in full.index if i not in asg]), ("assigned", [i for i in full.index if i in asg])]:
        ids = rng.choice(ids, size=min(20000, len(ids)), replace=False)
        w = full.loc[ids, "r1"].str[UMI_LEN:UMI_LEN + BC_LEN]
        d1 = half_distances(list(w.str[:SEG]), first)
        d2 = half_distances(list(w.str[SEG:]), second)
        rows.append({"layout": name, "reads": label, "n": len(ids),
                     "both halves ≤ 2 (%)": 100 * np.mean((d1 <= 2) & (d2 <= 2)),
                     "first only ≤ 2 (%)": 100 * np.mean((d1 <= 2) & (d2 > 2)),
                     "second only ≤ 2 (%)": 100 * np.mean((d1 > 2) & (d2 <= 2)),
                     "neither ≤ 2 (%)": 100 * np.mean((d1 > 2) & (d2 > 2)),
                     "median d first": np.median(d1), "median d second": np.median(d2)})
halves = pd.DataFrame(rows).set_index(["layout", "reads"])
halves.round(2).to_csv(OUT / "half_matches.csv")
halves.round(2)
