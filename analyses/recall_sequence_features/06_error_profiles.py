# %% [markdown]
# # Error profiles of the printed molecules
#
# Step 5. Every full-length read assigned to a real barcode at r = 7 is aligned to its called barcode, separately in Read 1 and in Read 2 (reverse complemented). Read 2 covers the barcode from the other end, so the two reads of a pair come from the same molecule but have independent sequencing errors.
#
# * An edit (position, type, base) found in **both** reads of a pair is in the molecule: a synthesis error (or an early PCR error).
# * An edit found in only one read is a sequencing error (or an alignment ambiguity).
#
# Alignment: the barcode is aligned as an infix (edlib, `HW` mode) to Read 1 positions 9-47 and to positions 10-47 of the reverse-complemented Read 2, so that indels in the UMI or barcode can shift the window by a few bases. Edit types are expressed relative to the designed barcode: **insertion** = extra base in the molecule, **deletion** = designed base missing, **substitution** = wrong base.
#
# Limitation: only reads assigned at r = 7 are included, i.e. molecules with at most about 7 edits; molecules with more errors are not represented.

# %%
import edlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import ARRAYS, BC_LEN, COLORS, OUT, load_assignments, load_raw, revcomp, setup_matplotlib

plt = setup_matplotlib()
data = {name: pd.read_parquet(OUT / f"analysis_table_{cfg['key']}.parquet") for name, cfg in ARRAYS.items()}


def edits(barcode, read):
    """Set of (barcode position, type, base) edits of read relative to barcode; position of an insertion = barcode position before it."""
    aln = edlib.align(barcode, read, mode="HW", task="path")
    nice = edlib.getNiceAlignment(aln, barcode, read)
    q, t = nice["query_aligned"], nice["target_aligned"]
    out, pos = set(), -1
    for a, b in zip(q, t):
        if a == "-":
            out.add((pos, "insertion", b))
        else:
            pos += 1
            if b == "-":
                out.add((pos, "deletion", a))
            elif a != b:
                out.add((pos, "substitution", b))
    return out, aln["editDistance"]


rows, per_read = [], []
for name in ARRAYS:
    raw = load_raw(name).set_index("read_id")
    asg = load_assignments(name, 7)
    asg = asg[asg.is_real]
    asg = asg[raw.loc[asg.read_id, "full_length"].values]
    r1 = raw.loc[asg.read_id, "r1"].values
    r2 = raw.loc[asg.read_id, "r2"].values
    for bc, a1, a2 in zip(asg.called_seq.values, r1, r2):
        e1, d1 = edits(bc, a1[8:])
        e2, d2 = edits(bc, revcomp(a2)[9:])
        both = e1 & e2
        per_read.append((name, bc, d1, d2, len(both),
                         sum(e[1] == "insertion" for e in both), sum(e[1] == "deletion" for e in both),
                         sum(e[1] == "substitution" for e in both)))
        for e in both:
            rows.append((name, e[0], e[1], e[2], "molecule"))
        for e in e1 - e2:
            rows.append((name, e[0], e[1], e[2], "read 1 only"))
        for e in e2 - e1:
            rows.append((name, e[0], e[1], e[2], "read 2 only"))
per_read = pd.DataFrame(per_read, columns=["layout", "barcode", "dist_r1", "dist_r2", "mol_edits", "mol_ins", "mol_del", "mol_sub"])
edits_df = pd.DataFrame(rows, columns=["layout", "position", "type", "base", "source"])
per_read.to_parquet(OUT / "error_profiles_per_read.parquet")
edits_df.to_parquet(OUT / "error_profiles_edits.parquet")
per_read.groupby("layout")[["dist_r1", "dist_r2", "mol_edits", "mol_ins", "mol_del", "mol_sub"]].mean().round(3)

# %% [markdown]
# ## Error rates per nucleotide, by source and type

# %%
n_reads = per_read.groupby("layout").size()
rates = (edits_df.groupby(["layout", "source", "type"]).size().unstack("type").div(n_reads * BC_LEN, axis=0, level="layout") * 100)
rates.round(3).to_csv(OUT / "error_rates_per_nt.csv")
rates.round(3)

# %% [markdown]
# ## Molecular error rate along the barcode (synthesis order)

# %%
fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), sharey=True)
for ax, t in zip(axes, ["insertion", "deletion", "substitution"]):
    for name in ARRAYS:
        e = edits_df[(edits_df.layout == name) & (edits_df.source == "molecule") & (edits_df.type == t)]
        prof = e.position.value_counts().reindex(range(-1, BC_LEN), fill_value=0) / n_reads[name] * 100
        ax.plot(prof.index + 1, prof.values, "-", color=COLORS[name], label=name)
    ax.set_title(f"Molecular {t}s")
    ax.set_xlabel("Barcode position (synthesis order)")
axes[0].set_ylabel("Per read, % at position")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "molecular_error_profile.pdf")
plt.show()

# %% [markdown]
# Bases involved in molecular errors: the base inserted, the designed base that was deleted, and the base substituted in.

# %%
mol = edits_df[edits_df.source == "molecule"]
bases = mol.groupby(["layout", "type"]).base.value_counts(normalize=True).unstack().reindex(columns=list("ACGT")) * 100
bases.round(1).to_csv(OUT / "molecular_error_bases.csv")
bases.round(1)

# %% [markdown]
# ## Molecular errors per feature, stray-light exposure and count group
#
# Per feature: mean molecular insertions, deletions and substitutions per read (features with at least 2 full-length reads). The stray-light hypothesis predicts more insertions where neighbours are illuminated while the feature is not (`stray_events`).

# %%
per_feature = per_read.groupby(["layout", "barcode"]).agg(n=("mol_edits", "size"), mol_ins=("mol_ins", "mean"),
                                                          mol_del=("mol_del", "mean"), mol_sub=("mol_sub", "mean")).reset_index()
rows = []
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
for ax, name in zip(axes, ARRAYS):
    df = data[name][data[name].analysis].merge(per_feature[per_feature.layout == name], on="barcode")
    df = df[df.n >= 2]
    for t in ["mol_ins", "mol_del", "mol_sub"]:
        for f in ["stray_events", "stray_cycles", "hamming_neighbours", "frac_G", "n_cycles"]:
            rho, p = spearmanr(df[f], df[t])
            rows.append({"layout": name, "error": t, "feature": f, "spearman_rho": rho, "p": p, "n_features": len(df)})
    df["bin"] = pd.qcut(df.stray_events.rank(method="first"), 10, labels=False)
    g = df.groupby("bin").agg(x=("stray_events", "mean"), ins=("mol_ins", "mean"), dele=("mol_del", "mean"), sub=("mol_sub", "mean"))
    ax.plot(g.x, g.ins, "o-", label="insertions")
    ax.plot(g.x, g.dele, "s-", label="deletions")
    ax.plot(g.x, g.sub, "^-", label="substitutions")
    ax.set_title(name)
    ax.set_xlabel("Stray-light events (decile mean)")
axes[0].set_ylabel("Molecular edits per read")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "molecular_errors_vs_stray_light.pdf")
plt.show()
stray = pd.DataFrame(rows)
stray.to_csv(OUT / "molecular_errors_vs_features.csv", index=False)
stray.pivot_table(index=["error", "feature"], columns="layout", values="spearman_rho").round(3)

# %% [markdown]
# Molecular errors by count group (features with at least 1 full-length read; zero-count features have no reads by definition).

# %%
rows = []
for name in ARRAYS:
    df = data[name][data[name].analysis].merge(per_feature[per_feature.layout == name], on="barcode")
    rows.append(df.groupby("group")[["mol_ins", "mol_del", "mol_sub"]].mean().assign(layout=name))
pd.concat(rows).reset_index().set_index(["layout", "group"]).round(3)

# %% [markdown]
# Stray-light exposure per layout: the quantity the layouts were designed to reduce.

# %%
pd.DataFrame({name: df[df.analysis][["stray_events", "stray_cycles", "hamming_neighbours"]].mean() for name, df in data.items()}).T.round(2)
