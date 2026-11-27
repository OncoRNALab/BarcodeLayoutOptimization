# %% [markdown]
# # Multivariate models of feature recovery
#
# Step 3 of `results_and_analysis_plan.md`. Response: UMI count per feature at r = 7, analysis area only (see `03_groups_and_tests`).
#
# 1. **Negative binomial regression** with standardized sequence features and the log local mean as covariate. Coefficients are reported as rate ratios per standard deviation.
# 2. **Gradient-boosted trees** (Poisson loss) with spatially blocked cross-validation (20 x 20 tiles held out together). Feature sets are compared by the fraction of Poisson deviance explained on held-out tiles ($D^2$).
# 3. **Shared-code units**: for Zipcodes, row and column effects from a two-way Poisson model, related to the sequence of the row and column codes; for Dual barcodes, block-level and within-block variation, related to the block identifier and to the spot code.

# %%
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import d2_tweedie_score
from sklearn.model_selection import GroupKFold

from common import ARRAYS, COLORS, OUT, SEG, setup_matplotlib, leftmost_schedule

plt = setup_matplotlib()
COUNT = "umis"
data = {name: pd.read_parquet(OUT / f"analysis_table_{cfg['key']}.parquet") for name, cfg in ARRAYS.items()}
SEQ = ["frac_A", "frac_C", "frac_G", "max_homopolymer", "n_GGG", "g4hunter", "mfe_barcode", "mfe_probe",
       "n_cycles", "hamming_neighbours", "stray_events", "nn_dist", "n_within_7", "kmer4_mean_log_freq"]

# %% [markdown]
# ## Negative binomial regression

# %%
rows = []
for name, df in data.items():
    a = df[df.analysis].copy()
    X = (a[SEQ] - a[SEQ].mean()) / a[SEQ].std()
    X = X.loc[:, a[SEQ].std() > 0]
    X["log_local_mean"] = np.log(a.local_mean + 0.1)
    X = sm.add_constant(X)
    fit = sm.NegativeBinomial(a[COUNT], X).fit(disp=0, maxiter=200)
    ci = fit.conf_int()
    for f in X.columns.drop("const"):
        rows.append({"layout": name, "feature": f, "rate_ratio_per_sd": np.exp(fit.params[f]),
                     "lo": np.exp(ci.loc[f, 0]), "hi": np.exp(ci.loc[f, 1]), "p": fit.pvalues[f]})
    print(f"{name}: alpha = {fit.params['alpha']:.3f}, pseudo-R2 = {fit.prsquared:.4f}")
nb = pd.DataFrame(rows)
nb.to_csv(OUT / "nb_regression.csv", index=False)
nb.pivot(index="feature", columns="layout", values="rate_ratio_per_sd").round(3)

# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
order = SEQ + ["log_local_mean"]
for ax, name in zip(axes, ARRAYS):
    d = nb[(nb.layout == name) & (nb.feature != "log_local_mean")].set_index("feature").reindex(SEQ).dropna()
    y = np.arange(len(d))
    ax.errorbar(d.rate_ratio_per_sd, y, xerr=[d.rate_ratio_per_sd - d.lo, d.hi - d.rate_ratio_per_sd],
                fmt="o", color=COLORS[name], ms=4)
    ax.axvline(1, color="grey", lw=0.8)
    ax.set_yticks(y, d.index)
    ax.set_title(name)
    ax.set_xlabel("Rate ratio per SD (95% CI)")
fig.tight_layout()
fig.savefig(OUT / "nb_rate_ratios.pdf")
plt.show()

# %% [markdown]
# ## Gradient boosting with spatially blocked cross-validation
#
# Feature sets:
# * `local`: log local mean only (smooth spatial variation);
# * `local + row/col`: adds row and column index. In Zipcodes, row and column index are identical to the row and column codes, so this set also captures code effects;
# * `sequence`: sequence features only;
# * `local + sequence`: the question of interest: does sequence add to the local expectation?
# * `all`.

# %%
SETS = {
    "local": ["log_local_mean"],
    "local + row/col": ["log_local_mean", "row", "col"],
    "sequence": SEQ,
    "local + sequence": ["log_local_mean"] + SEQ,
    "all": ["log_local_mean", "row", "col"] + SEQ,
}
rows, importances = [], []
for name, df in data.items():
    a = df[df.analysis].copy()
    a["log_local_mean"] = np.log(a.local_mean + 0.1)
    tiles = (a.row // 20).astype(str) + "_" + (a.col // 20).astype(str)
    y = a[COUNT].values
    for set_name, cols in SETS.items():
        pred = np.zeros(len(a))
        for train, test in GroupKFold(n_splits=5).split(a, y, tiles):
            m = HistGradientBoostingRegressor(loss="poisson", max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                              min_samples_leaf=100, random_state=0)
            m.fit(a[cols].values[train], y[train])
            pred[test] = m.predict(a[cols].values[test])
        rows.append({"layout": name, "features": set_name, "D2_heldout": d2_tweedie_score(y, pred, power=1)})
        if set_name == "local + sequence":
            train, test = next(GroupKFold(n_splits=5).split(a, y, tiles))
            m.fit(a[cols].values[train], y[train])
            imp = permutation_importance(m, a[cols].values[test], y[test], n_repeats=5, random_state=0,
                                         scoring=lambda est, X, yy: d2_tweedie_score(yy, est.predict(X), power=1))
            importances.append(pd.Series(imp.importances_mean, index=cols, name=name))
gbm = pd.DataFrame(rows).pivot(index="features", columns="layout", values="D2_heldout").reindex(SETS)
gbm.to_csv(OUT / "gbm_heldout_D2.csv")
importance = pd.concat(importances, axis=1)
importance.to_csv(OUT / "gbm_permutation_importance.csv")
display((gbm * 100).round(2))
(importance * 100).round(3).sort_values("Zipcodes", ascending=False)

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
(gbm * 100).T.plot.bar(ax=axes[0], rot=0)
axes[0].set_ylabel("Held-out Poisson deviance explained (%)")
axes[0].legend(frameon=False, fontsize=7)
(importance * 100).plot.barh(ax=axes[1], color=[COLORS[c] for c in importance.columns])
axes[1].set_xlabel("Permutation importance (loss of D², percentage points)")
axes[1].legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "gbm_feature_sets_and_importance.pdf")
plt.show()

# %% [markdown]
# ## Zipcodes: row-code and column-code effects
#
# Two-way Poisson model $\log \mu = \beta_0 + a_{\text{row}} + b_{\text{col}}$ fitted on the analysis area. Each row effect belongs to one row code (last 18 nt) and each column effect to one column code (first 18 nt). The effects are related to the sequence features of the codes.

# %%
def code_features(codes):
    out = pd.DataFrame({"code": codes})
    for b in "ACGT":
        out[f"frac_{b}"] = out.code.str.count(b) / SEG
    out["max_homopolymer"] = out.code.map(lambda s: max(len(r) for r in __import__("re").findall(r"(A+|C+|G+|T+)", s)))
    out["n_GG"] = out.code.map(lambda s: sum(s.startswith("GG", i) for i in range(len(s) - 1)))
    out["n_cycles"] = out.code.map(lambda s: leftmost_schedule(s)[-1] + 1)
    return out


z = data["Zipcodes"]
a = z[z.analysis]
X = pd.get_dummies(a[["row", "col"]].astype(str), drop_first=True, dtype=float)
fit = sm.GLM(a[COUNT].values, sm.add_constant(X), family=sm.families.Poisson()).fit()
row_eff = pd.Series({int(k.split("_")[1]): v for k, v in fit.params.items() if k.startswith("row_")})
col_eff = pd.Series({int(k.split("_")[1]): v for k, v in fit.params.items() if k.startswith("col_")})
row_eff.loc[int(sorted(a.row.astype(str).unique())[0])] = 0.0
col_eff.loc[int(sorted(a.col.astype(str).unique())[0])] = 0.0
units = []
for axis, eff, part in [("row", row_eff, slice(SEG, 2 * SEG)), ("col", col_eff, slice(0, SEG))]:
    codes = a.groupby(axis).barcode.first().str[part]
    cf = code_features(codes.values)
    cf[axis] = codes.index
    cf["effect"] = eff.reindex(codes.index).values
    cf["axis"] = axis
    units.append(cf)
units = pd.concat(units)
units.to_csv(OUT / "zipcodes_code_effects.csv", index=False)

rows = []
for axis in ["row", "col"]:
    u = units[units.axis == axis]
    for f in ["frac_A", "frac_C", "frac_G", "frac_T", "max_homopolymer", "n_GG", "n_cycles"]:
        rho, p = spearmanr(u[f], u.effect)
        rows.append({"axis": axis, "feature": f, "spearman_rho": rho, "p": p, "n_codes": len(u)})
zip_corr = pd.DataFrame(rows)
var = pd.DataFrame({"axis": ["row", "col"], "variance of log effects": [row_eff.var(), col_eff.var()]})
display(var.round(4))
zip_corr.pivot(index="feature", columns="axis", values="spearman_rho").round(3)

# %% [markdown]
# Partial association: does anything besides G content explain the code effects? Linear regression of the effect on G fraction, then Spearman correlation of the residual with the other features.

# %%
rows = []
for axis in ["row", "col"]:
    u = units[units.axis == axis]
    resid = sm.OLS(u.effect, sm.add_constant(u.frac_G)).fit().resid
    for f in ["frac_A", "frac_C", "frac_T", "max_homopolymer", "n_GG", "n_cycles"]:
        rho, p = spearmanr(u[f], resid)
        rows.append({"axis": axis, "feature": f, "rho_with_residual": rho, "p": p})
pd.DataFrame(rows).pivot(index="feature", columns="axis", values="rho_with_residual").round(3)

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, axis in zip(axes, ["col", "row"]):
    u = units[units.axis == axis]
    ax.scatter(u.frac_G * 100, u.effect, s=10, alpha=0.6, color=COLORS["Zipcodes"])
    rho = spearmanr(u.frac_G, u.effect)[0]
    ax.set_title(f"Zipcodes {axis} codes (rho = {rho:.2f})")
    ax.set_xlabel(f"G content of {axis} code (%)")
axes[0].set_ylabel("log effect (two-way Poisson model)")
fig.tight_layout()
fig.savefig(OUT / "zipcodes_code_effects_vs_G.pdf")
plt.show()

# %% [markdown]
# ## Dual barcodes: block identifier vs spot code
#
# * **Block level:** mean UMI count of each 4 x 4 block relative to the mean of its 8 neighbouring blocks (log2 ratio), related to the features of the block identifier (first 18 nt).
# * **Within block:** log2 ratio of each feature's count to its block mean, related to the features of its spot code (last 18 nt). Because the block identifier is shared within a block, this isolates the spot code.

# %%
d = data["Dual barcodes"]
a = d[d.analysis].copy()
a["brow"], a["bcol"] = a.row // 4, a.col // 4
blocks = a.groupby(["brow", "bcol"]).agg(mean=(COUNT, "mean"), n=(COUNT, "size"), block_id=("barcode", lambda s: s.iloc[0][:SEG]))
blocks = blocks[blocks.n >= 12]
bgrid = blocks["mean"].unstack()
neigh = {}
for (br, bc) in blocks.index:
    vals = [bgrid.at[br + i, bc + j] for i in (-1, 0, 1) for j in (-1, 0, 1)
            if (i or j) and (br + i) in bgrid.index and (bc + j) in bgrid.columns and not np.isnan(bgrid.at[br + i, bc + j])]
    neigh[(br, bc)] = np.mean(vals) if vals else np.nan
blocks["neigh_mean"] = pd.Series(neigh)
blocks["log2_ratio"] = np.log2((blocks["mean"] + 0.1) / (blocks.neigh_mean + 0.1))
bf = code_features(blocks.block_id.values)
bf.index = blocks.index
blocks = blocks.join(bf.drop(columns="code"))

a = a.merge(blocks[["mean"]].rename(columns={"mean": "block_mean"}), left_on=["brow", "bcol"], right_index=True)
a["within_log2_ratio"] = np.log2((a[COUNT] + 0.5) / (a.block_mean + 0.5))
sf = code_features(a.barcode.str[SEG:].values)
sf.index = a.index
a = a.join(sf.drop(columns="code"), rsuffix="_spot")

rows = []
for f in ["frac_A", "frac_C", "frac_G", "frac_T", "max_homopolymer", "n_GG", "n_cycles"]:
    rho_b, p_b = spearmanr(blocks[f], blocks.log2_ratio, nan_policy="omit")
    rho_s, p_s = spearmanr(a[f + "_spot"], a.within_log2_ratio)
    rows.append({"feature": f, "block identifier: rho": rho_b, "p_block": p_b, "spot code: rho": rho_s, "p_spot": p_s})
dual_units = pd.DataFrame(rows).set_index("feature")
dual_units.to_csv(OUT / "dual_block_vs_spot.csv")
dual_units.round(4)

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for ax, (x, y, title) in zip(axes, [(blocks.frac_G * 100, blocks.log2_ratio, "Block identifier"),
                                     (a.frac_G_spot * 100, a.within_log2_ratio, "Spot code")]):
    g = pd.DataFrame({"x": x.round(1), "y": y}).groupby("x").y.agg(["mean", "sem", "size"])
    g = g[g["size"] >= 30]
    ax.errorbar(g.index, g["mean"], yerr=1.96 * g["sem"], fmt="o-", color=COLORS["Dual barcodes"], ms=4)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlabel(f"G content of {title.lower()} (%)")
    ax.set_title(f"Dual barcodes: {title.lower()}")
axes[0].set_ylabel("Block log2 ratio to neighbouring blocks")
axes[1].set_ylabel("Feature log2 ratio to block mean")
fig.tight_layout()
fig.savefig(OUT / "dual_block_vs_spot_G.pdf")
plt.show()
