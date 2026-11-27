#!/usr/bin/env python3
"""Group tests, k-mer enrichment, negative-binomial models and sampling baseline.

Requires barcode_features.csv (and optionally simulation_recall.csv).
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

from io_utils import ARRAYS, OUT_DIR, SEG

FEATURES = [
    "G_frac", "G_frac_half1", "G_frac_half2", "GC_frac", "purine_frac",
    "longest_homopolymer", "GG_count", "GGG_count", "G4_score",
    "n_synthesis_cycles", "cycles_per_nt", "neighbour_schedule_hamming",
    "nn_hamming", "n_hamming_le8", "mean_4mer_share",
    "adapter_rc_k", "self_rc_k",
]
RNG = np.random.default_rng(2026)


def cliffs_delta(x, y, max_n=4000):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) > max_n:
        x = RNG.choice(x, max_n, replace=False)
    if len(y) > max_n:
        y = RNG.choice(y, max_n, replace=False)
    # P(x>y) - P(x<y)
    gt = 0
    lt = 0
    for xi in x:
        gt += np.sum(xi > y)
        lt += np.sum(xi < y)
    n = len(x) * len(y)
    return (gt - lt) / n if n else np.nan


def sampling_baseline(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        lam = a.reads_r7.mean()
        var = a.reads_r7.var()
        obs0 = (a.reads_r7 == 0).mean()
        pois0 = float(np.exp(-lam))
        # method of moments NB: var = mu + mu^2/k
        k = lam ** 2 / max(var - lam, 1e-9)
        nb0 = (k / (k + lam)) ** k
        loc = a.local_median_r7.clip(lower=0)
        exp0_local = np.exp(-loc).mean()
        rows.append({
            "layout": name,
            "n_analysis": int(len(a)),
            "mean_reads": lam,
            "variance": var,
            "var_over_mean": var / lam if lam else np.nan,
            "zero_observed": obs0,
            "zero_poisson": pois0,
            "zero_nb": float(nb0),
            "zero_poisson_local": float(exp0_local),
            "excess_zeros_poisson": obs0 - pois0,
            "excess_zeros_local": obs0 - exp0_local,
            "mean_umis": a.umis_r7.mean(),
        })
    return pd.DataFrame(rows)


def group_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        for contrast, left, right in (("zero_vs_high", "zero", "high"), ("low_vs_high", "low", "high")):
            A = a[a.count_group == left]
            B = a[a.count_group == right]
            for feat in FEATURES:
                x, y = A[feat], B[feat]
                if x.notna().sum() < 20 or y.notna().sum() < 20:
                    continue
                try:
                    u, p = stats.mannwhitneyu(x.dropna(), y.dropna(), alternative="two-sided")
                except ValueError:
                    u, p = np.nan, np.nan
                rows.append({
                    "layout": name,
                    "contrast": contrast,
                    "feature": feat,
                    "n_left": int(x.notna().sum()),
                    "n_right": int(y.notna().sum()),
                    "median_left": float(x.median()),
                    "median_right": float(y.median()),
                    "cliffs_delta": cliffs_delta(x, y),
                    "p_raw": p,
                })
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj"] = multipletests(out.p_raw.fillna(1), method="fdr_bh")[1]
    return out


def kmer_enrichment(df: pd.DataFrame, k: int = 4) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        z = a[a.count_group == "zero"]
        h = a[a.count_group == "high"]
        if len(z) < 50 or len(h) < 50:
            continue

        def present(seqs):
            bags = []
            for s in seqs:
                bags.append(set(s[i : i + k] for i in range(len(s) - k + 1)))
            return bags

        z_bags, h_bags = present(z.barcode), present(h.barcode)
        vocab = set().union(*z_bags, *h_bags)
        nz, nh = len(z_bags), len(h_bags)
        for mer in vocab:
            a_ = sum(mer in b for b in z_bags)
            c_ = sum(mer in b for b in h_bags)
            table = [[a_, nz - a_], [c_, nh - c_]]
            oddsratio, p = stats.fisher_exact(table)
            rows.append({
                "layout": name, "k": k, "kmer": mer,
                "frac_zero": a_ / nz, "frac_high": c_ / nh,
                "oddsratio": oddsratio, "p_raw": p,
            })
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj"] = multipletests(out.p_raw.fillna(1), method="fdr_bh")[1]
        out["log2_or"] = np.log2(out.oddsratio.replace(0, np.nan))
    return out


def position_composition(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        for grp, sub in a.groupby("count_group"):
            enc = np.array([list(s) for s in sub.barcode])
            for pos in range(enc.shape[1]):
                col = enc[:, pos]
                for b in "ACGT":
                    rows.append({
                        "layout": name, "group": grp, "position": pos + 1,
                        "base": b, "frequency": (col == b).mean(), "n": len(sub),
                    })
    return pd.DataFrame(rows)


def unit_level(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        if name == "Greedy + 2-OPT":
            continue
        for unit, seqcol, label in (("unit_1", "bar_1", "half1"), ("unit_2", "bar_2", "half2")):
            t = a.groupby(unit).agg(
                n=("reads_r7", "size"),
                mean_reads=("reads_r7", "mean"),
                seq=(seqcol, "first"),
            )
            t["G_frac"] = t.seq.str.count("G") / SEG
            rho, p = stats.spearmanr(t.G_frac, t.mean_reads)
            rows.append({
                "layout": name, "unit": label, "n_units": len(t),
                "between_unit_var": float(t.mean_reads.var()),
                "spearman_G": float(rho), "p_G": float(p),
            })
    return pd.DataFrame(rows)


def deviance(y, mu):
    y = np.asarray(y, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-8, None)
    return float(2 * np.sum(y * np.log((y + 1e-8) / mu) - (y - mu)))


def nb_models(df: pd.DataFrame) -> pd.DataFrame:
    model_feats = [
        "G_frac", "longest_homopolymer", "GG_count", "cycles_per_nt",
        "neighbour_schedule_hamming", "nn_hamming", "mean_4mer_share",
        "adapter_rc_k",
    ]
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set].dropna(subset=model_feats + ["reads_r7", "local_median_r7", "row", "col"])
        y = a.reads_r7.to_numpy(dtype=float)
        offset = np.log(a.local_median_r7.clip(lower=0) + 0.5)
        spat = sm.add_constant(pd.DataFrame({
            "row": a.row / 235, "col": a.col / 363,
        }))
        seq = sm.add_constant(a[model_feats].apply(lambda s: (s - s.mean()) / (s.std() + 1e-9)))
        both = sm.add_constant(pd.concat([
            spat.drop(columns="const"),
            a[model_feats].apply(lambda s: (s - s.mean()) / (s.std() + 1e-9)),
        ], axis=1))
        fam = sm.families.NegativeBinomial()
        fits = {}
        for label, X, off in (
            ("null", np.ones((len(a), 1)), None),
            ("spatial", spat, None),
            ("sequence", seq, None),
            ("spatial_offset", np.ones((len(a), 1)), offset),
            ("sequence_plus_offset", seq, offset),
            ("both", both, None),
        ):
            try:
                res = sm.GLM(y, X, family=fam, offset=off).fit(disp=False, maxiter=100)
                mu = res.fittedvalues
                fits[label] = res
                rows.append({
                    "layout": name, "model": label,
                    "pseudo_r2": 1 - res.deviance / res.null_deviance if res.null_deviance else np.nan,
                    "deviance": float(res.deviance),
                    "aic": float(res.aic),
                    "n": len(a),
                })
            except Exception as exc:
                rows.append({"layout": name, "model": label, "error": str(exc), "n": len(a)})
        # blocked CV for sequence_plus_offset vs offset-only
        tiles = (a.row.to_numpy() // 20) * 100 + (a.col.to_numpy() // 20)
        uniq = np.unique(tiles)
        RNG.shuffle(uniq)
        folds = np.array_split(uniq, 5)
        for model_name, X, off in (("offset", np.ones((len(a), 1)), offset), ("seq+offset", seq, offset)):
            preds = np.zeros(len(a))
            Xv = np.asarray(X)
            for held in folds:
                test = np.isin(tiles, held)
                train = ~test
                if train.sum() < 100 or test.sum() < 20:
                    continue
                off_tr = off[train] if off is not None else None
                res = sm.GLM(y[train], Xv[train], family=fam, offset=off_tr).fit(disp=False, maxiter=80)
                off_te = off[test] if off is not None else 0
                eta = Xv[test] @ res.params + (off_te if off is not None else 0)
                preds[test] = np.exp(eta)
            rows.append({
                "layout": name, "model": f"cv_{model_name}",
                "poisson_deviance": deviance(y, preds),
                "n": len(a),
            })
    return pd.DataFrame(rows)


def add_groups(df: pd.DataFrame) -> pd.DataFrame:
    def group(n):
        if n == 0:
            return "zero"
        if n == 1:
            return "low"
        if n <= 4:
            return "medium"
        return "high"

    df = df.copy()
    df["count_group"] = df.reads_r7.map(group)
    return df


def group_stability(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, g in df.groupby("layout"):
        a = g[g.analysis_set]
        for r in (5, 6, 7):
            col = f"reads_r{r}"
            gg = a[col].map(lambda n: "zero" if n == 0 else "low" if n == 1 else "medium" if n <= 4 else "high")
            rows.append({"layout": name, "threshold": r, **gg.value_counts(normalize=True).to_dict(),
                         "n": len(a)})
    return pd.DataFrame(rows)


def main():
    df = add_groups(pd.read_csv(OUT_DIR / "barcode_features.csv"))
    sim_path = OUT_DIR / "simulation_recall.csv"
    if sim_path.exists():
        sim = pd.read_csv(sim_path)
        df = df.merge(sim[["layout", "barcode", "sim_recall", "sim_n_correct"]], on=["layout", "barcode"], how="left")
    sampling = sampling_baseline(df)
    tests = group_tests(df)
    kmers = kmer_enrichment(df, k=3)
    pos = position_composition(df)
    units = unit_level(df)
    models = nb_models(df)
    stab = group_stability(df)
    sampling.to_csv(OUT_DIR / "sampling_baseline.csv", index=False)
    tests.to_csv(OUT_DIR / "group_tests.csv", index=False)
    kmers.to_csv(OUT_DIR / "kmer_enrichment_k3.csv", index=False)
    pos.to_csv(OUT_DIR / "position_composition.csv", index=False)
    units.to_csv(OUT_DIR / "unit_level.csv", index=False)
    models.to_csv(OUT_DIR / "nb_models.csv", index=False)
    stab.to_csv(OUT_DIR / "group_stability.csv", index=False)
    top = tests.sort_values("p_adj").groupby(["layout", "contrast"]).head(5)
    print("sampling")
    print(sampling.to_string(index=False))
    print("top tests")
    print(top.to_string(index=False))
    print("units")
    print(units.to_string(index=False))
    print("models")
    print(models.to_string(index=False))
    (OUT_DIR / "analyse_summary.json").write_text(json.dumps({
        "n_tests": int(len(tests)),
        "n_sig": int((tests.p_adj < 0.05).sum()) if len(tests) else 0,
    }))


if __name__ == "__main__":
    main()
