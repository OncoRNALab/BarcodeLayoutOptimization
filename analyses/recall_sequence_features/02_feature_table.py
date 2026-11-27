# %% [markdown]
# # Per-barcode sequence features
#
# One row per printed feature, for the three layouts. The table is cached in `output/features_<layout>.parquet`; delete the files to recompute (about 25 minutes on 8 cores, dominated by the nearest-neighbour search).
#
# | Group | Feature | Definition |
# |---|---|---|
# | Composition | `frac_A` ... `frac_T`, `GC`, `purine` | Fraction over the 36 nt |
# | | `G_first18`, `G_last18` | G fraction of each 18-nt half (Dual: block identifier / spot code; Zipcodes: column / row code) |
# | Motifs | `max_homopolymer`, `n_GG`, `n_GGG` | Longest run; overlapping GG and GGG counts |
# | | `g4hunter` | G4Hunter score (runs of G scored +1 to +4, runs of C -1 to -4, averaged) |
# | Synthesis | `n_cycles` | Cycles needed under left-most embedding in the A, C, G, T order |
# | | `hamming_neighbours` | Mean Hamming distance between the feature's schedule and those of its 8 neighbours |
# | | `stray_events` | Number of (cycle, neighbour) pairs in which a neighbour is illuminated while the feature is not, over all barcode cycles of the array |
# | | `stray_cycles` | Number of cycles in which at least one neighbour is illuminated while the feature is not |
# | Structure | `mfe_barcode`, `mfe_probe` | Minimum free energy (kcal/mol, ViennaRNA, DNA parameters of Mathews 2004) of the barcode alone, and of Read 1 handle + N10 + barcode + Read 2 handle |
# | Decoding | `nn_dist` | Levenshtein distance to the nearest other entry of the combined whitelist (real + decoys) |
# | | `n_within_7` | Number of whitelist entries within distance 7 |
# | | `kmer4_mean_log_freq`, `kmer4_min_freq` | Library frequency of the barcode's 4-mers (number of whitelist entries containing each 4-mer) |
#
# Assumption: the arrays were synthesized with left-most embedding of the barcodes in the repeating A, C, G, T order, all barcodes starting in the same cycle after the UMI.

# %%
from concurrent.futures import ProcessPoolExecutor
from itertools import product
import time

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist

from common import ARRAYS, OUT, ROWS, COLS, SEG, BC_LEN, fold_batch, load_layout, load_whitelist, leftmost_schedule

NN_CUTOFF = 20


def g4hunter(seq):
    scores = np.zeros(len(seq))
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if seq[i] in "GC":
            scores[i:j] = min(j - i, 4) * (1 if seq[i] == "G" else -1)
        i = j
    return scores.mean()


def max_run(seq):
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def count_overlapping(seq, motif):
    return sum(seq.startswith(motif, i) for i in range(len(seq) - len(motif) + 1))


def schedule_features(lay):
    cycles = [leftmost_schedule(s) for s in lay.barcode]
    T = max(c[-1] for c in cycles) + 1
    sched = np.zeros((ROWS, COLS, T), dtype=bool)
    for (r, c), cyc in zip(zip(lay.row, lay.col), cycles):
        sched[r, c, cyc] = True
    pad = np.zeros((ROWS + 2, COLS + 2, T), dtype=bool)
    pad[1:-1, 1:-1] = sched
    present = np.zeros((ROWS + 2, COLS + 2), dtype=bool)
    present[1:-1, 1:-1] = True
    ham = np.zeros((ROWS, COLS))
    n_nb = np.zeros((ROWS, COLS))
    events = np.zeros((ROWS, COLS))
    any_lit = np.zeros((ROWS, COLS, T), dtype=bool)
    for dr, dc in product((-1, 0, 1), repeat=2):
        if dr == dc == 0:
            continue
        nb = pad[1 + dr:ROWS + 1 + dr, 1 + dc:COLS + 1 + dc]
        ok = present[1 + dr:ROWS + 1 + dr, 1 + dc:COLS + 1 + dc]
        ham += ok * (nb ^ sched).sum(axis=2)
        n_nb += ok
        lit_while_dark = nb & ~sched
        events += lit_while_dark.sum(axis=2)
        any_lit |= lit_while_dark
    feats = pd.DataFrame({
        "n_cycles": [c[-1] + 1 for c in cycles],
        "hamming_neighbours": (ham / n_nb)[lay.row, lay.col],
        "stray_events": events[lay.row, lay.col],
        "stray_cycles": any_lit.sum(axis=2)[lay.row, lay.col],
    })
    return feats, T


def decoding_features(name, lay):
    real, decoys = load_whitelist(name)
    wl = real + decoys
    index = {s: i for i, s in enumerate(wl)}
    nn, within = np.empty(len(lay), dtype=int), np.empty(len(lay), dtype=int)
    chunk = 2000
    for start in range(0, len(lay), chunk):
        q = list(lay.barcode[start:start + chunk])
        d = cdist(q, wl, scorer=Levenshtein.distance, dtype=np.uint8, workers=-1, score_cutoff=NN_CUTOFF)
        d[np.arange(len(q)), [index[s] for s in q]] = 255
        nn[start:start + chunk] = d.min(axis=1)
        within[start:start + chunk] = (d <= 7).sum(axis=1)
    kfreq = {}
    for s in wl:
        for k in {s[i:i + 4] for i in range(len(s) - 3)}:
            kfreq[k] = kfreq.get(k, 0) + 1
    kmer_mean, kmer_min = [], []
    for s in lay.barcode:
        f = np.array([kfreq[s[i:i + 4]] for i in range(len(s) - 3)])
        kmer_mean.append(np.log10(f).mean())
        kmer_min.append(f.min())
    return pd.DataFrame({"nn_dist": nn, "n_within_7": within,
                         "kmer4_mean_log_freq": kmer_mean, "kmer4_min_freq": kmer_min})


def build(name):
    lay = load_layout(name)
    s = lay.barcode
    feats = pd.DataFrame({"row": lay.row, "col": lay.col, "barcode": s})
    for b in "ACGT":
        feats[f"frac_{b}"] = s.str.count(b) / BC_LEN
    feats["GC"] = feats.frac_G + feats.frac_C
    feats["purine"] = feats.frac_A + feats.frac_G
    feats["G_first18"] = s.str[:SEG].str.count("G") / SEG
    feats["G_last18"] = s.str[SEG:].str.count("G") / SEG
    feats["max_homopolymer"] = s.map(max_run)
    feats["n_GG"] = s.map(lambda x: count_overlapping(x, "GG"))
    feats["n_GGG"] = s.map(lambda x: count_overlapping(x, "GGG"))
    feats["g4hunter"] = s.map(g4hunter)
    sched, T = schedule_features(lay)
    feats = pd.concat([feats, sched], axis=1)
    t = time.time()
    batches = [list(s[i:i + 2000]) for i in range(0, len(s), 2000)]
    with ProcessPoolExecutor() as ex:
        folds = [x for batch in ex.map(fold_batch, batches) for x in batch]
    feats["mfe_barcode"] = [f[0] for f in folds]
    feats["mfe_probe"] = [f[1] for f in folds]
    print(f"{name}: folding {time.time() - t:.0f}s, {T} barcode cycles")
    t = time.time()
    feats = pd.concat([feats, decoding_features(name, lay)], axis=1)
    print(f"{name}: nearest-neighbour search {time.time() - t:.0f}s")
    feats["edge_dist"] = np.minimum.reduce([feats.row, feats.col, ROWS - 1 - feats.row, COLS - 1 - feats.col])
    return feats


features = {}
for name, cfg in ARRAYS.items():
    path = OUT / f"features_{cfg['key']}.parquet"
    if not path.exists():
        build(name).to_parquet(path)
    features[name] = pd.read_parquet(path)

# %% [markdown]
# ## Summary per layout

# %%
FEATURES = ["frac_A", "frac_C", "frac_G", "frac_T", "GC", "G_first18", "G_last18", "max_homopolymer", "n_GGG",
            "g4hunter", "n_cycles", "hamming_neighbours", "stray_events", "stray_cycles", "mfe_barcode", "mfe_probe",
            "nn_dist", "n_within_7", "kmer4_mean_log_freq", "kmer4_min_freq"]
summary = pd.concat({name: df[FEATURES].describe().T[["mean", "std", "min", "50%", "max"]]
                     for name, df in features.items()}, axis=1).round(3)
summary.to_csv(OUT / "feature_summary.csv")
summary

# %% [markdown]
# Correlation between features (Spearman), pooled over layouts, to identify redundant features before modelling.

# %%
pooled = pd.concat(features.values())
corr = pooled[FEATURES].corr(method="spearman").round(2)
corr.to_csv(OUT / "feature_correlation.csv")
corr
