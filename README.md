# barcode_design — `manuscript` branch

**Frozen snapshot** for the BMC Bioinformatics paper on spatial barcode layout optimization.

This branch contains the printed layouts, structure-matched decoys, and generator configs needed to reproduce all paper results. It does **not** contain the NN≈8 redesigns or simulation code — see the `main` branch for those.

For ongoing design work, switch to `main`:

```bash
git checkout main
```

## What is here

| Item | Path |
|------|------|
| Printed Dual layout | `data/1in1/layout_235x363_dual_4x4.csv` |
| Printed Zipcodes layout | `data/1in1/layout_235x363_zipcodes_nosep.csv` |
| Printed Greedy layout | `data/1in1/layout_235x363_randombarcodes_opt_coords.csv` |
| Dual half-density | `data/1in2/layout_1in2_dual_4x4.csv` |
| Manuscript decoys | `decoys/output/dual/`, `dual_1in2/`, `zipcodes_balanced_khalf/`, `greedy/` |
| Authoritative methods | `docs/methods.md` |
| Regenerate instructions | `REPRODUCE.md` |
| Decode FASTQs with QUIK | `docs/decode/README.md` |

## Reproduce

See [`REPRODUCE.md`](REPRODUCE.md).

### Decode reads (QUIK only)

Raw and simulated FASTQs are decoded with QUIK through the companion pipeline
[BarcodeCalling_benchmark](https://github.com/OncoRNALab/BarcodeCalling_benchmark)
(`--tool quik`). Step-by-step instructions and example params:

→ [`docs/decode/README.md`](docs/decode/README.md)

## Large files (FASTQs)

Real and simulated FASTQs are on Zenodo
([10.5281/zenodo.24366834](https://doi.org/10.5281/zenodo.24366834)).
See [`docs/ZENODO.md`](docs/ZENODO.md).
