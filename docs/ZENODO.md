# Zenodo deposit

Large files from this repository are deposited on Zenodo because they are too big for git
(FASTQs ~150 MB each; full set ~1.2 GB).

**DOI:** `TBD` *(update after upload)*

---

## What is deposited

### Simulated FASTQ reads (four NN≈8 arrays, `main`)

Unequal walk channel: substitution 0.04, insertion 0.05, deletion 0.14 per nucleotide.
1,000,000 paired-end reads per array. R1 = garbled barcode read; R2 = reverse_complement(R1).
Headers include `true_barcode=` for precision scoring.

| Preset | Files |
|--------|-------|
| `dual_shared_pool_g10_minH6` | `reads_R1.fastq`, `reads_R2.fastq`, `truth.csv`, `params.json` |
| `dual_fixed_bar2_8x8_g10_nn8` | same |
| `zipcodes_flip_g10_nn8` | same |
| `hsdb_g10_nn8_colored` | same |

Expected directory layout on Zenodo:

```
simulation/output/
  dual_shared_pool_g10_minH6/
    reads_R1.fastq
    reads_R2.fastq
    truth.csv
    params.json
  dual_fixed_bar2_8x8_g10_nn8/
    ...
  zipcodes_flip_g10_nn8/
    ...
  hsdb_g10_nn8_colored/
    ...
```

### Regenerate locally

```bash
cd simulation
for p in dual_shared_pool_g10_minH6 dual_fixed_bar2_8x8_g10_nn8 zipcodes_flip_g10_nn8 hsdb_g10_nn8_colored; do
  python make_reads.py --preset $p --out-dir output/$p --n-reads 1000000 --seed 1 \
    --srate 0.04 --irate 0.05 --drate 0.14
done
```

Outputs land in `simulation/output/` (gitignored).

---

## Checklist before upload

- [ ] All four preset outputs present under `simulation/output/`
- [ ] Each folder contains `reads_R1.fastq`, `reads_R2.fastq`, `truth.csv`, `params.json`
- [ ] `params.json` records seed, rates, layout path, n_reads
- [ ] MD5 sums computed and included in Zenodo upload metadata
- [ ] DOI updated in this file and in both branch READMEs / `REPRODUCE.md`
