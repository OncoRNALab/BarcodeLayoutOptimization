# Barcode read simulation

Article-style incremental noisy channel with **unequal** error rates. Used to build labeled datasets for decode-precision evaluation of the four NN≈8 layouts.

## Channel

Walk a cursor over the true 36-nt barcode. At each step draw:

| Action | Default probability | Emit | Advance |
|--------|---------------------|------|---------|
| match | 0.77 | reference base | yes |
| substitution | **0.04** | other base | yes |
| deletion | **0.14** | nothing | yes |
| insertion | **0.05** | random base | **no** (stay) |

No pad/truncate — read length is variable. Downstream decoders that need 36 nt should crop from the read.

## Paired FASTQ

- **R1** = garbled forward read  
- **R2** = `reverse_complement(R1)` (same erroneous molecule, both directions)  
- Phred+33 qualities drawn i.i.d. ~ Normal(mean_q=30, sd=3), clipped  

Header contract (same `true_barcode` on both mates):

```text
@read{id} true_barcode={TRUE36} row={r} col={c} mate=1
@read{id} true_barcode={TRUE36} row={r} col={c} mate=2
```

## Generate

```bash
cd simulation
for p in dual_shared_pool_g10_minH6 dual_fixed_bar2_8x8_g10_nn8 zipcodes_flip_g10_nn8 hsdb_g10_nn8_colored; do
  python make_reads.py --preset $p --out-dir output/$p --n-reads 1000000 --seed 1 \
    --srate 0.04 --irate 0.05 --drate 0.14
done
```


Outputs per run: `reads_R1.fastq`, `reads_R2.fastq`, `reads_R1.txt`, `reads_R2.txt`, `truth.csv`, `params.json`.
