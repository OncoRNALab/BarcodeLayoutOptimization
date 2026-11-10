# Methods (barcode design and decoys)

Excerpt of the Methods of the BMC Bioinformatics manuscript (`phd_thesis_revision/layoutopt_manuscript/manuscript.tex`, read-only copy, Oct 2026), lightly converted from LaTeX. This is the authoritative description of the three printed layouts and their decoys. Citations were removed.

## Array geometry

All layouts were designed for the full-density arrays of our maskless array synthesizer, which consist of $235 × 363 = 85,305$ square features of $13.7 × 13.7$ µm separated by 1 µm gaps, so that the centre-to-centre distance between edge-adjacent features is 14.7 µm and between diagonal neighbours 20.8 µm. Every position of the grid carries a 36-nt barcode.

## Greedy placement and 2-OPT refinement

A random candidate pool of unique 36-nt barcodes was generated as previously described: sequences were drawn uniformly at random and retained only if they contained no homopolymer longer than 3 nt and at most 25% G. The printed library requires on average 2.35 synthesis cycles per nucleotide under the A, C, G, T deposition order (maximum 3.03).

Barcodes were placed with BarLay, the GPU implementation of greedy placement and 2-OPT local search by Jatzkowski et al., with default settings. BarLay encodes each barcode by the left-most embedding of its synthesis schedule under the A, C, G, T deposition order: a binary vector with one entry per cycle that equals 1 if the barcode is extended, and therefore illuminated, in that cycle. It then minimizes the summed Hamming distance between the schedules of all pairs of neighbouring features in the 8-neighbourhood, that is, the number of cycles in which one feature is illuminated while its neighbour is not. The greedy step fills the grid one position at a time and, at each position, evaluates every barcode that has not yet been placed, choosing the one with the smallest summed schedule distance to its already placed neighbours. The 2-OPT step then repeatedly searches for exchanges between two placed barcodes, or between a placed barcode and an unplaced barcode of the pool, that reduce this cost, and applies them until no improving exchange remains. 85,305 of the candidates were placed on the array. Computations were run on NVIDIA A100-SXM4-80GB GPUs at the Vlaams Supercomputer Centrum (VSC).

## Co-designed layouts

**Dual barcodes.** The grid was partitioned into $4 × 4$ blocks, giving $59 × 91 = 5,369$ blocks, with truncated blocks at the lower and right edges. Each 36-nt barcode is the concatenation of an 18-nt block identifier, bar_1, shared by all features of a block, and an 18-nt spot code, bar_2, that distinguishes the features within the block. Both parts were drawn from a single library of 85,305 18-nt sequences, obtained as the 18-nt prefixes of a random 36-nt library filtered as described above. The library was shuffled, and each block received the least-used sequence as its bar_1, so that all 5,369 block identifiers are distinct. Each spot code was then drawn at random from the same library, excluding sequences already paired with the bar_1 of that block, so that every full barcode is unique. Spot codes are reused between blocks: the printed layout contains 53,810 distinct bar_2 sequences, and 3,391 sequences occur both as a block identifier and as a spot code. No composition filter was applied to the concatenated barcodes, which contain up to 42% G and, at the junction between the two parts, homopolymers of up to 6 nt. Because every feature of a block shares the first 18 synthesis steps of its barcode, blocks of 16 features are co-activated for about half of the synthesis.

**Zipcodes.** Zipcodes encode the column and the row of a feature in two separate 18-nt axis codes, concatenated without separator as [column code][row code]. Each axis code is a 36-bit binary codeword converted to DNA at 2 bits per nucleotide (00→A, 01→C, 10→T, 11→G). Codewords were generated one position at a time along each axis by a randomized search with backtracking. For each position, up to 10,000 random 36-bit candidates were tested against two constraints:

- soft targets: the Hamming distance to the codewords 1, 2, 3 and 4 positions back should be 2, 4, 6 and 8 bits, respectively. A candidate was accepted if the mean squared deviation from these targets was at most 16 (a tolerance of 0.5 times the largest target, squared);
- a hard constraint: at least 8 bits from every codeword 5 or more positions back.

If no candidate was accepted, the search backtracked by up to 5 positions. The targets make neighbouring columns, and neighbouring rows, share more of their code than distant ones. Neither co-design scheme uses synthesis schedules or an explicit model of stray light; co-activation of neighbouring features follows from the shared sequence. The targets are not reached exactly: in the printed layout, codewords of adjacent columns differ by 8.8 bits on average (range 4–12), rising to 10.5 bits at an offset of 4 and 17.9 bits beyond. Column and row codes were generated with the same random seed and settings, so the 235 row codes are identical to the first 235 column codes. As a consequence, 235 barcodes are tandem repeats of one axis code, and 27,495 pairs of barcodes are half-swapped versions of each other ($X_aX_b$ and $X_bX_a$). Sequence composition was reported but not enforced during generation; the printed barcodes contain up to 53% G and homopolymers of up to 11 nt.

## Decoy whitelists

Precision cannot be measured directly on real arrays, because the true origin of each read is unknown. We therefore estimated it by adding decoy barcodes, which are not printed, to the whitelist. Reads assigned to decoys are necessarily false, and their frequency indicates how many false assignments land on real barcodes. This estimate is unbiased only if a misassigned read is as likely to land on a decoy as on a real barcode, in proportion to their numbers. For a random library, where misassignments are spread over the whole whitelist, an independent set of random decoys satisfies this condition. In co-designed libraries, however, misassigned reads land mostly on the structural neighbours of the correct barcode, for example another spot of the same block, or the same column code combined with an adjacent row code. Decoys generated independently of the layout lie far from all real barcodes and would therefore collect almost no false assignments. We therefore generated decoys that reproduce the local structure of each layout. For Dual barcodes and Greedy + 2-OPT, the decoy-to-real ratio was 1:2, so decoys make up one third of the whitelist. For Zipcodes, a 1:1 ratio was used so that the local decoy share on each axis is also one third.

**Greedy + 2-OPT.** The 42,652 decoys were drawn from a second independently generated pool of 150,000 barcodes, produced with the same procedure and filters as the printed library and similarly unrelated to array position.

**Dual barcodes.** Each block kept its real block identifier and received additional spot codes drawn from the same 18-nt library, excluding sequences already paired with that identifier. The number of decoys per block was half its number of features (8 for complete blocks), so the local decoy share within every block equals the global share. Candidate decoys identical to a real barcode or to another decoy were rejected, giving 42,652 decoys.

**Zipcodes.** Virtual axis codes were generated in the gaps between adjacent real codes, including the ends of each axis, with Bresenham spacing: 181 virtual column codes and 117 virtual row codes. Each virtual code was obtained by flipping bits in one of the two flanking real codewords. Of 4,000 such candidates, we kept the one whose bit distances to the real codewords 1 to 4 positions away on either side best matched the mean distances observed between real codewords at the same offsets. Candidates were also required to be at least as far from all other codewords, real or virtual, as the closest pair of real codewords. Each virtual column code was combined with all 235 row codes, and each virtual row code with all 363 column codes, giving 85,006 decoys that differ from real barcodes in one axis code only (virtual×virtual combinations were not included). A 1:2 decoy set cannot place one third of the codes on both axes at once; this 1:1 design puts 117 decoys among the 351 whitelist entries that share a real column code and 181 among the 543 that share a real row code, so the local decoy share is $1/3$ on each axis.

For Dual barcodes and Zipcodes, all decoys are unique and none coincides with a real barcode, and their G content matches that of the printed barcodes. The Hamming distance from a decoy to its nearest real barcode matched the distance from a real barcode to its nearest other real barcode (Dual: mean 10.0 for both; Zipcodes: 6.6 and 5.6), whereas decoys generated independently of the layout were much further away (14.5 and 16.0). Among the ten whitelist entries closest to a real barcode, decoys made up 35% for Dual barcodes, close to their global share of 33%, and 28% for Zipcodes, close to the combinatorial local share of 33%.

## Precision estimation and calibration

Let $A$ be the number of reads assigned to any whitelist entry. When the origin of each read is known, the true precision is
$$
P = \frac{C}{A},

$$
where $C$ is the number of reads assigned to their correct barcode. On printed arrays $C$ is unknown. Let $D$ be the number of reads assigned to decoys and $\phi$ the expected fraction of misassigned reads that land on decoys. The quantity $D/A$ is the fraction of assignments that landed on decoys. Dividing by $\phi$ corrects that fraction by the expected decoy share, giving the estimated false-assignment rate among assigned reads. The decoy-based estimator of $P$ is
$$
\hat{P} = 1 - \frac{D/A}{\phi}.

$$
Equivalently, because every decoy assignment is false, the number of false assignments is estimated as $D/\phi$.
If misassignments are spread over the whitelist in proportion to its composition, $\phi$ equals the global decoy share
$$
\phi_{\text{global}} = \frac{n_D}{n_R + n_D},

$$
where $n_R$ and $n_D$ are the numbers of real and decoy barcodes. For the layouts printed here, $\phi_{\text{global}} = 1/3$ for Dual barcodes and Greedy + 2-OPT and $\phi_{\text{global}} \approx 1/2$ for Zipcodes. The same estimator, with $\phi = \phi_{\text{global}}$, was used previously on randomly placed barcodes.

To test this assumption for each layout, we simulated 1,000,000 reads from the printed barcodes of each layout, with no reads from decoys, using the error model of Press as implemented in the RandomBarcodes repository. We used the high-error regime used previously (3% substitutions, 15% insertions and 14% deletions per nucleotide), which reflects the insertion rates expected at full density. Reads were decoded with QUIK against the whitelist of real barcodes plus decoys, at the same thresholds as the sequenced reads. Because the origin of each simulated read is known, the true precision $P$ from the equation above was computed and compared with $\hat{P}$ from the equation above. The number of reads assigned to a wrong real barcode, $F$, is also known, so the effective decoy share can be measured directly as
$$
\phi_{\text{eff}} = \frac{D}{D + F}.

$$
For Dual barcodes and Greedy + 2-OPT, $\phi_{\text{eff}}$ matched $\phi_{\text{global}} = 1/3$, so printed-array precision used that value. For Zipcodes, misassignments stay on one axis (see Results), so $\phi$ is the local decoy share on that axis, which is $1/3$ by construction. Simulation confirmed that this design value recovers $P$ to within 0.6 percentage points, whereas $\phi_{\text{global}} \approx 1/2$ does not; printed Zipcodes precision was therefore estimated with $\phi = 1/3$, not $\phi_{\text{global}}$. Confidence intervals of 95% were computed from the exact Poisson interval for $D$.


## Half-density Dual barcode array

The full-density comparison identified Dual barcodes as the layout with the best balance of recall and precision (see Results). We therefore continued with Dual barcodes only, printed at half density. Alternate features were left empty in a checkerboard, so that each printed feature has four diagonal printed neighbours at 20.8 µm and no edge-adjacent printed neighbour. The 42,471 barcodes from the first 117 rows of the full-density Dual layout ($117 × 363$) were placed so that $x+y$ is even. Each $4 × 4$ block keeps its shared bar_1; on the printer the block occupies an $8 × 4$ window with empty sites in between.

The array was printed and processed with the same photochemistry, probe architecture and library preparation as the full-density arrays, and sequenced on a MiSeq Micro flow cell with the same 47-8-8-47 cycle recipe. Demultiplexing yielded 526,625 reads passing filter. Reads were decoded with QUIK at $r \in \{6, 7, 8\}$, because insertion rates are expected to be lower when printed neighbours are more distant. Decoded reads at $r = 7$ were counted per printed feature and mapped to printer coordinates for spatial heatmaps. Decoys were generated with the same per-block rule as for full-density Dual barcodes, giving 21,236 decoys (8 per complete 16-spot block; $\phi_{\text{global}} = 1/3$). One million reads were simulated from the printed barcodes with the medium-error Press model (4% substitutions, 5% insertions and 14% deletions per nucleotide) and decoded against the combined whitelist at the same thresholds.
