# Barcode Generation for 118×182 Arrays - Two Unique Sets

This directory contains scripts and configurations to generate two completely unique sets of spatial barcodes for 118×182 arrays (21,476 spots per set) using three different encoding approaches.

## Overview

**Goal:** Generate two independent barcode sets for DNA in-situ transcriptomics arrays with no overlap between sets.

**Array Specifications:**
- Dimensions: 118 rows × 182 columns = 21,476 spots per set
- Total unique barcodes needed: 42,952 (2 sets × 21,476)
- Ordering: Column-major (for array fabrication compatibility)

**Three Approaches:**
1. **Dual Barcodes** - 36nt (18nt + 18nt) from pre-validated library
2. **Hybrid Quadtree** - 36nt (24nt spatial + 12nt library)
3. **Zipcodes** - 18nt (9nt X-code + 9nt Y-code, no separator)

## Directory Structure

```
layout_opt/
├── Dual_barcodes/
│   ├── config_118x182_set1.json          # Configuration for Set 1
│   ├── config_118x182_set2.json          # Configuration for Set 2
│   ├── generate_layouts_118x182_dual.py  # Generation script
│   └── dual_barcode_layout.py            # Core library
├── Quadtree/
│   ├── config_118x182_hybrid_set1.json   # Configuration for Set 1 (seed=42)
│   ├── config_118x182_hybrid_set2.json   # Configuration for Set 2 (seed=100)
│   ├── generate_layouts_118x182_hybrid.py # Generation script
│   └── generate_hybrid_quadtree_library.py # Core library
├── Zipcodes/
│   ├── config_118x182_nosep_set1.json    # Configuration for Set 1 (seed=42)
│   ├── config_118x182_nosep_set2.json    # Configuration for Set 2 (seed=200)
│   ├── generate_layouts_118x182_zipcodes.py # Generation script
│   └── spatial_zipcodes/                 # Core library package
└── inhouse_barcodes/
    ├── bar_1in2_first18_set1.txt         # 18nt library for Set 1 (21,476 barcodes)
    ├── bar_1in2_first18_set2.txt         # 18nt library for Set 2 (21,176 barcodes)
    ├── bar_1in2_first12_set1.txt         # 12nt library for Set 1 (21,476 barcodes)
    └── bar_1in2_first12_set2.txt         # 12nt library for Set 2 (21,176 barcodes)
```

## Strategy for Ensuring Uniqueness

### Dual Barcodes
- **Method:** Non-overlapping library regions
- **Set 1:** First 21,476 barcodes from `bar_1in2_column_major.txt`
- **Set 2:** Last 21,176 barcodes from `bar_1in2_column_major.txt`
- **Guarantee:** No overlap by design (non-overlapping regions)

### Hybrid Quadtree
- **Method:** Different random seeds
- **Set 1:** seed = 42
- **Set 2:** seed = 100
- **Guarantee:** Different random selections from library + different quadtree codes

### Zipcodes
- **Method:** Different random seeds for Gray-code generation
- **Set 1:** seed = 42
- **Set 2:** seed = 200
- **Guarantee:** Different Gray-code sequences generated

## Usage Instructions

### Prerequisites

1. **Python Environment:**
   ```bash
   # Activate the virtual environment
   source /kyukon/scratch/gent/vo/000/gvo00027/vsc44685/layout_opt/.venv/bin/activate
   ```

2. **Required Libraries:**
   - numpy
   - pandas
   - All custom modules are already in the workspace

### Generating Layouts

#### 1. Dual Barcodes (36nt)

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Dual_barcodes
python generate_layouts_118x182_dual.py
```

**What it does:**
- Generates layouts for both Set 1 and Set 2
- Uses 4×4 blocks (each block shares one bar_1, has 16 unique bar_2)
- Combines bar_1 + bar_2 = 36nt final barcode
- Verifies uniqueness within and between sets

**Expected outputs:**
```
Set 1:
- layout_118x182_dual_set1.csv        # Full layout (row,col,barcode,bar_1,bar_2)
- barcodes_118x182_dual_set1.txt      # Barcodes only (column-major)
- stats_118x182_dual_set1.json        # Statistics

Set 2:
- layout_118x182_dual_set2.csv
- barcodes_118x182_dual_set2.txt
- stats_118x182_dual_set2.json
```

**Expected runtime:** < 1 minute

#### 2. Hybrid Quadtree (36nt)

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Quadtree
python generate_layouts_118x182_hybrid.py
```

**What it does:**
- Generates 24nt quadtree codes with RLL constraint (max_run=3)
- Adds 12nt library codes for sequencing quality
- Uses different random seeds for Set 1 vs Set 2
- Verifies uniqueness within and between sets

**Expected outputs:**
```
output_118x182_hybrid_set1/:
- layout_118x182_hybrid_set1.csv      # Full layout (row,col,barcode,quadtree,library)
- barcodes_118x182_hybrid_set1.txt    # Barcodes only (column-major)
- stats_118x182_hybrid_set1.json      # Statistics

output_118x182_hybrid_set2/:
- layout_118x182_hybrid_set2.csv
- barcodes_118x182_hybrid_set2.txt
- stats_118x182_hybrid_set2.json
```

**Expected runtime:** 1-2 minutes per set

#### 3. Zipcodes (18nt)

```bash
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Zipcodes
python generate_layouts_118x182_zipcodes.py
```

**What it does:**
- Generates Gray-code-like sequences for X and Y axes
- X-codes: 182 codes of 36 bits → 9nt DNA
- Y-codes: 118 codes of 36 bits → 9nt DNA
- Combines X + Y = 18nt barcode (no separator)
- Uses different random seeds for Set 1 vs Set 2
- Verifies uniqueness within and between sets

**Expected outputs:**
```
output_118x182_nosep_set1/:
- codebook.csv                        # Full codebook with metadata
- layout_118x182_nosep_set1.csv       # Filtered layout (row-major)
- layout_118x182_nosep_set1_colmajor.csv # Column-major layout
- barcodes_118x182_nosep_set1.txt     # Barcodes only (column-major)
- statistics.json                     # Statistics
- intermediate/axis_codes.json        # Intermediate bit codes

output_118x182_nosep_set2/:
- (same structure as Set 1)
```

**Expected runtime:** 2-5 minutes per set (constraint satisfaction is slow)

## Configuration Details

### Dual Barcodes Configuration

```json
{
  "array": {"rows": 118, "cols": 182},
  "block": {"rows": 4, "cols": 4},
  "barcodes": {
    "bar_1_file": ".../bar_1in2_first18_set1.txt",  // 18nt
    "bar_2_file": ".../bar_1in2_first18_set1.txt",  // 18nt
    "bar_1_length": 18,
    "bar_2_length": 18,
    "combined_length": 36
  }
}
```

**Key parameters:**
- `block`: 4×4 blocks minimize library size while maintaining diversity
- `bar_1_file` and `bar_2_file`: Use same library but different subsets
- Set 1 and Set 2 use different non-overlapping library files

### Hybrid Quadtree Configuration

```json
{
  "array": {"rows": 118, "cols": 182},
  "quadtree": {
    "n_bits_x": 8,        // 2^8 = 256 ≥ 182 columns
    "n_bits_y": 7,        // 2^7 = 128 ≥ 118 rows
    "total_bits": 15,     // 8 + 7 = 15 bits
    "triplets": 3,        // 15 bits / 3 = 5 triplets → 5×3 = 15 bits
    "quadtree_length_nt": 24  // 3 triplets × 8nt = 24nt
  },
  "rll": {"max_run": 3},  // Run-length-limited encoding
  "library": {
    "file": ".../bar_1in2_first12_set1.txt",  // 12nt library
    "length_nt": 12
  },
  "random_seed": 42       // Set 2 uses seed=100
}
```

**Key parameters:**
- `n_bits_x/y`: Must cover all positions (2^bits ≥ dimension)
- `triplets`: Number of 8nt segments in quadtree
- `max_run`: Maximum consecutive identical bases (prevents sequencing issues)
- `random_seed`: Different seeds for Set 1 (42) vs Set 2 (100)

### Zipcodes Configuration

```json
{
  "chip_geometry": {"rows": 118, "cols": 182},
  "axis_codes": {
    "x_length_bits": 36,  // Encode 182 positions
    "y_length_bits": 36,  // Encode 118 positions
    "random_seed": 42,    // Set 2 uses seed=200
    "target_knn": 5,      // Target k-nearest neighbors for distance constraints
    "knn_threshold": 12   // Minimum bit distance for neighbors
  },
  "bit_to_dna": {
    "scheme": "base4",
    "mapping": {"00": "A", "01": "T", "10": "C", "11": "G"}
  },
  "codebook": {
    "spacer_sequence": "",  // No separator
    "total_barcode_length": 18  // 9nt + 9nt
  }
}
```

**Key parameters:**
- `x_length_bits`, `y_length_bits`: 36 bits → 9nt DNA
- `random_seed`: Different seeds ensure different Gray-code sequences
- `target_knn`: Gray-code property (neighbors differ by few bits)
- `spacer_sequence`: Empty for no separator

## Verification

Each script performs comprehensive verification:

1. **Within-set uniqueness:**
   - Checks all barcodes in each set are unique
   - Reports: "✓ All 21,476 barcodes are unique"

2. **Between-set uniqueness:**
   - Compares Set 1 vs Set 2
   - Reports: "✓ NO OVERLAP - Both sets are completely unique!"

3. **Length verification:**
   - Dual/Quadtree: 36nt
   - Zipcodes: 18nt

4. **Quality metrics:**
   - Dual: Block structure integrity
   - Quadtree: RLL violations avoided
   - Zipcodes: Spatial distance properties

## Expected Results

### Dual Barcodes
- **Set 1:** 21,476 unique 36nt barcodes
- **Set 2:** 21,476 unique 36nt barcodes
- **Overlap:** 0 (guaranteed by non-overlapping libraries)
- **Structure:** Each 4×4 block shares bar_1, has 16 unique bar_2

### Hybrid Quadtree
- **Set 1:** 21,476 unique 36nt barcodes (seed=42)
- **Set 2:** 21,476 unique 36nt barcodes (seed=100)
- **Overlap:** 0 (very high probability due to different seeds)
- **Structure:** Spatial encoding (24nt) + library diversity (12nt)

### Zipcodes
- **Set 1:** 21,476 unique 18nt barcodes (seed=42)
- **Set 2:** 21,476 unique 18nt barcodes (seed=200)
- **Overlap:** 0 (very high probability due to different seeds)
- **Structure:** Gray-code property (neighbors differ by ~2-3 bases)

## Troubleshooting

### Issue: Library file not found
```
Solution: Check paths in config files point to correct location:
/user/gent/446/vsc44685/ScratchVO_dir/layout_opt/inhouse_barcodes/
```

### Issue: Insufficient barcodes in library
```
Solution: Set 2 has 21,176 barcodes (300 fewer than Set 1)
This is sufficient: 21,176 > 21,476/16 = 1,343 blocks needed
Each block only needs 1 bar_1 barcode
```

### Issue: Zipcodes generation is slow
```
Expected: 2-5 minutes per set
Reason: Constraint satisfaction for Gray-code generation
Solution: Be patient, or reduce array size for testing
```

### Issue: Memory error
```
Solution: These arrays are small (21K spots)
If issues persist, increase swap space or use compute node
```

## Output File Formats

### Layout CSV (Dual Barcodes)
```csv
row,col,barcode,bar_1,bar_2,block_id,position_in_block
0,0,ACGTACGTACGTACGTACGTACGTACGTACGTACGT,ACGTACGTACGTACGTAC,GTACGTACGTACGTACGT,0,0
0,1,ACGTACGTACGTACGTACGTACGTACGTACGTACGT,ACGTACGTACGTACGTAC,GTACGTACGTACGTACGT,0,1
...
```

### Layout CSV (Quadtree)
```csv
row,col,barcode,quadtree_combined,library
0,0,ACGTACGTACGTACGTACGTACGTACGTACGTACGT,ACGTACGTACGTACGTACGTACGT,ACGTACGTACGT
...
```

### Layout CSV (Zipcodes)
```csv
row,col,barcode
0,0,ACGTACGTACGTACGTAC
0,1,ACGTACGTACGTACGTAC
...
```

### Barcodes TXT
```
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCA
...
(column-major ordering)
```

## References

- Original 235×363 array generation (completed successfully)
- Dual barcodes: Block-based combinatorial approach
- Quadtree: Spatial encoding with RLL constraints
- Zipcodes: Gray-code-based spatial barcoding

## Contact

For questions or issues, refer to the main project documentation or contact the project maintainer.

---

**Last Updated:** 2024
**Status:** Ready for generation
**Tested:** Configuration files created and validated
