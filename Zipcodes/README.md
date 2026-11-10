> **Note (2026):** the printed Zipcodes are **not** Gray codes. Each axis code is a 36-bit codeword from a randomized search with backtracking, and the distance targets below (2/4/6/8, ≥8) are in **bits**, not bases. In the printed layout the bits map to DNA as 00→A, 01→C, 10→T, 11→G. See `../AGENTS.md` and `../docs/methods.md`.

# Spatial Zipcodes

A Python package for generating spatial barcodes ("zip codes") for photolithography-based spatial transcriptomics platforms. This tool implements a **two independent 1D Gray-code design** where each axis (X and Y) has its own error-tolerant barcode system, optimized for spatial shift tolerance rather than pure error correction.

## Design Philosophy

Unlike traditional error-correcting codes that maximize distances between all barcodes, spatial zipcodes are designed for **spatial shift tolerance**:

- **Nearby positions have similar barcodes**: Positions 1 step apart differ by ~2 bases
- **Distance increases gradually**: Positions 2-4 steps apart differ by 4-8 bases (smooth gradient)
- **Far positions are maximally different**: Positions ≥5 steps apart differ by ≥8 bases

This creates a "spatial Gray code" where sequencing errors that cause misassignment tend to assign to nearby (not distant) positions, preserving spatial structure in the data.

## Algorithm Overview

### Two Independent 1D Designs

Each barcode is composed of two independent components:
- **X-code**: One code per column (determined by column position only)
- **Y-code**: One code per row (determined by row position only)
- **Full barcode**: `X_code[col] + spacer + Y_code[row]`

### Distance Constraints (Per Axis)

For each axis (X and Y independently):
- **Positions 1 step apart** → Hamming distance ≈ 2 bases
- **Positions 2 steps apart** → Hamming distance ≈ 4 bases
- **Positions 3 steps apart** → Hamming distance ≈ 6 bases
- **Positions 4 steps apart** → Hamming distance ≈ 8 bases
- **Positions ≥5 steps apart** → Hamming distance ≥ 8 bases (minimum)

### Emergent 2D Behavior

The 2D spatial properties emerge naturally from combining two 1D designs:
- **Horizontal neighbors** (same row, adjacent columns): differ only in X_seq → ~2 bases
- **Vertical neighbors** (same column, adjacent rows): differ only in Y_seq → ~2 bases
- **Diagonal neighbors**: differ in both X_seq and Y_seq → ~4 bases total

## Features

- **Independent 1D axis designs**: Separate Gray-code-like sequences for X and Y dimensions
- **Configurable chip geometry**: Support for rectangular grids with arbitrary dimensions
- **Spatial Gray code**: Gradual distance increase with position difference (not maximum separation)
- **Randomized greedy algorithm**: Generates axis codes with backtracking to satisfy distance constraints
- **DNA sequence encoding**: Multiple encoding schemes including photolithography-realistic 4-bit one-hot encoding
- **Quality control**: Validates sequences for GC content, homopolymer runs, and problematic motifs
- **Distance analysis**: Computes and visualizes distance distributions to verify constraints
- **Error simulation**: Tests decoding accuracy under various error models
- **Modular design**: Easy to extend and customize for different platforms

## Installation

```bash
# Clone the repository
cd layout_opt

# Install in development mode
pip install -e .
```

### Requirements

- Python >= 3.9
- numpy
- (Optional) matplotlib for visualization

See `requirements.txt` for complete dependencies.

## Quick Start

### 1. Create a configuration file

Create `my_config.json`:

```json
{
  "chip_geometry": {
    "type": "rectangular",
    "rows": 118,
    "cols": 182,
    "feature_size_microns": 14,
    "feature_spacing_microns": 16
  },
  "error_correction": {
    "type": "edit_distance",
    "distances": [2, 4, 6, 8],
    "far_features_min_distance": 8
  },
  "bit_to_dna": {
    "scheme": "2-bits-per-base",
    "mapping": {
      "00": "A",
      "01": "C",
      "10": "G",
      "11": "T"
    }
  },
  "sequence_qc": {
    "gc_content_min": 0.3,
    "gc_content_max": 0.7,
    "max_homopolymer_length": 4,
    "reject_sequences_with_long_homopolymers": true
  },
  "codebook": {
    "spacer_sequence": "GGG",
    "feature_id_format": "row_col",
    "output_format": "csv"
  }
}
```

### 2. Run the pipeline

#### Using the example script:

```bash
python scripts/run_generate_codebook.py
```

#### Using the CLI:

```bash
# Full pipeline in one command
python -m spatial_zipcodes.cli full-pipeline \
  --config my_config.json \
  --output-dir output/

# Or run steps individually:

# Step 1: Generate axis codes
python -m spatial_zipcodes.cli generate-axis-codes \
  --config my_config.json \
  --output output/

# Step 2: Analyze distances
python -m spatial_zipcodes.cli analyze-distances \
  --config my_config.json \
  --codes output/axis_codes.json \
  --output output/

# Step 3: Build codebook
python -m spatial_zipcodes.cli build-codebook \
  --config my_config.json \
  --codes output/axis_codes.json \
  --output output/codebook.csv

# Step 4: Simulate errors
python -m spatial_zipcodes.cli simulate-errors \
  --config my_config.json \
  --codebook output/codebook.csv \
  --n-samples 1000 \
  --output output/simulation.json
```

### 3. Use as a Python library

```python
from spatial_zipcodes import load_config
from spatial_zipcodes.axis_code_generator import generate_axis_codes
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.codebook import build_codebook, export_codebook

# Load configuration
config = load_config("my_config.json")

# Generate X and Y axis codes
x_codes, x_verification = generate_axis_codes(
    n_positions=config.chip_geometry.cols,
    bit_length=config.x_bit_length,
    target_distances=config.error_correction.distances,
    far_features_min_distance=config.error_correction.far_features_min_distance,
    random_seed=42,
)

y_codes, y_verification = generate_axis_codes(
    n_positions=config.chip_geometry.rows,
    bit_length=config.y_bit_length,
    target_distances=config.error_correction.distances,
    far_features_min_distance=config.error_correction.far_features_min_distance,
    random_seed=42,
)

# Convert to DNA sequences
encoder = BitToDNAEncoder(scheme="2-bits-per-base")
x_sequences = [encoder.encode(code) for code in x_codes]
y_sequences = [encoder.encode(code) for code in y_codes]

# Build and export codebook
codebook = build_codebook(
    rows=config.chip_geometry.rows,
    cols=config.chip_geometry.cols,
    x_codes=x_codes,
    y_codes=y_codes,
    x_sequences=x_sequences,
    y_sequences=y_sequences,
    spacer="GGG",
)

export_codebook(codebook, "output/codebook.csv")
```

## Package Structure

```
spatial_zipcodes/
├── __init__.py              # Package initialization
├── config.py                # Configuration dataclasses and loading
├── axis_code_generator.py   # 1D axis code generation with distance constraints
├── bit_to_dna.py           # Bit-to-DNA encoding/decoding
├── sequence_qc.py          # DNA sequence quality control
├── codebook.py             # Codebook assembly and export
├── distance_analysis.py    # Distance distribution analysis
├── simulation.py           # Error simulation and decoding validation
└── cli.py                  # Command-line interface

scripts/
└── run_generate_codebook.py  # End-to-end example script

tests/
├── test_axis_code_generator.py
├── test_bit_to_dna.py
└── test_codebook.py
```

## Algorithm Details

### 1D Axis Code Generation

The package generates **two independent 1D Gray codes** (one for X-axis, one for Y-axis) using a **randomized greedy algorithm with backtracking**:

#### Algorithm Steps:

1. **Sequential Position Processing**: For each position along an axis (0, 1, 2, ..., n-1):
   - Sample random binary codes until one satisfies all constraints
   - Constraints are evaluated based on **position step difference** along that axis only

2. **Hard Constraints** (must be satisfied):
   - Positions ≥5 steps apart → Hamming distance ≥ `far_features_min_distance`
   - Example: column 10 and column 20 (10 steps apart) must differ by ≥8 bases

3. **Soft Constraints** (approximate target):
   - Position difference 1 → target distance = `distances[0]` (typically 2)
   - Position difference 2 → target distance = `distances[1]` (typically 4)
   - Position difference 3 → target distance = `distances[2]` (typically 6)
   - Position difference 4 → target distance = `distances[3]` (typically 8)
   - Example: column 10 and column 11 should differ by ~2 bases

4. **Backtracking**: If no valid code found after many attempts:
   - Revert `backtrack_depth` positions (default: 5)
   - Try alternative codes for those positions

5. **Verification**: Compute full pairwise distance distributions and report statistics

#### Key Properties:

- **1D design**: Distance constraints depend only on position difference along that single axis
- **Spatial Gray code**: Increasing target distances [2, 4, 6, 8] create smooth spatial gradient
- **Independent axes**: X and Y codes generated completely independently

### Bit-to-DNA Encoding

Two encoding schemes are supported:

#### 2-bits-per-base (Standard):
- `00` → A
- `01` → C
- `10` → G
- `11` → T

Preserves Hamming distances: each bit flip = single base change.

#### 4-bits-per-base (Photolithography-realistic, One-hot):
- `1000` → A
- `0100` → C
- `0010` → G
- `0001` → T

Models actual photolithography synthesis where each base requires one cycle. Input codes are 2-bit patterns that get converted to 4-bit one-hot during encoding:
- 36-bit input → 18-base DNA sequence
- 60-bit input → 30-base DNA sequence

### Quality Control

Sequences are validated for:
- **GC content**: Within specified range (default: 30-70%)
- **Homopolymer runs**: Maximum length (default: 4)
- **Problematic motifs**: G-quadruplexes, long runs, etc.

### Error Simulation

Two error models are supported:
- **Substitution**: Per-base substitution probability
- **Indel**: Insertions, deletions, and substitutions

Decoding uses **nearest-neighbor** with configurable maximum edit distance.

### Why This Design?

**Traditional error-correcting codes** maximize minimum distance between all codewords, which is optimal for random errors but not for spatial transcriptomics:
- Random misassignment to any position
- Doesn't preserve spatial structure

**Spatial Gray codes** (this package) minimize distance between nearby positions:
- Errors preferentially assign to nearby positions
- Preserves spatial structure even under errors
- "Graceful degradation": small errors → small spatial shifts

**Example behavior:**
- A barcode at position (50, 90) with sequencing errors is more likely to decode as (50, 91) or (51, 90) than as (10, 10)
- Spatial patterns (expression gradients, tissue domains) are preserved even when individual spots are misassigned

## Configuration Reference

### Chip Geometry
- `type`: Chip type (default: "rectangular")
- `rows`: Number of rows
- `cols`: Number of columns
- `feature_size_microns`: Feature size (for documentation)
- `feature_spacing_microns`: Feature spacing (for documentation)

### Error Correction
- `type`: Error correction type (default: "edit_distance")
- `distances`: Target Hamming distances for position step differences [1, 2, 3, 4] along each axis
  - Example: `[2, 4, 6, 8]` means positions 1 step apart differ by ~2 bases, 2 steps → ~4 bases, etc.
- `far_features_min_distance`: Minimum distance for positions ≥5 steps apart on the same axis

### Bit-to-DNA
- `scheme`: Encoding scheme ("2-bits-per-base" or "4-bits-per-base")
- `mapping`: Custom bit-to-base mapping
  - For "2-bits-per-base": 4 entries (00, 01, 10, 11)
  - For "4-bits-per-base": 4 entries with one-hot patterns (1000, 0100, 0010, 0001)

### Sequence QC
- `gc_content_min`: Minimum GC content (0-1)
- `gc_content_max`: Maximum GC content (0-1)
- `max_homopolymer_length`: Maximum homopolymer run length
- `reject_sequences_with_long_homopolymers`: Reject sequences with long homopolymers

### Codebook
- `spacer_sequence`: Sequence between X and Y segments (default: "GGG")
- `feature_id_format`: "row_col" or "linear"
- `output_format`: "csv" or "tsv"

### Axis Code Algorithm
- `max_attempts_per_position`: Max random samples before backtracking (default: 10000)
- `backtrack_depth`: Positions to backtrack when stuck (default: 5)
- `tolerance`: Tolerance for soft constraints (default: 0.5)
- `random_seed`: Random seed for reproducibility (optional)

### Simulation
- `error_model`: "substitution" or "indel"
- `per_base_error_rate`: Per-base error probability (default: 0.01)
- `distance_metric`: "hamming" or "levenshtein"
- `max_edit_distance_for_decoding`: Maximum distance for successful decoding (default: 3)

## Output Files

### Codebook CSV/TSV
Columns:
- `feature_id`: Unique feature identifier
- `row`: Row position
- `col`: Column position
- `X_bits`: Binary X code
- `Y_bits`: Binary Y code
- `X_seq`: DNA sequence for X
- `Y_seq`: DNA sequence for Y
- `zip_seq`: Full zip code (X + spacer + Y)

### Distance Analysis JSON
Contains distance distributions and constraint satisfaction statistics for X and Y axes.

### Simulation Results JSON
Contains decoding accuracy metrics:
- Overall accuracy
- Near error rate (misassignment to nearby positions)
- Far error rate (misassignment to distant positions)
- Ambiguous decoding rate
- Failed decoding rate

## Running Tests

```bash
# Run all tests
python -m unittest discover tests/

# Run specific test module
python -m unittest tests.test_axis_code_generator

# Run with verbose output
python -m unittest discover tests/ -v
```

## Examples

See `scripts/run_generate_codebook.py` for a complete end-to-end example.

Example configurations are provided in:
- `example_config.json`: Standard configuration
- `example_config_small.json`: Small test configuration

## Citation

This package implements spatial barcoding principles inspired by photolithography-based spatial transcriptomics platforms. If you use this tool, please cite relevant papers in the field.

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Contact

For questions or issues, please open an issue on the GitHub repository.
