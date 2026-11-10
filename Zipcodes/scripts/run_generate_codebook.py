#!/usr/bin/env python3
"""
Example script demonstrating end-to-end spatial zipcode generation.

This script shows how to:
1. Load configuration
2. Generate X and Y axis codes
3. Convert to DNA sequences
4. Perform QC checks
5. Build and export codebook
6. Analyze distance distributions
7. Optionally simulate errors
"""

import logging
from pathlib import Path
import json

from spatial_zipcodes.config import load_config
from spatial_zipcodes.axis_code_generator import generate_axis_codes
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.sequence_qc import validate_axis_sequences
from spatial_zipcodes.codebook import build_codebook, export_codebook
from spatial_zipcodes.distance_analysis import analyze_axis_codes
from spatial_zipcodes.simulation import (
    SubstitutionErrorModel,
    ZipcodeDecoder,
    ZipcodeSimulator,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)


def main():
    """Main execution function."""
    
    # Configuration
    config_path = Path("example_config.json")
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)
    
    logger.info("=" * 70)
    logger.info("Spatial Zipcode Generation Pipeline")
    logger.info("=" * 70)
    
    # Step 1: Load configuration
    logger.info("\n[Step 1] Loading configuration...")
    config = load_config(config_path)
    
    logger.info(f"  Chip: {config.chip_geometry.rows}x{config.chip_geometry.cols}")
    logger.info(f"  X bit length: {config.x_bit_length}")
    logger.info(f"  Y bit length: {config.y_bit_length}")
    logger.info(f"  Target distances: {config.error_correction.distances}")
    logger.info(f"  Far features min distance: {config.error_correction.far_features_min_distance}")
    
    # Step 2: Generate X axis codes
    logger.info("\n[Step 2] Generating X axis codes...")
    x_codes, x_verification = generate_axis_codes(
        n_positions=config.chip_geometry.cols,
        bit_length=config.x_bit_length,
        target_distances=config.error_correction.distances,
        far_features_min_distance=config.error_correction.far_features_min_distance,
        max_attempts_per_position=config.axis_code.max_attempts_per_position,
        backtrack_depth=config.axis_code.backtrack_depth,
        tolerance=config.axis_code.tolerance,
        random_seed=config.axis_code.random_seed,
        candidate_mode=config.axis_code.candidate_mode,
        flip_counts=tuple(config.axis_code.flip_counts),
        require_tolerance=config.axis_code.require_tolerance,
    )
    
    logger.info(f"  Generated {len(x_codes)} X codes")
    logger.info(f"  Hard constraint violations: {x_verification['hard_constraint_violations']}")
    logger.info(f"  Positions outside tolerance: {x_verification['positions_outside_tolerance']}")
    
    # Step 3: Generate Y axis codes
    logger.info("\n[Step 3] Generating Y axis codes...")
    y_codes, y_verification = generate_axis_codes(
        n_positions=config.chip_geometry.rows,
        bit_length=config.y_bit_length,
        target_distances=config.error_correction.distances,
        far_features_min_distance=config.error_correction.far_features_min_distance,
        max_attempts_per_position=config.axis_code.max_attempts_per_position,
        backtrack_depth=config.axis_code.backtrack_depth,
        tolerance=config.axis_code.tolerance,
        random_seed=config.axis_code.random_seed,
        candidate_mode=config.axis_code.candidate_mode,
        flip_counts=tuple(config.axis_code.flip_counts),
        require_tolerance=config.axis_code.require_tolerance,
    )
    
    logger.info(f"  Generated {len(y_codes)} Y codes")
    logger.info(f"  Hard constraint violations: {y_verification['hard_constraint_violations']}")
    logger.info(f"  Positions outside tolerance: {y_verification['positions_outside_tolerance']}")
    
    # Step 4: Convert to DNA sequences
    logger.info("\n[Step 4] Converting codes to DNA sequences...")
    encoder = BitToDNAEncoder(
        scheme=config.bit_to_dna.scheme,
        mapping=config.bit_to_dna.mapping,
    )
    
    x_sequences = [encoder.encode(code) for code in x_codes]
    y_sequences = [encoder.encode(code) for code in y_codes]
    
    logger.info(f"  X sequences: {len(x_sequences[0])} bp each")
    logger.info(f"  Y sequences: {len(y_sequences[0])} bp each")
    
    # Step 5: Perform QC checks
    logger.info("\n[Step 5] Performing sequence QC...")
    from dataclasses import asdict
    qc_config = asdict(config.sequence_qc)
    
    x_passed, x_qc = validate_axis_sequences(x_sequences, "X", qc_config)
    y_passed, y_qc = validate_axis_sequences(y_sequences, "Y", qc_config)
    
    logger.info(f"  X QC: {x_qc['n_passed']}/{x_qc['n_sequences']} passed")
    logger.info(f"  Y QC: {y_qc['n_passed']}/{y_qc['n_sequences']} passed")
    
    # Step 6: Build codebook
    logger.info("\n[Step 6] Building codebook...")
    codebook = build_codebook(
        rows=config.chip_geometry.rows,
        cols=config.chip_geometry.cols,
        x_codes=x_codes,
        y_codes=y_codes,
        x_sequences=x_sequences,
        y_sequences=y_sequences,
        spacer=config.codebook.spacer_sequence,
        feature_id_format=config.codebook.feature_id_format,
    )
    
    stats = codebook.get_statistics()
    logger.info(f"  Total features: {stats['n_features']}")
    logger.info(f"  Zip code length: {stats['zip_length']} bp")
    
    # Step 7: Export codebook
    logger.info("\n[Step 7] Exporting codebook...")
    codebook_path = output_dir / "codebook.csv"
    export_codebook(codebook, codebook_path, file_format="csv")
    logger.info(f"  Codebook exported to {codebook_path}")
    
    # Step 8: Analyze distance distributions
    logger.info("\n[Step 8] Analyzing distance distributions...")
    x_dist_results, y_dist_results = analyze_axis_codes(
        x_codes=x_codes,
        y_codes=y_codes,
        target_distances=config.error_correction.distances,
        far_features_min_distance=config.error_correction.far_features_min_distance,
    )
    
    # Save analysis results
    analysis_path = output_dir / "distance_analysis.json"
    with open(analysis_path, 'w') as f:
        json.dump({
            'x_analysis': x_dist_results,
            'y_analysis': y_dist_results,
        }, f, indent=2, default=str)
    logger.info(f"  Distance analysis saved to {analysis_path}")
    
    # Step 9: Simulate errors (optional)
    logger.info("\n[Step 9] Simulating sequencing errors...")
    
    error_model = SubstitutionErrorModel(
        error_rate=config.simulation.per_base_error_rate,
        random_seed=42,
    )
    
    decoder = ZipcodeDecoder(
        codebook=codebook,
        distance_metric=config.simulation.distance_metric,
        max_edit_distance=config.simulation.max_edit_distance_for_decoding,
        spacer=config.codebook.spacer_sequence,
    )
    
    simulator = ZipcodeSimulator(
        codebook=codebook,
        error_model=error_model,
        decoder=decoder,
        random_seed=42,
    )
    
    sim_results = simulator.simulate_batch(n_samples=500, sample_all=False)
    simulator.print_summary(sim_results)
    
    # Save simulation results
    sim_path = output_dir / "simulation_results.json"
    with open(sim_path, 'w') as f:
        json.dump({
            k: v for k, v in sim_results.items()
            if k != 'individual_results'
        }, f, indent=2)
    logger.info(f"  Simulation results saved to {sim_path}")
    
    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("Pipeline Complete!")
    logger.info("=" * 70)
    logger.info(f"Outputs:")
    logger.info(f"  - Codebook: {codebook_path}")
    logger.info(f"  - Distance analysis: {analysis_path}")
    logger.info(f"  - Simulation results: {sim_path}")
    logger.info(f"\nKey metrics:")
    logger.info(f"  - Total features: {stats['n_features']}")
    logger.info(f"  - Zip code length: {stats['zip_length']} bp")
    logger.info(f"  - Decoding accuracy: {sim_results['accuracy']:.2%}")
    logger.info(f"  - Near error rate: {sim_results['near_error_rate']:.2%}")
    logger.info(f"  - Far error rate: {sim_results['far_error_rate']:.2%}")


if __name__ == "__main__":
    main()
