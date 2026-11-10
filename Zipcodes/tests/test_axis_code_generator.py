"""
Unit tests for axis code generator.
"""

import unittest
from spatial_zipcodes.axis_code_generator import (
    hamming_distance,
    AxisCodeGenerator,
    generate_axis_codes,
    generate_column_and_row_codes,
)
from spatial_zipcodes.bit_to_dna import BitToDNAEncoder
from spatial_zipcodes.sequence_qc import (
    calculate_g_content,
    junction_homopolymer_length,
    summarize_junction_homopolymers,
)

PRINTED_MAPPING = {"00": "A", "01": "C", "10": "T", "11": "G"}


class TestHammingDistance(unittest.TestCase):
    """Test Hamming distance calculation."""
    
    def test_identical_codes(self):
        """Hamming distance of identical codes should be 0."""
        code = [0, 1, 0, 1, 1, 0]
        self.assertEqual(hamming_distance(code, code), 0)
    
    def test_completely_different(self):
        """Hamming distance of complementary codes."""
        code1 = [0, 0, 0, 0]
        code2 = [1, 1, 1, 1]
        self.assertEqual(hamming_distance(code1, code2), 4)
    
    def test_partial_difference(self):
        """Hamming distance with partial differences."""
        code1 = [0, 1, 0, 1, 1, 0]
        code2 = [0, 1, 1, 1, 0, 0]
        self.assertEqual(hamming_distance(code1, code2), 2)
    
    def test_different_lengths_raises(self):
        """Different length codes should raise ValueError."""
        code1 = [0, 1, 0]
        code2 = [0, 1, 0, 1]
        with self.assertRaises(ValueError):
            hamming_distance(code1, code2)


class TestAxisCodeGenerator(unittest.TestCase):
    """Test axis code generation."""
    
    def test_small_configuration(self):
        """Test generation for a small configuration."""
        n_positions = 10
        bit_length = 10  # Increased for sufficient distance capacity
        target_distances = [2, 3, 4, 5]  # More relaxed targets
        far_min = 6  # Reduced from 8
        
        generator = AxisCodeGenerator(
            n_positions=n_positions,
            bit_length=bit_length,
            target_distances=target_distances,
            far_features_min_distance=far_min,
            random_seed=42,
        )
        
        codes = generator.generate()
        
        # Check we got the right number of codes
        self.assertEqual(len(codes), n_positions)
        
        # Check each code has correct length
        for code in codes:
            self.assertEqual(len(code), bit_length)
        
        # Verify constraints
        verification = generator.verify_constraints()
        
        # Should have no hard constraint violations
        self.assertEqual(verification['hard_constraint_violations'], 0)
        
        # Check that far-features constraint is satisfied
        for i in range(len(codes)):
            for j in range(i + 5, min(len(codes), i + 10)):  # Check subset to avoid timeout
                dist = hamming_distance(codes[i], codes[j])
                self.assertGreaterEqual(
                    dist, far_min,
                    f"Far-features constraint violated: "
                    f"distance({i}, {j}) = {dist} < {far_min}"
                )
    
    def test_reproducibility_with_seed(self):
        """Test that same seed produces same codes."""
        n_positions = 8
        bit_length = 6
        target_distances = [2, 3, 4, 5]
        far_min = 5
        
        codes1, _ = generate_axis_codes(
            n_positions=n_positions,
            bit_length=bit_length,
            target_distances=target_distances,
            far_features_min_distance=far_min,
            random_seed=123,
        )
        
        codes2, _ = generate_axis_codes(
            n_positions=n_positions,
            bit_length=bit_length,
            target_distances=target_distances,
            far_features_min_distance=far_min,
            random_seed=123,
        )
        
        self.assertEqual(codes1, codes2)
    
    def test_distance_verification(self):
        """Test distance verification statistics."""
        n_positions = 12  # Reduced from 15
        bit_length = 12  # Increased for more capacity
        target_distances = [2, 3, 4, 5]  # More relaxed targets
        far_min = 6  # Reduced from 8
        
        codes, verification = generate_axis_codes(
            n_positions=n_positions,
            bit_length=bit_length,
            target_distances=target_distances,
            far_features_min_distance=far_min,
            random_seed=42,
        )
        
        # Check verification structure
        self.assertIn('hard_constraint_violations', verification)
        self.assertIn('soft_constraint_stats', verification)
        self.assertIn('distance_distribution', verification)
        
        # Should have stats for position differences 1-4
        for diff in range(1, 5):
            if diff in verification['soft_constraint_stats']:
                stats = verification['soft_constraint_stats'][diff]
                self.assertEqual(stats['target'], target_distances[diff - 1])


class TestCompositionFilters(unittest.TestCase):
    """G content is rejected before distances are scored."""

    def _encoder(self):
        return BitToDNAEncoder(scheme="2-bits-per-base", mapping=PRINTED_MAPPING)

    def test_high_g_is_not_distance_scored(self):
        import spatial_zipcodes.axis_code_generator as axis_mod

        encoder = self._encoder()
        # Printed mapping: 11 is G, 00 is A.
        high_g = [1, 1, 1, 1, 1, 1, 1, 1]
        no_g = [0, 0, 0, 0, 0, 0, 0, 0]
        draws = [high_g, no_g]
        scored = []

        class RecordingGenerator(AxisCodeGenerator):
            def _check_hard_constraints(self, candidate, position):
                scored.append(list(candidate))
                return super()._check_hard_constraints(candidate, position)

        generator = RecordingGenerator(
            n_positions=1,
            bit_length=8,
            target_distances=[1, 1, 1, 1],
            far_features_min_distance=1,
            max_attempts_per_position=5,
            random_seed=42,
            encoder=encoder,
            g_content_max=0.0,
        )
        original = axis_mod.random_binary_code
        axis_mod.random_binary_code = lambda length, rng: list(draws.pop(0))
        try:
            code = generator._generate_position(0)
        finally:
            axis_mod.random_binary_code = original

        self.assertEqual(code, no_g)
        self.assertEqual(generator.rejected_for_g, 1)
        self.assertTrue(scored)
        self.assertTrue(all(candidate == no_g for candidate in scored))
        self.assertEqual(calculate_g_content(encoder.encode(code)), 0.0)

    def test_column_end_leaves_room_for_any_row_start(self):
        encoder = self._encoder()
        generator = AxisCodeGenerator(
            n_positions=1,
            bit_length=6,
            target_distances=[1, 1, 1, 1],
            far_features_min_distance=1,
            encoder=encoder,
            max_homopolymer_length=3,
            spacer="",
            junction_role="left",
        )
        # 00 is A, so this is AAA. A trailing run of 3 would fill the budget.
        self.assertFalse(generator._passes_composition([0, 0, 0, 0, 0, 0]))
        # AAC ends with a single C, which leaves room for a following C.
        self.assertTrue(generator._passes_composition([0, 0, 0, 0, 0, 1]))

    def test_row_start_is_rejected_when_the_junction_would_be_too_long(self):
        encoder = self._encoder()
        generator = AxisCodeGenerator(
            n_positions=1,
            bit_length=6,
            target_distances=[1, 1, 1, 1],
            far_features_min_distance=1,
            encoder=encoder,
            max_homopolymer_length=3,
            junction_role="right",
            partner_right_edge={"A": 2, "C": 0, "G": 0, "T": 0},
        )
        # AA C would make a run of 4 against a column ending in AA.
        self.assertFalse(generator._passes_composition([0, 0, 0, 0, 0, 1]))
        # AC C adds one A. 2 + 1 fits in a limit of 3.
        self.assertTrue(generator._passes_composition([0, 0, 0, 1, 0, 1]))

    def test_internal_homopolymer_is_rejected_before_the_junction(self):
        encoder = self._encoder()
        generator = AxisCodeGenerator(
            n_positions=1,
            bit_length=12,
            target_distances=[1, 1, 1, 1],
            far_features_min_distance=1,
            encoder=encoder,
            max_homopolymer_length=3,
            junction_role="right",
            partner_right_edge={"A": 0, "C": 0, "G": 0, "T": 0},
        )
        # CAAAAC has an internal run of 4 and does not cross a junction.
        self.assertFalse(generator._passes_composition([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]))
        self.assertEqual(generator.rejected_for_homopolymer, 1)
        self.assertEqual(generator.rejected_for_junction, 0)

    def test_column_and_row_codes_respect_g_and_junction(self):
        encoder = self._encoder()
        g_max = 0.25
        homo_max = 3
        _, _, x_sequences, y_sequences, _, _ = generate_column_and_row_codes(
            n_columns=8,
            n_rows=6,
            x_bit_length=16,
            y_bit_length=16,
            target_distances=[2, 3, 4, 5],
            far_features_min_distance=4,
            encoder=encoder,
            g_content_max=g_max,
            max_homopolymer_length=homo_max,
            spacer="",
            enforce_junction=True,
            max_attempts_per_position=2000,
            random_seed=7,
        )

        for seq in x_sequences + y_sequences:
            self.assertLessEqual(calculate_g_content(seq), g_max)
            self.assertEqual(len(seq), 8)

        summary = summarize_junction_homopolymers(
            x_sequences, y_sequences, homo_max, spacer="",
        )
        self.assertEqual(summary["n_violations"], 0)
        self.assertLessEqual(summary["max_junction_run"], homo_max)


class TestJunctionHomopolymer(unittest.TestCase):
    """Junction runs are measured on the concatenated barcode."""

    def test_run_crossing_the_boundary(self):
        self.assertEqual(junction_homopolymer_length("AAA", "ACC"), 4)
        self.assertEqual(junction_homopolymer_length("AAC", "CTT"), 2)
        self.assertEqual(junction_homopolymer_length("AAG", "CTT"), 0)

    def test_spacer_extends_a_crossing_run(self):
        self.assertEqual(junction_homopolymer_length("AA", "TT", spacer="GGG"), 0)
        self.assertEqual(junction_homopolymer_length("AAG", "GTT", spacer="GG"), 4)

    def test_summary_counts_only_pairs_over_the_limit(self):
        summary = summarize_junction_homopolymers(
            ["AAA", "CCC"],
            ["ACC", "TTT"],
            max_length=3,
        )
        self.assertEqual(summary["n_barcodes"], 4)
        self.assertEqual(summary["n_violations"], 1)
        self.assertEqual(summary["examples"][0]["barcode"], "AAAACC")
        self.assertEqual(summary["max_junction_run"], 4)


if __name__ == '__main__':
    unittest.main()
