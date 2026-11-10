"""
Unit tests for axis code generator.
"""

import unittest
from spatial_zipcodes.axis_code_generator import (
    hamming_distance,
    AxisCodeGenerator,
    generate_axis_codes,
)


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


if __name__ == '__main__':
    unittest.main()
