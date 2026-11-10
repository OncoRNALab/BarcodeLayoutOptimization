"""
Unit tests for codebook generation.
"""

import unittest
from spatial_zipcodes.codebook import build_codebook, Feature


class TestCodebook(unittest.TestCase):
    """Test codebook building and operations."""
    
    def setUp(self):
        """Set up test data."""
        # Create simple test codes
        self.rows = 3
        self.cols = 4
        
        self.x_codes = [
            [0, 0], [0, 1], [1, 0], [1, 1]
        ]
        self.y_codes = [
            [0, 0], [0, 1], [1, 0]
        ]
        
        self.x_sequences = ["AA", "AC", "GA", "GC"]
        self.y_sequences = ["AA", "AC", "GA"]
        
        self.spacer = "GGG"
    
    def test_codebook_creation(self):
        """Test basic codebook creation."""
        codebook = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
        )
        
        # Check total number of features
        self.assertEqual(len(codebook.features), self.rows * self.cols)
    
    def test_feature_indexing(self):
        """Test feature access by row and column."""
        codebook = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
        )
        
        # Get feature at row=1, col=2
        feature = codebook.get_feature(1, 2)
        
        self.assertIsNotNone(feature)
        self.assertEqual(feature.row, 1)
        self.assertEqual(feature.col, 2)
        
        # Check X and Y sequences
        self.assertEqual(feature.x_seq, self.x_sequences[2])
        self.assertEqual(feature.y_seq, self.y_sequences[1])
        
        # Check zip sequence
        expected_zip = self.x_sequences[2] + self.spacer + self.y_sequences[1]
        self.assertEqual(feature.zip_seq, expected_zip)
    
    def test_zip_sequence_structure(self):
        """Test that zip sequences are correctly assembled."""
        codebook = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
        )
        
        for feature in codebook.features:
            expected_zip = feature.x_seq + self.spacer + feature.y_seq
            self.assertEqual(feature.zip_seq, expected_zip)
    
    def test_unique_zip_sequences(self):
        """Test that all zip sequences are unique."""
        codebook = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
        )
        
        zip_sequences = [f.zip_seq for f in codebook.features]
        
        # All sequences should be unique
        self.assertEqual(len(zip_sequences), len(set(zip_sequences)))
    
    def test_feature_id_formats(self):
        """Test different feature ID formats."""
        # Row-col format
        codebook_rc = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
            feature_id_format="row_col",
        )
        
        feature = codebook_rc.get_feature(1, 2)
        self.assertEqual(feature.feature_id, "R001_C002")
        
        # Linear format
        codebook_lin = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
            feature_id_format="linear",
        )
        
        feature = codebook_lin.get_feature(1, 2)
        expected_id = str(1 * self.cols + 2)
        self.assertEqual(feature.feature_id, expected_id)
    
    def test_statistics(self):
        """Test codebook statistics."""
        codebook = build_codebook(
            rows=self.rows,
            cols=self.cols,
            x_codes=self.x_codes,
            y_codes=self.y_codes,
            x_sequences=self.x_sequences,
            y_sequences=self.y_sequences,
            spacer=self.spacer,
        )
        
        stats = codebook.get_statistics()
        
        self.assertEqual(stats['n_features'], self.rows * self.cols)
        self.assertEqual(stats['rows'], self.rows)
        self.assertEqual(stats['cols'], self.cols)
        self.assertEqual(stats['x_bit_length'], 2)
        self.assertEqual(stats['y_bit_length'], 2)
        self.assertEqual(stats['x_dna_length'], 2)
        self.assertEqual(stats['y_dna_length'], 2)
        self.assertEqual(stats['spacer'], self.spacer)
        self.assertEqual(stats['spacer_length'], len(self.spacer))
        self.assertEqual(stats['zip_length'], 2 + len(self.spacer) + 2)


if __name__ == '__main__':
    unittest.main()
