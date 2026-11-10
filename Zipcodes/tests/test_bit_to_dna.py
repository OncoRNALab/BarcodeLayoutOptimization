"""
Unit tests for bit-to-DNA encoding.
"""

import unittest
from spatial_zipcodes.bit_to_dna import (
    BitToDNAEncoder,
    encode_codes_to_dna,
    decode_dna_to_codes,
    hamming_distance_dna,
)


class TestBitToDNAEncoder(unittest.TestCase):
    """Test bit-to-DNA encoding and decoding."""
    
    def test_default_2bit_encoding(self):
        """Test default 2-bits-per-base encoding."""
        encoder = BitToDNAEncoder(scheme="2-bits-per-base")
        
        # Test encoding
        bits = [0, 0, 0, 1, 1, 0, 1, 1]
        dna = encoder.encode(bits)
        
        # 00->A, 01->C, 10->G, 11->T
        expected = "ACGT"
        self.assertEqual(dna, expected)
    
    def test_encoding_decoding_invertibility(self):
        """Test that encoding and decoding are inverse operations."""
        encoder = BitToDNAEncoder(scheme="2-bits-per-base")
        
        bits = [0, 1, 1, 0, 1, 1, 0, 0, 0, 1]
        dna = encoder.encode(bits)
        decoded_bits = encoder.decode(dna)
        
        self.assertEqual(bits, decoded_bits)
    
    def test_1bit_encoding(self):
        """Test 1-bit-per-base encoding."""
        encoder = BitToDNAEncoder(scheme="1-bit-per-base")
        
        bits = [0, 1, 0, 1]
        dna = encoder.encode(bits)
        
        # Default: 0->A, 1->T
        expected = "ATAT"
        self.assertEqual(dna, expected)
    
    def test_custom_mapping(self):
        """Test custom bit-to-base mapping."""
        custom_mapping = {
            "00": "T",
            "01": "G",
            "10": "C",
            "11": "A",
        }
        
        encoder = BitToDNAEncoder(scheme="2-bits-per-base", mapping=custom_mapping)
        
        bits = [0, 0, 0, 1, 1, 0, 1, 1]
        dna = encoder.encode(bits)
        expected = "TGCA"
        self.assertEqual(dna, expected)
        
        # Test decoding
        decoded = encoder.decode(dna)
        self.assertEqual(bits, decoded)
    
    def test_invalid_bit_length(self):
        """Test that invalid bit length raises error."""
        encoder = BitToDNAEncoder(scheme="2-bits-per-base")
        
        # Odd number of bits should raise error
        bits = [0, 1, 0]
        with self.assertRaises(ValueError):
            encoder.encode(bits)
    
    def test_batch_encoding(self):
        """Test encoding multiple codes."""
        codes = [
            [0, 0, 1, 1],
            [0, 1, 1, 0],
            [1, 0, 0, 1],
        ]
        
        sequences = encode_codes_to_dna(codes, scheme="2-bits-per-base")
        
        # 00->A, 01->C, 10->G, 11->T
        expected = ["AT", "CG", "GC"]
        self.assertEqual(sequences, expected)
        
        # Test round-trip
        decoded = decode_dna_to_codes(sequences, scheme="2-bits-per-base")
        self.assertEqual(codes, decoded)
    
    def test_dna_length_calculations(self):
        """Test DNA and bit length calculations."""
        encoder = BitToDNAEncoder(scheme="2-bits-per-base")
        
        self.assertEqual(encoder.get_dna_length(8), 4)
        self.assertEqual(encoder.get_dna_length(12), 6)
        self.assertEqual(encoder.get_bit_length(4), 8)
        self.assertEqual(encoder.get_bit_length(6), 12)


class TestHammingDistanceDNA(unittest.TestCase):
    """Test Hamming distance for DNA sequences."""
    
    def test_identical_sequences(self):
        """Hamming distance of identical sequences."""
        seq = "ACGTACGT"
        self.assertEqual(hamming_distance_dna(seq, seq), 0)
    
    def test_different_sequences(self):
        """Hamming distance of different sequences."""
        seq1 = "ACGT"
        seq2 = "TGCA"
        self.assertEqual(hamming_distance_dna(seq1, seq2), 4)
    
    def test_partial_difference(self):
        """Hamming distance with partial differences."""
        seq1 = "ACGTACGT"
        seq2 = "ACGTGCGT"
        self.assertEqual(hamming_distance_dna(seq1, seq2), 1)
    
    def test_different_lengths_raises(self):
        """Different length sequences should raise error."""
        seq1 = "ACGT"
        seq2 = "ACGTAC"
        with self.assertRaises(ValueError):
            hamming_distance_dna(seq1, seq2)


if __name__ == '__main__':
    unittest.main()
