"""
Parse and validate a seed string of Pauli operators.

This module can reuse the original implementation since seed parsing
doesn't require polynomial arithmetic.
"""

# Reuse the original implementation
from stabilizer.generator import parse_seed, VALID_PAULI

__all__ = ['parse_seed', 'VALID_PAULI'] 