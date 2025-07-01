"""
Replicate a single generator around a 1D ring with periodic boundary conditions.

This module can reuse the original implementation since lattice construction
works with Pauli strings and doesn't require polynomial arithmetic.
"""

# Reuse the original implementation
from stabilizer.lattice import build_lattice

__all__ = ['build_lattice'] 