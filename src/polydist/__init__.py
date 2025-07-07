"""
Polynomial distance computation for stabilizer codes.

This package implements stabilizer code analysis using Laurent polynomials
over F₂[x,x⁻¹] modulo (x^N - 1).
"""

from .orthogonality import check_seed_orthogonality
from .gcd import compute_logical_qubits

__all__ = ['check_seed_orthogonality', 'compute_logical_qubits'] 