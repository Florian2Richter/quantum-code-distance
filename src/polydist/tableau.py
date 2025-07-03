"""
Polynomial tableau operations for stabilizer codes.

This module implements tableau construction and linear algebra operations
using Laurent polynomials over F₂[x,x⁻¹].
"""

import numpy as np
from typing import List
from .polynomial import LaurentPolynomial


def build_tableau_poly(stab_ops: List[List[str]]):
    """
    Build polynomial stabilizer tableau from Pauli operators.
    
    Args:
        stab_ops: List of stabilizer operators, each a list of Pauli strings
        
    Returns:
        Polynomial tableau (matrix of Laurent polynomials)
    """
    from .polynomial import pauli_to_polynomial_vector
    
    # Convert each stabilizer operator to polynomial symplectic vector
    poly_tableau = []
    for stab_op in stab_ops:
        poly_vec = pauli_to_polynomial_vector(stab_op)
        poly_tableau.append(poly_vec)
    
    return poly_tableau


def compute_rank_poly(tableau) -> int:
    """
    Compute rank of polynomial tableau.
    
    Args:
        tableau: Polynomial stabilizer tableau
        
    Returns:
        Rank of the tableau
    """
    from .polynomial import rank_poly
    return rank_poly(tableau)


# For backward compatibility, also provide the function signatures that match stabilizer.tableau
def build_tableau(stab_ops: List[List[str]]):
    """Wrapper for build_tableau_poly to match original signature."""
    return build_tableau_poly(stab_ops)


def compute_rank(tableau) -> int:
    """Wrapper for compute_rank_poly to match original signature.""" 
    return compute_rank_poly(tableau) 