"""
Polynomial distance computation for stabilizer codes.

This module implements logical operator extraction and distance computation
using Laurent polynomials over F₂[x,x⁻¹].
"""

import numpy as np
from typing import List, Tuple, Union
from .polynomial import LaurentPolynomial


def find_logical_operators(tableau) -> List:
    """
    Compute a basis (2k vectors) for the logical operators modulo the stabilizer.
    
    Uses polynomial nullspace computation to find logical operators.
    
    Args:
        tableau: Polynomial stabilizer tableau
        
    Returns:
        List of polynomial logical operator vectors
    """
    from .polynomial import (build_symplectic_form_poly, matrix_multiply_poly, 
                           nullspace_poly, rank_poly)
    
    if not tableau or not tableau[0]:
        return []
    
    n = len(tableau)  # number of stabilizers
    twoL = len(tableau[0])  # 2L where L is number of qubits
    L = twoL // 2
    
    # Build symplectic form J = [[0, I], [I, 0]]
    J = build_symplectic_form_poly(L)
    
    # Centralizer = nullspace of (tableau · J)
    tableau_J = matrix_multiply_poly(tableau, J)
    central_basis = nullspace_poly(tableau_J)
    
    # Pick those that extend the stabilizer span to get exactly 2k vectors
    stab_rank = rank_poly(tableau)
    k = L - stab_rank
    logical_basis = []
    
    # Start with tableau (stabilizer generators)
    aug = [row[:] for row in tableau]  # Deep copy
    
    for v in central_basis:
        if len(logical_basis) >= 2 * k:
            break
        
        # Check if v extends the span
        test_aug = aug + [v]
        if rank_poly(test_aug) > rank_poly(aug):
            logical_basis.append(v[:])  # Copy vector
            aug.append(v[:])  # Add to augmented matrix
    
    return logical_basis


def find_distance(tableau, *, return_logical_ops: bool = False) -> Union[int, Tuple[int, List]]:
    """
    Compute the code distance by brute-forcing all nonzero combos of the 2k logical generators.
    
    Uses polynomial weight (number of nonzero monomials) instead of Hamming weight.
    
    Args:
        tableau: Polynomial stabilizer tableau
        return_logical_ops: Whether to also return logical operators
        
    Returns:
        Code distance, or (distance, logical_ops) if return_logical_ops=True
    """
    from .polynomial import rank_poly
    from typing import Union, Tuple, List
    import itertools
    
    if not tableau or not tableau[0]:
        if return_logical_ops:
            return 0, []
        return 0
    
    twoL = len(tableau[0])
    L = twoL // 2
    stab_rank = rank_poly(tableau)
    k = L - stab_rank
    
    if k == 0:
        if return_logical_ops:
            return 0, []
        return 0   # no logical qubits ⇒ distance 0

    logical_ops = find_logical_operators(tableau)
    best_distance = L  # Initialize to maximum possible distance
    
    def poly_weight(poly_vec):
        """Compute polynomial weight: number of nonzero monomials in each qubit."""
        x_part = poly_vec[:L]
        z_part = poly_vec[L:]
        
        weight = 0
        for i in range(L):
            if not x_part[i].is_zero() or not z_part[i].is_zero():
                weight += 1
        return weight
    
    # Enumerate all nonempty subsets of logical operators
    n = len(logical_ops)
    for r in range(1, n + 1):
        for combo in itertools.combinations(range(n), r):
            # Linear combination of logical operators
            vec = [logical_ops[combo[0]][i] for i in range(twoL)]  # Start with first operator
            
            # Add remaining operators
            for j in range(1, len(combo)):
                for i in range(twoL):
                    vec[i] = vec[i] + logical_ops[combo[j]][i]
            
            w = poly_weight(vec)
            if w < best_distance:
                best_distance = w
                if best_distance == 1:
                    if return_logical_ops:
                        return 1, logical_ops
                    else:
                        return 1
    
    if return_logical_ops:
        return best_distance, logical_ops
    return best_distance


def format_polynomial_vector(v: List[LaurentPolynomial]) -> str:
    """
    Convert a polynomial symplectic vector into its Pauli-string representation.
    
    Args:
        v: Polynomial symplectic vector
        
    Returns:
        Pauli string representation
    """
    raise NotImplementedError("Polynomial vector formatting not yet implemented") 