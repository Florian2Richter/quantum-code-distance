"""
Brute-force search (or heuristic) for minimal logical operator weight.
"""

import itertools
import numpy as np
from .utils import is_logical_op


def find_distance(tableau: np.ndarray) -> int:
    """
    Find the minimum weight of a logical operator (code distance).
    WARNING: NP-complete, exponential in L for brute force.
    """
    L = tableau.shape[1] // 2
    
    # Early exit if no logical qubits
    if tableau.shape[0] >= L:
        return 0
    
    # Start from weight 1 and increase until we find a logical operator
    for w in range(1, L + 1):
        print(f"  Searching weight {w} operators...")
        for positions in itertools.combinations(range(L), w):
            if is_logical_op(positions, tableau):
                return w
    
    # If no logical operator found (shouldn't happen for valid codes)
    return L


def find_logical_operators(tableau: np.ndarray, max_weight: int = None) -> list[tuple]:
    """
    Find all logical operators up to a given weight.
    Useful for debugging and understanding the code structure.
    """
    L = tableau.shape[1] // 2
    if max_weight is None:
        max_weight = min(L, 5)  # Default limit to avoid exponential explosion
    
    logical_ops = []
    for w in range(1, max_weight + 1):
        for positions in itertools.combinations(range(L), w):
            if is_logical_op(positions, tableau):
                logical_ops.append(positions)
    
    return logical_ops 