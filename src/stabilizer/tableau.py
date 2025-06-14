import numpy as np
from .utils import pauli_to_symplectic

"""
Build binary symplectic stabilizer tableau and compute its rank.
"""

def build_tableau(stab_ops: list[list[str]]) -> np.ndarray:
    """
    Build binary symplectic tableau from stabilizer operators.
    Each operator → 2L-bit vector in symplectic form [x|z].
    """
    if not stab_ops:
        raise ValueError("No stabilizer operators provided")
    
    L = len(stab_ops[0])
    M = len(stab_ops)
    
    # Validate all operators have same length
    for i, op in enumerate(stab_ops):
        if len(op) != L:
            raise ValueError(f"Operator {i} has length {len(op)}, expected {L}")
    
    tab = np.zeros((M, 2*L), dtype=int)
    for i, op in enumerate(stab_ops):
        tab[i, :] = pauli_to_symplectic(op)
    return tab


def compute_rank(tableau: np.ndarray) -> int:
    """
    Compute rank of binary matrix over GF(2) using Gaussian elimination.
    """
    if tableau.size == 0:
        return 0
        
    # Work over GF(2)
    A = tableau.copy() % 2
    rank = 0
    rows, cols = A.shape
    col = 0
    
    for row in range(rows):
        # Find pivot in current column
        while col < cols and not any(A[r, col] for r in range(row, rows)):
            col += 1
        if col == cols:
            break
            
        # Find first row with 1 in current column
        pivot_row = next(r for r in range(row, rows) if A[r, col])
        
        # Swap rows to bring pivot to current row
        if pivot_row != row:
            A[[row, pivot_row]] = A[[pivot_row, row]]
        
        # Eliminate column in all other rows
        for r in range(rows):
            if r != row and A[r, col]:
                A[r, :] ^= A[row, :]
        
        rank += 1
        col += 1
    
    return rank 