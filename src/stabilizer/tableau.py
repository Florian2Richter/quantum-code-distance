import numpy as np
from .utils import pauli_to_symplectic

"""
Build binary symplectic stabilizer tableau and compute its rank.
"""

def build_tableau(stab_ops: list[list[str]]) -> np.ndarray:
    # Each op → 2L-bit vector
    L = len(stab_ops[0])
    M = len(stab_ops)
    tab = np.zeros((M, 2*L), dtype=int)
    for i, op in enumerate(stab_ops):
        tab[i, :] = pauli_to_symplectic(op)
    return tab


def compute_rank(tableau: np.ndarray) -> int:
    # Over GF(2)
    A = tableau.copy() % 2
    rank = 0
    rows, cols = A.shape
    c = 0
    for r in range(rows):
        while c < cols and not any(A[i, c] for i in range(r, rows)):
            c += 1
        if c == cols:
            break
        pivot = next(i for i in range(r, rows) if A[i, c])
        A[[r, pivot]] = A[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        rank += 1
        c += 1
    return rank 