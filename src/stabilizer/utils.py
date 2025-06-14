import numpy as np

"""
Helper functions: Pauli ↔ symplectic, weight, logical checks.
"""

PAULI_MAP = {
    'I': (0, 0),
    'X': (1, 0),
    'Z': (0, 1),
    'Y': (1, 1),
}

def pauli_to_symplectic(op: list[str]) -> np.ndarray:
    """Convert list of Paulis to 2L-bit vector [x|z]"""
    L = len(op)
    x = np.zeros(L, dtype=int)
    z = np.zeros(L, dtype=int)
    for i, p in enumerate(op):
        xi, zi = PAULI_MAP[p]
        x[i], z[i] = xi, zi
    return np.concatenate([x, z])


def symplectic_product(v1: np.ndarray, v2: np.ndarray) -> int:
    """
    Compute symplectic inner product between two 2L-dimensional vectors.
    Returns 0 if operators commute, 1 if they anticommute.
    """
    L = len(v1) // 2
    x1, z1 = v1[:L], v1[L:]
    x2, z2 = v2[:L], v2[L:]
    return (np.dot(x1, z2) + np.dot(z1, x2)) % 2


def positions_to_symplectic(positions: tuple[int, ...], L: int, pauli_type: str = 'X') -> np.ndarray:
    """
    Convert qubit positions to symplectic vector.
    pauli_type: 'X' for X operators, 'Z' for Z operators
    """
    vec = np.zeros(2 * L, dtype=int)
    if pauli_type == 'X':
        vec[list(positions)] = 1
    elif pauli_type == 'Z':
        vec[L + np.array(list(positions))] = 1
    return vec


def is_in_stabilizer_group(op_vector: np.ndarray, tableau: np.ndarray) -> bool:
    """
    Check if an operator (symplectic vector) is in the stabilizer group.
    This uses Gaussian elimination over GF(2).
    """
    # Augment tableau with the operator vector
    augmented = np.vstack([tableau, op_vector.reshape(1, -1)]) % 2
    
    # Perform Gaussian elimination
    rows, cols = augmented.shape
    for c in range(min(rows - 1, cols)):  # Don't eliminate the last row (our test vector)
        # Find pivot in current column
        pivot_row = None
        for r in range(c, rows - 1):  # Exclude last row from pivot search
            if augmented[r, c] == 1:
                pivot_row = r
                break
        
        if pivot_row is None:
            continue
            
        # Swap pivot row to current position
        if pivot_row != c:
            augmented[[c, pivot_row]] = augmented[[pivot_row, c]]
        
        # Eliminate column in other rows (including last row)
        for r in range(rows):
            if r != c and augmented[r, c] == 1:
                augmented[r] = (augmented[r] + augmented[c]) % 2
    
    # Check if last row (our test vector) became zero
    return np.all(augmented[-1] == 0)


def is_logical_vec(v: np.ndarray, tableau: np.ndarray) -> bool:
    """
    Check if symplectic vector v (shape (2L,)) is a non-trivial logical:
      1) it commutes with every stabilizer (symplectic product = 0)
      2) it is not itself in the span of the stabilizer rows
    """
    # 1) Commutation check
    for stab in tableau:
        if symplectic_product(v, stab) != 0:
            return False

    # 2) Exclude pure stabilizers
    return not is_in_stabilizer_group(v, tableau) 