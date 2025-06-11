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


def weight_of_op(op: list[str]) -> int:
    """Count non-Identity Paulis"""
    return sum(1 for p in op if p != 'I')


def is_logical_op(positions: tuple[int, ...], tableau: np.ndarray) -> bool:
    """
    Placeholder: implement commutation with all stabilizers
    and ensure op is not in the generated stabilizer group.
    """
    return False 