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


def seed_is_valid(seed: str) -> bool:
    """
    Check if a Pauli seed string generates valid commuting stabilizer operators.
    
    Args:
        seed: String of Pauli operators (e.g., 'XZY')
        
    Returns:
        True if all cyclic translations of the seed commute with each other
    """
    N = len(seed)
    
    # Generate all N cyclic translations of the seed
    translations = []
    for i in range(N):
        # Cyclic shift: move first i characters to the end
        shifted = seed[i:] + seed[:i]
        translations.append(list(shifted))
    
    # Convert each translation to binary symplectic vector
    symplectic_vectors = []
    for translation in translations:
        vec = pauli_to_symplectic(translation)
        symplectic_vectors.append(vec)
    
    # Stack vectors into N × 2N matrix S
    S = np.array(symplectic_vectors, dtype=int)
    
    # Construct symplectic form J = [[0, I], [-I, 0]] of size 2N × 2N
    I_N = np.eye(N, dtype=int)
    zero_N = np.zeros((N, N), dtype=int)
    J = np.block([[zero_N, I_N], [-I_N, zero_N]])
    
    # Compute C = S · J · S^T (mod 2)
    C = np.dot(np.dot(S, J), S.T) % 2
    
    # Return True if C is the zero matrix
    return np.all(C == 0)


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





