import itertools
import numpy as np
import logging

def _rref_mod2(M: np.ndarray):
    """
    Reduced row-echelon form over GF(2).
    Returns (R, pivots) where pivots is a list of pivot columns.
    """
    A = M.copy() % 2
    rows, cols = A.shape
    pivots = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        # find a pivot in column c at or below row r
        for i in range(r, rows):
            if A[i, c] == 1:
                A[[r, i]] = A[[i, r]]
                break
        else:
            continue
        # eliminate all other 1’s in column c
        for j in range(rows):
            if j != r and A[j, c] == 1:
                A[j] ^= A[r]
        pivots.append(c)
        r += 1
    return A, pivots

def _rank_mod2(M: np.ndarray) -> int:
    """Rank over GF(2)."""
    _, pivots = _rref_mod2(M)
    return len(pivots)

def _nullspace_mod2(M: np.ndarray) -> list[np.ndarray]:
    """
    Basis for the null-space of M over GF(2).
    Returns a list of column vectors v with M·v = 0 (mod 2).
    """
    R, pivots = _rref_mod2(M)
    rows, cols = R.shape
    free_vars = [c for c in range(cols) if c not in pivots]
    basis = []
    for fv in free_vars:
        v = np.zeros(cols, dtype=int)
        v[fv] = 1
        # back-substitute for pivot rows
        for i, pc in enumerate(pivots):
            if R[i, fv] == 1:
                v[pc] = 1
        basis.append(v)
    return basis

def find_logical_operators(tableau: np.ndarray) -> list[np.ndarray]:
    """
    Compute a basis (2k vectors) for the logical operators modulo the stabilizer.
    Same signature as your old function, so you can swap it in directly.
    """
    logger = logging.getLogger("qca")
    n, twoL = tableau.shape
    L = twoL // 2

    # build symplectic form J = [[0, I],[I, 0]]
    J = np.block([
        [np.zeros((L, L), dtype=int), np.eye(L, dtype=int)],
        [np.eye(L, dtype=int),       np.zeros((L, L), dtype=int)]
    ]) % 2

    # centralizer = null-space of (S·J)
    C = (tableau.dot(J)) % 2
    central_basis = _nullspace_mod2(C)

    logger.info("Computing logical operator basis:")
    logger.info("  Centralizer basis size: %d", len(central_basis))

    # pick those that extend the stabilizer span to get exactly 2k vectors
    stab_rank = _rank_mod2(tableau)
    k = L - stab_rank
    logical_basis = []
    aug = tableau.copy()
    
    logger.info("  Target logical qubits (k): %d", k)
    logger.info("  Need %d logical operators (2k)", 2 * k)
    
    for i, v in enumerate(central_basis):
        if len(logical_basis) >= 2 * k:
            break
        new_aug = np.vstack([aug, v.reshape(1, -1)]) % 2
        if _rank_mod2(new_aug) > _rank_mod2(aug):
            pauli_str = format_symplectic_vector(v)
            weight = sum(1 for char in pauli_str if char != 'I')
            logger.info("  Logical basis[%d]: %s (weight: %d)", len(logical_basis), pauli_str, weight)
            logical_basis.append(v.copy())
            aug = new_aug

    logger.info("  Final logical basis size: %d", len(logical_basis))
    return logical_basis

def find_distance(tableau: np.ndarray, *, return_logical_ops: bool = False):
    """
    Compute the code distance by brute-forcing all nonzero combos of the 2k logical generators.
    """
    logger = logging.getLogger("qca")
    L = tableau.shape[1] // 2
    stab_rank = _rank_mod2(tableau)
    if L - stab_rank == 0:
        if return_logical_ops:
            return 0, []
        return 0   # no logical qubits ⇒ distance 0

    log_ops = find_logical_operators(tableau)
    best = L

    def _weight(v: np.ndarray) -> int:
        x, z = v[:L], v[L:]
        return int(np.count_nonzero(x | z))

    logger.info("Finding code distance:")
    logger.info("  Testing combinations of %d logical operators", len(log_ops))

    # enumerate all nonempty subsets
    n = len(log_ops)
    combinations_tested = 0
    best_combination = None
    
    for r in range(1, n + 1):
        for combo in itertools.combinations(range(n), r):
            vec = sum(log_ops[i] for i in combo) % 2
            w = _weight(vec)
            combinations_tested += 1
            
            if w < best:
                best = w
                best_combination = combo
                combined_pauli = format_symplectic_vector(vec)
                logger.info("  New minimum: weight %d from combination %s → %s", w, combo, combined_pauli)
                
                if best == 1:
                    logger.info("  Found distance 1 - stopping search")
                    if not return_logical_ops:
                        return 1
                    else:
                        return (1, log_ops)
    
    logger.info("  Distance search complete: tested %d combinations", combinations_tested)
    logger.info("  Final distance: %d (from combination %s)", best, best_combination)
    
    if return_logical_ops:
        return best, log_ops
    return best

def format_symplectic_vector(v: np.ndarray) -> str:
    """
    Convert a 2L-length symplectic vector v into its Pauli-string representation.
    """
    L = len(v) // 2
    x_bits, z_bits = v[:L], v[L:]
    pauli_chars = []
    for i in range(L):
        xi, zi = int(x_bits[i]), int(z_bits[i])
        if   xi == 0 and zi == 0: pauli_chars.append('I')
        elif xi == 1 and zi == 0: pauli_chars.append('X')
        elif xi == 0 and zi == 1: pauli_chars.append('Z')
        else:                     pauli_chars.append('Y')
    return ''.join(pauli_chars)
