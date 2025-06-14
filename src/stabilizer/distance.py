"""
Brute-force search (or heuristic) for minimal logical operator weight.
"""

import itertools
import numpy as np
from tqdm import tqdm
from .utils import (
    symplectic_product,
    is_in_stabilizer_group,
    is_logical_vec,
)


def find_distance(tableau: np.ndarray) -> int:
    """
    Find the minimal qubit-weight of a logical operator.
    We now enumerate binary vectors v of length 2L, 
    grouped by qubit-weight.
    """
    L = tableau.shape[1] // 2
    
    # Early exit: if rank >= L, no logicals at all
    if tableau.shape[0] >= L:
        return 0

    # Precompute all non-zero single-qubit symplectic patterns
    # (X=(1,0), Z=(0,1), Y=(1,1)) for convenience
    single_qubit_patterns = [(1, 0), (0, 1), (1, 1)]

    # For each target weight w = 1..L
    for w in range(1, L + 1):
        print(f"  Searching symplectic vectors of qubit-weight {w}...")
        
        # Calculate total number of combinations for progress bar
        num_positions = len(list(itertools.combinations(range(L), w)))
        num_patterns = len(single_qubit_patterns) ** w
        total_combinations = num_positions * num_patterns
        
        # Use progress bar for this weight
        with tqdm(total=total_combinations, desc=f"Weight {w}", leave=False) as pbar:
            # 1) Choose which w qubits are non-identity
            for qubit_positions in itertools.combinations(range(L), w):

                # 2) For each assignment of X/Z/Y on those positions
                for labels in itertools.product(single_qubit_patterns, repeat=w):
                    # build the 2L-vector
                    v = np.zeros(2 * L, dtype=int)
                    for i, (qx, qz) in zip(qubit_positions, labels):
                        v[i] = qx      # X-bit
                        v[i + L] = qz  # Z-bit

                    # 3) test if logical
                    if is_logical_vec(v, tableau):
                        pbar.close()  # Close progress bar early
                        return w
                    
                    pbar.update(1)

    return L


def find_logical_operators(tableau: np.ndarray, max_weight: int = None) -> list[np.ndarray]:
    """
    Find all logical operators up to a given weight.
    Returns list of symplectic vectors instead of position tuples.
    """
    L = tableau.shape[1] // 2
    if max_weight is None:
        max_weight = min(L, 5)  # Default limit to avoid exponential explosion

    logical_ops = []
    single_qubit_patterns = [(1, 0), (0, 1), (1, 1)]
    
    # Calculate total combinations across all weights for progress bar
    total_combinations = 0
    for w in range(1, max_weight + 1):
        num_positions = len(list(itertools.combinations(range(L), w)))
        num_patterns = len(single_qubit_patterns) ** w
        total_combinations += num_positions * num_patterns

    with tqdm(total=total_combinations, desc="Finding logical operators") as pbar:
        for w in range(1, max_weight + 1):
            for qubit_positions in itertools.combinations(range(L), w):
                for labels in itertools.product(single_qubit_patterns, repeat=w):
                    v = np.zeros(2 * L, dtype=int)
                    for i, (qx, qz) in zip(qubit_positions, labels):
                        v[i] = qx
                        v[i + L] = qz

                    if is_logical_vec(v, tableau):
                        logical_ops.append(v.copy())
                    
                    pbar.update(1)

    return logical_ops


def find_distance_with_early_stop(tableau: np.ndarray, target_distance: int = None) -> int:
    """
    Optimized distance finder that stops early if target distance is reached.
    Useful when you have an estimate of the expected distance.
    """
    L = tableau.shape[1] // 2
    
    if tableau.shape[0] >= L:
        return 0

    single_qubit_patterns = [(1, 0), (0, 1), (1, 1)]
    max_weight = target_distance if target_distance else L

    for w in range(1, max_weight + 1):
        print(f"  Searching weight {w} (target: {target_distance})...")
        
        num_positions = len(list(itertools.combinations(range(L), w)))
        num_patterns = len(single_qubit_patterns) ** w
        total_combinations = num_positions * num_patterns
        
        with tqdm(total=total_combinations, desc=f"Weight {w}", leave=False) as pbar:
            for qubit_positions in itertools.combinations(range(L), w):
                for labels in itertools.product(single_qubit_patterns, repeat=w):
                    v = np.zeros(2 * L, dtype=int)
                    for i, (qx, qz) in zip(qubit_positions, labels):
                        v[i] = qx
                        v[i + L] = qz

                    if is_logical_vec(v, tableau):
                        pbar.close()
                        return w
                    
                    pbar.update(1)

    return L


def format_symplectic_vector(v: np.ndarray) -> str:
    """
    Convert symplectic vector back to Pauli string for display.
    """
    L = len(v) // 2
    x_bits, z_bits = v[:L], v[L:]
    
    pauli_chars = []
    for i in range(L):
        if x_bits[i] == 0 and z_bits[i] == 0:
            pauli_chars.append('I')
        elif x_bits[i] == 1 and z_bits[i] == 0:
            pauli_chars.append('X')
        elif x_bits[i] == 0 and z_bits[i] == 1:
            pauli_chars.append('Z')
        elif x_bits[i] == 1 and z_bits[i] == 1:
            pauli_chars.append('Y')
    
    return ''.join(pauli_chars) 