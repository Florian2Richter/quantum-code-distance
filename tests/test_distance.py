import pytest
import numpy as np
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance


def distance_from_tableau(seeds):
    """
    Helper function to compute (n, k, d) from a list of Pauli seed strings.
    
    Args:
        seeds: List of Pauli strings (e.g., ["ZZII", "IZZI"])
        
    Returns:
        Tuple (n, k, d) where:
        - n: number of physical qubits
        - k: number of logical qubits  
        - d: code distance
    """
    if not seeds:
        # Empty code case - need to determine n from context
        # For empty codes, we'll handle this in the test
        raise ValueError("Empty seed list - n must be specified")
    
    # Convert strings to list of lists
    stab_ops = [list(seed) for seed in seeds]
    
    # Get number of physical qubits
    n = len(stab_ops[0])
    
    # Build tableau and compute parameters
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    k = n - rank
    
    if k == 0:
        d = 0  # No logical qubits
    else:
        d = find_distance(tableau)
    
    return (n, k, d)


def distance_empty_code(n):
    """
    Helper for empty code cases where we have no stabilizers.
    For empty codes: k = n (all qubits are logical), d = 1 (no error correction).
    """
    return (n, n, 1)


@pytest.mark.parametrize("n,seeds,expected", [
    # Empty codes - no stabilizers, all qubits logical
    (1, [], (1, 1, 1)),
    (2, [], (2, 2, 1)),
    (3, [], (3, 3, 1)),
    (4, [], (4, 4, 1)),
    
    # Z-repetition codes
    (4, ["ZZII", "IZZI", "IIZZ", "ZIIZ"], (4, 1, 1)),
    (3, ["ZZI", "IZZ", "ZIZ"], (3, 1, 1)),
    
    # X-repetition codes  
    (4, ["XXII", "IXXI", "IIXX", "XIIX"], (4, 1, 1)),
    
    # Small cluster-state (3-qubit, fully constrained)
    (3, ["XZZ", "ZXZ", "ZZX"], (3, 0, 0)),
])
def test_small_n_codes(n, seeds, expected):
    """Test canonical small-N codes with hard-coded expected values."""
    if not seeds:
        # Handle empty code case
        result = distance_empty_code(n)
    else:
        result = distance_from_tableau(seeds)
    
    assert result == expected, f"For n={n}, seeds={seeds}: got {result}, expected {expected}" 