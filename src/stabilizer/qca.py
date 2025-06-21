"""
Quantum Cellular Automaton (QCA) evolution for Pauli strings.

Implements the symplectic formulation with matrix M(x) = (x^-1+1+x, 1; 1, 0) mod 2
for evolving Pauli operators under the QCA unitary transformation.
"""

import numpy as np
from typing import Tuple, List


def pauli_to_symplectic_vectors(pauli_string: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert a Pauli string to symplectic binary vectors (f, g).
    
    Args:
        pauli_string: List of Pauli operators ['X', 'Z', 'I', 'Y', ...]
        
    Returns:
        Tuple of (f, g) where f, g are binary vectors of length N
        Encoding: (0,0)→I, (1,0)→X, (0,1)→Z, (1,1)→Y
    """
    N = len(pauli_string)
    f = np.zeros(N, dtype=int)
    g = np.zeros(N, dtype=int)
    
    pauli_map = {
        'I': (0, 0),
        'X': (1, 0), 
        'Z': (0, 1),
        'Y': (1, 1)
    }
    
    for i, pauli in enumerate(pauli_string):
        f[i], g[i] = pauli_map[pauli]
    
    return f, g


def symplectic_vectors_to_pauli(f: np.ndarray, g: np.ndarray) -> List[str]:
    """
    Convert symplectic binary vectors back to Pauli string.
    
    Args:
        f, g: Binary vectors of length N
        
    Returns:
        List of Pauli operators corresponding to (f, g)
    """
    pauli_map = {
        (0, 0): 'I',
        (1, 0): 'X',
        (0, 1): 'Z', 
        (1, 1): 'Y'
    }
    
    N = len(f)
    pauli_string = []
    
    for i in range(N):
        pauli_string.append(pauli_map[(f[i], g[i])])
    
    return pauli_string


def apply_qca_matrix(f: np.ndarray, g: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply one step of QCA evolution using matrix M(x).
    
    The transformation is:
    f'_i = f_{i-1} + f_i + f_{i+1} + g_i (mod 2)
    g'_i = f_i (mod 2)
    
    Args:
        f, g: Current symplectic vectors
        
    Returns:  
        f', g': Updated symplectic vectors after one QCA step
    """
    N = len(f)
    f_new = np.zeros(N, dtype=int)
    g_new = np.zeros(N, dtype=int)
    
    for i in range(N):
        # Periodic boundary conditions
        i_prev = (i - 1) % N
        i_next = (i + 1) % N
        
        # Apply M(x) transformation
        f_new[i] = (f[i_prev] + f[i] + f[i_next] + g[i]) % 2
        g_new[i] = f[i] % 2
    
    return f_new, g_new


def evolve_pauli_string(pauli_string: List[str], num_steps: int = 1, debug: bool = False) -> List[List[str]]:
    """
    Evolve a Pauli string through multiple QCA time steps.
    
    Args:
        pauli_string: Initial Pauli string
        num_steps: Number of QCA evolution steps
        debug: Whether to print debug information
        
    Returns:
        List of evolved Pauli strings, one for each time step
    """
    if debug:
        print(f"    QCA DEBUG: Starting evolution for {num_steps} steps")
        print(f"    QCA DEBUG: Initial string: {''.join(pauli_string)}")
    
    # Convert to symplectic form
    f, g = pauli_to_symplectic_vectors(pauli_string)
    
    if debug:
        print(f"    QCA DEBUG: Initial f = {f}")
        print(f"    QCA DEBUG: Initial g = {g}")
    
    evolved_strings = []
    
    for step in range(num_steps):
        # Apply one QCA evolution step
        f, g = apply_qca_matrix(f, g)
        
        # Convert back to Pauli string
        evolved_pauli = symplectic_vectors_to_pauli(f, g)
        evolved_strings.append(evolved_pauli)
        
        if debug:
            print(f"    QCA DEBUG: Step {step + 1}:")
            print(f"    QCA DEBUG:   f = {f}")
            print(f"    QCA DEBUG:   g = {g}")
            print(f"    QCA DEBUG:   Pauli = {''.join(evolved_pauli)}")
    
    return evolved_strings


def qca_evolution_step(pauli_string: List[str]) -> List[str]:
    """
    Apply a single QCA evolution step to a Pauli string.
    
    Args:
        pauli_string: Input Pauli string
        
    Returns:
        Evolved Pauli string after one time step
    """
    f, g = pauli_to_symplectic_vectors(pauli_string)
    f_new, g_new = apply_qca_matrix(f, g)
    return symplectic_vectors_to_pauli(f_new, g_new) 