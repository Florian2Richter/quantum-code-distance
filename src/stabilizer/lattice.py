"""
Replicate a single generator around a 1D ring with periodic boundary conditions.
"""

def build_lattice(pauli_seed: list[str]) -> list[list[str]]:
    """
    Return a list of length-L stabilizer generators, each a list of L Paulis.
    """
    L = len(pauli_seed)
    ops = []
    for shift in range(L):
        ops.append([pauli_seed[(i - shift) % L] for i in range(L)])
    return ops 