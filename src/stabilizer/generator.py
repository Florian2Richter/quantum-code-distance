"""
Parse and validate a seed string of Pauli operators.
"""

VALID_PAULI = {'I', 'X', 'Y', 'Z'}

def parse_seed(seed: str) -> list[str]:
    seed = seed.upper()
    if not all(ch in VALID_PAULI for ch in seed):
        raise ValueError(f"Invalid Pauli in seed: {seed}")
    return list(seed) 