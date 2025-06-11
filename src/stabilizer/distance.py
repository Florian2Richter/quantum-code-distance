"""
Brute-force search (or heuristic) for minimal logical operator weight.
"""

import itertools
from .utils import is_logical_op


def find_distance(tableau) -> int:
    # WARNING: NP-complete, exponential in L
    L = tableau.shape[1] // 2
    for w in range(1, L+1):
        for positions in itertools.combinations(range(L), w):
            if is_logical_op(positions, tableau):
                return w
    return L 