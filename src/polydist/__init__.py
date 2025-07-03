"""
Laurent polynomial formalism for stabilizer codes.

This module implements stabilizer operations using Laurent polynomials over F₂[x,x⁻¹]
instead of GF(2) bit-vectors, maintaining exact equivalence with the original implementation.
"""

__version__ = "0.1.0"

# Import sub-modules
from . import lattice
from . import tableau
from . import distance
from . import utils
from . import polynomial
from . import qca 