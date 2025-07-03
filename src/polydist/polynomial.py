"""
Laurent polynomial arithmetic over F₂[x,x⁻¹].

This module implements the core polynomial data structures and operations
needed for the polynomial formalism of stabilizer codes.
"""

from typing import Dict, Union, List, Tuple
import numpy as np


class LaurentPolynomial:
    """
    A Laurent polynomial over F₂[x,x⁻¹].
    
    Represents polynomials of the form: ∑ₖ aₖ xᵏ where aₖ ∈ {0,1} and k ∈ ℤ.
    Stored as a dictionary mapping exponents to coefficients.
    """
    
    def __init__(self, coeffs: Dict[int, int] = None):
        """
        Initialize Laurent polynomial.
        
        Args:
            coeffs: Dictionary mapping exponents to coefficients (mod 2)
        """
        self.coeffs = {}
        if coeffs:
            # Only store non-zero coefficients
            for exp, coeff in coeffs.items():
                if coeff % 2 != 0:
                    self.coeffs[exp] = 1
    
    @classmethod
    def zero(cls):
        """Return the zero polynomial."""
        return cls({})
    
    @classmethod
    def one(cls):
        """Return the constant polynomial 1."""
        return cls({0: 1})
    
    @classmethod
    def x(cls, power: int = 1):
        """Return the monomial x^power."""
        return cls({power: 1})
    
    def is_zero(self) -> bool:
        """Check if polynomial is zero."""
        return len(self.coeffs) == 0
    
    def degree(self) -> Union[int, float]:
        """Return the degree (highest exponent) or -∞ for zero polynomial."""
        if self.is_zero():
            return float('-inf')
        return max(self.coeffs.keys())
    
    def min_degree(self) -> Union[int, float]:
        """Return the minimum exponent or +∞ for zero polynomial."""
        if self.is_zero():
            return float('+inf')
        return min(self.coeffs.keys())
    
    def weight(self) -> int:
        """Return the number of nonzero monomials."""
        return len(self.coeffs)
    
    def __add__(self, other):
        """Add two polynomials (XOR of coefficients)."""
        if not isinstance(other, LaurentPolynomial):
            return NotImplemented
        
        result_coeffs = self.coeffs.copy()
        for exp, coeff in other.coeffs.items():
            if exp in result_coeffs:
                result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                if result_coeffs[exp] == 0:
                    del result_coeffs[exp]
            else:
                result_coeffs[exp] = coeff % 2
        
        return LaurentPolynomial(result_coeffs)
    
    def __mul__(self, other):
        """Multiply two polynomials."""
        if not isinstance(other, LaurentPolynomial):
            return NotImplemented
        
        if self.is_zero() or other.is_zero():
            return LaurentPolynomial.zero()
        
        result_coeffs = {}
        for exp1, coeff1 in self.coeffs.items():
            for exp2, coeff2 in other.coeffs.items():
                exp = exp1 + exp2
                coeff = (coeff1 * coeff2) % 2
                if exp in result_coeffs:
                    result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                else:
                    result_coeffs[exp] = coeff
        
        # Remove zero coefficients
        result_coeffs = {exp: coeff for exp, coeff in result_coeffs.items() if coeff != 0}
        
        return LaurentPolynomial(result_coeffs)
    
    def __eq__(self, other) -> bool:
        """Check equality of polynomials."""
        if not isinstance(other, LaurentPolynomial):
            return False
        return self.coeffs == other.coeffs
    
    def __repr__(self) -> str:
        """String representation of polynomial."""
        if self.is_zero():
            return "0"
        
        terms = []
        for exp in sorted(self.coeffs.keys()):
            if exp == 0:
                terms.append("1")
            elif exp == 1:
                terms.append("x")
            elif exp > 1:
                terms.append(f"x^{exp}")
            elif exp == -1:
                terms.append("x^(-1)")
            else:
                terms.append(f"x^({exp})")
        
        return " + ".join(terms)
    
    def __hash__(self) -> int:
        """Make polynomial hashable."""
        return hash(frozenset(self.coeffs.items()))


def rref_poly(matrix: List[List[LaurentPolynomial]], max_degree: int = 10) -> List[List[LaurentPolynomial]]:
    """
    Reduced row echelon form for polynomial matrices over F₂[x,x⁻¹].
    
    For now, implements a simple version that works for constant polynomial matrices
    (i.e., matrices where all entries are just 0 or 1).
    
    Args:
        matrix: Matrix of Laurent polynomials
        max_degree: Maximum degree bound for polynomial arithmetic (unused for now)
        
    Returns:
        Matrix in reduced row echelon form
    """
    if not matrix or not matrix[0]:
        return matrix
    
    # Make a deep copy to avoid modifying the original
    result = []
    for row in matrix:
        result.append([p for p in row])  # Copy each polynomial
    
    rows = len(result)
    cols = len(result[0])
    
    current_row = 0
    for col in range(cols):
        # Find pivot row (first non-zero entry in current column)
        pivot_row = -1
        for row in range(current_row, rows):
            if not result[row][col].is_zero():
                pivot_row = row
                break
        
        if pivot_row == -1:
            continue  # No pivot in this column
        
        # Swap rows if needed
        if pivot_row != current_row:
            result[current_row], result[pivot_row] = result[pivot_row], result[current_row]
        
        # The pivot should be 1 (for F₂, any non-zero element is 1)
        pivot = result[current_row][col]
        if pivot.is_zero():
            continue
        
        # For now, we assume pivot is always 1 (since we're working over F₂)
        # In the future, we'd need polynomial division here
        
        # Eliminate other entries in this column
        for row in range(rows):
            if row != current_row and not result[row][col].is_zero():
                # Add current_row to row (which subtracts in F₂)
                for c in range(cols):
                    result[row][c] = result[row][c] + result[current_row][c]
        
        current_row += 1
        if current_row >= rows:
            break
    
    return result


def nullspace_poly(matrix: List[List[LaurentPolynomial]]) -> List[List[LaurentPolynomial]]:
    """
    Compute nullspace of polynomial matrix over F₂[x,x⁻¹].
    
    For now, implements a simple version for constant polynomial matrices.
    
    Args:
        matrix: Matrix of Laurent polynomials
        
    Returns:
        Basis for the nullspace
    """
    if not matrix or not matrix[0]:
        # Empty matrix has full nullspace
        return []
    
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Get RREF to find pivot columns
    rref_matrix = rref_poly(matrix)
    
    # Find pivot columns
    pivot_cols = []
    pivot_rows = []
    
    for row_idx in range(rows):
        for col_idx in range(cols):
            if not rref_matrix[row_idx][col_idx].is_zero():
                pivot_cols.append(col_idx)
                pivot_rows.append(row_idx)
                break
    
    # Free columns are those not in pivot_cols
    free_cols = [i for i in range(cols) if i not in pivot_cols]
    
    # Build basis vectors for nullspace
    nullspace_basis = []
    
    for free_col in free_cols:
        # Create basis vector for this free variable
        basis_vector = [LaurentPolynomial.zero() for _ in range(cols)]
        basis_vector[free_col] = LaurentPolynomial.one()
        
        # Set dependent variables based on RREF
        for i, pivot_col in enumerate(pivot_cols):
            if i < len(pivot_rows):
                pivot_row = pivot_rows[i]
                # basis_vector[pivot_col] = -rref_matrix[pivot_row][free_col]
                # In F₂, -x = x, so:
                basis_vector[pivot_col] = rref_matrix[pivot_row][free_col]
        
        nullspace_basis.append(basis_vector)
    
    return nullspace_basis


def rank_poly(matrix: List[List[LaurentPolynomial]]) -> int:
    """
    Compute rank of polynomial matrix over F₂[x,x⁻¹].
    
    Args:
        matrix: Matrix of Laurent polynomials
        
    Returns:
        Rank of the matrix
    """
    if not matrix or not matrix[0]:
        return 0
    
    # Get RREF and count non-zero rows
    rref_matrix = rref_poly(matrix)
    
    rank = 0
    for row in rref_matrix:
        # Check if row is non-zero
        row_is_zero = True
        for entry in row:
            if not entry.is_zero():
                row_is_zero = False
                break
        if not row_is_zero:
            rank += 1
    
    return rank


def pauli_to_polynomial_vector(pauli_op: List[str]) -> List[LaurentPolynomial]:
    """
    Convert a Pauli operator to polynomial symplectic vector.
    
    For now, this just converts the GF(2) representation to constant polynomials.
    In the future, this will be enhanced to use actual polynomial encoding.
    
    Args:
        pauli_op: List of Pauli operators ['X', 'Z', 'I', 'Y']
        
    Returns:
        Polynomial symplectic vector [x_polys|z_polys]
    """
    L = len(pauli_op)
    x_polys = []
    z_polys = []
    
    pauli_map = {
        'I': (0, 0),
        'X': (1, 0),
        'Z': (0, 1),
        'Y': (1, 1),
    }
    
    for pauli in pauli_op:
        x_bit, z_bit = pauli_map[pauli]
        
        if x_bit == 0:
            x_polys.append(LaurentPolynomial.zero())
        else:
            x_polys.append(LaurentPolynomial.one())
            
        if z_bit == 0:
            z_polys.append(LaurentPolynomial.zero())
        else:
            z_polys.append(LaurentPolynomial.one())
    
    return x_polys + z_polys


def polynomial_symplectic_inner_product(v1: List[LaurentPolynomial], 
                                       v2: List[LaurentPolynomial]) -> LaurentPolynomial:
    """
    Compute symplectic inner product of two polynomial vectors.
    
    For vectors v1 = [x1|z1] and v2 = [x2|z2], the symplectic inner product is:
    <v1, v2> = ∑ᵢ (x1ᵢ * z2ᵢ + z1ᵢ * x2ᵢ)
    
    Args:
        v1, v2: Polynomial symplectic vectors of length 2L
        
    Returns:
        Polynomial inner product
    """
    if len(v1) != len(v2) or len(v1) % 2 != 0:
        raise ValueError("Vectors must have equal even length")
    
    L = len(v1) // 2
    
    # Split into X and Z parts
    x1, z1 = v1[:L], v1[L:]
    x2, z2 = v2[:L], v2[L:]
    
    result = LaurentPolynomial.zero()
    
    for i in range(L):
        # x1[i] * z2[i] + z1[i] * x2[i]
        term1 = x1[i] * z2[i]
        term2 = z1[i] * x2[i]
        result = result + term1 + term2
    
    return result


def build_symplectic_form_poly(L: int) -> List[List[LaurentPolynomial]]:
    """
    Build symplectic form matrix J = [[0, I], [I, 0]] for polynomial vectors.
    
    Args:
        L: Half the size of symplectic vectors
        
    Returns:
        2L × 2L symplectic form matrix
    """
    # Create 2L × 2L matrix
    J = []
    for i in range(2 * L):
        row = []
        for j in range(2 * L):
            if i < L and j >= L and i == j - L:
                # Upper right block: I
                row.append(LaurentPolynomial.one())
            elif i >= L and j < L and i - L == j:
                # Lower left block: I
                row.append(LaurentPolynomial.one())
            else:
                # Zero elsewhere
                row.append(LaurentPolynomial.zero())
        J.append(row)
    
    return J


def matrix_multiply_poly(A: List[List[LaurentPolynomial]], 
                        B: List[List[LaurentPolynomial]]) -> List[List[LaurentPolynomial]]:
    """
    Multiply two polynomial matrices A × B.
    
    Args:
        A: m × n matrix
        B: n × p matrix
        
    Returns:
        m × p product matrix
    """
    if not A or not A[0] or not B or not B[0]:
        return []
    
    m, n = len(A), len(A[0])
    n2, p = len(B), len(B[0])
    
    if n != n2:
        raise ValueError(f"Matrix dimension mismatch: {n} != {n2}")
    
    result = []
    for i in range(m):
        row = []
        for j in range(p):
            # Compute A[i,:] · B[:,j]
            element = LaurentPolynomial.zero()
            for k in range(n):
                element = element + A[i][k] * B[k][j]
            row.append(element)
        result.append(row)
    
    return result 