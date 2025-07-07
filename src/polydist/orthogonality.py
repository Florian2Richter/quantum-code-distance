"""
Orthogonality check for stabilizer code seeds using Laurent polynomials.

This module implements the self-orthogonality condition:
s(x)^T J s(x) = s_X(x) * s_Z(x^-1) + s_Z(x) * s_X(x^-1) ≡ 0 (mod x^N - 1)

Where s(x) = (s_X(x), s_Z(x)) is the seed written as Laurent polynomials.
"""

from typing import List, Dict


class LaurentPolynomialGF2:
    """
    Laurent polynomial over GF(2) modulo (x^N - 1).
    
    Represents polynomials of the form sum_i c_i * x^i where c_i ∈ {0, 1}
    and all operations are performed modulo (x^N - 1).
    """
    
    def __init__(self, coeffs: Dict[int, int], N: int):
        """
        Initialize Laurent polynomial.
        
        Args:
            coeffs: Dictionary mapping exponents to coefficients {exp: coeff}
            N: Ring size for modular reduction (x^N - 1)
        """
        self.N = N
        # Normalize coefficients modulo 2 and reduce exponents modulo N
        self.coeffs = {}
        for exp, coeff in coeffs.items():
            if coeff % 2 != 0:  # Only keep non-zero coefficients in GF(2)
                self.coeffs[exp % N] = 1
    
    @classmethod
    def from_pauli_part(cls, pauli_list: List[str], pauli_type: str, N: int):
        """
        Create Laurent polynomial from X or Z part of Pauli string.
        
        Args:
            pauli_list: List of Pauli operators ['X', 'Z', 'I', 'Y', ...]
            pauli_type: Either 'X' or 'Z' to extract that component
            N: Ring size
            
        Returns:
            Laurent polynomial representing the X or Z part
        """
        coeffs = {}
        for i, pauli in enumerate(pauli_list):
            if pauli_type == 'X' and pauli in ['X', 'Y']:
                coeffs[i] = 1
            elif pauli_type == 'Z' and pauli in ['Z', 'Y']:
                coeffs[i] = 1
        
        return cls(coeffs, N)
    
    def __add__(self, other):
        """Add two Laurent polynomials (XOR in GF(2))."""
        if self.N != other.N:
            raise ValueError("Cannot add polynomials with different ring sizes")
        
        result_coeffs = self.coeffs.copy()
        for exp, coeff in other.coeffs.items():
            if exp in result_coeffs:
                # XOR in GF(2): 1 + 1 = 0, 1 + 0 = 1, 0 + 1 = 1
                result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                if result_coeffs[exp] == 0:
                    del result_coeffs[exp]
            else:
                result_coeffs[exp] = coeff
        
        return LaurentPolynomialGF2(result_coeffs, self.N)
    
    def __mul__(self, other):
        """Multiply two Laurent polynomials modulo (x^N - 1)."""
        if self.N != other.N:
            raise ValueError("Cannot multiply polynomials with different ring sizes")
        
        result_coeffs = {}
        for exp1, coeff1 in self.coeffs.items():
            for exp2, coeff2 in other.coeffs.items():
                exp = (exp1 + exp2) % self.N  # Reduce modulo x^N - 1
                coeff = (coeff1 * coeff2) % 2  # Multiply in GF(2)
                
                if exp in result_coeffs:
                    result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                    if result_coeffs[exp] == 0:
                        del result_coeffs[exp]
                else:
                    if coeff != 0:
                        result_coeffs[exp] = coeff
        
        return LaurentPolynomialGF2(result_coeffs, self.N)
    
    def inverse_powers(self):
        """
        Apply x^-1 transformation: x^k -> x^(N-k) mod (x^N - 1).
        
        Since x^N ≡ 1, we have x^-1 ≡ x^(N-1).
        """
        inv_coeffs = {}
        for exp, coeff in self.coeffs.items():
            # x^exp -> x^(-exp) = x^(N-exp) mod (x^N - 1)
            inv_exp = (-exp) % self.N
            inv_coeffs[inv_exp] = coeff
        
        return LaurentPolynomialGF2(inv_coeffs, self.N)
    
    def is_zero(self):
        """Check if polynomial is zero."""
        return len(self.coeffs) == 0
    
    def __str__(self):
        """String representation of the polynomial."""
        if self.is_zero():
            return "0"
        
        terms = []
        for exp in sorted(self.coeffs.keys()):
            if exp == 0:
                terms.append("1")
            elif exp == 1:
                terms.append("x")
            else:
                terms.append(f"x^{exp}")
        
        return " + ".join(terms)
    
    def __repr__(self):
        return f"LaurentPolynomialGF2({self.coeffs}, N={self.N})"


def check_seed_orthogonality(seed: str) -> bool:
    """
    Check if a Pauli seed satisfies the self-orthogonality condition.
    
    The condition is: s_X(x) * s_Z(x^-1) + s_Z(x) * s_X(x^-1) ≡ 0 (mod x^N - 1)
    
    Args:
        seed: String of Pauli operators (e.g., "XZIY")
        
    Returns:
        True if seed is self-orthogonal, False otherwise
    """
    # Convert seed string to list
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract X and Z parts as Laurent polynomials
    s_X = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    # Apply inverse transformation: s_Z(x^-1) and s_X(x^-1)
    s_Z_inv = s_Z.inverse_powers()
    s_X_inv = s_X.inverse_powers()
    
    # Compute orthogonality condition: s_X(x) * s_Z(x^-1) + s_Z(x) * s_X(x^-1)
    term1 = s_X * s_Z_inv
    term2 = s_Z * s_X_inv
    orthogonality_poly = term1 + term2
    
    # Check if result is zero polynomial
    return orthogonality_poly.is_zero()


def analyze_seed_orthogonality(seed: str, verbose: bool = False) -> Dict[str, any]:
    """
    Detailed analysis of seed orthogonality with intermediate results.
    
    Args:
        seed: String of Pauli operators
        verbose: Whether to include detailed intermediate polynomials
        
    Returns:
        Dictionary with analysis results
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract X and Z parts
    s_X = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    # Apply inverse transformation
    s_Z_inv = s_Z.inverse_powers()
    s_X_inv = s_X.inverse_powers()
    
    # Compute orthogonality terms
    term1 = s_X * s_Z_inv
    term2 = s_Z * s_X_inv
    orthogonality_poly = term1 + term2
    
    result = {
        'seed': seed,
        'N': N,
        'is_orthogonal': orthogonality_poly.is_zero(),
        'orthogonality_polynomial': str(orthogonality_poly)
    }
    
    if verbose:
        result.update({
            's_X': str(s_X),
            's_Z': str(s_Z),
            's_X_inv': str(s_X_inv),
            's_Z_inv': str(s_Z_inv),
            'term1_sX_sZ_inv': str(term1),
            'term2_sZ_sX_inv': str(term2)
        })
    
    return result 