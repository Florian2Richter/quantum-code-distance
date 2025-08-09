"""
GCD computation for Laurent polynomials over GF(2).

This module implements polynomial GCD to compute the number of logical qubits k
using the formula: k = deg(gcd(s_X(x), s_Z(x), x^N - 1))
"""

from typing import List, Dict
from .orthogonality import LaurentPolynomialGF2


class PolynomialGF2:
    """
    Polynomial over GF(2) WITHOUT automatic modular reduction.
    Used for GCD computation in the full polynomial ring.
    """
    
    def __init__(self, coeffs: Dict[int, int]):
        """
        Initialize polynomial.
        
        Args:
            coeffs: Dictionary mapping exponents to coefficients {exp: coeff}
        """
        # Normalize coefficients modulo 2 but DON'T reduce exponents
        self.coeffs = {}
        for exp, coeff in coeffs.items():
            if coeff % 2 != 0:  # Only keep non-zero coefficients in GF(2)
                self.coeffs[exp] = 1
    
    @classmethod
    def from_laurent(cls, laurent_poly: LaurentPolynomialGF2):
        """Convert from LaurentPolynomialGF2 to PolynomialGF2."""
        return cls(laurent_poly.coeffs.copy())
    
    def to_laurent(self, N: int) -> LaurentPolynomialGF2:
        """Convert back to LaurentPolynomialGF2 with given ring size."""
        return LaurentPolynomialGF2(self.coeffs, N)
    
    def __add__(self, other):
        """Add two polynomials (XOR in GF(2))."""
        result_coeffs = self.coeffs.copy()
        for exp, coeff in other.coeffs.items():
            if exp in result_coeffs:
                # XOR in GF(2): 1 + 1 = 0, 1 + 0 = 1, 0 + 1 = 1
                result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                if result_coeffs[exp] == 0:
                    del result_coeffs[exp]
            else:
                result_coeffs[exp] = coeff
        
        return PolynomialGF2(result_coeffs)
    
    def __mul__(self, other):
        """Multiply two polynomials WITHOUT modular reduction."""
        result_coeffs = {}
        for exp1, coeff1 in self.coeffs.items():
            for exp2, coeff2 in other.coeffs.items():
                exp = exp1 + exp2  # NO modular reduction
                coeff = (coeff1 * coeff2) % 2  # Multiply in GF(2)
                
                if exp in result_coeffs:
                    result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                    if result_coeffs[exp] == 0:
                        del result_coeffs[exp]
                else:
                    if coeff != 0:
                        result_coeffs[exp] = coeff
        
        return PolynomialGF2(result_coeffs)
    
    def is_zero(self):
        """Check if polynomial is zero."""
        return len(self.coeffs) == 0
    
    def degree(self):
        """Get the degree of the polynomial."""
        if self.is_zero():
            return -1
        return max(self.coeffs.keys())
    
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


def polynomial_gcd_full_ring(poly1: PolynomialGF2, poly2: PolynomialGF2) -> PolynomialGF2:
    """
    Compute GCD of two polynomials over GF(2) using Euclidean algorithm.
    Works in the full polynomial ring, not the quotient ring.
    
    Args:
        poly1, poly2: Polynomials in the full ring
        
    Returns:
        GCD polynomial
    """"""
GCD computation for Laurent polynomials over GF(2).

This module implements polynomial GCD to compute the number of logical qubits k
using the formula: k = deg(gcd(s_X(x), s_Z(x), x^N - 1))
"""

from typing import List, Dict
from .orthogonality import LaurentPolynomialGF2


class PolynomialGF2:
    """
    Polynomial over GF(2) WITHOUT automatic modular reduction.
    Used for GCD computation in the full polynomial ring.
    """
    
    def __init__(self, coeffs: Dict[int, int]):
        """
        Initialize polynomial.
        
        Args:
            coeffs: Dictionary mapping exponents to coefficients {exp: coeff}
        """
        # Normalize coefficients modulo 2 but DON'T reduce exponents
        self.coeffs = {}
        for exp, coeff in coeffs.items():
            if coeff % 2 != 0:  # Only keep non-zero coefficients in GF(2)
                self.coeffs[exp] = 1
    
    @classmethod
    def from_laurent(cls, laurent_poly: LaurentPolynomialGF2):
        """Convert from LaurentPolynomialGF2 to PolynomialGF2."""
        return cls(laurent_poly.coeffs.copy())
    
    def to_laurent(self, N: int) -> LaurentPolynomialGF2:
        """Convert back to LaurentPolynomialGF2 with given ring size."""
        return LaurentPolynomialGF2(self.coeffs, N)
    
    def __add__(self, other):
        """Add two polynomials (XOR in GF(2))."""
        result_coeffs = self.coeffs.copy()
        for exp, coeff in other.coeffs.items():
            if exp in result_coeffs:
                # XOR in GF(2): 1 + 1 = 0, 1 + 0 = 1, 0 + 1 = 1
                result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                if result_coeffs[exp] == 0:
                    del result_coeffs[exp]
            else:
                result_coeffs[exp] = coeff
        
        return PolynomialGF2(result_coeffs)
    
    def __mul__(self, other):
        """Multiply two polynomials WITHOUT modular reduction."""
        result_coeffs = {}
        for exp1, coeff1 in self.coeffs.items():
            for exp2, coeff2 in other.coeffs.items():
                exp = exp1 + exp2  # NO modular reduction
                coeff = (coeff1 * coeff2) % 2  # Multiply in GF(2)
                
                if exp in result_coeffs:
                    result_coeffs[exp] = (result_coeffs[exp] + coeff) % 2
                    if result_coeffs[exp] == 0:
                        del result_coeffs[exp]
                else:
                    if coeff != 0:
                        result_coeffs[exp] = coeff
        
        return PolynomialGF2(result_coeffs)
    
    def is_zero(self):
        """Check if polynomial is zero."""
        return len(self.coeffs) == 0
    
    def degree(self):
        """Get the degree of the polynomial."""
        if self.is_zero():
            return -1
        return max(self.coeffs.keys())
    
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


def polynomial_gcd_full_ring(poly1: PolynomialGF2, poly2: PolynomialGF2) -> PolynomialGF2:
    """
    Compute GCD of two polynomials over GF(2) using Euclidean algorithm.
    Works in the full polynomial ring, not the quotient ring.
    
    Args:
        poly1, poly2: Polynomials in the full ring
        
    Returns:
        GCD polynomial
    """
    # Handle zero polynomials
    if poly1.is_zero():
        return poly2
    if poly2.is_zero():
        return poly1
    
    # Ensure we work with copies
    a = PolynomialGF2(poly1.coeffs.copy())
    b = PolynomialGF2(poly2.coeffs.copy())
    
    # Euclidean algorithm: while b != 0, (a, b) = (b, a mod b)
    while not b.is_zero():
        remainder = polynomial_mod_full_ring(a, b)
        a = b
        b = remainder
    
    return a


def polynomial_mod_full_ring(dividend: PolynomialGF2, divisor: PolynomialGF2) -> PolynomialGF2:
    """
    Compute polynomial division remainder: dividend mod divisor.
    Works in the full polynomial ring GF(2)[x], not the quotient ring.
    
    Args:
        dividend: Polynomial to be divided
        divisor: Polynomial to divide by
        
    Returns:
        Remainder polynomial
    """
    if divisor.is_zero():
        raise ValueError("Cannot divide by zero polynomial")
    
    if dividend.is_zero():
        return PolynomialGF2({})
    
    # Work with copies
    remainder = PolynomialGF2(dividend.coeffs.copy())
    
    while not remainder.is_zero():
        # Get leading terms
        remainder_deg = remainder.degree()
        divisor_deg = divisor.degree()
        
        # If remainder degree < divisor degree, we're done
        if remainder_deg < divisor_deg:
            break
        
        # Compute quotient term: x^(remainder_deg - divisor_deg)
        quotient_exp = remainder_deg - divisor_deg
        quotient_term = PolynomialGF2({quotient_exp: 1})
        
        # Subtract divisor * quotient_term from remainder
        to_subtract = divisor * quotient_term
        remainder = remainder + to_subtract  # Addition = subtraction in GF(2)
    
    return remainder


# Old helper functions removed - now using methods on PolynomialGF2 class


def create_x_power_minus_one(N: int) -> PolynomialGF2:
    """
    Create the polynomial x^N - 1 = x^N + 1 (since -1 = 1 in GF(2)).
    
    Args:
        N: Ring size
        
    Returns:
        Polynomial representing x^N + 1 with degree N
    """
    # In GF(2) we have x^N - 1 = x^N + 1
    # We need the actual polynomial of degree N with coefficients {N: 1, 0: 1}
    coeffs = {N: 1, 0: 1}
    
    return PolynomialGF2(coeffs)


def gcd_three_polynomials(s_X: PolynomialGF2, s_Z: PolynomialGF2, 
                         x_n_minus_1: PolynomialGF2) -> PolynomialGF2:
    """
    Compute GCD of three polynomials: gcd(s_X(x), s_Z(x), x^N - 1).
    
    Args:
        s_X: X part polynomial
        s_Z: Z part polynomial  
        x_n_minus_1: The polynomial x^N - 1
        
    Returns:
        GCD polynomial
    """
    # gcd(a, b, c) = gcd(gcd(a, b), c)
    gcd_xy = polynomial_gcd_full_ring(s_X, s_Z)
    return polynomial_gcd_full_ring(gcd_xy, x_n_minus_1)


def compute_logical_qubits(seed: str) -> int:
    """
    Compute the number of logical qubits k using polynomial GCD.
    
    The formula is: k = deg(gcd(s_X(x), s_Z(x), x^N - 1))
    
    Args:
        seed: Pauli seed string (e.g., "XZIY")
        
    Returns:
        Number of logical qubits k
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract X and Z parts as Laurent polynomials, then convert to full ring polynomials
    s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    # Convert to full ring polynomials for GCD computation
    s_X = PolynomialGF2.from_laurent(s_X_laurent)
    s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
    
    # Create x^N - 1 polynomial in the full ring
    x_n_minus_1 = create_x_power_minus_one(N)
    
    # Compute GCD in the full ring
    gcd_poly = gcd_three_polynomials(s_X, s_Z, x_n_minus_1)
    
    # Return degree
    return gcd_poly.degree()


def analyze_gcd_computation(seed: str, verbose: bool = False) -> Dict[str, any]:
    """
    Detailed analysis of GCD computation with intermediate results.
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed intermediate polynomials
        
    Returns:
        Dictionary with analysis results
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract Laurent polynomials and convert to full ring
    s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    s_X = PolynomialGF2.from_laurent(s_X_laurent)
    s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
    x_n_minus_1 = create_x_power_minus_one(N)
    
    # Compute intermediate GCDs in full ring
    gcd_xy = polynomial_gcd_full_ring(s_X, s_Z)
    gcd_final = polynomial_gcd_full_ring(gcd_xy, x_n_minus_1)
    
    k = gcd_final.degree()
    
    result = {
        'seed': seed,
        'N': N,
        'k_logical': k,
        'gcd_polynomial': str(gcd_final),
        'gcd_degree': k
    }
    
    if verbose:
        result.update({
            's_X': str(s_X),
            's_Z': str(s_Z),
            'x_n_minus_1': str(x_n_minus_1),
            'gcd_sX_sZ': str(gcd_xy),
            's_X_degree': s_X.degree(),
            's_Z_degree': s_Z.degree(),
            'gcd_xy_degree': gcd_xy.degree()
        })
    
    return result 
    # Handle zero polynomials
    if poly1.is_zero():
        return poly2
    if poly2.is_zero():
        return poly1
    
    # Ensure we work with copies
    a = PolynomialGF2(poly1.coeffs.copy())
    b = PolynomialGF2(poly2.coeffs.copy())
    
    # Euclidean algorithm: while b != 0, (a, b) = (b, a mod b)
    while not b.is_zero():
        remainder = polynomial_mod_full_ring(a, b)
        a = b
        b = remainder
    
    return a


def polynomial_mod_full_ring(dividend: PolynomialGF2, divisor: PolynomialGF2) -> PolynomialGF2:
    """
    Compute polynomial division remainder: dividend mod divisor.
    Works in the full polynomial ring GF(2)[x], not the quotient ring.
    
    Args:
        dividend: Polynomial to be divided
        divisor: Polynomial to divide by
        
    Returns:
        Remainder polynomial
    """
    if divisor.is_zero():
        raise ValueError("Cannot divide by zero polynomial")
    
    if dividend.is_zero():
        return PolynomialGF2({})
    
    # Work with copies
    remainder = PolynomialGF2(dividend.coeffs.copy())
    
    while not remainder.is_zero():
        # Get leading terms
        remainder_deg = remainder.degree()
        divisor_deg = divisor.degree()
        
        # If remainder degree < divisor degree, we're done
        if remainder_deg < divisor_deg:
            break
        
        # Compute quotient term: x^(remainder_deg - divisor_deg)
        quotient_exp = remainder_deg - divisor_deg
        quotient_term = PolynomialGF2({quotient_exp: 1})
        
        # Subtract divisor * quotient_term from remainder
        to_subtract = divisor * quotient_term
        remainder = remainder + to_subtract  # Addition = subtraction in GF(2)
    
    return remainder


# Old helper functions removed - now using methods on PolynomialGF2 class


def create_x_power_minus_one(N: int) -> PolynomialGF2:
    """
    Create the polynomial x^N - 1 = x^N + 1 (since -1 = 1 in GF(2)).
    
    Args:
        N: Ring size
        
    Returns:
        Polynomial representing x^N + 1 with degree N
    """
    # In GF(2) we have x^N - 1 = x^N + 1
    # We need the actual polynomial of degree N with coefficients {N: 1, 0: 1}
    coeffs = {N: 1, 0: 1}
    
    return PolynomialGF2(coeffs)


def gcd_three_polynomials(s_X: PolynomialGF2, s_Z: PolynomialGF2, 
                         x_n_minus_1: PolynomialGF2) -> PolynomialGF2:
    """
    Compute GCD of three polynomials: gcd(s_X(x), s_Z(x), x^N - 1).
    
    Args:
        s_X: X part polynomial
        s_Z: Z part polynomial  
        x_n_minus_1: The polynomial x^N - 1
        
    Returns:
        GCD polynomial
    """
    # gcd(a, b, c) = gcd(gcd(a, b), c)
    gcd_xy = polynomial_gcd_full_ring(s_X, s_Z)
    return polynomial_gcd_full_ring(gcd_xy, x_n_minus_1)


def compute_logical_qubits(seed: str) -> int:
    """
    Compute the number of logical qubits k using polynomial GCD.
    
    The formula is: k = deg(gcd(s_X(x), s_Z(x), x^N - 1))
    
    Args:
        seed: Pauli seed string (e.g., "XZIY")
        
    Returns:
        Number of logical qubits k
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract X and Z parts as Laurent polynomials, then convert to full ring polynomials
    s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    # Convert to full ring polynomials for GCD computation
    s_X = PolynomialGF2.from_laurent(s_X_laurent)
    s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
    
    # Create x^N - 1 polynomial in the full ring
    x_n_minus_1 = create_x_power_minus_one(N)
    
    # Compute GCD in the full ring
    gcd_poly = gcd_three_polynomials(s_X, s_Z, x_n_minus_1)
    
    # Return degree
    return gcd_poly.degree()


def analyze_gcd_computation(seed: str, verbose: bool = False) -> Dict[str, any]:
    """
    Detailed analysis of GCD computation with intermediate results.
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed intermediate polynomials
        
    Returns:
        Dictionary with analysis results
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract Laurent polynomials and convert to full ring
    s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    s_X = PolynomialGF2.from_laurent(s_X_laurent)
    s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
    x_n_minus_1 = create_x_power_minus_one(N)
    
    # Compute intermediate GCDs in full ring
    gcd_xy = polynomial_gcd_full_ring(s_X, s_Z)
    gcd_final = polynomial_gcd_full_ring(gcd_xy, x_n_minus_1)
    
    k = gcd_final.degree()
    
    result = {
        'seed': seed,
        'N': N,
        'k_logical': k,
        'gcd_polynomial': str(gcd_final),
        'gcd_degree': k
    }
    
    if verbose:
        result.update({
            's_X': str(s_X),
            's_Z': str(s_Z),
            'x_n_minus_1': str(x_n_minus_1),
            'gcd_sX_sZ': str(gcd_xy),
            's_X_degree': s_X.degree(),
            's_Z_degree': s_Z.degree(),
            'gcd_xy_degree': gcd_xy.degree()
        })
    
    return result 