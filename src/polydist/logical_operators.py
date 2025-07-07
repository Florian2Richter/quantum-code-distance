"""
Logical operator extraction for stabilizer codes using polynomial formalism.

This module implements the extraction of logical Z and X operators using
the extended Euclidean algorithm on Laurent polynomials over GF(2).
"""

from typing import Tuple, Dict, List
from .gcd import PolynomialGF2, create_x_power_minus_one
from .orthogonality import LaurentPolynomialGF2


def extended_gcd_polynomials(poly1: PolynomialGF2, poly2: PolynomialGF2) -> Tuple[PolynomialGF2, PolynomialGF2, PolynomialGF2]:
    """
    Extended Euclidean algorithm for polynomials over GF(2).
    
    Returns (gcd, u, v) such that u*poly1 + v*poly2 = gcd.
    
    Args:
        poly1, poly2: Input polynomials
        
    Returns:
        Tuple of (gcd, u, v) where u*poly1 + v*poly2 = gcd
    """
    if poly1.is_zero():
        return poly2, PolynomialGF2({}), PolynomialGF2({0: 1})
    if poly2.is_zero():
        return poly1, PolynomialGF2({0: 1}), PolynomialGF2({})
    
    # Initialize: 
    # old_r, r = poly1, poly2
    # old_s, s = 1, 0  (coefficients for poly1)
    # old_t, t = 0, 1  (coefficients for poly2)
    
    old_r, r = PolynomialGF2(poly1.coeffs.copy()), PolynomialGF2(poly2.coeffs.copy())
    old_s, s = PolynomialGF2({0: 1}), PolynomialGF2({})
    old_t, t = PolynomialGF2({}), PolynomialGF2({0: 1})
    
    while not r.is_zero():
        # Division: old_r = quotient * r + remainder
        quotient = polynomial_divide(old_r, r)
        
        # Update remainders: old_r, r = r, old_r - quotient * r
        new_r = old_r + quotient * r  # Addition = subtraction in GF(2)
        old_r, r = r, new_r
        
        # Update s coefficients: old_s, s = s, old_s - quotient * s  
        new_s = old_s + quotient * s
        old_s, s = s, new_s
        
        # Update t coefficients: old_t, t = t, old_t - quotient * t
        new_t = old_t + quotient * t
        old_t, t = t, new_t
    
    # old_r is the GCD, old_s and old_t are the Bézout coefficients
    return old_r, old_s, old_t


def polynomial_divide(dividend: PolynomialGF2, divisor: PolynomialGF2) -> PolynomialGF2:
    """
    Compute polynomial division quotient: dividend // divisor.
    
    Args:
        dividend: Polynomial to be divided
        divisor: Polynomial to divide by
        
    Returns:
        Quotient polynomial
    """
    if divisor.is_zero():
        raise ValueError("Cannot divide by zero polynomial")
    
    if dividend.is_zero():
        return PolynomialGF2({})
    
    quotient_coeffs = {}
    remainder = PolynomialGF2(dividend.coeffs.copy())
    
    while not remainder.is_zero():
        remainder_deg = remainder.degree()
        divisor_deg = divisor.degree()
        
        # If remainder degree < divisor degree, we're done
        if remainder_deg < divisor_deg:
            break
        
        # Quotient term: x^(remainder_deg - divisor_deg)
        quotient_exp = remainder_deg - divisor_deg
        quotient_coeffs[quotient_exp] = 1
        
        # Create quotient term and subtract divisor * quotient_term from remainder
        quotient_term = PolynomialGF2({quotient_exp: 1})
        to_subtract = divisor * quotient_term
        remainder = remainder + to_subtract  # Addition = subtraction in GF(2)
    
    return PolynomialGF2(quotient_coeffs)


def extract_logical_z_operator(seed: str) -> Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2]:
    """
    Extract logical Z operator using polynomial formalism.
    
    The algorithm:
    1. Compute (g1, u1, v1) = egcd(s_X, s_Z) such that u1*s_X + v1*s_Z = g1
    2. Compute (f, w, t) = egcd(g1, x^N - 1) such that w*g1 + t*(x^N - 1) = f  
    3. Set u = w*u1, v = w*v1, so u*s_X + v*s_Z ≡ f (mod x^N - 1)
    4. Return L_Z = (v, u) as the logical Z operator
    
    Args:
        seed: Pauli seed string (e.g., "XZIY")
        
    Returns:
        Tuple of (v(x), u(x)) representing the logical Z operator L_Z = (v, u)
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract X and Z parts as Laurent polynomials, then convert to full ring
    s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    s_X = PolynomialGF2.from_laurent(s_X_laurent)
    s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
    
    # Stage 1: Extended GCD of s_X and s_Z
    g1, u1, v1 = extended_gcd_polynomials(s_X, s_Z)
    
    # Stage 2: Extended GCD of g1 and x^N - 1
    x_n_minus_1 = create_x_power_minus_one(N)
    f, w, t = extended_gcd_polynomials(g1, x_n_minus_1)
    
    # Stage 3: Combine Bézout coefficients
    u = w * u1  # Coefficient for s_X
    v = w * v1  # Coefficient for s_Z
    
    # Convert back to Laurent polynomials in the quotient ring
    v_laurent = v.to_laurent(N)
    u_laurent = u.to_laurent(N)
    
    return v_laurent, u_laurent


def verify_logical_operator(seed: str, v: LaurentPolynomialGF2, u: LaurentPolynomialGF2) -> Dict[str, any]:
    """
    Verify that the extracted logical operator is correct.
    
    Checks that u*s_X + v*s_Z ≡ f (mod x^N - 1) where f = gcd(s_X, s_Z, x^N - 1).
    
    Args:
        seed: Original Pauli seed
        v, u: Extracted logical operator coefficients
        
    Returns:
        Dictionary with verification results
    """
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Extract original polynomials
    s_X = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
    s_Z = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
    
    # Compute u*s_X + v*s_Z
    term1 = u * s_X
    term2 = v * s_Z
    combination = term1 + term2
    
    # Compute expected f = gcd(s_X, s_Z, x^N - 1)
    from .gcd import compute_logical_qubits, gcd_three_polynomials
    
    s_X_full = PolynomialGF2.from_laurent(s_X)
    s_Z_full = PolynomialGF2.from_laurent(s_Z)
    x_n_minus_1 = create_x_power_minus_one(N)
    
    f_full = gcd_three_polynomials(s_X_full, s_Z_full, x_n_minus_1)
    f_laurent = f_full.to_laurent(N)
    
    # Check if combination equals f
    difference = combination + f_laurent  # Should be zero if they're equal
    is_correct = difference.is_zero()
    
    result = {
        'seed': seed,
        'N': N,
        'is_correct': is_correct,
        'u_s_X_plus_v_s_Z': str(combination),
        'expected_f': str(f_laurent),
        'difference': str(difference),
        'logical_qubits_k': f_full.degree(),
        'v_polynomial': str(v),
        'u_polynomial': str(u)
    }
    
    return result


def analyze_logical_operator_extraction(seed: str, verbose: bool = False) -> Dict[str, any]:
    """
    Complete analysis of logical operator extraction process.
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed intermediate results
        
    Returns:
        Dictionary with analysis results including verification
    """
    try:
        # Extract logical operator
        v, u = extract_logical_z_operator(seed)
        
        # Verify the result
        verification = verify_logical_operator(seed, v, u)
        
        result = {
            'seed': seed,
            'N': len(seed),
            'extraction_successful': True,
            'logical_Z_operator': (str(v), str(u)),
            'verification': verification,
            'is_verified': verification['is_correct']
        }
        
        if verbose:
            # Add intermediate computation details
            pauli_list = list(seed)
            N = len(pauli_list)
            
            s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
            s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
            s_X = PolynomialGF2.from_laurent(s_X_laurent)
            s_Z = PolynomialGF2.from_laurent(s_Z_laurent)
            
            g1, u1, v1 = extended_gcd_polynomials(s_X, s_Z)
            x_n_minus_1 = create_x_power_minus_one(N)
            f, w, t = extended_gcd_polynomials(g1, x_n_minus_1)
            
            result['verbose_details'] = {
                's_X': str(s_X),
                's_Z': str(s_Z),
                'stage1_g1': str(g1),
                'stage1_u1': str(u1),
                'stage1_v1': str(v1),
                'x_n_minus_1': str(x_n_minus_1),
                'stage2_f': str(f),
                'stage2_w': str(w),
                'stage2_t': str(t),
                'final_u': str(u),
                'final_v': str(v)
            }
        
        return result
        
    except Exception as e:
        return {
            'seed': seed,
            'N': len(seed),
            'extraction_successful': False,
            'error': str(e),
            'verification': None,
            'is_verified': False
        } 