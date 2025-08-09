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


def reverse_polynomial(poly: LaurentPolynomialGF2, N: int) -> LaurentPolynomialGF2:
    """
    Compute polynomial reversal: p_rev(x) = x^(deg p) * p(x^-1).
    
    In the quotient ring modulo x^N - 1, this flips the exponents around
    the center and adjusts by the degree.
    
    Args:
        poly: Input Laurent polynomial
        N: Ring size for quotient ring Z[x]/(x^N - 1)
        
    Returns:
        Reversed polynomial
    """
    if poly.is_zero():
        return LaurentPolynomialGF2({}, N)
    
    # Find the actual degree considering the quotient ring
    max_exp = max(poly.coeffs.keys())
    
    # For reversal, we map exponent i to (max_exp - i) mod N
    reversed_coeffs = {}
    for exp, coeff in poly.coeffs.items():
        if coeff == 1:  # Only keep non-zero coefficients
            new_exp = (max_exp - exp) % N
            reversed_coeffs[new_exp] = 1
    
    return LaurentPolynomialGF2(reversed_coeffs, N)


def extract_logical_x_operator(v: LaurentPolynomialGF2, u: LaurentPolynomialGF2, N: int) -> Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2]:
    """
    Extract logical X operator that satisfies symplectic orthogonality with logical Z.
    
    Algorithm:
    1. Reverse v and u: vrev(x) = x^(deg v) * v(x^-1), urev(x) = x^(deg u) * u(x^-1)
    2. Solve r'(x) * vrev(x) + t'(x) * urev(x) = 1 (mod x^N - 1) using extended GCD
    3. Reverse back: r(x) = t'(x^-1), t(x) = r'(x^-1)
    4. Return LX(x) = (r(x), t(x))
    
    This ensures LZ^T * J * LX ≡ 1 while both commute with stabilizers.
    
    Args:
        v, u: Logical Z operator components (v, u)
        N: Ring size for quotient ring Z[x]/(x^N - 1)
        
    Returns:
        Tuple (r, t) representing logical X operator LX = (r, t)
    """
    # Step 1: Reverse v and u polynomials
    vrev = reverse_polynomial(v, N)
    urev = reverse_polynomial(u, N)
    
    # Step 2: Convert to full polynomials for extended GCD
    vrev_full = PolynomialGF2.from_laurent(vrev)
    urev_full = PolynomialGF2.from_laurent(urev)
    
    # Create target polynomial: 1
    target = PolynomialGF2({0: 1})
    
    # We need to solve r' * vrev + t' * urev = 1 (mod x^N - 1)
    # This is equivalent to finding the extended GCD coefficients
    # where gcd(vrev, urev) should divide 1, so gcd should be 1
    
    # First, check if vrev and urev are coprime modulo x^N - 1
    x_n_minus_1 = create_x_power_minus_one(N)
    
    # Use extended GCD in the quotient ring
    # We solve: r' * vrev + t' * urev ≡ gcd(vrev, urev) (mod x^N - 1)
    gcd_poly, r_prime_poly, t_prime_poly = extended_gcd_polynomials(vrev_full, urev_full)
    
    # If gcd is not 1, we need to work in the quotient ring
    # Compute gcd with x^N - 1 to ensure we can find inverse
    final_gcd, alpha, beta = extended_gcd_polynomials(gcd_poly, x_n_minus_1)
    
    # Combine to get final coefficients: multiply by alpha to get unit
    r_prime_final = (alpha * r_prime_poly)
    t_prime_final = (alpha * t_prime_poly)
    
    # Convert back to Laurent polynomials in quotient ring
    r_prime_laurent = r_prime_final.to_laurent(N)
    t_prime_laurent = t_prime_final.to_laurent(N)
    
    # Step 3: Reverse back to get final r and t
    # r(x) = t'(x^-1), t(x) = r'(x^-1)
    r = reverse_polynomial(t_prime_laurent, N)
    t = reverse_polynomial(r_prime_laurent, N)
    
    return r, t


def generate_logical_qubit_operators(LX: Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2], 
                                   LZ: Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2], 
                                   k: int, N: int) -> Dict[str, List[Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2]]]:
    """
    Generate k independent logical qubit operators using polynomial shifts.
    
    The shifts Xi = LX(x) * x^i and Zi = LZ(x) * x^i for i = 0, ..., k-1
    give k independent logical qubit operators.
    
    Args:
        LX: Logical X operator (r, t)
        LZ: Logical Z operator (v, u) 
        k: Number of logical qubits
        N: Ring size for quotient ring Z[x]/(x^N - 1)
        
    Returns:
        Dictionary with 'X_operators' and 'Z_operators' lists of k operators each
    """
    r, t = LX
    v, u = LZ
    
    X_operators = []
    Z_operators = []
    
    for i in range(k):
        # Create shift polynomial x^i
        shift_coeffs = {i % N: 1}
        shift = LaurentPolynomialGF2(shift_coeffs, N)
        
        # Apply shift: Xi = LX * x^i = (r * x^i, t * x^i)
        Xi_r = r * shift
        Xi_t = t * shift
        X_operators.append((Xi_r, Xi_t))
        
        # Apply shift: Zi = LZ * x^i = (v * x^i, u * x^i)  
        Zi_v = v * shift
        Zi_u = u * shift
        Z_operators.append((Zi_v, Zi_u))
    
    return {
        'X_operators': X_operators,
        'Z_operators': Z_operators
    }


def verify_symplectic_orthogonality(LX: Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2],
                                   LZ: Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2],
                                   N: int) -> bool:
    """
    Verify that LZ^T * J * LX ≡ 1 (mod x^N - 1) where J is the symplectic form.
    
    For Pauli operators, the symplectic form gives:
    LZ^T * J * LX = v*t + u*r (in the quotient ring)
    
    Args:
        LX: Logical X operator (r, t)
        LZ: Logical Z operator (v, u)
        N: Ring size
        
    Returns:
        True if symplectic orthogonality holds
    """
    r, t = LX
    v, u = LZ
    
    # Compute symplectic product: v*t + u*r
    term1 = v * t
    term2 = u * r
    symplectic_product = term1 + term2
    
    # Should equal 1 in the quotient ring
    expected = LaurentPolynomialGF2({0: 1}, N)
    difference = symplectic_product + expected  # Addition = subtraction in GF(2)
    
    return difference.is_zero()


def extract_complete_logical_operators(seed: str, verbose: bool = False) -> Dict[str, any]:
    """
    Complete logical operator extraction: Z operators, matching X operators, and k logical qubits.
    
    Implements the full algorithm:
    1. Extract logical Z operator using extended Euclidean algorithm
    2. Build matching logical X operator using symplectic form
    3. Spread into k logical qubits using polynomial shifts
    4. Verify all commutation relations
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed intermediate results
        
    Returns:
        Dictionary with complete logical operator analysis
    """
    N = len(seed)
    
    try:
        # Step 1: Extract logical Z operator
        v, u = extract_logical_z_operator(seed)
        LZ = (v, u)
        
        # Step 2: Extract matching logical X operator
        r, t = extract_logical_x_operator(v, u, N)
        LX = (r, t)
        
        # Step 3: Verify symplectic orthogonality
        symplectic_correct = verify_symplectic_orthogonality(LX, LZ, N)
        
        # Step 4: Determine number of logical qubits
        from .gcd import compute_logical_qubits
        k = compute_logical_qubits(seed)
        
        # Step 5: Generate k logical qubit operators
        logical_qubits = generate_logical_qubit_operators(LX, LZ, k, N)
        
        # Step 6: Verify original logical Z operator
        z_verification = verify_logical_operator(seed, v, u)
        
        result = {
            'seed': seed,
            'N': N,
            'k_logical_qubits': k,
            'extraction_successful': True,
            'logical_Z_operator': (str(v), str(u)),
            'logical_X_operator': (str(r), str(t)),
            'symplectic_orthogonality_verified': symplectic_correct,
            'z_operator_verified': z_verification['is_correct'],
            'logical_qubit_operators': {
                'X_ops': [(str(Xi[0]), str(Xi[1])) for Xi in logical_qubits['X_operators']],
                'Z_ops': [(str(Zi[0]), str(Zi[1])) for Zi in logical_qubits['Z_operators']]
            }
        }
        
        if verbose:
            # Add detailed intermediate results
            pauli_list = list(seed)
            s_X_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'X', N)
            s_Z_laurent = LaurentPolynomialGF2.from_pauli_part(pauli_list, 'Z', N)
            
            vrev = reverse_polynomial(v, N)
            urev = reverse_polynomial(u, N)
            
            result['verbose_details'] = {
                's_X_polynomial': str(s_X_laurent),
                's_Z_polynomial': str(s_Z_laurent),
                'v_reversed': str(vrev),
                'u_reversed': str(urev),
                'symplectic_product': str((v * t) + (u * r)),
                'all_verification_details': z_verification
            }
        
        return result
        
    except Exception as e:
        return {
            'seed': seed,
            'N': N,
            'k_logical_qubits': 0,
            'extraction_successful': False,
            'error': str(e),
            'symplectic_orthogonality_verified': False,
            'z_operator_verified': False
        }


def extract_logical_operators_for_distance(seed: str) -> Tuple[bool, int, List[Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2]]]:
    """
    Extract logical operators as polynomial objects for distance calculation.
    
    Args:
        seed: Pauli seed string
        
    Returns:
        Tuple of (success, k, logical_operators) where logical_operators is a list
        of (v, u) polynomial pairs representing both X and Z logical operators
    """
    N = len(seed)
    
    try:
        # Step 1: Extract logical Z operator
        v, u = extract_logical_z_operator(seed)
        LZ = (v, u)
        
        # Step 2: Extract matching logical X operator  
        r, t = extract_logical_x_operator(v, u, N)
        LX = (r, t)
        
        # Step 3: Determine number of logical qubits
        from .gcd import compute_logical_qubits
        k = compute_logical_qubits(seed)
        
        if k == 0:
            return True, 0, []
        
        # Step 4: Generate k logical qubit operators
        logical_qubits = generate_logical_qubit_operators(LX, LZ, k, N)
        
        # Combine X and Z operators into single list
        all_operators = []
        all_operators.extend(logical_qubits['Z_operators'])  # Add Z operators first
        all_operators.extend(logical_qubits['X_operators'])  # Add X operators  
        
        return True, k, all_operators
        
    except Exception as e:
        return False, 0, [] 