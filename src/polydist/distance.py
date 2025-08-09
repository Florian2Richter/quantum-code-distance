"""
Distance calculation for stabilizer codes using polynomial formalism.

This module implements code distance calculation using logical operators
extracted via the Laurent polynomial extended Euclidean algorithm.
"""

import itertools
from typing import List, Tuple, Dict, Any, Optional
from .logical_operators import (
    extract_complete_logical_operators, 
    analyze_logical_operator_extraction,
    extract_logical_operators_for_distance,
    LaurentPolynomialGF2
)
from .gcd import compute_logical_qubits


def laurent_polynomial_to_pauli_string(poly: LaurentPolynomialGF2, N: int, pauli_type: str = 'Z') -> str:
    """
    Convert a Laurent polynomial to its Pauli string representation.
    
    The polynomial represents positions where the specified Pauli operator appears.
    For example, if poly = 1 + x^2, then the Pauli string has the operator at
    positions 0 and 2.
    
    Args:
        poly: Laurent polynomial in quotient ring Z[x]/(x^N - 1)
        N: Length of the Pauli string (number of qubits)
        pauli_type: Type of Pauli operator ('X', 'Y', 'Z', or 'I')
        
    Returns:
        Pauli string representation
    """
    pauli_chars = ['I'] * N
    
    for exp, coeff in poly.coeffs.items():
        if coeff == 1:  # Non-zero coefficient in GF(2)
            pos = exp % N  # Ensure position is within [0, N-1]
            pauli_chars[pos] = pauli_type
    
    return ''.join(pauli_chars)


def logical_operator_to_pauli_string(logical_op: Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2], 
                                   N: int) -> str:
    """
    Convert a logical operator (v, u) to its full Pauli string representation.
    
    The logical operator is represented as (v, u) where:
    - v corresponds to Z positions 
    - u corresponds to X positions
    
    Args:
        logical_op: Tuple (v, u) of Laurent polynomials
        N: Number of qubits
        
    Returns:
        Pauli string where each position shows the combined X/Z operator
    """
    v_poly, u_poly = logical_op
    
    pauli_chars = []
    for i in range(N):
        # Check if X operator (u polynomial) has coefficient at position i
        has_x = i in u_poly.coeffs and u_poly.coeffs[i] == 1
        # Check if Z operator (v polynomial) has coefficient at position i  
        has_z = i in v_poly.coeffs and v_poly.coeffs[i] == 1
        
        if has_x and has_z:
            pauli_chars.append('Y')  # X and Z = Y
        elif has_x:
            pauli_chars.append('X')
        elif has_z:
            pauli_chars.append('Z')
        else:
            pauli_chars.append('I')
    
    return ''.join(pauli_chars)


def compute_pauli_weight(pauli_string: str) -> int:
    """
    Compute the weight of a Pauli operator (number of non-identity positions).
    
    Args:
        pauli_string: Pauli string (e.g., "XZYI")
        
    Returns:
        Number of non-I positions
    """
    return sum(1 for char in pauli_string if char != 'I')


def combine_logical_operators(ops: List[Tuple[LaurentPolynomialGF2, LaurentPolynomialGF2]], 
                            indices: List[int], N: int) -> str:
    """
    Combine multiple logical operators by addition in GF(2).
    
    Args:
        ops: List of logical operators (v, u) pairs
        indices: Indices of operators to combine
        N: Number of qubits
        
    Returns:
        Combined Pauli string
    """
    if not indices:
        return 'I' * N
    
    # Start with the first operator
    combined_v = LaurentPolynomialGF2(ops[indices[0]][0].coeffs.copy(), N)
    combined_u = LaurentPolynomialGF2(ops[indices[0]][1].coeffs.copy(), N)
    
    # Add remaining operators
    for idx in indices[1:]:
        v_poly, u_poly = ops[idx]
        combined_v = combined_v + v_poly  # Addition in GF(2)
        combined_u = combined_u + u_poly
    
    return logical_operator_to_pauli_string((combined_v, combined_u), N)


def find_distance_polynomial(seed: str, verbose: bool = False) -> Tuple[int, List[str]]:
    """
    Compute code distance using polynomial-extracted logical operators.
    
    Algorithm:
    1. Extract logical operators using polynomial extended Euclidean algorithm
    2. Generate all k logical qubit operators via polynomial shifts
    3. Enumerate all non-empty combinations of logical operators
    4. Find the combination with minimum Pauli weight
    
    Args:
        seed: Pauli seed string (e.g., "XZIY")
        verbose: Whether to include detailed computation info
        
    Returns:
        Tuple of (distance, logical_operators_used) where logical_operators_used 
        is a list of Pauli strings representing the logical operators
    """
    N = len(seed)
    
    # Extract logical operators as polynomial objects
    success, k, logical_operators = extract_logical_operators_for_distance(seed)
    
    if not success:
        if verbose:
            print(f"Failed to extract logical operators for seed: {seed}")
        return 0, []
    
    if k == 0:
        return 0, []  # No logical qubits means distance is 0
    
    # Find minimum weight combination
    min_distance = N  # Upper bound
    best_operators = []
    
    num_ops = len(logical_operators)
    
    if verbose:
        print(f"Found {num_ops} logical operators for {k} logical qubits")
        for i, (v, u) in enumerate(logical_operators):
            pauli_str = logical_operator_to_pauli_string((v, u), N)
            print(f"  Operator {i}: {pauli_str} (weight: {compute_pauli_weight(pauli_str)})")
    
    # Enumerate all non-empty subsets of logical operators
    for r in range(1, num_ops + 1):
        for combo in itertools.combinations(range(num_ops), r):
            # Combine the selected logical operators
            combined_pauli = combine_logical_operators(logical_operators, list(combo), N)
            weight = compute_pauli_weight(combined_pauli)
            
            if verbose and weight <= min_distance:
                print(f"  Combination {combo}: {combined_pauli} (weight: {weight})")
            
            if weight < min_distance:
                min_distance = weight
                best_operators = [logical_operator_to_pauli_string(logical_operators[i], N) for i in combo]
                
                # Early termination if distance 1 found
                if min_distance == 1:
                    break
        
        if min_distance == 1:
            break
    
    if verbose:
        print(f"Minimum distance found: {min_distance}")
    
    return min_distance, best_operators


def verify_distance_consistency(seed: str, verbose: bool = False) -> Dict[str, Any]:
    """
    Verify consistency between polynomial and GF(2) distance calculations.
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed comparison
        
    Returns:
        Dictionary with comparison results
    """
    try:
        # Polynomial-based calculation
        poly_distance, poly_ops = find_distance_polynomial(seed, verbose=verbose)
        
        # GF(2) calculation for comparison (if available)
        try:
            import sys
            import os
            sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
            from stabilizer.lattice import build_lattice
            from stabilizer.tableau import build_tableau
            from stabilizer.distance import find_distance
            
            pauli_list = list(seed)
            stab_ops = build_lattice(pauli_list)
            tableau = build_tableau(stab_ops)
            gf2_distance = find_distance(tableau)
            
            consistency_check = (poly_distance == gf2_distance)
            
        except ImportError:
            gf2_distance = None
            consistency_check = None
        
        result = {
            'seed': seed,
            'polynomial_distance': poly_distance,
            'gf2_distance': gf2_distance,
            'distances_match': consistency_check,
            'polynomial_logical_ops': poly_ops,
            'comparison_successful': True
        }
        
        if verbose:
            result['polynomial_extraction_details'] = extract_complete_logical_operators(seed, verbose=True)
        
        return result
        
    except Exception as e:
        return {
            'seed': seed,
            'polynomial_distance': None,
            'gf2_distance': None,
            'distances_match': False,
            'comparison_successful': False,
            'error': str(e)
        }


def analyze_distance_calculation(seed: str, verbose: bool = False) -> Dict[str, Any]:
    """
    Complete analysis of distance calculation using polynomial formalism.
    
    Args:
        seed: Pauli seed string
        verbose: Whether to include detailed intermediate results
        
    Returns:
        Dictionary with complete distance analysis
    """
    N = len(seed)
    k = compute_logical_qubits(seed)
    
    result = {
        'seed': seed,
        'N': N,
        'k_logical_qubits': k,
        'analysis_successful': False,
        'distance': None,
        'logical_operators': [],
        'extraction_verified': False,
        'consistency_verified': False
    }
    
    try:
        # Step 1: Compute distance using polynomial method
        distance, logical_ops = find_distance_polynomial(seed, verbose=verbose)
        
        result.update({
            'analysis_successful': True,
            'distance': distance,
            'logical_operators': logical_ops,
            'num_logical_operators': len(logical_ops)
        })
        
        # Step 2: Verify logical operator extraction
        extraction_analysis = analyze_logical_operator_extraction(seed, verbose=verbose)
        result['extraction_verified'] = extraction_analysis['is_verified']
        
        # Step 3: Verify consistency with GF(2) if available
        consistency_check = verify_distance_consistency(seed, verbose=verbose)
        result['consistency_verified'] = consistency_check['distances_match']
        result['gf2_distance'] = consistency_check['gf2_distance']
        
        if verbose:
            result['verbose_details'] = {
                'extraction_analysis': extraction_analysis,
                'consistency_check': consistency_check,
                'logical_operator_weights': [compute_pauli_weight(op) for op in logical_ops]
            }
        
        return result
        
    except Exception as e:
        result.update({
            'analysis_successful': False,
            'error': str(e)
        })
        return result


def batch_distance_analysis(seeds: List[str], verbose: bool = False) -> Dict[str, Any]:
    """
    Analyze distance calculation for multiple seeds.
    
    Args:
        seeds: List of Pauli seed strings
        verbose: Whether to include detailed results for each seed
        
    Returns:
        Dictionary with batch analysis results
    """
    results = []
    summary = {
        'total_seeds': len(seeds),
        'successful_analyses': 0,
        'verified_extractions': 0,
        'consistent_with_gf2': 0,
        'distance_distribution': {}
    }
    
    for seed in seeds:
        analysis = analyze_distance_calculation(seed, verbose=verbose)
        results.append(analysis)
        
        if analysis['analysis_successful']:
            summary['successful_analyses'] += 1
            
            distance = analysis['distance']
            if distance is not None:
                summary['distance_distribution'][distance] = summary['distance_distribution'].get(distance, 0) + 1
        
        if analysis['extraction_verified']:
            summary['verified_extractions'] += 1
            
        if analysis['consistency_verified']:
            summary['consistent_with_gf2'] += 1
    
    return {
        'individual_results': results,
        'summary': summary,
        'batch_successful': summary['successful_analyses'] == summary['total_seeds']
    } 