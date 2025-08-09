#!/usr/bin/env python3
"""
Comprehensive comparison of logical operators between GF(2) and polynomial approaches.

This script extracts all logical operators from both methods and shows their
combinations to understand discrepancies in distance calculations.
"""

import sys
import itertools
sys.path.append('src')

import numpy as np
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_logical_operators, format_symplectic_vector

from polydist.distance import extract_logical_operators_for_distance, logical_operator_to_pauli_string, compute_pauli_weight, combine_logical_operators


def analyze_gf2_logical_operators(seed: str):
    """Extract and analyze all logical operators using GF(2) stabilizer tableau."""
    print(f"🔬 GF(2) STABILIZER APPROACH")
    print("=" * 50)
    
    pauli_list = list(seed)
    N = len(pauli_list)
    
    # Build lattice and tableau
    stab_ops = build_lattice(pauli_list)
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    k = N - rank
    
    print(f"Physical qubits (n): {N}")
    print(f"Stabilizer rank: {rank}")
    print(f"Logical qubits (k): {k}")
    
    if k == 0:
        print("No logical qubits found.")
        return None, []
    
    # Find all logical operators
    logical_ops = find_logical_operators(tableau)
    
    print(f"\nFound {len(logical_ops)} logical operators:")
    pauli_operators = []
    for i, vec in enumerate(logical_ops):
        pauli_str = format_symplectic_vector(vec)
        weight = sum(1 for char in pauli_str if char != 'I')
        pauli_operators.append(pauli_str)
        print(f"  L{i:2d}: {pauli_str} (weight: {weight})")
    
    return k, pauli_operators


def analyze_polynomial_logical_operators(seed: str):
    """Extract and analyze all logical operators using polynomial approach."""
    print(f"\n🧮 POLYNOMIAL APPROACH") 
    print("=" * 50)
    
    N = len(seed)
    
    # Extract logical operators
    success, k, logical_operators = extract_logical_operators_for_distance(seed)
    
    print(f"Physical qubits (n): {N}")
    print(f"Logical qubits (k): {k}")
    print(f"Extraction successful: {success}")
    
    if not success or k == 0:
        print("No logical qubits found.")
        return None, []
    
    print(f"\nFound {len(logical_operators)} logical operators:")
    pauli_operators = []
    for i, (v_poly, u_poly) in enumerate(logical_operators):
        pauli_str = logical_operator_to_pauli_string((v_poly, u_poly), N)
        weight = compute_pauli_weight(pauli_str)
        pauli_operators.append(pauli_str)
        print(f"  L{i:2d}: {pauli_str} (weight: {weight})")
    
    return k, pauli_operators


def find_minimum_weight_combinations(operators: list, method_name: str, max_combinations: int = 20):
    """Find all combinations of logical operators and their weights."""
    print(f"\n🔍 {method_name} - OPERATOR COMBINATIONS")
    print("=" * 60)
    
    if not operators:
        print("No operators to combine.")
        return None, []
    
    N = len(operators[0])
    num_ops = len(operators)
    
    all_combinations = []
    
    # Generate all non-empty combinations
    for r in range(1, num_ops + 1):
        for combo_indices in itertools.combinations(range(num_ops), r):
            # Combine operators using XOR (GF(2) addition)
            if method_name == "POLYNOMIAL":
                # For polynomial, we need special combination logic
                # This is a simplified version - we'll combine Pauli strings directly
                combined_pauli = combine_pauli_strings([operators[i] for i in combo_indices])
            else:
                # For GF(2), combine Pauli strings using XOR
                combined_pauli = combine_pauli_strings([operators[i] for i in combo_indices])
            
            weight = sum(1 for char in combined_pauli if char != 'I')
            
            all_combinations.append({
                'indices': combo_indices,
                'operators': [operators[i] for i in combo_indices],
                'combined': combined_pauli,
                'weight': weight
            })
    
    # Sort by weight
    all_combinations.sort(key=lambda x: x['weight'])
    
    print(f"Total combinations: {len(all_combinations)}")
    print(f"Showing first {min(max_combinations, len(all_combinations))} combinations (by weight):")
    
    min_weight = all_combinations[0]['weight'] if all_combinations else 0
    
    for i, combo in enumerate(all_combinations[:max_combinations]):
        indices_str = f"({', '.join(map(str, combo['indices']))})"
        print(f"  {i+1:3d}. {indices_str:15} → {combo['combined']} (weight: {combo['weight']})")
        if combo['weight'] > min_weight and i > 5:  # Stop showing after a few minimum weight examples
            break
    
    if len(all_combinations) > max_combinations:
        print(f"  ... and {len(all_combinations) - max_combinations} more combinations")
    
    return min_weight, all_combinations


def combine_pauli_strings(pauli_list: list) -> str:
    """Combine multiple Pauli strings using GF(2) addition (XOR)."""
    if not pauli_list:
        return ""
    
    N = len(pauli_list[0])
    result = ['I'] * N
    
    for pauli_str in pauli_list:
        for i, char in enumerate(pauli_str):
            if result[i] == 'I':
                result[i] = char
            elif char == 'I':
                continue  # result[i] stays the same
            elif result[i] == char:
                result[i] = 'I'  # X+X=I, Y+Y=I, Z+Z=I
            else:
                # Handle X+Z=Y, X+Y=Z, Z+Y=X cases
                if (result[i] == 'X' and char == 'Z') or (result[i] == 'Z' and char == 'X'):
                    result[i] = 'Y'
                elif (result[i] == 'X' and char == 'Y') or (result[i] == 'Y' and char == 'X'):
                    result[i] = 'Z'
                elif (result[i] == 'Z' and char == 'Y') or (result[i] == 'Y' and char == 'Z'):
                    result[i] = 'X'
    
    return ''.join(result)


def compare_methods(seed: str):
    """Complete comparison between GF(2) and polynomial approaches."""
    print(f"🔬 LOGICAL OPERATOR COMPARISON")
    print(f"Seed: {seed}")
    print("=" * 80)
    
    # Analyze GF(2) approach
    gf2_k, gf2_operators = analyze_gf2_logical_operators(seed)
    
    # Analyze polynomial approach  
    poly_k, poly_operators = analyze_polynomial_logical_operators(seed)
    
    # Compare basic parameters
    print(f"\n📊 COMPARISON SUMMARY")
    print("=" * 50)
    print(f"GF(2) logical qubits (k):      {gf2_k}")
    print(f"Polynomial logical qubits (k): {poly_k}")
    print(f"Parameters match: {gf2_k == poly_k}")
    
    if gf2_operators and poly_operators:
        print(f"GF(2) operators found:      {len(gf2_operators)}")
        print(f"Polynomial operators found: {len(poly_operators)}")
        print(f"Operator count match: {len(gf2_operators) == len(poly_operators)}")
        
        # Check if any operators are identical
        common_operators = set(gf2_operators) & set(poly_operators)
        print(f"Identical operators: {len(common_operators)}")
        if common_operators:
            for op in common_operators:
                print(f"  Common: {op}")
    
    # Find minimum weight combinations for both methods
    if gf2_operators:
        gf2_min_weight, gf2_combinations = find_minimum_weight_combinations(gf2_operators, "GF(2)")
    else:
        gf2_min_weight = None
        
    if poly_operators:
        poly_min_weight, poly_combinations = find_minimum_weight_combinations(poly_operators, "POLYNOMIAL")
    else:
        poly_min_weight = None
    
    # Compare distances
    print(f"\n🎯 DISTANCE COMPARISON")
    print("=" * 50)
    print(f"GF(2) minimum weight (distance):      {gf2_min_weight}")
    print(f"Polynomial minimum weight (distance): {poly_min_weight}")
    
    if gf2_min_weight is not None and poly_min_weight is not None:
        print(f"Distance match: {gf2_min_weight == poly_min_weight}")
        if gf2_min_weight != poly_min_weight:
            print(f"⚠️  DISCREPANCY: GF(2) found {gf2_min_weight}, Polynomial found {poly_min_weight}")
    
    return {
        'gf2_k': gf2_k,
        'poly_k': poly_k,
        'gf2_operators': gf2_operators,
        'poly_operators': poly_operators,
        'gf2_distance': gf2_min_weight,
        'poly_distance': poly_min_weight
    }


def main():
    """Main function to test different seeds."""
    test_cases = [
        ("XZIY", "Simple 4-qubit case"),
        ("IIIIIIIIXZZYZZZZYZZXIIIIIIII", "Complex 28-qubit QCA case (known discrepancy)"),
        ("ZZII", "Simple Z repetition"),
        ("XZZX", "4-qubit mixed case"),
    ]
    
    for seed, description in test_cases:
        print(f"\n" + "="*100)
        print(f"TESTING: {description}")
        print(f"SEED: {seed}")
        print("="*100)
        
        try:
            result = compare_methods(seed)
            
            # Summary for this test case
            print(f"\n📋 TEST CASE SUMMARY")
            print("-" * 30)
            if result['gf2_distance'] == result['poly_distance']:
                print(f"✅ MATCH: Both methods found distance {result['gf2_distance']}")
            else:
                print(f"❌ MISMATCH: GF(2)={result['gf2_distance']}, Polynomial={result['poly_distance']}")
                
        except Exception as e:
            print(f"❌ ERROR in test case: {e}")
            import traceback
            traceback.print_exc()
        
        print("\n" + "="*100)


if __name__ == "__main__":
    main() 