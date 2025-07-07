#!/usr/bin/env python3
"""
Comparison script for rank calculation and logical qubit computation.

This script compares the results of:
1. Stabilizer tableau method: rank = rank(tableau), k = L - rank
2. Polynomial GCD method: k = deg(gcd(s_X(x), s_Z(x), x^N - 1)), rank = L - k

Both should give the same results for the same seeds.
"""

import sys
import time
import numpy as np
sys.path.insert(0, 'src')

from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.utils import seed_is_valid
from polydist.gcd import compute_logical_qubits, analyze_gcd_computation
from polydist.orthogonality import check_seed_orthogonality


def compute_stabilizer_rank_and_logical(seed: str):
    """Compute rank and logical qubits using stabilizer tableau method."""
    try:
        start_time = time.time()
        
        # Build lattice from seed
        pauli_list = list(seed)
        stab_ops = build_lattice(pauli_list)
        
        # Build tableau
        tableau = build_tableau(stab_ops)
        
        # Compute rank
        rank = compute_rank(tableau)
        
        L = len(seed)
        n_logical = L - rank
        
        elapsed_time = time.time() - start_time
        
        return {
            'rank': rank,
            'n_logical': n_logical,
            'L': L,
            'success': True,
            'time': elapsed_time
        }
        
    except Exception as e:
        return {
            'rank': None,
            'n_logical': None,
            'L': len(seed),
            'success': False,
            'error': str(e),
            'time': 0.0
        }


def compute_polynomial_rank_and_logical(seed: str):
    """Compute rank and logical qubits using polynomial GCD method."""
    try:
        start_time = time.time()
        
        # Compute logical qubits via GCD
        k = compute_logical_qubits(seed)
        
        L = len(seed)
        rank = L - k  # In polynomial formalism: rank = L - k
        
        elapsed_time = time.time() - start_time
        
        return {
            'rank': rank,
            'n_logical': k,
            'L': L,
            'success': True,
            'time': elapsed_time
        }
        
    except Exception as e:
        return {
            'rank': None,
            'n_logical': None,
            'L': len(seed),
            'success': False,
            'error': str(e),
            'time': 0.0
        }


def compare_methods(seed: str):
    """Compare both rank calculation methods for a given seed."""
    stab_result = compute_stabilizer_rank_and_logical(seed)
    poly_result = compute_polynomial_rank_and_logical(seed)
    
    return stab_result, poly_result


def test_seed(seed: str):
    """Test a single seed and display comparison results."""
    print(f"\nTesting seed: '{seed}' (length: {len(seed)})")
    print("-" * 70)
    
    # First check if seed generates valid stabilizers (orthogonality check)
    try:
        stab_orthogonal = seed_is_valid(seed)
        poly_orthogonal = check_seed_orthogonality(seed)
        
        print(f"Orthogonality check:  Stabilizer={stab_orthogonal}, Polynomial={poly_orthogonal}")
        
        if not stab_orthogonal or not poly_orthogonal:
            if stab_orthogonal != poly_orthogonal:
                print(f"⚠️  ORTHOGONALITY MISMATCH! Cannot proceed with rank comparison.")
                return "orthogonality_mismatch", 0.0, 0.0
            else:
                print(f"❌ INVALID SEED: Does not generate valid stabilizers (not orthogonal).")
                return "invalid_seed", 0.0, 0.0
        else:
            print(f"✓ Valid seed: Generates orthogonal stabilizers")
            
    except Exception as e:
        print(f"❌ ERROR in orthogonality check: {e}")
        return "orthogonality_error", 0.0, 0.0
    
    stab_result, poly_result = compare_methods(seed)
    
    # Display results
    if stab_result['success']:
        print(f"Stabilizer method:   rank={stab_result['rank']:<3} k={stab_result['n_logical']:<3} (time: {stab_result['time']:.6f}s)")
    else:
        print(f"Stabilizer method:   ERROR: {stab_result['error']:<30} (time: {stab_result['time']:.6f}s)")
    
    if poly_result['success']:
        print(f"Polynomial method:   rank={poly_result['rank']:<3} k={poly_result['n_logical']:<3} (time: {poly_result['time']:.6f}s)")
    else:
        print(f"Polynomial method:   ERROR: {poly_result['error']:<30} (time: {poly_result['time']:.6f}s)")
    
    # Show speedup if both methods succeeded
    if stab_result['success'] and poly_result['success'] and stab_result['time'] > 0 and poly_result['time'] > 0:
        if poly_result['time'] < stab_result['time']:
            speedup = stab_result['time'] / poly_result['time']
            print(f"Speedup:             Polynomial is {speedup:.2f}x faster")
        else:
            speedup = poly_result['time'] / stab_result['time']
            print(f"Speedup:             Stabilizer is {speedup:.2f}x faster")
    
    # Check if results match
    if stab_result['success'] and poly_result['success']:
        rank_match = stab_result['rank'] == poly_result['rank']
        logical_match = stab_result['n_logical'] == poly_result['n_logical']
        overall_match = rank_match and logical_match
        
        print(f"Results match:")
        print(f"  Rank:              {'✓' if rank_match else '✗'} ({stab_result['rank']} vs {poly_result['rank']})")
        print(f"  Logical qubits:    {'✓' if logical_match else '✗'} ({stab_result['n_logical']} vs {poly_result['n_logical']})")
        print(f"  Overall:           {'✓ MATCH' if overall_match else '✗ MISMATCH'}")
        
        return overall_match, stab_result['time'], poly_result['time']
    else:
        print(f"Results match:       ⚠ ERROR (cannot compare)")
        return None, stab_result['time'], poly_result['time']


def main():
    """Run rank comparison tests on various seeds."""
    print("="*80)
    print("RANK CALCULATION METHOD COMPARISON")
    print("="*80)
    print("Comparing stabilizer tableau vs polynomial GCD methods")
    
    # Test cases covering various scenarios
    test_cases = [
        # Simple cases
        "X",
        "Z", 
        "Y",
        "I",
        
        # Two-qubit cases
        "XY",
        "XZ",
        "YZ",
        "XX",
        "ZZ",
        "YY",
        "II",
        "XI",
        "ZI",
        "YI",
        
        # Four-qubit cases
        "XYXY", 
        "XXZZ",
        "XZXZ",
        "YXYX",
        "ZYZY",
        "XXXX",
        "ZZZZ",
        "YYYY",
        "IIII",
        
        # Mixed cases
        "XYZ",
        "XZY", 
        "YXZ",
        "YZX",
        "ZXY",
        "ZYX",
        "XIYI",
        "XIZI",
        "XYZI",
        "XXYY",
        "YYZZ",
        
        # Larger cases
        "XYXYXY",
        "XZXZXZ",
        "YZYZYZ",
        "XXYYZZII",
        "XIYZXIYZ",
        
        # Large test case: N=20 with X at position 10
        "I" * 10 + "X" + "I" * 9,
        
        # Large test case: N=50 with pattern
        "XY" * 25,
        
        # Very large test case: N=100 with X at position 50  
        "I" * 50 + "X" + "I" * 49,
    ]
    
    results = []
    match_count = 0
    error_count = 0
    total_stab_time = 0.0
    total_poly_time = 0.0
    
    invalid_seed_count = 0
    orthogonality_mismatch_count = 0
    orthogonality_error_count = 0
    
    for seed in test_cases:
        try:
            match, stab_time, poly_time = test_seed(seed)
            results.append((seed, match, stab_time, poly_time))
            
            if match == "invalid_seed":
                invalid_seed_count += 1
            elif match == "orthogonality_mismatch":
                orthogonality_mismatch_count += 1
            elif match == "orthogonality_error":
                orthogonality_error_count += 1
            elif match is True:
                match_count += 1
            elif match is None:
                error_count += 1
            
            # Accumulate timing statistics (only for actual computations)
            if match not in ["invalid_seed", "orthogonality_mismatch", "orthogonality_error"]:
                total_stab_time += stab_time
                total_poly_time += poly_time
            
        except Exception as e:
            print(f"Unexpected error testing '{seed}': {e}")
            results.append((seed, None, 0.0, 0.0))
            error_count += 1
    
    # Summary
    print("\n" + "="*90)
    print("SUMMARY")
    print("="*90)
    print(f"{'Seed':<15} {'Length':<8} {'Match':<8} {'Stab Time':<12} {'Poly Time':<12} {'Speedup':<10} {'Status'}")
    print("-" * 90)
    
    for seed, match, stab_time, poly_time in results:
        if match is True:
            status = "MATCH"
        elif match is False:
            status = "MISMATCH"
        elif match == "invalid_seed":
            status = "INVALID"
        elif match == "orthogonality_mismatch":
            status = "ORTH_MISMATCH"
        elif match == "orthogonality_error":
            status = "ORTH_ERROR"
        else:
            status = "ERROR"
        
        match_str = str(match) if match not in ["invalid_seed", "orthogonality_mismatch", "orthogonality_error", None] else "N/A"
        seed_display = seed if len(seed) <= 12 else f"{seed[:5]}...{seed[-4:]}"
        
        # Calculate speedup
        if stab_time > 0 and poly_time > 0:
            if poly_time < stab_time:
                speedup = f"{stab_time/poly_time:.2f}x (P)"
            else:
                speedup = f"{poly_time/stab_time:.2f}x (S)"
        else:
            speedup = "N/A"
        
        print(f"{seed_display:<15} {len(seed):<8} {match_str:<8} {stab_time:<12.6f} {poly_time:<12.6f} {speedup:<10} {status}")
    
    total_tests = len(results)
    rank_mismatch_count = sum(1 for _, match, _, _ in results if match is False)
    
    print(f"\nTest Results:")
    print(f"  Total tests:             {total_tests}")
    print(f"  Valid seeds tested:      {match_count + rank_mismatch_count}")
    print(f"  Rank matches:            {match_count}")
    print(f"  Rank mismatches:         {rank_mismatch_count}")
    print(f"  Invalid seeds:           {invalid_seed_count}")
    print(f"  Orthogonality mismatches: {orthogonality_mismatch_count}")
    print(f"  Orthogonality errors:    {orthogonality_error_count}")
    print(f"  Other errors:            {error_count}")
    
    print(f"\nTiming Results:")
    print(f"  Total stabilizer time: {total_stab_time:.6f}s")
    print(f"  Total polynomial time: {total_poly_time:.6f}s")
    if total_stab_time > 0 and total_poly_time > 0:
        if total_poly_time < total_stab_time:
            overall_speedup = total_stab_time / total_poly_time
            print(f"  Overall speedup:       Polynomial is {overall_speedup:.2f}x faster")
        else:
            overall_speedup = total_poly_time / total_stab_time
            print(f"  Overall speedup:       Stabilizer is {overall_speedup:.2f}x faster")
    
    # Performance analysis by size
    print(f"\nPerformance by Size:")
    small_cases = [(seed, stab_t, poly_t) for seed, _, stab_t, poly_t in results if len(seed) <= 4]
    medium_cases = [(seed, stab_t, poly_t) for seed, _, stab_t, poly_t in results if 5 <= len(seed) <= 20]
    large_cases = [(seed, stab_t, poly_t) for seed, _, stab_t, poly_t in results if len(seed) > 20]
    
    for category, cases in [("Small (≤4)", small_cases), ("Medium (5-20)", medium_cases), ("Large (>20)", large_cases)]:
        if cases:
            avg_stab = sum(stab_t for _, stab_t, _ in cases) / len(cases)
            avg_poly = sum(poly_t for _, _, poly_t in cases) / len(cases)
            print(f"  {category:<12}: Stabilizer {avg_stab:.6f}s, Polynomial {avg_poly:.6f}s")
    
    valid_seeds_tested = match_count + rank_mismatch_count
    
    if valid_seeds_tested == 0:
        print(f"\n❌ No valid seeds found for testing!")
    elif rank_mismatch_count == 0 and orthogonality_mismatch_count == 0:
        print(f"\n🎉 All methods agree! Both implementations are consistent for valid seeds.")
        if invalid_seed_count > 0:
            print(f"    Note: {invalid_seed_count} invalid seeds were correctly identified by both methods.")
    elif rank_mismatch_count > 0:
        print(f"\n⚠️  Found {rank_mismatch_count} rank mismatches among valid seeds - implementations may differ!")
    elif orthogonality_mismatch_count > 0:
        print(f"\n⚠️  Found {orthogonality_mismatch_count} orthogonality method mismatches - check orthogonality implementations!")
    else:
        print(f"\n⚠️  Some tests had errors - check implementation details.")


if __name__ == "__main__":
    main() 