#!/usr/bin/env python3
"""
Comparison script for orthogonality checks.

This script compares the results of:
1. polydist.orthogonality.check_seed_orthogonality (polynomial method)
2. stabilizer.utils.seed_is_valid (symplectic tableau method)

Both should give the same results for the same seeds.
"""

import sys
import time
sys.path.insert(0, 'src')

from polydist.orthogonality import check_seed_orthogonality
from stabilizer.utils import seed_is_valid


def compare_methods(seed: str):
    """Compare both orthogonality check methods for a given seed with timing."""
    # Time polynomial method
    try:
        start_time = time.time()
        poly_result = check_seed_orthogonality(seed)
        poly_time = time.time() - start_time
    except Exception as e:
        poly_result = f"ERROR: {e}"
        poly_time = 0.0
    
    # Time symplectic method
    try:
        start_time = time.time()
        symplectic_result = seed_is_valid(seed)
        symplectic_time = time.time() - start_time
    except Exception as e:
        symplectic_result = f"ERROR: {e}"
        symplectic_time = 0.0
    
    return poly_result, poly_time, symplectic_result, symplectic_time


def test_seed(seed: str):
    """Test a single seed and display comparison results."""
    print(f"\nTesting seed: '{seed}' (length: {len(seed)})")
    print("-" * 60)
    
    poly_result, poly_time, symplectic_result, symplectic_time = compare_methods(seed)
    
    print(f"Polynomial method:  {poly_result:<20} (time: {poly_time:.6f}s)")
    print(f"Symplectic method:  {symplectic_result:<20} (time: {symplectic_time:.6f}s)")
    
    # Show speedup if both methods succeeded
    if poly_time > 0 and symplectic_time > 0:
        if poly_time < symplectic_time:
            speedup = symplectic_time / poly_time
            print(f"Speedup:            Polynomial is {speedup:.2f}x faster")
        else:
            speedup = poly_time / symplectic_time
            print(f"Speedup:            Symplectic is {speedup:.2f}x faster")
    
    # Check if results match (handle both bool and numpy.bool_)
    import numpy as np
    if (isinstance(poly_result, (bool, np.bool_)) and isinstance(symplectic_result, (bool, np.bool_))):
        match = bool(poly_result) == bool(symplectic_result)
        status = "✓ MATCH" if match else "✗ MISMATCH"
        print(f"Results match:      {status}")
        return match, poly_time, symplectic_time
    else:
        print(f"Results match:      ⚠ ERROR (cannot compare)")
        print(f"Debug: poly_result type: {type(poly_result)}, symplectic_result type: {type(symplectic_result)}")
        return None, poly_time, symplectic_time


def main():
    """Run comparison tests on various seeds."""
    print("="*60)
    print("ORTHOGONALITY METHOD COMPARISON")
    print("="*60)
    print("Comparing polynomial vs symplectic tableau methods")
    
    # Test cases covering various scenarios
    test_cases = [
        # Known valid cases
        "XY",
        "XYXY", 
        "XXZZ",
        "XZXZ",
        "YXYX",
        "ZYZY",
        
        # Edge cases
        "X",
        "Z", 
        "Y",
        "I",
        "II",
        "XI",
        "ZI",
        "YI",
        
        # All same operator cases
        "XX",
        "ZZ", 
        "YY",
        "XXX",
        "ZZZ",
        "YYY",
        "XXXX",
        "ZZZZ",
        "YYYY",
        
        # Mixed cases
        "XYZ",
        "XZY", 
        "YXZ",
        "YZX",
        "ZXY",
        "ZYX",
                 # "XYZW",  # Invalid Pauli - removed
        "XIYI",
        "XIZI",
        "XYZI",
        "XXYY",
        "XXZZ",
        "YYZZ",
        
        # Larger cases
        "XYXYXY",
        "XZXZXZ",
        "YZYZYZ",
        "XXYYZZII",
        "XIYZXIYZ",
        
        # Large test case: N=100 with X at position 50
        "I" * 50 + "X" + "I" * 49,
    ]
    
    results = []
    match_count = 0
    error_count = 0
    total_poly_time = 0.0
    total_symplectic_time = 0.0
    
    for seed in test_cases:
        try:
            match, poly_time, symplectic_time = test_seed(seed)
            results.append((seed, match, poly_time, symplectic_time))
            if match is True:
                match_count += 1
            elif match is None:
                error_count += 1
            
            # Accumulate timing statistics
            total_poly_time += poly_time
            total_symplectic_time += symplectic_time
            
        except Exception as e:
            print(f"Unexpected error testing '{seed}': {e}")
            results.append((seed, None, 0.0, 0.0))
            error_count += 1
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"{'Seed':<15} {'Length':<8} {'Match':<8} {'Poly Time':<12} {'Sympl Time':<12} {'Status'}")
    print("-" * 80)
    
    for seed, match, poly_time, symplectic_time in results:
        if match is True:
            status = "MATCH"
        elif match is False:
            status = "MISMATCH"
        else:
            status = "ERROR"
        
        match_str = str(match) if match is not None else "ERROR"
        seed_display = seed if len(seed) <= 12 else f"{seed[:5]}...{seed[-4:]}"
        
        print(f"{seed_display:<15} {len(seed):<8} {match_str:<8} {poly_time:<12.6f} {symplectic_time:<12.6f} {status}")
    
    total_tests = len(results)
    mismatch_count = total_tests - match_count - error_count
    
    print(f"\nTest Results:")
    print(f"  Total tests:        {total_tests}")
    print(f"  Matches:            {match_count}")
    print(f"  Mismatches:         {mismatch_count}")
    print(f"  Errors:             {error_count}")
    
    print(f"\nTiming Results:")
    print(f"  Total polynomial time:  {total_poly_time:.6f}s")
    print(f"  Total symplectic time:  {total_symplectic_time:.6f}s")
    if total_poly_time > 0 and total_symplectic_time > 0:
        if total_poly_time < total_symplectic_time:
            overall_speedup = total_symplectic_time / total_poly_time
            print(f"  Overall speedup:        Polynomial is {overall_speedup:.2f}x faster")
        else:
            overall_speedup = total_poly_time / total_symplectic_time
            print(f"  Overall speedup:        Symplectic is {overall_speedup:.2f}x faster")
    
    if mismatch_count == 0 and error_count == 0:
        print(f"\n🎉 All methods agree! Both implementations are consistent.")
    elif mismatch_count > 0:
        print(f"\n⚠️  Found {mismatch_count} mismatches - implementations may differ!")
    else:
        print(f"\n⚠️  Some tests had errors - check implementation details.")


if __name__ == "__main__":
    main() 