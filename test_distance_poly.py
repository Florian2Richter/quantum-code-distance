#!/usr/bin/env python3
"""
Test script for polynomial-based distance calculation.
"""

import sys
sys.path.append('src')

from polydist.distance import (
    find_distance_polynomial,
    analyze_distance_calculation,
    verify_distance_consistency
)

def test_simple_seeds():
    """Test distance calculation on simple known seeds."""
    
    test_cases = [
        ("XZIY", "Small test case"),
        ("ZZII", "Simple Z pattern"),
        ("XZZX", "Mixed pattern"),
        ("IIII", "All identity"),
        ("IIIIIIIIXZZYZZZZYZZXIIIIIIII", "Complex QCA evolution case - Expected distance 9"),
    ]
    
    print("Testing polynomial distance calculation:")
    print("=" * 50)
    
    for seed, description in test_cases:
        print(f"\nTesting seed: {seed} ({description})")
        print("-" * 30)
        
        try:
            # Test basic distance calculation
            distance, logical_ops = find_distance_polynomial(seed, verbose=True)
            print(f"Distance: {distance}")
            print(f"Logical operators used: {logical_ops}")
            
            # Test consistency verification
            consistency = verify_distance_consistency(seed, verbose=False)
            print(f"Consistency check: {consistency}")
            
            # For the complex case, show detailed comparison
            if "IIIIIIIIXZZYZZZZYZZXIIIIIIII" in seed:
                print(f"\n🔍 DETAILED ANALYSIS for complex case:")
                print(f"  Polynomial distance: {distance}")
                print(f"  Expected GF(2) distance: 9")
                print(f"  Difference: {distance - 9 if distance else 'N/A'}")
                print(f"  Best logical operator found: {logical_ops[0] if logical_ops else 'None'}")
                print(f"  Operator weight: {distance}")
                if distance != 9:
                    print(f"  ⚠️  MISMATCH: Polynomial found {distance}, expected 9")
            
            # Test complete analysis
            analysis = analyze_distance_calculation(seed, verbose=False)
            print(f"Analysis successful: {analysis['analysis_successful']}")
            if analysis['analysis_successful']:
                print(f"Verification: extraction={analysis['extraction_verified']}, consistency={analysis['consistency_verified']}")
            
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    test_simple_seeds() 