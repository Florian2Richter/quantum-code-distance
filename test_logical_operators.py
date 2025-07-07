#!/usr/bin/env python3
"""
Test script for logical operator extraction using polynomial formalism.

This script demonstrates the extraction of logical Z operators from stabilizer
code seeds using the extended Euclidean algorithm approach.
"""

import sys
import time
from src.polydist.logical_operators import (
    extract_logical_z_operator,
    verify_logical_operator, 
    analyze_logical_operator_extraction
)
from src.polydist.orthogonality import check_seed_orthogonality


def test_logical_operator_extraction():
    """Test logical operator extraction on various seeds."""
    
    # Test seeds with different properties
    test_seeds = [
        # Basic seeds
        "I",          # k=1 logical qubit
        "X",          # k=0 logical qubits  
        "Z",          # k=0 logical qubits
        "Y",          # k=0 logical qubits
        
        # Two-qubit seeds
        "II",         # k=2 logical qubits
        "XX",         # k=1 logical qubit
        "ZZ",         # k=1 logical qubit  
        "XZ",         # k=0 logical qubits
        "YY",         # k=1 logical qubit
        
        # Four-qubit seeds
        "IIII",       # k=4 logical qubits
        "XXXX",       # k=3 logical qubits
        "ZZZZ",       # k=3 logical qubits
        "XYXY",       # k=2 logical qubits
        "XZXZ",       # k=2 logical qubits
        "XXZZ",       # k=1 logical qubit
        
        # Larger seeds
        "XYXYXY",     # k=4 logical qubits
        "XZXZXZ",     # k=4 logical qubits
        
        # Invalid seeds (non-orthogonal)
        "XYZ",        # Should fail orthogonality check
    ]
    
    print("="*80)
    print("LOGICAL OPERATOR EXTRACTION TEST")
    print("="*80)
    print("Testing polynomial-based extraction of logical Z operators\n")
    
    successful_extractions = 0
    verified_operators = 0
    total_valid_seeds = 0
    
    for seed in test_seeds:
        print(f"Testing seed: '{seed}' (length: {len(seed)})")
        print("-" * 70)
        
        # Check orthogonality first
        is_orthogonal = check_seed_orthogonality(seed)
        print(f"Orthogonality check: {'✓ Valid' if is_orthogonal else '❌ Invalid'}")
        
        if not is_orthogonal:
            print("❌ SKIPPING: Seed does not generate valid stabilizers (not orthogonal).\n")
            continue
            
        total_valid_seeds += 1
        
        try:
            # Time the extraction
            start_time = time.time()
            analysis = analyze_logical_operator_extraction(seed, verbose=True)
            extraction_time = time.time() - start_time
            
            print(f"Extraction time: {extraction_time:.6f}s")
            
            if analysis['extraction_successful']:
                successful_extractions += 1
                v_str, u_str = analysis['logical_Z_operator']
                print(f"✓ Logical Z operator extracted:")
                print(f"  L_Z = (v(x), u(x)) = ({v_str}, {u_str})")
                
                verification = analysis['verification']
                print(f"Verification: {'✓ PASSED' if verification['is_correct'] else '❌ FAILED'}")
                
                if verification['is_correct']:
                    verified_operators += 1
                    
                print(f"Logical qubits k: {verification['logical_qubits_k']}")
                print(f"Relation: u*s_X + v*s_Z = {verification['u_s_X_plus_v_s_Z']}")
                print(f"Expected f(x): {verification['expected_f']}")
                
                if 'verbose_details' in analysis:
                    details = analysis['verbose_details']
                    print(f"\nIntermediate results:")
                    print(f"  s_X(x) = {details['s_X']}")
                    print(f"  s_Z(x) = {details['s_Z']}")
                    print(f"  Stage 1 - g1 = gcd(s_X, s_Z) = {details['stage1_g1']}")
                    print(f"  Stage 1 - u1*s_X + v1*s_Z = g1: u1={details['stage1_u1']}, v1={details['stage1_v1']}")
                    print(f"  Stage 2 - f = gcd(g1, x^N-1) = {details['stage2_f']}")
                    print(f"  Stage 2 - w*g1 + t*(x^N-1) = f: w={details['stage2_w']}, t={details['stage2_t']}")
                    print(f"  Final - u = w*u1 = {details['final_u']}")
                    print(f"  Final - v = w*v1 = {details['final_v']}")
                
            else:
                print(f"❌ Extraction failed: {analysis.get('error', 'Unknown error')}")
                
        except Exception as e:
            print(f"❌ Exception during extraction: {e}")
            
        print()
    
    # Summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total seeds tested: {len(test_seeds)}")
    print(f"Valid orthogonal seeds: {total_valid_seeds}")
    print(f"Successful extractions: {successful_extractions}")
    print(f"Verified operators: {verified_operators}")
    print(f"Success rate: {successful_extractions/total_valid_seeds*100:.1f}% extraction")
    print(f"Verification rate: {verified_operators/successful_extractions*100:.1f}% verification" if successful_extractions > 0 else "Verification rate: N/A")
    
    if verified_operators == successful_extractions == total_valid_seeds:
        print("\n🎉 All valid seeds successfully processed with verified logical operators!")
    else:
        print(f"\n⚠️  Some extractions failed or verification issues found.")


def test_specific_example():
    """Test a specific example in detail."""
    print("\n" + "="*80)
    print("DETAILED EXAMPLE: SEED 'XXZZ'")
    print("="*80)
    
    seed = "XXZZ"
    print(f"Analyzing seed: '{seed}'")
    print(f"Expected: k=1 logical qubit (from previous GCD analysis)")
    
    # Extract logical operator
    analysis = analyze_logical_operator_extraction(seed, verbose=True)
    
    if analysis['extraction_successful']:
        verification = analysis['verification']
        print(f"\nLogical Z operator: L_Z = {analysis['logical_Z_operator']}")
        print(f"Verification: {'PASSED' if verification['is_correct'] else 'FAILED'}")
        print(f"Logical qubits k: {verification['logical_qubits_k']}")
        
        # Show the mathematical relationship
        print(f"\nMathematical verification:")
        print(f"u(x)*s_X(x) + v(x)*s_Z(x) = {verification['u_s_X_plus_v_s_Z']}")
        print(f"f(x) = gcd(s_X, s_Z, x^N-1) = {verification['expected_f']}")
        print(f"Difference = {verification['difference']} (should be 0)")
        
        # Show commutation property
        print(f"\nCommutation check:")
        print(f"The logical operator L_Z = (v, u) should commute with all")
        print(f"cyclic shifts of the stabilizer seed s = (s_X, s_Z).")
        print(f"This means [L_Z, x^j * s] = 0 for all j = 0, 1, ..., N-1.")
        
    else:
        print(f"❌ Extraction failed: {analysis.get('error', 'Unknown error')}")


if __name__ == "__main__":
    try:
        test_logical_operator_extraction()
        test_specific_example()
        
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nUnexpected error: {e}")
        sys.exit(1) 