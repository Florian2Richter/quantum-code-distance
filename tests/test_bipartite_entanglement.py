"""
Comprehensive tests for the existing bipartite entanglement calculation.

Tests the compute_entanglement() function in stabilizer.utils without modifying
any source code - only testing the existing functionality thoroughly.
"""

import pytest
import numpy as np
import json
from hypothesis import given, assume, strategies as st
from hypothesis import settings, HealthCheck

from stabilizer.utils import compute_entanglement, seed_is_valid, pauli_to_symplectic
from stabilizer.generator import parse_seed
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance
from stabilizer.qca import qca_evolution_step


def build_code_tableau(seed: str):
    """Helper to build tableau and get logical operators from a seed."""
    pauli_list = list(seed)
    stab_ops = build_lattice(pauli_list)
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    
    n = len(seed)
    k = n - rank
    
    if k == 0:
        logical_ops = []
    else:
        _, logical_ops = find_distance(tableau, return_logical_ops=True)
    
    return tableau, logical_ops, k


class TestBasicEntanglementFunctionality:
    """Test basic functionality of compute_entanglement()."""
    
    def test_golden_reference_cases(self):
        """Test entanglement matches known golden reference values."""
        with open("tests/data/golden.json") as f:
            golden_data = json.load(f)
        
        for seed, expected_results in golden_data.items():
            if not seed_is_valid(seed):
                continue
                
            # Test step 0 (initial state)
            initial_result = expected_results[0]
            tableau, logical_ops, k = build_code_tableau(seed)
            
            entanglement = compute_entanglement(tableau, logical_ops if logical_ops else None)
            expected_ent = initial_result["entanglement"]
            
            assert entanglement == expected_ent, \
                f"Seed {seed}: expected entanglement {expected_ent}, got {entanglement}"
    
    def test_entanglement_without_logical_operators(self):
        """Test entanglement calculation without logical operators."""
        test_cases = [
            "XZIY",  # Test case from golden data
            "ZZII",  # Simple Z-repetition case
        ]
        
        for seed in test_cases:
            if not seed_is_valid(seed):
                continue
                
            tableau, logical_ops, k = build_code_tableau(seed)
            
            # Test without logical operators
            ent_without = compute_entanglement(tableau, None)
            assert ent_without >= 0, f"Negative entanglement for {seed}: {ent_without}"
            
            # Test with logical operators  
            ent_with = compute_entanglement(tableau, logical_ops if logical_ops else None)
            assert ent_with >= 0, f"Negative entanglement with logical ops for {seed}: {ent_with}"
    
    def test_entanglement_deterministic(self):
        """Test that entanglement calculation is deterministic."""
        seed = "XZIY"
        if not seed_is_valid(seed):
            pytest.skip("Seed not valid")
            
        tableau, logical_ops, k = build_code_tableau(seed)
        
        # Multiple calls should give same result
        ent1 = compute_entanglement(tableau, logical_ops)
        ent2 = compute_entanglement(tableau, logical_ops)
        ent3 = compute_entanglement(tableau, logical_ops)
        
        assert ent1 == ent2 == ent3, f"Non-deterministic results: {ent1}, {ent2}, {ent3}"


class TestEntanglementBounds:
    """Test fundamental bounds and constraints on entanglement values."""
    
    def test_non_negativity(self):
        """Test that entanglement is always non-negative."""
        test_seeds = ["XZIY", "ZZII", "XZZX", "IIII"]
        
        for seed in test_seeds:
            if not seed_is_valid(seed):
                continue
                
            tableau, logical_ops, k = build_code_tableau(seed)
            
            # Test both with and without logical operators
            ent_without = compute_entanglement(tableau, None)
            ent_with = compute_entanglement(tableau, logical_ops if logical_ops else None)
            
            assert ent_without >= 0, f"Negative entanglement without logical ops for {seed}: {ent_without}"
            assert ent_with >= 0, f"Negative entanglement with logical ops for {seed}: {ent_with}"
    
    def test_system_size_bounds(self):
        """Test that entanglement doesn't exceed reasonable bounds."""
        test_seeds = ["XZIY", "ZZII", "XZZX"]
        
        for seed in test_seeds:
            if not seed_is_valid(seed):
                continue
                
            tableau, logical_ops, k = build_code_tableau(seed)
            n = len(seed)
            
            ent = compute_entanglement(tableau, logical_ops if logical_ops else None)
            
            # Entanglement should not exceed system size
            assert ent <= n, f"Entanglement {ent} exceeds system size {n} for {seed}"


class TestLogicalOperatorEffects:
    """Test how logical operators affect entanglement calculation."""
    
    def test_z_logical_filtering(self):
        """Test that only Z-logical operators are included."""
        seed = "XZIY"
        if not seed_is_valid(seed):
            pytest.skip("Seed not valid")
            
        tableau, logical_ops, k = build_code_tableau(seed)
        
        if logical_ops:
            # Manually filter to only Z-logicals (X part all zeros)
            L = tableau.shape[1] // 2
            z_logicals = [op for op in logical_ops if not np.any(op[:L])]
            
            # Test with all logical ops vs only Z-logicals
            ent_all = compute_entanglement(tableau, logical_ops)
            ent_z_only = compute_entanglement(tableau, z_logicals)
            
            # Should be the same since function filters internally
            assert ent_all == ent_z_only, \
                f"Different entanglement: all logical ops {ent_all} vs Z-only {ent_z_only}"
    
    def test_empty_logical_operators(self):
        """Test behavior with empty logical operators."""
        seed = "XZIY"
        if not seed_is_valid(seed):
            pytest.skip("Seed not valid")
            
        tableau, logical_ops, k = build_code_tableau(seed)
        
        # Test with empty list vs None
        ent_none = compute_entanglement(tableau, None)
        ent_empty = compute_entanglement(tableau, [])
        
        assert ent_none == ent_empty, f"Different results for None vs empty list: {ent_none} vs {ent_empty}"


class TestEntanglementDynamics:
    """Test entanglement under quantum cellular automaton evolution."""
    
    def test_qca_evolution_entanglement(self):
        """Test how entanglement changes under QCA evolution."""
        initial_seeds = ["XZIY", "ZZII", "XZZX"]
        
        for seed in initial_seeds:
            if not seed_is_valid(seed):
                continue
                
            # Initial entanglement
            tableau0, logical_ops0, k0 = build_code_tableau(seed)
            ent0 = compute_entanglement(tableau0, logical_ops0 if logical_ops0 else None)
            
            # Evolve one step
            pauli_list = list(seed)
            evolved_pauli = qca_evolution_step(pauli_list)
            evolved_seed = ''.join(evolved_pauli)
            
            if seed_is_valid(evolved_seed):
                tableau1, logical_ops1, k1 = build_code_tableau(evolved_seed)
                ent1 = compute_entanglement(tableau1, logical_ops1 if logical_ops1 else None)
                
                # Both entanglements should be non-negative
                assert ent0 >= 0, f"Initial entanglement negative for {seed}: {ent0}"
                assert ent1 >= 0, f"Evolved entanglement negative for {evolved_seed}: {ent1}"


class TestPropertyBasedEntanglement:
    """Property-based tests using Hypothesis for robust coverage."""
    
    @given(seed=st.text(alphabet='IXYZ', min_size=2, max_size=6))
    @settings(max_examples=30, deadline=30000, suppress_health_check=[HealthCheck.filter_too_much])
    def test_entanglement_always_non_negative(self, seed):
        """Property test: entanglement is always non-negative for valid seeds."""
        assume(seed_is_valid(seed))
        
        tableau, logical_ops, k = build_code_tableau(seed)
        
        # Test both with and without logical operators
        ent_without = compute_entanglement(tableau, None)
        ent_with = compute_entanglement(tableau, logical_ops if logical_ops else None)
        
        assert ent_without >= 0, f"Negative entanglement without logical ops for {seed}: {ent_without}"
        assert ent_with >= 0, f"Negative entanglement with logical ops for {seed}: {ent_with}"


class TestEntanglementConsistency:
    """Test consistency of entanglement calculation across different scenarios."""
    
    def test_cross_reference_with_golden_data(self):
        """Cross-reference entanglement values with golden data."""
        with open("tests/data/golden.json") as f:
            golden_data = json.load(f)
        
        for seed, expected_results in golden_data.items():
            if not seed_is_valid(seed):
                continue
                
            # Test that our calculation matches the expected value
            initial_result = expected_results[0]
            tableau, logical_ops, k = build_code_tableau(seed)
            our_ent = compute_entanglement(tableau, logical_ops if logical_ops else None)
            expected_ent = initial_result["entanglement"]
            
            assert our_ent == expected_ent, \
                f"Entanglement mismatch for {seed}: got {our_ent}, expected {expected_ent}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"]) 