"""
Property-based tests for quantum error correction invariants using Hypothesis.

These tests validate fundamental mathematical properties that should hold
for all valid stabilizer codes, going beyond specific test cases.
"""

import pytest
from hypothesis import given, assume, strategies as st
from hypothesis import settings, Verbosity, HealthCheck
import numpy as np

from stabilizer.utils import seed_is_valid
from stabilizer.generator import parse_seed
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance
from stabilizer.qca import qca_evolution_step


def lattice_based_analysis(seed: str):
    """
    Analyze a lattice-based stabilizer code (like the main CLI does).
    This is different from distance_from_tableau which treats seed as single stabilizer.
    """
    pauli_list = list(seed)
    stab_ops = build_lattice(pauli_list)  # Create all cyclic translations
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    
    n = len(seed)
    k = n - rank
    
    if k == 0:
        d = 0
    else:
        d = find_distance(tableau)
    
    return (n, k, d)


# Strategy for generating valid commuting seeds
@st.composite
def valid_seed(draw, min_N=3, max_N=11):
    """Generate valid Pauli seeds that produce commuting stabilizers."""
    N = draw(st.integers(min_value=min_N, max_value=max_N))
    
    # Generate random Pauli string
    pauli_chars = ['I', 'X', 'Y', 'Z']
    seed = ''.join(draw(st.sampled_from(pauli_chars)) for _ in range(N))
    
    # Filter to only valid (commuting) seeds
    assume(seed_is_valid(seed))
    return seed


# Strategy for single Pauli characters
pauli_char = st.sampled_from(['I', 'X', 'Y', 'Z'])


class TestStabilizerInvariants:
    """Test fundamental invariants of stabilizer codes."""
    
    @given(seed=valid_seed(min_N=3, max_N=8))
    @settings(max_examples=50, deadline=30000)  # Allow more time for complex computations
    def test_translation_invariance(self, seed):
        """Code parameters should be invariant under cyclic translations."""
        n0, k0, d0 = lattice_based_analysis(seed)
        
        # Test all cyclic shifts
        for shift in range(1, len(seed)):
            shifted = seed[shift:] + seed[:shift]
            n1, k1, d1 = lattice_based_analysis(shifted)
            assert (n1, k1, d1) == (n0, k0, d0), \
                f"Translation broke invariance: {seed} -> {shifted}, {(n0,k0,d0)} -> {(n1,k1,d1)}"

    @given(css_type=st.sampled_from(['X', 'Z']), n=st.integers(min_value=3, max_value=8))
    @settings(max_examples=30, deadline=30000, suppress_health_check=[HealthCheck.filter_too_much])
    def test_css_duality_pure_codes(self, css_type, n):
        """X↔Z swap should preserve n and k for pure CSS codes."""
        # Create a pure CSS code (only X's or only Z's, plus I's)
        base_seed = css_type * min(n//2 + 1, n-1) + 'I' * (n - min(n//2 + 1, n-1))
        seed = base_seed[:n]  # Truncate to exact length
        
        # Only test if it's a valid commuting seed
        assume(seed_is_valid(seed))
        
        n0, k0, d0 = lattice_based_analysis(seed)
        
        # Swap X↔Z
        css_dual = seed.translate(str.maketrans("XZ", "ZX"))
        n1, k1, d1 = lattice_based_analysis(css_dual)
        
        # Physical qubits and logical qubits should be preserved for pure CSS codes
        assert n1 == n0, f"Physical qubits changed: {n0} -> {n1}"
        assert k1 == k0, f"Logical qubits changed: {k0} -> {k1}"

    @given(seed=valid_seed(min_N=3, max_N=10))
    @settings(max_examples=100, deadline=20000)
    def test_distance_bounds(self, seed):
        """Test fundamental bounds on code distance."""
        n, k, d = lattice_based_analysis(seed)
        
        # Basic sanity checks
        assert n == len(seed), f"Physical qubits mismatch: {n} != {len(seed)}"
        assert k >= 0, f"Negative logical qubits: {k}"
        assert 0 <= k <= n, f"Logical qubits out of bounds: {k} not in [0, {n}]"
        
        if k == 0:
            assert d == 0, f"Distance should be 0 when k=0, got {d}"
        else:
            assert d >= 1, f"Distance should be ≥1 when k>0, got {d}"
            # Singleton bound: d ≤ n (can't have distance larger than code length)
            assert d <= n, f"Distance {d} exceeds code length {n}"

    @given(seed=valid_seed(min_N=3, max_N=8))
    @settings(max_examples=30, deadline=30000)
    def test_idempotence(self, seed):
        """Multiple calls should give identical results."""
        result1 = lattice_based_analysis(seed)
        result2 = lattice_based_analysis(seed)
        result3 = lattice_based_analysis(seed)
        
        assert result1 == result2 == result3, \
            f"Non-idempotent results: {result1}, {result2}, {result3}"

    @given(seed=valid_seed(min_N=4, max_N=7))
    @settings(max_examples=20, deadline=60000)
    def test_qca_evolution_consistency(self, seed):
        """QCA evolution should maintain certain properties."""
        n0, k0, d0 = lattice_based_analysis(seed)
        
        # Apply one QCA step
        pauli_list = list(seed)
        evolved_pauli = qca_evolution_step(pauli_list)
        evolved_seed = ''.join(evolved_pauli)
        
        # Check if evolved seed is still valid (commuting)
        if seed_is_valid(evolved_seed):
            n1, k1, d1 = lattice_based_analysis(evolved_seed)
            
            # Physical qubits should be preserved
            assert n1 == n0, f"QCA changed physical qubits: {n0} -> {n1}"
            
            # Distance should be non-negative
            assert d1 >= 0, f"QCA produced negative distance: {d1}"
            
            # If no logical qubits, distance should be 0
            if k1 == 0:
                assert d1 == 0, f"QCA: k=0 but d={d1}"


class TestSpecialCases:
    """Test known special cases and edge conditions."""
    
    @given(n=st.integers(min_value=2, max_value=8))
    def test_all_identity_seed(self, n):
        """All-I seed gives trivial stabilizers: rank=0, so k=n."""
        seed = 'I' * n
        assert seed_is_valid(seed), "All-I should be valid"
        
        result = lattice_based_analysis(seed)
        n_actual, k_actual, d_actual = result
        
        # All-I creates n identical I...I stabilizers (all zeros in symplectic form)
        # Rank of zero matrix is 0, so k = n - 0 = n (all qubits are logical)
        expected_k = n
        assert k_actual == expected_k, f"All-I: expected k={expected_k}, got {k_actual}"
        
        # With all qubits logical and no stabilizers, distance should be 1
        if k_actual > 0:
            assert d_actual >= 1, f"All-I: distance should be ≥1 when k>0, got {d_actual}"

    @given(n=st.integers(min_value=3, max_value=8))
    def test_repetition_codes(self, n):
        """Test X and Z repetition codes."""
        # For lattice-based construction, repetition codes behave differently
        # than single-generator codes
        
        # X-repetition lattice
        x_seed = 'X' * n
        if seed_is_valid(x_seed):
            n_x, k_x, d_x = lattice_based_analysis(x_seed)
            # All X's: should give rank=1 (all generators are identical)
            # So k = n - 1
            assert k_x == n - 1, f"X-repetition lattice should have k={n-1}, got {k_x}"
        
        # Z-repetition lattice  
        z_seed = 'Z' * n
        if seed_is_valid(z_seed):
            n_z, k_z, d_z = lattice_based_analysis(z_seed)
            # All Z's: should give rank=1 (all generators are identical) 
            # So k = n - 1
            assert k_z == n - 1, f"Z-repetition lattice should have k={n-1}, got {k_z}"


# Run with verbose output for debugging
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--hypothesis-show-statistics"]) 