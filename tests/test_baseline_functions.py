"""
Baseline tests for stabilizer functions using canonical fixtures.

These parametrized tests will run against both the original GF(2) implementation
and the future Laurent polynomial implementation to ensure bit-for-bit parity.
"""

import pytest
import numpy as np
from fixtures import CANONICAL_FIXTURES, get_valid_fixtures, get_invalid_fixtures

# Import the original GF(2) implementations
from stabilizer.utils import seed_is_valid, compute_entanglement
from stabilizer.distance import find_distance, find_logical_operators
from stabilizer.generator import parse_seed
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank


def build_stabilizer_tableau(seed: str):
    """Helper to build stabilizer tableau from seed string."""
    pauli_list = list(seed)
    stab_ops = build_lattice(pauli_list)
    tableau = build_tableau(stab_ops)
    return tableau


# Placeholder imports for polynomial implementations (will be added in Step 3)
# These will raise NotImplementedError initially
try:
    import polydist.utils
    import polydist.distance
    import polydist.tableau
    POLYDIST_AVAILABLE = True
except ImportError:
    POLYDIST_AVAILABLE = False


class TestSeedValidation:
    """Test seed_is_valid() function with canonical fixtures."""
    
    @pytest.mark.parametrize("fixture", CANONICAL_FIXTURES, ids=lambda f: f.name)
    def test_seed_validity_gf2(self, fixture):
        """Test seed validity using original GF(2) implementation."""
        result = seed_is_valid(fixture.seed)
        assert result == fixture.is_valid, \
            f"GF(2): Expected {fixture.is_valid}, got {result} for {fixture.seed}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", CANONICAL_FIXTURES, ids=lambda f: f.name)
    def test_seed_validity_poly(self, fixture):
        """Test seed validity using polynomial implementation."""
        result = polydist.utils.seed_is_valid(fixture.seed)
        assert result == fixture.is_valid, \
            f"Poly: Expected {fixture.is_valid}, got {result} for {fixture.seed}"


class TestLogicalOperators:
    """Test find_logical_operators() function with canonical fixtures."""
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_logical_operators_gf2(self, fixture):
        """Test logical operator extraction using original GF(2) implementation."""
        if fixture.k == 0:
            pytest.skip("No logical operators for k=0")
            
        tableau = build_stabilizer_tableau(fixture.seed)
        logical_ops = find_logical_operators(tableau)
        
        # Should return 2k logical operators
        assert len(logical_ops) == 2 * fixture.k, \
            f"GF(2): Expected {2 * fixture.k} logical ops, got {len(logical_ops)}"
        
        # Each logical operator should be a 2L-length binary vector
        L = len(fixture.seed)
        for i, op in enumerate(logical_ops):
            assert len(op) == 2 * L, \
                f"GF(2): Logical op {i} has wrong length {len(op)}, expected {2 * L}"
            assert np.all((op == 0) | (op == 1)), \
                f"GF(2): Logical op {i} has non-binary entries"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_logical_operators_poly(self, fixture):
        """Test logical operator extraction using polynomial implementation."""
        if fixture.k == 0:
            pytest.skip("No logical operators for k=0")

        # Build lattice from seed, then tableau from lattice
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        tableau = polydist.tableau.build_tableau_poly(stab_ops)
        logical_ops = polydist.distance.find_logical_operators(tableau)

        # Should return 2k logical operators
        assert len(logical_ops) == 2 * fixture.k, \
            f"Poly: Expected {2 * fixture.k} logical ops, got {len(logical_ops)}"


class TestCodeDistance:
    """Test find_distance() function with canonical fixtures."""
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_distance_gf2(self, fixture):
        """Test distance calculation using original GF(2) implementation."""
        tableau = build_stabilizer_tableau(fixture.seed)
        distance = find_distance(tableau)
        
        assert distance == fixture.d, \
            f"GF(2): Expected distance {fixture.d}, got {distance}"
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_distance_with_logical_ops_gf2(self, fixture):
        """Test distance calculation with logical operators using GF(2)."""
        tableau = build_stabilizer_tableau(fixture.seed)
        distance, logical_ops = find_distance(tableau, return_logical_ops=True)
        
        assert distance == fixture.d, \
            f"GF(2): Expected distance {fixture.d}, got {distance}"
        
        if fixture.k > 0:
            assert len(logical_ops) == 2 * fixture.k, \
                f"GF(2): Expected {2 * fixture.k} logical ops, got {len(logical_ops)}"
        else:
            assert len(logical_ops) == 0, \
                f"GF(2): Expected 0 logical ops for k=0, got {len(logical_ops)}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_distance_poly(self, fixture):
        """Test distance calculation using polynomial implementation."""
        # Build lattice from seed, then tableau from lattice
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        tableau = polydist.tableau.build_tableau_poly(stab_ops)
        distance = polydist.distance.find_distance(tableau)
        
        assert distance == fixture.d, \
            f"Poly: Expected distance {fixture.d}, got {distance}"


class TestBipartiteEntanglement:
    """Test compute_entanglement() function with canonical fixtures."""
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_entanglement_gf2(self, fixture):
        """Test entanglement calculation using original GF(2) implementation."""
        tableau = build_stabilizer_tableau(fixture.seed)
        
        if fixture.k > 0:
            logical_ops = find_logical_operators(tableau)
            entanglement = compute_entanglement(tableau, logical_ops)
        else:
            entanglement = compute_entanglement(tableau, None)
        
        assert entanglement == fixture.entanglement, \
            f"GF(2): Expected entanglement {fixture.entanglement}, got {entanglement}"
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_entanglement_without_logical_ops_gf2(self, fixture):
        """Test entanglement calculation without logical operators using GF(2)."""
        tableau = build_stabilizer_tableau(fixture.seed)
        entanglement = compute_entanglement(tableau, None)
        
        # Should be non-negative
        assert entanglement >= 0, \
            f"GF(2): Negative entanglement {entanglement}"
        
        # For codes with no logical qubits, should match expected
        if fixture.k == 0:
            assert entanglement == fixture.entanglement, \
                f"GF(2): Expected entanglement {fixture.entanglement}, got {entanglement}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_entanglement_poly(self, fixture):
        """Test entanglement calculation using polynomial implementation."""
        # Build lattice from seed, then tableau from lattice
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        tableau = polydist.tableau.build_tableau_poly(stab_ops)
        
        if fixture.k > 0:
            logical_ops = polydist.distance.find_logical_operators(tableau)
            entanglement = polydist.qca.compute_entanglement(tableau, logical_ops)
        else:
            entanglement = polydist.qca.compute_entanglement(tableau, None)
        
        assert entanglement == fixture.entanglement, \
            f"Poly: Expected entanglement {fixture.entanglement}, got {entanglement}"


class TestCodeParameters:
    """Test complete code parameter computation (n, k, d) with canonical fixtures."""
    
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_code_parameters_gf2(self, fixture):
        """Test complete code parameter calculation using GF(2)."""
        tableau = build_stabilizer_tableau(fixture.seed)
        rank = compute_rank(tableau)
        
        n = len(fixture.seed)
        k = n - rank
        
        if k > 0:
            d = find_distance(tableau)
        else:
            d = 0
        
        assert n == fixture.n, f"GF(2): Expected n={fixture.n}, got {n}"
        assert k == fixture.k, f"GF(2): Expected k={fixture.k}, got {k}"
        assert d == fixture.d, f"GF(2): Expected d={fixture.d}, got {d}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_code_parameters_poly(self, fixture):
        """Test complete code parameter calculation using polynomials."""
        # Build lattice from seed, then tableau from lattice
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        tableau = polydist.tableau.build_tableau_poly(stab_ops)
        rank = polydist.tableau.compute_rank_poly(tableau)
        
        n = len(fixture.seed)
        k = n - rank
        
        if k > 0:
            d = polydist.distance.find_distance(tableau)
            entanglement = polydist.qca.compute_entanglement(tableau)
        else:
            d = 0
            entanglement = 0
        
        assert n == fixture.n, f"Poly: Expected n={fixture.n}, got {n}"
        assert k == fixture.k, f"Poly: Expected k={fixture.k}, got {k}"
        assert d == fixture.d, f"Poly: Expected d={fixture.d}, got {d}"


class TestCrossValidation:
    """Cross-validation tests between GF(2) and polynomial implementations."""
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", CANONICAL_FIXTURES, ids=lambda f: f.name)
    def test_seed_validity_equivalence(self, fixture):
        """Test that both implementations agree on seed validity."""
        gf2_result = seed_is_valid(fixture.seed)
        poly_result = polydist.utils.seed_is_valid(fixture.seed)
        
        assert gf2_result == poly_result, \
            f"Validity mismatch for {fixture.seed}: GF(2)={gf2_result}, Poly={poly_result}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_distance_equivalence(self, fixture):
        """Test that both implementations give the same distance."""
        gf2_tableau = build_stabilizer_tableau(fixture.seed)
        gf2_distance = find_distance(gf2_tableau)
        
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        poly_tableau = polydist.tableau.build_tableau_poly(stab_ops)
        poly_distance = polydist.distance.find_distance(poly_tableau)
        
        assert gf2_distance == poly_distance, \
            f"Distance mismatch for {fixture.seed}: GF(2)={gf2_distance}, Poly={poly_distance}"
    
    @pytest.mark.skipif(not POLYDIST_AVAILABLE, reason="polydist module not available")
    @pytest.mark.parametrize("fixture", get_valid_fixtures(), ids=lambda f: f.name)
    def test_entanglement_equivalence(self, fixture):
        """Test that both implementations give the same entanglement."""
        gf2_tableau = build_stabilizer_tableau(fixture.seed)
        if fixture.k > 0:
            gf2_logical = find_logical_operators(gf2_tableau)
            gf2_ent = compute_entanglement(gf2_tableau, gf2_logical)
        else:
            gf2_ent = compute_entanglement(gf2_tableau, None)
        
        stab_ops = polydist.lattice.build_lattice(list(fixture.seed))
        poly_tableau = polydist.tableau.build_tableau_poly(stab_ops)
        if fixture.k > 0:
            poly_logical = polydist.distance.find_logical_operators(poly_tableau)
            poly_ent = polydist.qca.compute_entanglement(poly_tableau, poly_logical)
        else:
            poly_ent = polydist.qca.compute_entanglement(poly_tableau, None)
        
        assert gf2_ent == poly_ent, \
            f"Entanglement mismatch for {fixture.seed}: GF(2)={gf2_ent}, Poly={poly_ent}" 