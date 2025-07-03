"""
Unit tests for polynomial linear algebra operations.

These tests verify the core polynomial matrix operations before using them
in the main stabilizer code functions.
"""

import pytest
import numpy as np
from polydist.polynomial import (
    LaurentPolynomial, rref_poly, nullspace_poly, rank_poly,
    pauli_to_polynomial_vector, polynomial_symplectic_inner_product
)


class TestLaurentPolynomial:
    """Test basic Laurent polynomial operations."""
    
    def test_polynomial_creation(self):
        """Test creating polynomials."""
        # Zero polynomial
        p0 = LaurentPolynomial.zero()
        assert p0.is_zero()
        assert p0.weight() == 0
        
        # Constant polynomial 1
        p1 = LaurentPolynomial.one()
        assert not p1.is_zero()
        assert p1.weight() == 1
        assert p1.degree() == 0
        
        # Monomial x^2
        p2 = LaurentPolynomial.x(2)
        assert p2.weight() == 1
        assert p2.degree() == 2
        assert p2.min_degree() == 2
    
    def test_polynomial_arithmetic(self):
        """Test polynomial addition and multiplication."""
        # (1 + x) + (x + x^2) = 1 + x^2 (mod 2)
        p1 = LaurentPolynomial({0: 1, 1: 1})  # 1 + x
        p2 = LaurentPolynomial({1: 1, 2: 1})  # x + x^2
        p3 = p1 + p2
        
        expected = LaurentPolynomial({0: 1, 2: 1})  # 1 + x^2
        assert p3 == expected
        assert p3.weight() == 2
        
        # (1 + x) * (x + x^2) = x + x^2 + x^2 + x^3 = x + x^3 (mod 2)
        p4 = p1 * p2
        expected_mult = LaurentPolynomial({1: 1, 3: 1})  # x + x^3
        assert p4 == expected_mult
        assert p4.weight() == 2
    
    def test_negative_exponents(self):
        """Test polynomials with negative exponents."""
        # x^(-1) + 1 + x
        p = LaurentPolynomial({-1: 1, 0: 1, 1: 1})
        assert p.degree() == 1
        assert p.min_degree() == -1
        assert p.weight() == 3
        
        # Multiplication with negative exponents
        p1 = LaurentPolynomial({-1: 1})  # x^(-1)
        p2 = LaurentPolynomial({1: 1})   # x
        p3 = p1 * p2  # Should be 1
        assert p3 == LaurentPolynomial.one()


class TestPolynomialMatrixOperations:
    """Test matrix operations with polynomials."""
    
    def test_simple_rref(self):
        """Test RREF on simple polynomial matrices."""
        # Test with 1x1 matrix [1]
        matrix = [[LaurentPolynomial.one()]]
        result = rref_poly(matrix)
        expected = [[LaurentPolynomial.one()]]
        assert result == expected
    
    @pytest.mark.skip(reason="Will implement after basic functions")
    def test_2x2_rref(self):
        """Test RREF on 2x2 polynomial matrix."""
        # Test with [[1, x], [x, 1]]
        matrix = [
            [LaurentPolynomial.one(), LaurentPolynomial.x()],
            [LaurentPolynomial.x(), LaurentPolynomial.one()]
        ]
        result = rref_poly(matrix)
        # Expected result depends on implementation details
        assert len(result) == 2
        assert len(result[0]) == 2
    
    def test_simple_rank(self):
        """Test rank computation on simple matrices."""
        # Rank of 1x1 identity matrix should be 1
        matrix = [[LaurentPolynomial.one()]]
        rank = rank_poly(matrix)
        assert rank == 1
        
        # Rank of 1x1 zero matrix should be 0
        matrix = [[LaurentPolynomial.zero()]]
        rank = rank_poly(matrix)
        assert rank == 0
    
    @pytest.mark.skip(reason="Will implement after basic functions")
    def test_simple_nullspace(self):
        """Test nullspace computation."""
        # For 1x1 zero matrix, nullspace should have dimension 1
        matrix = [[LaurentPolynomial.zero()]]
        null_basis = nullspace_poly(matrix)
        assert len(null_basis) == 1  # Dimension of nullspace
        assert len(null_basis[0]) == 1  # Vector length


class TestPauliPolynomialConversion:
    """Test conversion between Pauli operators and polynomial vectors."""
    
    def test_pauli_to_polynomial(self):
        """Test converting Pauli operators to polynomial vectors."""
        # Single qubit cases
        pauli_i = ['I']
        pauli_x = ['X']
        pauli_z = ['Z']
        pauli_y = ['Y']
        
        vec_i = pauli_to_polynomial_vector(pauli_i)
        vec_x = pauli_to_polynomial_vector(pauli_x)
        vec_z = pauli_to_polynomial_vector(pauli_z)
        vec_y = pauli_to_polynomial_vector(pauli_y)
        
        # Check that vectors have correct length (2L for L qubits)
        assert len(vec_i) == 2
        assert len(vec_x) == 2
        assert len(vec_z) == 2
        assert len(vec_y) == 2
    
    def test_polynomial_symplectic_product(self):
        """Test polynomial symplectic inner product."""
        # Test commutation relations
        # [X, Z] = -2iY ≠ 0, so X and Z should have non-zero symplectic product
        vec_x = pauli_to_polynomial_vector(['X'])
        vec_z = pauli_to_polynomial_vector(['Z'])
        
        product = polynomial_symplectic_inner_product(vec_x, vec_z)
        assert not product.is_zero()  # X and Z don't commute
        
        # [X, X] = 0, so X should commute with itself
        product_xx = polynomial_symplectic_inner_product(vec_x, vec_x)
        assert product_xx.is_zero()  # X commutes with itself 