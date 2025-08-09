#!/usr/bin/env python3
"""
orthogonality.py

Check self-orthogonality of a stabilizer seed string using Laurent polynomials over GF(2).

The condition is:
    s_X(x) * s_Z(x^-1) + s_Z(x) * s_X(x^-1) ≡ 0 (mod x^N - 1)

Usage:
    python orthogonality.py <seed_string>

Example:
    python orthogonality.py IIIIIIIIIXZZYZZZZYZZXIIIIIIII  # prints: Orthogonal? True
"""
import argparse
from typing import Dict, Any


def _reduce_exponents(coeffs: Dict[int, int], N: int) -> Dict[int, int]:
    """
    Normalize coefficients to GF(2) and reduce exponents modulo N.
    """
    return {exp % N: 1 for exp, c in coeffs.items() if c % 2}


class LaurentPolynomialGF2:
    """
    Laurent polynomial over GF(2) modulo (x^N - 1).
    Represents polynomials sum_i c_i x^i with c_i in {0,1}.
    """

    def __init__(self, coeffs: Dict[int, int], N: int) -> None:
        self.N = N
        self.coeffs = _reduce_exponents(coeffs, N)

    @classmethod
    def from_pauli_part(cls, paulis: str, part: str) -> "LaurentPolynomialGF2":
        """
        Build polynomial from X or Z part of a Pauli string.

        Args:
            paulis: Seed string of Pauli letters (e.g. 'IXZY').
            part: 'X' to extract X or Y positions, 'Z' to extract Z or Y.
        """
        N = len(paulis)
        coeffs: Dict[int, int] = {}
        for i, p in enumerate(paulis):
            if (part == 'X' and p in 'XY') or (part == 'Z' and p in 'ZY'):
                coeffs[i] = 1
        return cls(coeffs, N)

    def __add__(self, other: "LaurentPolynomialGF2") -> "LaurentPolynomialGF2":
        if self.N != other.N:
            raise ValueError("Mismatched ring sizes for addition")
        result = self.coeffs.copy()
        for exp in other.coeffs:
            result[exp] = result.get(exp, 0) ^ 1
            if result[exp] == 0:
                del result[exp]
        return LaurentPolynomialGF2(result, self.N)

    def __mul__(self, other: "LaurentPolynomialGF2") -> "LaurentPolynomialGF2":
        if self.N != other.N:
            raise ValueError("Mismatched ring sizes for multiplication")
        result: Dict[int,int] = {}
        for a in self.coeffs:
            for b in other.coeffs:
                exp = (a + b) % self.N
                result[exp] = result.get(exp, 0) ^ 1
                if result[exp] == 0:
                    del result[exp]
        return LaurentPolynomialGF2(result, self.N)

    def inverse_powers(self) -> "LaurentPolynomialGF2":
        """
        Compute the inverse-power polynomial: x^k -> x^(N-k) mod (x^N - 1).
        """
        inv = {(-exp) % self.N: 1 for exp in self.coeffs}
        return LaurentPolynomialGF2(inv, self.N)

    def is_zero(self) -> bool:
        return not self.coeffs

    def __str__(self) -> str:
        if not self.coeffs:
            return "0"
        terms = []
        for exp in sorted(self.coeffs):
            if exp == 0:
                terms.append('1')
            elif exp == 1:
                terms.append('x')
            else:
                terms.append(f'x^{exp}')
        return ' + '.join(terms)


def is_self_orthogonal(seed: str, verbose: bool=False) -> bool:
    """
    Return True if the stabilizer seed is self-orthogonal.
    Optionally prints intermediate Laurent polynomials when verbose.
    """
    sX = LaurentPolynomialGF2.from_pauli_part(seed, 'X')
    sZ = LaurentPolynomialGF2.from_pauli_part(seed, 'Z')
    sZ_inv = sZ.inverse_powers()
    sX_inv = sX.inverse_powers()
    term1 = sX * sZ_inv
    term2 = sZ * sX_inv
    orth_poly = term1 + term2
    if verbose:
        print(f"sX(x)      = {sX}")
        print(f"sZ(x)      = {sZ}")
        print(f"sZ(x^-1)   = {sZ_inv}")
        print(f"sX(x^-1)   = {sX_inv}")
        print(f"term1      = {term1}")
        print(f"term2      = {term2}")
        print(f"orth_poly  = {orth_poly}")
    return orth_poly.is_zero()


def main():
    parser = argparse.ArgumentParser(
        description="Check self-orthogonality of a stabilizer seed string."
    )
    parser.add_argument('seed', help='Pauli seed string (e.g., "XZIY")')
    parser.add_argument('-v', '--verbose', action='store_true', help='show intermediate polynomials')
    args = parser.parse_args()

    ortho = is_self_orthogonal(args.seed, args.verbose)
    print(f"Orthogonal? {ortho}")


if __name__ == '__main__':
    main()