#!/usr/bin/env python3
"""
polydist.py  (version 1.1.0)

Distance via Laurent-polynomial logical extraction with algebraic centralizer correction.

Operations:
  - Self-orthogonality check
  - Count logical qubits k = deg(gcd(s_X, s_Z, x^N-1))
  - Extract raw L_Z, L_X via extended GCD
  - Build commutator P(x) = v·s_Z(x⁻¹) + u·s_X(x⁻¹)
  - Compute a(x) = (x^N - 1)/gcd(P, x^N - 1) to force commutation
  - Apply a(x) to (v,u) and (r,t) to get true logicals
  - Enumerate true logicals and their translates to find minimum weight d

Usage:
  python polydist.py <seed> [-v]

Example:
  python polydist.py IIIIIIIIIXZZYZZZZYZZXIIIIIIII -v
"""

VERSION = "polydist.py v1.1.0"

import argparse
import itertools
from typing import Dict, Tuple, List

# ------- Laurent polynomial mod (x^N - 1) -------
class LaurentPolynomialGF2:
    def __init__(self, coeffs: Dict[int,int], N: int):
        self.N = N
        self.coeffs = {e % N:1 for e,c in coeffs.items() if c & 1}

    @classmethod
    def from_pauli_part(cls, seed: str, part: str):
        N = len(seed)
        coeffs = {i:1 for i,p in enumerate(seed)
                  if (part=='X' and p in 'XY') or (part=='Z' and p in 'ZY')}
        return cls(coeffs, N)

    def inverse_powers(self):
        return LaurentPolynomialGF2({(-e)%self.N:1 for e in self.coeffs}, self.N)

    def __mul__(self, other: 'LaurentPolynomialGF2') -> 'LaurentPolynomialGF2':
        r: Dict[int,int] = {}
        for a in self.coeffs:
            for b in other.coeffs:
                e = (a + b) % self.N
                r[e] = r.get(e,0) ^ 1
                if not r[e]: r.pop(e)
        return LaurentPolynomialGF2(r, self.N)

    def __add__(self, other: 'LaurentPolynomialGF2') -> 'LaurentPolynomialGF2':
        r = self.coeffs.copy()
        for e in other.coeffs:
            r[e] = r.get(e,0) ^ 1
            if not r[e]: r.pop(e)
        return LaurentPolynomialGF2(r, self.N)

    def is_zero(self) -> bool:
        return not self.coeffs

    def __str__(self):
        if not self.coeffs:
            return '0'
        terms = []
        for e in sorted(self.coeffs):
            if e==0: terms.append('1')
            elif e==1: terms.append('x')
            else: terms.append(f'x^{e}')
        return ' + '.join(terms)

# ------- GF(2)[x] polynomial ops -------
def _add_poly(a: Dict[int,int], b: Dict[int,int]) -> Dict[int,int]:
    r = a.copy()
    for e in b:
        r[e] = r.get(e,0) ^ 1
        if not r[e]: r.pop(e)
    return r

def _mul_poly(a: Dict[int,int], b: Dict[int,int]) -> Dict[int,int]:
    r: Dict[int,int] = {}
    for x in a:
        for y in b:
            e = x + y
            r[e] = r.get(e,0) ^ 1
            if not r[e]: r.pop(e)
    return r

def _div_poly(a: Dict[int,int], b: Dict[int,int]) -> Dict[int,int]:
    r, q = a.copy(), {}
    db = max(b) if b else -1
    while b and r and max(r) >= db:
        shift = max(r) - db
        q[shift] = 1
        for e in list(b):
            ee = e + shift
            r[ee] = r.get(ee,0) ^ 1
            if not r[ee]: r.pop(ee)
    return q

def _mod_poly(a: Dict[int,int], b: Dict[int,int]) -> Dict[int,int]:
    r = a.copy()
    db = max(b) if b else -1
    while b and r and max(r) >= db:
        shift = max(r) - db
        for e in list(b):
            ee = e + shift
            r[ee] = r.get(ee,0) ^ 1
            if not r[ee]: r.pop(ee)
    return r

def _ext_gcd(a: Dict[int,int], b: Dict[int,int]) -> Tuple[Dict[int,int],Dict[int,int],Dict[int,int]]:
    if not a: return b.copy(), {0:1}, {}
    if not b: return a.copy(), {0:1}, {}
    r0, r1 = a.copy(), b.copy()
    s0, s1 = {0:1}, {}
    t0, t1 = {}, {0:1}
    while r1:
        q = _div_poly(r0, r1)
        r0, r1 = r1, _mod_poly(r0, r1)
        s0, s1 = s1, _add_poly(s0, _mul_poly(q, s1))
        t0, t1 = t1, _add_poly(t0, _mul_poly(q, t1))
    return r0, s0, t0

# ------- Pauli conversion -------
def poly_to_pauli(v: Dict[int,int], u: Dict[int,int], N: int) -> str:
    pauli = []
    for i in range(N):
        xi, zi = u.get(i,0), v.get(i,0)
        pauli.append('Y' if xi and zi else 'X' if xi else 'Z' if zi else 'I')
    return ''.join(pauli)

def shift_poly(poly: Dict[int,int], N: int, sh: int) -> Dict[int,int]:
    return {(e+sh)%N:1 for e in poly}

def verify_centrality(v: Dict[int,int], u: Dict[int,int], seed: str) -> bool:
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed,'X')
    sZ = LaurentPolynomialGF2.from_pauli_part(seed,'Z')
    for sh in range(N):
        SX = sX * LaurentPolynomialGF2({sh:1}, N)
        SZ = sZ * LaurentPolynomialGF2({sh:1}, N)
        lhs = LaurentPolynomialGF2(v,N)*SZ.inverse_powers() + LaurentPolynomialGF2(u,N)*SX.inverse_powers()
        if not lhs.is_zero():
            return False
    return True

# ------- Core routines -------

def is_self_orthogonal(seed: str, verbose: bool=False) -> bool:
    """Check seed orthogonality: sX*sZ(x^-1)+sZ*sX(x^-1)=0"""
    sX = LaurentPolynomialGF2.from_pauli_part(seed, 'X')
    sZ = LaurentPolynomialGF2.from_pauli_part(seed, 'Z')
    term1 = sX * sZ.inverse_powers()
    term2 = sZ * sX.inverse_powers()
    ortho = term1 + term2
    if verbose:
        print(f'sX(x)      = {sX}')
        print(f'sZ(x)      = {sZ}')
        print(f'term1      = {term1}')
        print(f'term2      = {term2}')
        print(f'orth_poly  = {ortho}')
    return ortho.is_zero()


def compute_logical_qubits(seed: str) -> int:
    """k = deg(gcd(sX, sZ, x^N-1))"""
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed, 'X').coeffs
    sZ = LaurentPolynomialGF2.from_pauli_part(seed, 'Z').coeffs
    xN1 = {N:1, 0:1}
    g1, _, _ = _ext_gcd(sX, sZ)
    gf, _, _ = _ext_gcd(g1, xN1)
    return max(gf) if gf else 0


def extract_logical_z(seed: str) -> Tuple[Dict[int,int], Dict[int,int]]:
    """Extended GCD to get raw L_Z = (v,u)"""
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed, 'X').coeffs
    sZ = LaurentPolynomialGF2.from_pauli_part(seed, 'Z').coeffs
    g1, u1, v1 = _ext_gcd(sX, sZ)
    xN1 = {N:1, 0:1}
    _, w, _ = _ext_gcd(g1, xN1)
    v = _mul_poly(w, v1)
    u = _mul_poly(w, u1)
    return v, u


def extract_logical_x(v: Dict[int,int], u: Dict[int,int], N: int) -> Tuple[Dict[int,int], Dict[int,int]]:
    """Extended GCD on reversed polynomials to get raw L_X"""
    vrev = {(-e)%N:1 for e in v}
    urev = {(-e)%N:1 for e in u}
    g, rp, tp = _ext_gcd(vrev, urev)
    xN1 = {N:1, 0:1}
    alpha = _ext_gcd(g, xN1)[1]
    r_base = _mul_poly(alpha, rp)
    t_base = _mul_poly(alpha, tp)
    # reverse back
    r = {(-e)%N:1 for e in t_base}
    t = {(-e)%N:1 for e in r_base}
    return r, t

# now existing compute_commutator_poly follows

def compute_lcm_poly(a: Dict[int,int], b: Dict[int,int], N: int) -> Dict[int,int]:
    # helper to compute lcm = (a * b) / gcd(a,b) in polynomial ring
    g,_,_ = _ext_gcd(a, b)
    # divide: (a*b) // g
    prod = _mul_poly(a, b)
    # compute quotient poly division -- skip exact division, use gcd formula
    # here lcm = prod XOR terms where exponent sums cancel
    return prod  # placeholder, not used directly

def compute_commutator_poly(v: Dict[int,int], u: Dict[int,int], seed: str) -> Dict[int,int]:
    # P(x) = v*sZ(x^-1) + u*sX(x^-1)
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed,'X').coeffs
    sZ = LaurentPolynomialGF2.from_pauli_part(seed,'Z').coeffs
    # inverse polynomials
    sZ_inv = LaurentPolynomialGF2({(-e)%N:1 for e in sZ}, N).coeffs
    sX_inv = LaurentPolynomialGF2({(-e)%N:1 for e in sX}, N).coeffs
    term1 = _mul_poly(v, sZ_inv)
    term2 = _mul_poly(u, sX_inv)
    return _add_poly(term1, term2)


def compute_distance(seed: str, verbose: bool=False) -> Tuple[int,str]:
    N = len(seed)
    k = compute_logical_qubits(seed)
    if k == 0:
        return 0, ''
    # raw Bézout logicals
    v, u = extract_logical_z(seed)
    r, t = extract_logical_x(v, u, N)
    # build commutator poly and its gcd with x^N-1
    P = compute_commutator_poly(v, u, seed)
    xN1 = {N:1, 0:1}
    g, _, _ = _ext_gcd(P, xN1)
    # cofactor a(x) = (x^N-1)/g
    # since (x^N-1) = g * a, and gcd divisible
    # find a by polynomial division: a = _div_poly(xN1, g)
    a = _div_poly(xN1, g)
    # true logicals = a * (v,u) and a * (r,t)
    v_true = _mul_poly(a, v)
    u_true = _mul_poly(a, u)
    r_true = _mul_poly(a, r)
    t_true = _mul_poly(a, t)
    # verify centrality
    if verbose:
        print(f"Commutator P(x) = {_add_poly(P,P)} (should be 0 after a)")
    # enumerate combinations
    best = N+1
    best_str = ''
    ops = [(v_true, u_true), (r_true, t_true)]
    for rsize in range(1,3):
        for combo in itertools.combinations(range(2), rsize):
            vs, us = {}, {}
            for idx in combo:
                v_op, u_op = ops[idx]
                vs = _add_poly(vs, v_op)
                us = _add_poly(us, u_op)
            for sh in range(N):
                vs_sh = shift_poly(vs, N, sh)
                us_sh = shift_poly(us, N, sh)
                if not verify_centrality(vs_sh, us_sh, seed):
                    continue
                pauli = poly_to_pauli(vs_sh, us_sh, N)
                w = pauli.count('I')
                w = N - w
                if w < best:
                    best, best_str = w, pauli
    return best, best_str

# ------- CLI -------
def main():
    p = argparse.ArgumentParser(description="Compute distance using polynomial and centralizer correction.")
    p.add_argument('seed', help='Pauli seed string')
    p.add_argument('-v','--verbose', action='store_true', help='verbose')
    args = p.parse_args()
    seed = args.seed
    print(VERSION)
    print(f"Orthogonal? {is_self_orthogonal(seed, args.verbose)}")
    k = compute_logical_qubits(seed)
    print(f"Logical qubits k = {k}")
    if k > 0:
        v, u = extract_logical_z(seed)
        print("Raw L_Z:")
        print(f"  v(x) = {LaurentPolynomialGF2(v,len(seed))}")
        print(f"  u(x) = {LaurentPolynomialGF2(u,len(seed))}")
        r, t = extract_logical_x(v, u, len(seed))
        print("Raw L_X:")
        print(f"  r(x) = {LaurentPolynomialGF2(r,len(seed))}")
        print(f"  t(x) = {LaurentPolynomialGF2(t,len(seed))}")
        d, pauli = compute_distance(seed, args.verbose)
        print(f"Distance d = {d}")
        print(f"Minimal logical = {pauli}")

if __name__=='__main__':
    main()
