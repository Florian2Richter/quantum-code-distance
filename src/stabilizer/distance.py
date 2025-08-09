#!/usr/bin/env python3
"""
polydist.py

Distance via Laurent-polynomial logical extraction with algebraic centralizer correction.

Operations:
  - Self-orthogonality check
  - Count logical qubits k = deg(gcd(s_X, s_Z, x^N-1))
  - Extract raw L_Z, L_X via extended GCD
  - Build commutator P(x) = v*s_Z(x^-1) + u*s_X(x^-1)
  - Compute a(x) = (x^N - 1)/gcd(P, x^N - 1) to force commutation
  - Apply a(x) to (v,u) and (r,t) to get true logicals
  - Enumerate combinations and shifts to find minimum weight d

Usage:
  python polydist.py <seed> [-v]

Example:
  python polydist.py IIIIIIIIIXZZYZZZZYZZXIIIIIIII -v
"""
import argparse
import itertools
from typing import Dict, Tuple

# ------- Laurent polynomial mod (x^N-1) -------
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
        r = {}
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
    r = {}
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

# ------- Main routines -------
def compute_logical_qubits(seed: str) -> int:
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed,'X').coeffs
    sZ = LaurentPolynomialGF2.from_pauli_part(seed,'Z').coeffs
    xN1 = {N:1,0:1}
    g1,_,_ = _ext_gcd(sX, sZ)
    gcd_full,_,_ = _ext_gcd(g1, xN1)
    return max(gcd_full) if gcd_full else 0

def extract_logical_z(seed: str) -> Tuple[Dict[int,int],Dict[int,int]]:
    N = len(seed)
    sX = LaurentPolynomialGF2.from_pauli_part(seed,'X').coeffs
    sZ = LaurentPolynomialGF2.from_pauli_part(seed,'Z').coeffs
    g1, u1, v1 = _ext_gcd(sX, sZ)
    xN1 = {N:1,0:1}
    _, w, _ = _ext_gcd(g1, xN1)
    return _mul_poly(w, v1), _mul_poly(w, u1)

def extract_logical_x(v: Dict[int,int], u: Dict[int,int], N: int) -> Tuple[Dict[int,int],Dict[int,int]]:
    vrev = {(-e)%N:1 for e in v}
    urev = {(-e)%N:1 for e in u}
    g, rp, tp = _ext_gcd(vrev, urev)
    xN1 = {N:1,0:1}
    alpha = _ext_gcd(g, xN1)[1]
    r_base = _mul_poly(alpha, rp)
    t_base = _mul_poly(alpha, tp)
    return {(-e)%N:1 for e in t_base}, {(-e)%N:1 for e in r_base}

# ------- CLI -------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('seed')
    p.add_argument('-v','--verbose', action='store_true')
    args = p.parse_args()
    seed = args.seed
    N = len(seed)

    # Orthogonality check
    sX = LaurentPolynomialGF2.from_pauli_part(seed,'X')
    sZ = LaurentPolynomialGF2.from_pauli_part(seed,'Z')
    ortho = (sX * sZ.inverse_powers()) + (sZ * sX.inverse_powers())
    print(f'Orthogonal? {ortho.is_zero()}')

    # Logical qubits
    k = compute_logical_qubits(seed)
    print(f'Logical qubits k = {k}')
    if k == 0:
        return

    # Raw Bézout logicals
    v, u = extract_logical_z(seed)
    r, t = extract_logical_x(v, u, N)
    if args.verbose:
        print('raw L_Z:', poly_to_pauli(v,u,N))
        print('raw L_X:', poly_to_pauli(r,t,N))

    # Commutator polynomial P(x)
    P = _add_poly(_mul_poly(v, LaurentPolynomialGF2.from_pauli_part(seed,'Z').inverse_powers().coeffs),
                  _mul_poly(u, LaurentPolynomialGF2.from_pauli_part(seed,'X').inverse_powers().coeffs))
    xN1 = {N:1,0:1}
    g = _ext_gcd(P, xN1)[0]
    a = _div_poly(xN1, g)
    if args.verbose:
        print('P(x) =', LaurentPolynomialGF2(P, N))
        print('g(x) =', LaurentPolynomialGF2(g, N))
        print('a(x) =', LaurentPolynomialGF2(a, N))

    # True logicals
    v_t = _mul_poly(a, v)
    u_t = _mul_poly(a, u)
    r_t = _mul_poly(a, r)
    t_t = _mul_poly(a, t)
    print('true L_Z:', poly_to_pauli(v_t,u_t,N))
    print('true L_X:', poly_to_pauli(r_t,t_t,N))

    # Verify centrality
    assert verify_centrality(v_t,u_t,seed)
    assert verify_centrality(r_t,t_t,seed)

    # Distance search
    best, best_str = N+1, ''
    for subset in itertools.chain.from_iterable(itertools.combinations([('Z',v_t,u_t),('X',r_t,t_t)], r) for r in range(1,3)):
        vs, us = {}, {}
        for _, vv, uu in subset:
            vs = _add_poly(vs, vv)
            us = _add_poly(us, uu)
        for sh in range(N):
            vs_sh = shift_poly(vs,N,sh)
            us_sh = shift_poly(us,N,sh)
            if not verify_centrality(vs_sh,us_sh,seed):
                continue
            pstr = poly_to_pauli(vs_sh,us_sh,N)
            w = N - pstr.count('I')
            if w < best:
                best, best_str = w, pstr
    print(f'Distance d = {best}')
    print(f'Minimal logical operator = {best_str}')

if __name__=='__main__':
    main()
