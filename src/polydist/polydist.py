#!/usr/bin/env python3
"""
polydist.py

Orthogonality, logical-qubit counting, and logical operator extraction
for a translationally-invariant stabilizer seed.

Operations:
  - Check self-orthogonality via Laurent polynomials modulo (x^N - 1).
  - Compute number of logical qubits k = deg(gcd(s_X, s_Z, x^N - 1)).
  - Extract logical Z and X operators via extended Euclidean algorithm.

Usage:
  python polydist.py <seed_string> [-v]

Example:
  python polydist.py IIIIIIIIIXZZYZZZZYZZXIIIIIIII -v
"""
import argparse
from typing import Dict, Tuple, List

# ------- Laurent polynomial class over GF(2) mod (x^N - 1) -------
class LaurentPolynomialGF2:
    def __init__(self, coeffs: Dict[int,int], N: int) -> None:
        self.N = N
        self.coeffs = {e % N:1 for e,c in coeffs.items() if c & 1}
    @classmethod
    def from_pauli_part(cls, paulis: str, part: str) -> 'LaurentPolynomialGF2':
        N = len(paulis)
        coeffs: Dict[int,int] = {}
        for i,p in enumerate(paulis):
            if (part=='X' and p in 'XY') or (part=='Z' and p in 'ZY'):
                coeffs[i] = 1
        return cls(coeffs, N)
    def __add__(self, other:'LaurentPolynomialGF2') -> 'LaurentPolynomialGF2':
        if self.N != other.N: raise ValueError('Mismatched N')
        r = self.coeffs.copy()
        for e in other.coeffs:
            r[e] = r.get(e,0) ^ 1
            if not r[e]: r.pop(e)
        return LaurentPolynomialGF2(r, self.N)
    def __mul__(self, other:'LaurentPolynomialGF2') -> 'LaurentPolynomialGF2':
        if self.N != other.N: raise ValueError('Mismatched N')
        r:Dict[int,int]={}
        for a in self.coeffs:
            for b in other.coeffs:
                e=(a+b)%self.N
                r[e]=r.get(e,0)^1
                if not r[e]: r.pop(e)
        return LaurentPolynomialGF2(r, self.N)
    def inverse_powers(self) -> 'LaurentPolynomialGF2':
        inv={(-e)%self.N:1 for e in self.coeffs}
        return LaurentPolynomialGF2(inv, self.N)
    def is_zero(self)->bool:
        return not self.coeffs
    def __str__(self)->str:
        if not self.coeffs: return '0'
        terms=[]
        for e in sorted(self.coeffs):
            terms.append('1' if e==0 else 'x' if e==1 else f'x^{e}')
        return ' + '.join(terms)

# ------- Dict-based polynomial utils over GF(2) -------
def _add_poly(a:Dict[int,int], b:Dict[int,int]) -> Dict[int,int]:
    r=a.copy()
    for e in b:
        r[e]=r.get(e,0)^1
        if not r[e]: r.pop(e)
    return r

def _mul_poly(a:Dict[int,int], b:Dict[int,int]) -> Dict[int,int]:
    r:Dict[int,int]={}
    for e1 in a:
        for e2 in b:
            e=e1+e2
            r[e]=r.get(e,0)^1
            if not r[e]: r.pop(e)
    return r

def _div_poly(a:Dict[int,int], b:Dict[int,int]) -> Dict[int,int]:
    # polynomial quotient a//b in GF(2)[x]
    if not b: raise ValueError('divide by zero')
    r=a.copy(); q:Dict[int,int]={}
    d_b=max(b)
    while r and max(r)>=d_b:
        shift=max(r)-d_b
        q[shift]=1
        # subtract b*x^shift
        for e in list(b):
            ee=e+shift
            r[ee]=r.get(ee,0)^1
            if not r[ee]: r.pop(ee)
    return q

def _mod_poly(a:Dict[int,int], b:Dict[int,int]) -> Dict[int,int]:
    # polynomial remainder a mod b
    if not b: return a.copy()
    r=a.copy(); d_b=max(b)
    while r and max(r)>=d_b:
        shift=max(r)-d_b
        for e in list(b):
            ee=e+shift
            r[ee]=r.get(ee,0)^1
            if not r[ee]: r.pop(ee)
    return r

def _ext_gcd(a:Dict[int,int], b:Dict[int,int]) -> Tuple[Dict[int,int],Dict[int,int],Dict[int,int]]:
    # returns (g, s, t) with s*a + t*b = g
    if not a: return b.copy(), {0:1}, {}
    if not b: return a.copy(), {0:1}, {}
    r0, r1 = a.copy(), b.copy()
    s0, s1 = {0:1}, {}
    t0, t1 = {}, {0:1}
    while r1:
        q=_div_poly(r0, r1)
        r0, r1 = r1, _mod_poly(r0, r1)
        s0, s1 = s1, _add_poly(s0, _mul_poly(q, s1))
        t0, t1 = t1, _add_poly(t0, _mul_poly(q, t1))
    return r0, s0, t0

# ------- Core operations -------
def is_self_orthogonal(seed:str, verbose:bool=False)->bool:
    sX=LaurentPolynomialGF2.from_pauli_part(seed,'X')
    sZ=LaurentPolynomialGF2.from_pauli_part(seed,'Z')
    term1=sX * sZ.inverse_powers()
    term2=sZ * sX.inverse_powers()
    ortho=term1+term2
    if verbose:
        print(f'sX(x)      = {sX}')
        print(f'sZ(x)      = {sZ}')
        print(f'term1      = {term1}')
        print(f'term2      = {term2}')
        print(f'orth_poly  = {ortho}')
    return ortho.is_zero()

def compute_logical_qubits(seed:str)->int:
    N=len(seed)
    sX=LaurentPolynomialGF2.from_pauli_part(seed,'X').coeffs
    sZ=LaurentPolynomialGF2.from_pauli_part(seed,'Z').coeffs
    xN1={N:1,0:1}
    g1=_ext_gcd(sX,sZ)[0]
    gf=_ext_gcd(g1,xN1)[0]
    return max(gf) if gf else 0

def extract_logical_z(seed:str)->Tuple[Dict[int,int],Dict[int,int]]:
    # returns (v,u) dicts for Z operator v(x), u(x)
    N=len(seed)
    sX=LaurentPolynomialGF2.from_pauli_part(seed,'X').coeffs
    sZ=LaurentPolynomialGF2.from_pauli_part(seed,'Z').coeffs
    # stage1: egcd(sX,sZ)
    g1,u1,v1=_ext_gcd(sX,sZ)
    # stage2: egcd(g1,x^N-1)
    xN1={N:1,0:1}
    f,w,_=_ext_gcd(g1,xN1)
    # combine
    u=_mul_poly(w,u1)
    v=_mul_poly(w,v1)
    return v,u

def extract_logical_x(v:Dict[int,int],u:Dict[int,int],N:int)->Tuple[Dict[int,int],Dict[int,int]]:
    # LZ=(v,u), returns LX=(r,t)
    # reverse polynomials
    vrev={(-e)%N:1 for e in v}
    urev={(-e)%N:1 for e in u}
    # egcd on full ring
    g, rp,tp = _ext_gcd(vrev,urev)
    # ensure g=1 by gcd with x^N-1
    xN1={N:1,0:1}
    alpha,_ ,_ = _ext_gcd(g,xN1)
    r_prime=_mul_poly(alpha,rp)
    t_prime=_mul_poly(alpha,tp)
    # reverse back
    r={(-e)%N:1 for e in t_prime}
    t={(-e)%N:1 for e in r_prime}
    return r,t

# ------- CLI -------
def main():
    p=argparse.ArgumentParser()
    p.add_argument('seed', help='Pauli seed string')
    p.add_argument('-v','--verbose',action='store_true',help='verbose')
    args=p.parse_args()
    seed=args.seed
    print(f'Orthogonal? {is_self_orthogonal(seed,args.verbose)}')
    k=compute_logical_qubits(seed)
    print(f'Logical qubits k = {k}')
    if k>0:
        v,u=extract_logical_z(seed)
        print('Logical Z operator:')
        print(f'  v(x) = {LaurentPolynomialGF2(v,len(seed))}')
        print(f'  u(x) = {LaurentPolynomialGF2(u,len(seed))}')
        # extract X
        r,t=extract_logical_x(v,u,len(seed))
        print('Logical X operator:')
        print(f'  r(x) = {LaurentPolynomialGF2(r,len(seed))}')
        print(f'  t(x) = {LaurentPolynomialGF2(t,len(seed))}')

if __name__=='__main__':
    main()
