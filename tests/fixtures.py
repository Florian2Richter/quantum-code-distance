"""
Canonical fixtures for test-driven development of Laurent polynomial formalism.

Contains small stabilizer codes with known properties (validity, logical qubits, 
distance, logical operators, entanglement) to ensure bit-for-bit parity between
GF(2) and polynomial implementations.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Optional, Any
from stabilizer.utils import seed_is_valid, compute_entanglement
from stabilizer.generator import parse_seed
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance, find_logical_operators


@dataclass
class CodeFixture:
    """A canonical stabilizer code example with known properties."""
    name: str
    seed: str
    is_valid: bool
    n: int  # physical qubits
    k: int  # logical qubits  
    d: int  # code distance
    entanglement: int  # bipartite entanglement for half-chain cut
    logical_ops: Optional[List[np.ndarray]] = None
    description: str = ""


def _compute_properties(seed: str) -> tuple:
    """Helper to compute all properties of a seed for fixture validation."""
    if not seed_is_valid(seed):
        return False, len(seed), 0, 0, 0, None
    
    pauli_list = list(seed)
    stab_ops = build_lattice(pauli_list)
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    
    n = len(seed)
    k = n - rank
    
    if k == 0:
        d = 0
        logical_ops = []
    else:
        d, logical_ops = find_distance(tableau, return_logical_ops=True)
    
    ent = compute_entanglement(tableau, logical_ops if logical_ops else None)
    
    return True, n, k, d, ent, logical_ops


# Canonical test fixtures
CANONICAL_FIXTURES = [
    # 1. Simple repetition codes
    CodeFixture(
        name="z_repetition_4",
        seed="ZZII", 
        is_valid=True,
        n=4, k=1, d=1, entanglement=1,
        description="Z-repetition code on 4 qubits, distance 1"
    ),
    
    CodeFixture(
        name="x_repetition_3", 
        seed="XXI",
        is_valid=True,
        n=3, k=1, d=1, entanglement=1,
        description="X-repetition code on 3 qubits, distance 1"
    ),
    
    # 2. Higher distance codes
    CodeFixture(
        name="distance_2_code",
        seed="XZIY",
        is_valid=True, 
        n=4, k=1, d=2, entanglement=2,
        description="Distance-2 code from golden data"
    ),
    
    CodeFixture(
        name="alternating_xz",
        seed="XZZX", 
        is_valid=True,
        n=4, k=1, d=2, entanglement=2,
        description="Alternating XZ pattern, distance 2"
    ),
    
    # 3. Invalid seed (non-commuting)
    CodeFixture(
        name="invalid_xzy",
        seed="XZY",
        is_valid=False,
        n=3, k=0, d=0, entanglement=0,
        description="XZY seed produces non-commuting stabilizers"
    ),
    
    # 4. All identity (trivial stabilizers)
    CodeFixture(
        name="all_identity",
        seed="IIII",
        is_valid=True,
        n=4, k=4, d=1, entanglement=0,
        description="All-I seed gives n logical qubits"
    ),
    
    # 5. Two-logical code
    CodeFixture(
        name="two_logical",
        seed="IIIII",
        is_valid=True,
        n=5, k=5, d=1, entanglement=0,
        description="5-qubit all-I gives 5 logical qubits"
    ),
    
    # 6. Invalid seed (non-commuting)
    CodeFixture(
        name="invalid_xyz",
        seed="XYZ",
        is_valid=False,
        n=3, k=0, d=0, entanglement=0,
        description="XYZ seed produces non-commuting stabilizers"
    ),
    
    # 7. No logical qubits (fully constrained)
    CodeFixture(
        name="no_logical_7qubit",
        seed="IIIXIII",
        is_valid=True,
        n=7, k=0, d=0, entanglement=0,
        description="Single X in sea of I's produces fully constrained code"
    ),
    
    # 8. CSS-type code
    CodeFixture(
        name="css_small",
        seed="ZZXX",
        is_valid=True, 
        n=4, k=1, d=2, entanglement=2,
        description="CSS-like code with Z and X blocks, distance 2"
    ),
]


def validate_fixtures():
    """Validate that all fixtures match their computed properties."""
    print("Validating canonical fixtures...")
    
    for fixture in CANONICAL_FIXTURES:
        if not fixture.is_valid:
            # Just check that it's indeed invalid
            assert not seed_is_valid(fixture.seed), f"Fixture {fixture.name} should be invalid"
            print(f"✓ {fixture.name}: correctly invalid")
            continue
            
        # Compute actual properties
        is_valid, n, k, d, ent, logical_ops = _compute_properties(fixture.seed)
        
        # Check each property
        assert is_valid == fixture.is_valid, f"{fixture.name}: validity mismatch"
        assert n == fixture.n, f"{fixture.name}: n mismatch {n} != {fixture.n}"
        assert k == fixture.k, f"{fixture.name}: k mismatch {k} != {fixture.k}"
        assert d == fixture.d, f"{fixture.name}: d mismatch {d} != {fixture.d}"
        assert ent == fixture.entanglement, f"{fixture.name}: entanglement mismatch {ent} != {fixture.entanglement}"
        
        # Store computed logical operators for later use
        fixture.logical_ops = logical_ops
        
        print(f"✓ {fixture.name}: n={n}, k={k}, d={d}, ent={ent}")
    
    print(f"All {len(CANONICAL_FIXTURES)} fixtures validated!")


def get_valid_fixtures():
    """Return only valid fixtures for testing."""
    return [f for f in CANONICAL_FIXTURES if f.is_valid]


def get_invalid_fixtures():
    """Return only invalid fixtures for testing."""
    return [f for f in CANONICAL_FIXTURES if not f.is_valid]


if __name__ == "__main__":
    # Validate fixtures when run directly
    validate_fixtures() 