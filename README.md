# My Stabilizer Code Project

Lean Python CLI to build stabilizer tableaus on a 1D periodic lattice and compute code parameters.

## Installation

```bash
pip install .
```

## Usage

```bash
stabcode --seed XZIY
```

This takes `--seed` (a string of Pauli operators), infers `L = len(seed)`, applies periodic boundary conditions, builds the full stabilizer tableau, computes rank & logical qubits, then (NP-complete) finds code distance.

The CLI also reports the bipartite entanglement across a half cut of the ring.

Example output:

```
Code parameters:
  Physical qubits (n): 4
  Stabilizer rank: 3
  Logical qubits (k): 1
  Code distance (d): 2
  Bipartite entanglement: 1
```
