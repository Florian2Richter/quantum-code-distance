import json
from stabilizer.generator import parse_seed
from stabilizer.qca import qca_evolution_step
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance
from stabilizer.utils import compute_entanglement

MAX_K = 10


def compute_results(seed: str, steps: int):
    pauli = parse_seed(seed)
    L = len(pauli)
    results = []
    current = pauli
    for step in range(steps):
        stab_ops = build_lattice(current)
        tab = build_tableau(stab_ops)
        rank = compute_rank(tab)
        k = L - rank
        if k > MAX_K:
            break
        if k == 0:
            dist, log_ops = 0, []
        else:
            dist, log_ops = find_distance(tab, return_logical_ops=True)
        ent = compute_entanglement(tab, log_ops if log_ops else None)
        results.append({
            "seed": seed,
            "step": step,
            "k": k,
            "distance": dist,
            "entanglement": ent,
        })
        current = qca_evolution_step(current)
    return results


def test_golden_data():
    with open("tests/data/golden.json") as f:
        golden = json.load(f)

    for seed, expected in golden.items():
        steps = max(item["step"] for item in expected) + 1
        result = compute_results(seed, steps)
        assert result == expected

