import json
from stabilizer.generator import parse_seed
from stabilizer.qca import qca_evolution_step
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance
from stabilizer.utils import compute_entanglement

seeds = ["XZIY", "ZZII", "XZZX"]
MAX_STEPS = 2

def compute(seed: str, max_steps: int = MAX_STEPS):
    pauli = parse_seed(seed)
    L = len(pauli)
    data = []
    current = pauli
    for step in range(max_steps):
        stab_ops = build_lattice(current)
        tab = build_tableau(stab_ops)
        rank = compute_rank(tab)
        k = L - rank
        if k > 10:
            break
        if k == 0:
            dist, log_ops = 0, []
        else:
            dist, log_ops = find_distance(tab, return_logical_ops=True)
        ent = compute_entanglement(tab, log_ops if log_ops else None)
        data.append({'seed': seed, 'step': step, 'k': k, 'distance': dist, 'entanglement': ent})
        current = qca_evolution_step(current)
    return data

def compute_all(seeds: list[str], max_steps: int = MAX_STEPS):
    return {s: compute(s, max_steps) for s in seeds}


def main() -> None:
    out = compute_all(seeds, MAX_STEPS)
    with open("tests/data/golden.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
