import json
from stabilizer.generator import parse_seed
from stabilizer.qca import qca_evolution_step
from stabilizer.lattice import build_lattice
from stabilizer.tableau import build_tableau, compute_rank
from stabilizer.distance import find_distance
from stabilizer.utils import compute_entanglement

seeds = ['XZIY', 'ZZII', 'XZZX']
max_steps = 2

def compute(seed):
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

out = {s: compute(s) for s in seeds}
print(json.dumps(out, indent=2))
