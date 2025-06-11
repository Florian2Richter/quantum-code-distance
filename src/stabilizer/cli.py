import click
from .generator import parse_seed
from .lattice import build_lattice
from .tableau import build_tableau, compute_rank
from .distance import find_distance

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--seed",
    required=True,
    help="String of Pauli ops (e.g. 'XZIY') to seed the 1D ring."
)
def main(seed: str):
    """
    CLI entry point for building a stabilizer code on a 1D ring.
    """
    # 1. Parse and validate the seed string
    pauli = parse_seed(seed)

    # 2. Build the list of shifted generators around the ring
    L = len(pauli)
    stab_ops = build_lattice(pauli)

    # 3. Construct the binary symplectic tableau & compute its rank
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    n_logical = L - rank
    click.echo(f"Logical qubits: {n_logical} (L={L}, rank={rank})")

    # 4. (NP-complete) Search for the code distance
    distance = find_distance(tableau)
    click.echo(f"Code distance: {distance}")

if __name__ == "__main__":
    main() 