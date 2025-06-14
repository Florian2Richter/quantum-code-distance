import click
from .generator import parse_seed
from .lattice import build_lattice
from .tableau import build_tableau, compute_rank
from .distance import find_distance, find_logical_operators, format_symplectic_vector

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--seed",
    required=True,
    help="String of Pauli ops (e.g. 'XZIY') to seed the 1D ring."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Show detailed information about the computation."
)
@click.option(
    "--show-tableau",
    is_flag=True,
    help="Display the binary symplectic tableau."
)

@click.option(
    "--no-progress",
    is_flag=True,
    help="Disable progress bars."
)
def main(seed: str, verbose: bool, show_tableau: bool, no_progress: bool):
    """
    CLI entry point for building a stabilizer code on a 1D ring.
    
    Takes a seed string of Pauli operators and computes the quantum error 
    correcting code parameters on a 1D lattice with periodic boundary conditions.
    """
    click.echo(f"Building stabilizer code from seed: {seed}")
    
    # 1. Parse and validate the seed string
    try:
        pauli = parse_seed(seed)
        if verbose:
            click.echo(f"Parsed seed: {pauli}")
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        return
    
    # 2. Build the list of shifted generators around the ring
    L = len(pauli)
    stab_ops = build_lattice(pauli)
    
    if verbose:
        click.echo(f"\nGenerated {len(stab_ops)} stabilizer operators on L={L} qubits:")
        for i, op in enumerate(stab_ops):
            click.echo(f"  S{i}: {''.join(op)}")
    
    # 3. Construct the binary symplectic tableau & compute its rank
    tableau = build_tableau(stab_ops)
    rank = compute_rank(tableau)
    n_logical = L - rank
    
    if show_tableau:
        click.echo(f"\nBinary symplectic tableau ({tableau.shape[0]}×{tableau.shape[1]}):")
        click.echo(tableau)
    
    click.echo(f"\nCode parameters:")
    click.echo(f"  Physical qubits (n): {L}")
    click.echo(f"  Stabilizer rank: {rank}")
    click.echo(f"  Logical qubits (k): {n_logical}")
    
    if n_logical == 0:
        click.echo(f"  Code distance (d): N/A (no logical qubits)")
        return
    
    # 4. (NP-complete) Search for the code distance
    click.echo(f"\nSearching for code distance...")
    
    if no_progress:
        # Disable progress bars by using a simple version
        import os
        os.environ['TQDM_DISABLE'] = '1'
    
    distance = find_distance(tableau)
    click.echo(f"  Code distance (d): {distance}")
    
    # Summary
    click.echo(f"\nQuantum Error Correcting Code: [[{L}, {n_logical}, {distance}]]")
    
    if verbose and n_logical > 0:
        click.echo(f"\nFinding logical operators (up to weight 3)...")
        logical_ops = find_logical_operators(tableau, max_weight=3)
        if logical_ops:
            click.echo(f"Found {len(logical_ops)} logical operators:")
            for i, vec in enumerate(logical_ops[:5]):  # Show first 5
                pauli_str = format_symplectic_vector(vec)
                click.echo(f"  Logical op {i+1}: {pauli_str}")
        else:
            click.echo("  No logical operators found in low-weight search.")

if __name__ == "__main__":
    main() 