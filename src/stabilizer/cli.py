import click
import time
from contextlib import contextmanager
from typing import Tuple, List
from .generator import parse_seed
from .lattice import build_lattice
from .tableau import build_tableau, compute_rank
from .distance import find_distance, find_logical_operators, format_symplectic_vector
from .qca import evolve_pauli_string, qca_evolution_step

@contextmanager
def timer(description: str, verbose: bool = False):
    """
    Context manager for timing operations with optional verbose output.
    
    Args:
        description: Description of the operation being timed
        verbose: Whether to display timing information
    """
    if verbose:
        start_time = time.time()
        yield
        elapsed = time.time() - start_time
        click.echo(f"  {description} time: {elapsed:.6f} seconds")
    else:
        yield

def compute_code_distance(pauli_string: List[str], L: int, verbose: bool) -> Tuple[int, int, str]:
    """
    Compute stabilizer code parameters for a given Pauli string.
    
    Args:
        pauli_string: List of Pauli operators
        L: Number of physical qubits
        verbose: Whether to show timing information
        
    Returns:
        Tuple of (n_logical, rank, distance) where distance is str or int
    """
    with timer("Building lattice", verbose):
        stab_ops = build_lattice(pauli_string)
    
    with timer("Building tableau", verbose):
        tableau = build_tableau(stab_ops)
    
    with timer("Computing rank", verbose):
        rank = compute_rank(tableau)
        n_logical = L - rank
    
    if n_logical == 0:
        distance = "N/A (no logical qubits)"
    else:
        with timer("Computing distance", verbose):
            distance = find_distance(tableau)
    
    return n_logical, rank, distance

def run_qca_evolution(initial_pauli: List[str], L: int, time_steps: int, verbose: bool):
    """
    Run QCA time evolution for specified number of steps.
    
    Args:
        initial_pauli: Initial Pauli string
        L: Number of physical qubits  
        time_steps: Number of evolution steps to perform
        verbose: Whether to show detailed timing
    """
    click.echo(f"\n" + "="*60)
    click.echo(f"QCA TIME EVOLUTION ({time_steps} steps)")
    click.echo(f"="*60)
    
    current_pauli = initial_pauli.copy()
    
    # Unified loop: step 0 = initial, steps 1+ = evolved
    for step in range(time_steps + 1):
        if step > 0:
            # Evolve from previous step
            with timer(f"QCA step {step} evolution", verbose):
                current_pauli = qca_evolution_step(current_pauli)
        
        # Compute distance for current state (same logic for all steps)
        if verbose and step > 0:
            click.echo(f"\nTime step {step}:")
        
        n_logical, rank, distance = compute_code_distance(current_pauli, L, verbose)
        
        # Display result in unified format
        click.echo(f"time step {step}: {''.join(current_pauli)} → distance = {distance}")
        
        if verbose and n_logical > 0:
            click.echo(f"  [[{L}, {n_logical}, {distance}]]")

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
@click.option(
    "--time-steps",
    type=int,
    default=0,
    help="Number of QCA time evolution steps to perform (default: 0, no evolution)."
)
def main(seed: str, verbose: bool, show_tableau: bool, no_progress: bool, time_steps: int):
    """
    CLI entry point for building a stabilizer code on a 1D ring.
    
    Takes a seed string of Pauli operators and computes the quantum error 
    correcting code parameters on a 1D lattice with periodic boundary conditions.
    """
    click.echo(f"Building stabilizer code from seed: {seed}")
    
    if verbose:
        total_start_time = time.time()
    
    # 1. Parse and validate the seed string
    try:
        with timer("Parsing seed", verbose):
            pauli = parse_seed(seed)
        if verbose:
            click.echo(f"Parsed seed: {pauli}")
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        return
    
    L = len(pauli)
    
    # 2. Build initial stabilizer code (only for display and logical operators)
    if verbose:
        click.echo(f"\nGenerating stabilizer operators on L={L} qubits...")
    
    with timer("Building initial lattice", verbose):
        stab_ops = build_lattice(pauli)
    
    if verbose:
        click.echo(f"Generated {len(stab_ops)} stabilizer operators:")
        for i, op in enumerate(stab_ops):
            click.echo(f"  S{i}: {''.join(op)}")
    
    with timer("Building initial tableau", verbose):
        tableau = build_tableau(stab_ops)
    
    with timer("Computing initial rank", verbose):
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
        if time_steps == 0:
            return
    else:
        # Compute distance for initial state
        if no_progress:
            import os
            os.environ['TQDM_DISABLE'] = '1'
        
        click.echo(f"\nSearching for code distance...")
        with timer("Computing initial distance", verbose):
            distance = find_distance(tableau)
        
        click.echo(f"  Code distance (d): {distance}")
        click.echo(f"\nQuantum Error Correcting Code: [[{L}, {n_logical}, {distance}]]")
        
        # Find logical operators if verbose
        if verbose and n_logical > 0:
            click.echo(f"\nFinding logical operators...")
            
            with timer("Finding logical operators", verbose):
                logical_ops = find_logical_operators(tableau)
            
            if logical_ops:
                click.echo(f"Found {len(logical_ops)} logical operators:")
                for i, vec in enumerate(logical_ops):  
                    pauli_str = format_symplectic_vector(vec)
                    click.echo(f"  Logical op {i+1}: {pauli_str}")
            else:
                click.echo("  No logical operators found.")
    
    # QCA Time Evolution (if requested)
    if time_steps > 0:
        run_qca_evolution(pauli, L, time_steps, verbose)

    if verbose:
        total_time = time.time() - total_start_time
        click.echo(f"\nTotal execution time: {total_time:.6f} seconds")

if __name__ == "__main__":
    main() 