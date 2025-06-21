import click
import time
from .generator import parse_seed
from .lattice import build_lattice
from .tableau import build_tableau, compute_rank
from .distance import find_distance, find_logical_operators, format_symplectic_vector
from .qca import evolve_pauli_string, qca_evolution_step

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
    if verbose:
        start_time = time.time()
    
    try:
        pauli = parse_seed(seed)
        if verbose:
            parse_time = time.time() - start_time
            click.echo(f"Parsed seed: {pauli}")
            click.echo(f"  Parsing time: {parse_time:.6f} seconds")
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        return
    
    # 2. Build the list of shifted generators around the ring
    L = len(pauli)
    
    if verbose:
        start_time = time.time()
    
    stab_ops = build_lattice(pauli)
    
    if verbose:
        lattice_time = time.time() - start_time
        click.echo(f"\nGenerated {len(stab_ops)} stabilizer operators on L={L} qubits:")
        for i, op in enumerate(stab_ops):
            click.echo(f"  S{i}: {''.join(op)}")
        click.echo(f"  Lattice building time: {lattice_time:.6f} seconds")
    
    # 3. Construct the binary symplectic tableau & compute its rank
    if verbose:
        start_time = time.time()
    
    tableau = build_tableau(stab_ops)
    
    if verbose:
        tableau_time = time.time() - start_time
        start_time = time.time()
    
    rank = compute_rank(tableau)
    n_logical = L - rank
    
    if verbose:
        rank_time = time.time() - start_time
        click.echo(f"  Tableau construction time: {tableau_time:.6f} seconds")
        click.echo(f"  Rank computation time: {rank_time:.6f} seconds")
    
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
    
    if verbose:
        start_time = time.time()
    
    distance = find_distance(tableau)
    
    if verbose:
        distance_time = time.time() - start_time
        click.echo(f"  Code distance (d): {distance}")
        click.echo(f"  Distance computation time: {distance_time:.6f} seconds")
    else:
        click.echo(f"  Code distance (d): {distance}")
    
    # Summary
    click.echo(f"\nQuantum Error Correcting Code: [[{L}, {n_logical}, {distance}]]")
    
    if verbose and n_logical > 0:
        click.echo(f"\nFinding logical operators...")
        
        start_time = time.time()
        logical_ops = find_logical_operators(tableau)
        logical_ops_time = time.time() - start_time
        
        if logical_ops:
            click.echo(f"Found {len(logical_ops)} logical operators:")
            for i, vec in enumerate(logical_ops):  
                pauli_str = format_symplectic_vector(vec)
                click.echo(f"  Logical op {i+1}: {pauli_str}")
            
        else:
            click.echo("  No logical operators found in low-weight search.")
        
        click.echo(f"  Logical operators search time: {logical_ops_time:.6f} seconds")
    
    # QCA Time Evolution (if requested)
    if time_steps > 0:
        click.echo(f"\n" + "="*60)
        click.echo(f"QCA TIME EVOLUTION ({time_steps} steps)")
        click.echo(f"="*60)
        
        # Show initial state
        click.echo(f"time step 0: {''.join(pauli)} → distance = {distance}")
        
        current_pauli = pauli.copy()
        
        for step in range(1, time_steps + 1):
            if verbose:
                qca_start_time = time.time()
            
            # Apply one QCA evolution step
            current_pauli = qca_evolution_step(current_pauli)
            
            if verbose:
                qca_step_time = time.time() - qca_start_time
                click.echo(f"\n  QCA step {step} evolution time: {qca_step_time:.6f} seconds")
            
            # Build stabilizer code for evolved string
            if verbose:
                start_time = time.time()
            
            evolved_stab_ops = build_lattice(current_pauli)
            evolved_tableau = build_tableau(evolved_stab_ops)
            evolved_rank = compute_rank(evolved_tableau)
            evolved_n_logical = L - evolved_rank
            
            if verbose:
                evolved_setup_time = time.time() - start_time
                click.echo(f"  Evolved code setup time: {evolved_setup_time:.6f} seconds")
            
            if evolved_n_logical == 0:
                evolved_distance = "N/A (no logical qubits)"
            else:
                if verbose:
                    start_time = time.time()
                
                evolved_distance = find_distance(evolved_tableau)
                
                if verbose:
                    evolved_distance_time = time.time() - start_time
                    click.echo(f"  Evolved distance computation time: {evolved_distance_time:.6f} seconds")
            
            # Display result for this time step
            click.echo(f"time step {step}: {''.join(current_pauli)} → distance = {evolved_distance}")
            
            if verbose and evolved_n_logical > 0:
                click.echo(f"  [[{L}, {evolved_n_logical}, {evolved_distance}]]")

    if verbose:
        total_time = time.time() - total_start_time
        click.echo(f"\nTotal execution time: {total_time:.6f} seconds")

if __name__ == "__main__":
    main() 