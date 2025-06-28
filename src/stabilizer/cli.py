import click
import logging
import numpy as np
import time
from contextlib import contextmanager
from typing import Tuple, List
from .generator import parse_seed
from .lattice import build_lattice
from .tableau import build_tableau, compute_rank
from .distance import find_distance, find_logical_operators, format_symplectic_vector
from .qca import qca_evolution_step
from .utils import seed_is_valid


def configure_logging(verbose: bool):
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-5s %(message)s",
        datefmt="%H:%M:%S"
    )


class Display:
    """Handles all display logic using Python's standard logging module."""
    
    def __init__(self):
        self.logger = logging.getLogger("qca")
    
    def info(self, message: str, *args):
        """Always display important information."""
        self.logger.info(message, *args)
    
    def debug(self, message: str, *args):
        """Display debug information only if debug level is enabled."""
        self.logger.debug(message, *args)
    
    def section(self, title: str):
        """Display a section header."""
        self.logger.info("\n%s", title)
    
    def separator(self, text: str, char: str = "=", width: int = 60):
        """Display a separator line with text."""
        line = char * width
        self.logger.info("\n%s\n%s\n%s", line, text, line)
    
    @contextmanager
    def timer(self, description: str):
        """Context manager for timing operations with optional debug output."""
        start_time = time.time()
        yield
        elapsed = time.time() - start_time
        self.logger.debug("  %s time: %.6f seconds", description, elapsed)
    
    def show_operators(self, stab_ops: List[List[str]], L: int):
        """Display stabilizer operators if in debug mode."""
        self.logger.debug("Generated %d stabilizer operators:", len(stab_ops))
        for i, op in enumerate(stab_ops):
            self.logger.debug("  S%d: %s", i, ''.join(op))
    
    def show_logical_operators(self, logical_ops: List, display_all: bool = False):
        """Display logical operators if found."""
        if logical_ops:
            self.logger.info("Found %d logical operators:", len(logical_ops))
            ops_to_show = logical_ops if display_all else logical_ops[:10]
            for i, vec in enumerate(ops_to_show):  
                pauli_str = format_symplectic_vector(vec)
                self.logger.info("  Logical op %d: %s", i+1, pauli_str)
            if not display_all and len(logical_ops) > 10:
                self.logger.info("  ... and %d more operators", len(logical_ops) - 10)
        else:
            self.logger.info("  No logical operators found.")


def analyze_stabilizer_code(pauli_string: List[str], L: int, step: int, 
                          display: Display, verbose: bool) -> Tuple[int, int, str]:
    """
    Complete analysis of a stabilizer code for a given Pauli string.
    
    Args:
        pauli_string: List of Pauli operators
        L: Number of physical qubits
        step: Current time step (for display purposes)
        display: Display handler for output
        verbose: Whether to show debug info (includes tableau and disables progress)
        
    Returns:
        Tuple of (n_logical, rank, distance) where distance is str or int
    """
    step_label = "Time step %d" % step if step > 0 else "Initial state"
    if step > 0:
        display.separator("%s: %s" % (step_label, ''.join(pauli_string)))
    
    # Build stabilizer operators
    display.debug("Generating stabilizer operators on L=%d qubits...", L)
    
    with display.timer("Building lattice"):
        stab_ops = build_lattice(pauli_string)
    
    display.show_operators(stab_ops, L)
    
    # Build tableau
    with display.timer("Building tableau"):
        tableau = build_tableau(stab_ops)
    
    # Compute rank
    with display.timer("Computing rank"):
        rank = compute_rank(tableau)
        n_logical = L - rank
    
    if verbose:
        display.section("Binary symplectic tableau (%dx%d):" % (tableau.shape[0], tableau.shape[1]))
        # Temporarily set numpy to show full array without truncation
        with np.printoptions(threshold=np.inf, linewidth=np.inf):
            display.debug("Tableau:\n%s", tableau)
    
    # Display code parameters
    display.section("Code parameters:")
    display.info("  Physical qubits (n): %d", L)
    display.info("  Stabilizer rank: %d", rank)
    display.info("  Logical qubits (k): %d", n_logical)
    
    if n_logical == 0:
        distance = "N/A (no logical qubits)"
        display.info("  Code distance (d): %s", distance)
    else:
        # Compute distance
        if verbose:
            import os
            os.environ['TQDM_DISABLE'] = '1'
        
        display.section("Searching for code distance...")
        with display.timer("Computing distance"):
            distance = find_distance(tableau)
        
        display.info("  Code distance (d): %s", distance)
        display.section("Quantum Error Correcting Code: [[%d, %d, %s]]" % (L, n_logical, distance))
        
        # Find logical operators if in debug mode
        if display.logger.isEnabledFor(logging.DEBUG) and n_logical > 0:
            display.section("Finding logical operators...")
            
            with display.timer("Finding logical operators"):
                logical_ops = find_logical_operators(tableau)
            
            display.show_logical_operators(logical_ops)
    
    return n_logical, rank, distance


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--seed",
    required=True,
    help="String of Pauli ops (e.g. 'XZIY') to seed the 1D ring."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Show detailed debug information, including timing, tableau display, and disable progress bars."
)
@click.option(
    "--time-steps",
    type=int,
    default=0,
    help="Number of QCA time evolution steps to perform (default: 0, no evolution)."
)
def main(seed: str, verbose: bool, time_steps: int):
    """
    CLI entry point for building a stabilizer code on a 1D ring.
    
    Takes a seed string of Pauli operators and computes the quantum error 
    correcting code parameters on a 1D lattice with periodic boundary conditions.
    """
    # Configure logging based on verbosity
    configure_logging(verbose)
    display = Display()
    
    display.info("Building stabilizer code from seed: %s", seed)
    
    total_start_time = time.time() if verbose else None
    
    # 1. Parse and validate the seed string
    try:
        with display.timer("Parsing seed"):
            pauli = parse_seed(seed)
        display.debug("Parsed seed: %s", pauli)
        
        # Check if seed generates commuting stabilizers
        with display.timer("Validating seed commutativity"):
            is_valid = seed_is_valid(seed)
        display.debug("Seed commutativity check: %s", "PASSED" if is_valid else "FAILED")
        if not is_valid:
            display.info("Error: Cyclic translations of this seed do not all commute")
            return
            
    except ValueError as e:
        display.info("Error: %s", e)
        return
    
    L = len(pauli)
    current_pauli = pauli.copy()
    
    # Time evolution loop: step 0 = initial, steps 1+ = evolved
    for step in range(time_steps + 1):
        if step > 0:
            # Evolve from previous step
            with display.timer("QCA step %d evolution" % step):
                current_pauli = qca_evolution_step(current_pauli)
        
        # Run full analysis for current state  
        n_logical, rank, distance = analyze_stabilizer_code(
            current_pauli, L, step, display, verbose
        )
        
        # Early exit if no logical qubits and no time evolution requested
        if n_logical == 0 and time_steps == 0:
            break

    if verbose and total_start_time:
        total_time = time.time() - total_start_time
        display.section("Total execution time: %.6f seconds" % total_time)


if __name__ == "__main__":
    main() 