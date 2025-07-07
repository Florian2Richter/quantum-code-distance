#!/usr/bin/env python3
"""
Simple runner for profiling comparison with command line options.
"""

import argparse
import sys

# Add src to path
sys.path.insert(0, 'src')

from profile_comparison import main as profile_main, create_initial_pauli_string, ProfilerResults, profile_tableau_approach, profile_polynomial_approach


def main():
    parser = argparse.ArgumentParser(description='Profile polynomial vs tableau stabilizer code approaches')
    parser.add_argument('--size', '-L', type=int, default=200, help='Total number of qubits (default: 200)')
    parser.add_argument('--x-position', type=int, default=100, help='Position of X gate (0-indexed, default: 100)')
    parser.add_argument('--time-steps', '-t', type=int, default=10, help='Number of time evolution steps (default: 10)')
    parser.add_argument('--approach', choices=['tableau', 'polynomial', 'both'], default='both', 
                       help='Which approach to profile (default: both)')
    parser.add_argument('--csv', type=str, help='Save results to CSV file')
    
    args = parser.parse_args()
    
    print("="*80)
    print("STABILIZER CODE PROFILING COMPARISON")
    print("="*80)
    print(f"System size: {args.size} qubits")
    print(f"X position: {args.x_position}")
    print(f"Time steps: {args.time_steps}")
    print(f"Approach: {args.approach}")
    print("="*80)
    
    # Create initial Pauli string
    initial_pauli = create_initial_pauli_string(args.size, args.x_position)
    print(f"Initial Pauli string: {''.join(initial_pauli[:10])}...{''.join(initial_pauli[-10:])}")
    
    # Create profiler
    profiler = ProfilerResults()
    
    try:
        if args.approach in ['tableau', 'both']:
            profile_tableau_approach(initial_pauli, args.time_steps, profiler)
        
        if args.approach in ['polynomial', 'both']:
            profile_polynomial_approach(initial_pauli, args.time_steps, profiler)
            
    except Exception as e:
        print(f"Error during profiling: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Print results
    profiler.print_summary()
    
    # Save to CSV if requested
    if args.csv:
        save_to_csv(profiler, args.csv)
        print(f"\nResults saved to {args.csv}")
    
    return 0


def save_to_csv(profiler: ProfilerResults, filename: str):
    """Save profiling results to CSV file."""
    import csv
    
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Step', 'Operation', 'Tableau_Time', 'Polynomial_Time', 'Speedup'])
        
        # Get all operations and steps
        all_ops = set()
        all_steps = set()
        for approach in profiler.results:
            for step in profiler.results[approach]:
                all_steps.add(step)
                for op in profiler.results[approach][step]:
                    all_ops.add(op)
        
        for step in sorted(all_steps):
            for op in sorted(all_ops):
                tab_time = profiler.results['tableau'].get(step, {}).get(op, 0)
                poly_time = profiler.results['polynomial'].get(step, {}).get(op, 0)
                
                if tab_time > 0 or poly_time > 0:
                    speedup = tab_time / poly_time if poly_time > 0 else 'N/A'
                    writer.writerow([step, op, tab_time, poly_time, speedup])


if __name__ == "__main__":
    sys.exit(main()) 