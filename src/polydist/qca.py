"""
Bipartite entanglement computation using polynomial formalism.

This module implements entanglement calculation using Laurent polynomials
over F₂[x,x⁻¹] instead of GF(2) matrices.
"""


def compute_entanglement(tableau, logical_ops=None) -> int:
    """
    Compute bipartite entanglement across a half-ring cut.
    
    The input tableau specifies stabilizer generators for a pure state.
    Optional logical_ops may be supplied, but only those whose X part
    is all zeros (Z logicals) are appended as additional stabilizers fixing
    the logical qubits. The entanglement is then half the number of rows
    that have support on both sides of the cut when the chain is split in
    half.
    
    Args:
        tableau: Polynomial stabilizer tableau
        logical_ops: Optional list of logical operators
        
    Returns:
        Bipartite entanglement value
    """
    
    if not tableau or not tableau[0]:
        return 0
    
    twoL = len(tableau[0])
    L = twoL // 2
    
    # Start with original tableau
    extended_tableau = [row[:] for row in tableau]  # Deep copy
    
    # If logical operators are provided, only add Z-logical operators
    if logical_ops:
        z_logicals = []
        for op in logical_ops:
            # Check if X part (first L entries) is all zeros
            x_part = op[:L]
            is_z_logical = True
            for poly in x_part:
                if not poly.is_zero():
                    is_z_logical = False
                    break
            if is_z_logical:
                z_logicals.append(op)
        
        # Add Z-logical operators to tableau
        for z_logical in z_logicals:
            extended_tableau.append(z_logical[:])
    
    # Cut position
    cut = L // 2
    
    # Count rows with support on both sides of the cut
    count = 0
    for row in extended_tableau:
        x = row[:L]  # X part
        z = row[L:]  # Z part
        
        # Check if there's support on the left side (positions 0 to cut-1)
        left_support = False
        for i in range(cut):
            if not x[i].is_zero() or not z[i].is_zero():
                left_support = True
                break
        
        # Check if there's support on the right side (positions cut to L-1)
        right_support = False
        for i in range(cut, L):
            if not x[i].is_zero() or not z[i].is_zero():
                right_support = True
                break
        
        # Count if both sides have support
        if left_support and right_support:
            count += 1
    
    return count // 2 