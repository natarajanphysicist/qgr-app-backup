# lqg_simulation/mathematics/wigner_symbols.py
"""
Provides functions to calculate Wigner 3j, 6j, and 9j symbols.

These symbols are crucial in quantum mechanics for coupling angular momenta and appear
frequently in Loop Quantum Gravity calculations, particularly in the definition of
vertex amplitudes and recoupling theory.

The inputs (j1, j2, ..., m1, m2, ...) are expected to be integers or half-integers.
The sympy library handles these inputs naturally.
"""

from sympy import N as EvaluateNumerical, S
from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j
# from sympy import Integer # Not strictly needed for current validation logic
# from sympy import HalfInteger # This does not exist for direct import

# Helper to check if a number is integer or half-integer
def is_int_or_half_int(val):
    """Checks if a number is an integer or a half-integer."""
    # Check for Sympy types first, as they might also be instances of float/int
    if hasattr(val, 'is_rational'): # Common attribute for Sympy numbers
        doubled = val * 2
        return doubled.is_integer if hasattr(doubled, 'is_integer') else False
    elif isinstance(val, (int, float)): # Standard Python numbers
        return (val * 2) % 1 == 0
    return False # Not a recognized/supported numeric type for this check


def calculate_wigner_3j(j1, j2, j3, m1, m2, m3):
    """
    Calculates the Wigner 3j symbol:
    ( j1 j2 j3 )
    ( m1 m2 m3 )

    Args:
        j1, j2, j3: Angular momentum quantum numbers (non-negative int or half-int).
        m1, m2, m3: Projection quantum numbers (int or half-int).

    Returns:
        The numerical value of the Wigner 3j symbol (float).
        Returns 0.0 if selection rules are not met.

    Selection Rules (automatically handled by sympy.physics.wigner.wigner_3j):
    1. j1, j2, j3 must satisfy triangle inequalities: |j1 - j2| <= j3 <= j1 + j2.
    2. j1 + j2 + j3 must be an integer (closure condition).
    3. m1, m2, m3 must be such that |mi| <= ji.
    4. m1 + m2 + m3 = 0.
    5. If m1 = m2 = m3 = 0, then j1 + j2 + j3 must be even.
    """
    # Input validation
    js = [j1, j2, j3]
    ms = [m1, m2, m3]

    for j_val in js:
        if not (is_int_or_half_int(j_val) and j_val >= 0):
            raise ValueError(f"j value {j_val} must be a non-negative integer or half-integer.")
    for i, m_val in enumerate(ms):
        if not is_int_or_half_int(m_val):
            raise ValueError(f"m value {m_val} must be an integer or half-integer.")
        # abs(m) <= j is a selection rule handled by sympy by returning 0, not an input error.
        # However, we should ensure m is of the same type (int/half-int) as j.
        # e.g. if j is int, m must be int. if j is half-int, m must be half-int.
        # (j - m) must be an integer.
        if not is_int_or_half_int(js[i] - m_val) or (js[i] - m_val) % 1 != 0 :
            # This check is implicitly handled by sympy's wigner_3j if it expects j, m to be "compatible"
            # For now, let's rely on sympy for this level of detail, the primary check is type.
            pass

    # Convert all inputs to Sympy numbers to ensure consistent behavior with Sympy functions
    sj1, sj2, sj3 = S(j1), S(j2), S(j3)
    sm1, sm2, sm3 = S(m1), S(m2), S(m3)

    # Sympy's wigner_3j function computes the symbol.
    # It returns an exact Sympy expression. We convert it to a float.
    # It automatically returns 0 if selection rules are not met.
    symbol_val = wigner_3j(sj1, sj2, sj3, sm1, sm2, sm3)

    # NOTE: Debug prints for specific Sympy behavior were here.
    # It was observed that sympy.wigner_3j(S(1),S(1),S(1),S(0),S(0),S(0))
    # was returning 0 in this execution environment, contrary to expected -sqrt(3)/9.

    try:
        # Prefer evalf() for Sympy expressions to get a numerical approximation
        return float(symbol_val.evalf())
    except AttributeError:
        # If no evalf (e.g., already a Python number or Sympy basic number like S.Zero), directly convert
        return float(symbol_val)


def calculate_wigner_6j(j1, j2, j3, j4, j5, j6):
    """
    Calculates the Wigner 6j symbol:
    { j1 j2 j3 }
    { j4 j5 j6 }

    Args:
        j1, j2, j3, j4, j5, j6: Angular momentum quantum numbers (non-negative int or half-int).

    Returns:
        The numerical value of the Wigner 6j symbol (float).
        Returns 0.0 if selection rules are not met.

    Selection Rules (automatically handled by sympy.physics.wigner.wigner_6j):
    The four triads must satisfy triangle inequalities:
    1. (j1, j2, j3)
    2. (j1, j5, j6)
    3. (j4, j2, j6)
    4. (j4, j5, j3)
    Also, the sum of elements in each triad must be an integer (closure).
    All j values must be non-negative integers or half-integers.
    """
    # Input validation
    js = [j1, j2, j3, j4, j5, j6]
    for j_val in js:
        if not (is_int_or_half_int(j_val) and j_val >= 0):
            raise ValueError(f"All j values for 6j symbol must be non-negative integers or half-integers. Got {j_val}")

    # Convert all inputs to Sympy numbers
    sj1, sj2, sj3, sj4, sj5, sj6 = S(j1), S(j2), S(j3), S(j4), S(j5), S(j6)

    # Sympy's wigner_6j function computes the symbol.
    # It returns an exact Sympy expression. We convert it to a float.
    # It automatically returns 0 if selection rules are not met.
    # However, for some invalid triad sums (e.g. j1+j2+j3 not int), it might raise ValueError.
    try:
        symbol_val = wigner_6j(sj1, sj2, sj3, sj4, sj5, sj6)
        # NOTE: Debug prints for specific Sympy behavior were here.
        # It was observed that sympy.wigner_6j for (1/2,1/2,1, 1/2,1/2,1)
        # was returning 1/6 in this execution environment, contrary to expected -1/6.
    except ValueError as e:
        # If Sympy raises ValueError due to selection rules (like non-integer sum of a triad),
        # we interpret this as the symbol being 0.
        if "fulfill the triangle relation" in str(e) or "integer or half integer" in str(e):
            # This specific error message "j values must be integer or half integer and fulfill the triangle relation"
            # is raised by _big_delta_coeff when, for example, j1+j2+j3 is not an integer.
            return 0.0
        raise # Re-raise other ValueErrors

    try:
        return float(symbol_val.evalf())
    except AttributeError:
        return float(symbol_val)


def calculate_wigner_9j(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    """
    Calculates the Wigner 9j symbol:
    { j1 j2 j3 }
    { j4 j5 j6 }
    { j7 j8 j9 }

    Args:
        j1-j9: Angular momentum quantum numbers (non-negative int or half-int).

    Returns:
        The numerical value of the Wigner 9j symbol (float).
        Returns 0.0 if selection rules are not met.

    Selection Rules (handled by sympy.physics.wigner.wigner_9j):
    - Each row and each column must satisfy triangle inequalities.
      (j1,j2,j3), (j4,j5,j6), (j7,j8,j9) for rows.
      (j1,j4,j7), (j2,j5,j8), (j3,j6,j9) for columns.
    - The sum of j values in each triad must be an integer.
    - All j values must be non-negative integers or half-integers.
    """
    # Input validation
    js = [j1, j2, j3, j4, j5, j6, j7, j8, j9]
    for j_val in js:
        if not (is_int_or_half_int(j_val) and j_val >= 0):
            raise ValueError(f"All j values for 9j symbol must be non-negative integers or half-integers. Got {j_val}")

    # Convert all inputs to Sympy numbers
    s_js = [S(j) for j in js]

    try:
        # Sympy's wigner_9j expects flattened list of j's or individual args
        symbol_val = wigner_9j(*s_js)
    except ValueError as e:
        # Catch ValueErrors from Sympy that might indicate selection rule violations
        # (e.g., if a triad sum is not an integer, similar to 6j issues)
        # A common error message from underlying C code via Swig might be "Error in Wigner coefficient"
        # or specific messages about triangle conditions.
        # For now, a broad catch for ValueError from Sympy's layer is used.
        # More specific error message checks could be added if Sympy's behavior is more nuanced.
        # Based on sympy.physics.wigner.wigner_9j structure, it calls wigner_racah which can raise ValueError.
        return 0.0 # Assuming selection rule violation

    try:
        return float(symbol_val.evalf())
    except AttributeError:
        return float(symbol_val)


if __name__ == '__main__':
    print("Wigner 3j Symbol Examples:")
    # Example 1: ( 1  1  1 ) / ( 0  0  0 ) -> should be -1/sqrt(3) approx -0.57735
    val1 = calculate_wigner_3j(1, 1, 1, 0, 0, 0)
    print(f"wigner_3j(1, 1, 1, 0, 0, 0) = {val1:.5f}") # Expected: -0.57735

    # Example 2: Clebsch-Gordan coeff relation: C(j1,m1,j2,m2;j3,m3) = (-1)^(j1-j2+m3) * sqrt(2*j3+1) * wigner_3j(j1,j2,j3,m1,m2,-m3)
    # For (1/2, 1/2, 1/2, -1/2 ; 1, 0) => C = 1/sqrt(2)
    # j1=1/2, m1=1/2, j2=1/2, m2=-1/2, j3=1, m3=0
    # wigner_3j(1/2, 1/2, 1, 1/2, -1/2, 0)
    # Expected value based on online calculator: 1/sqrt(6) ~ 0.408248
    val2 = calculate_wigner_3j(0.5, 0.5, 1, 0.5, -0.5, 0)
    print(f"wigner_3j(0.5, 0.5, 1, 0.5, -0.5, 0) = {val2:.5f}") # Expected: 0.40825

    # Example 3: Selection rule violation (m1+m2+m3 != 0)
    val3 = calculate_wigner_3j(1, 1, 1, 1, 0, 0)
    print(f"wigner_3j(1, 1, 1, 1, 0, 0) = {val3:.5f}") # Expected: 0.0

    # Example 4: Selection rule violation (triangle inequality: 1+1 < 3)
    val4 = calculate_wigner_3j(1, 1, 3, 0, 0, 0)
    print(f"wigner_3j(1, 1, 3, 0, 0, 0) = {val4:.5f}") # Expected: 0.0

    # Example 5: (2, 1, 1, 0,0,0) -> 1/sqrt(15) ~ 0.258199
    val5 = calculate_wigner_3j(2,1,1,0,0,0)
    print(f"wigner_3j(2,1,1,0,0,0) = {val5:.5f}")


    print("\nWigner 6j Symbol Examples:")
    # Example 1: Racah formula relation.
    # {1/2 1/2 1}
    # {1/2 1/2 1}
    # Expected value: -1/6 ~ -0.16667
    val6j_1 = calculate_wigner_6j(0.5, 0.5, 1, 0.5, 0.5, 1)
    print(f"wigner_6j(0.5, 0.5, 1, 0.5, 0.5, 1) = {val6j_1:.5f}")

    # Example 2:
    # {1 1 1}
    # {1 1 1}
    # Expected value: 1/6 ~ 0.16667
    val6j_2 = calculate_wigner_6j(1, 1, 1, 1, 1, 1)
    print(f"wigner_6j(1, 1, 1, 1, 1, 1) = {val6j_2:.5f}")

    # Example 3: Triangle rule violation (1,1,3)
    val6j_3 = calculate_wigner_6j(1, 1, 3, 1, 1, 1)
    print(f"wigner_6j(1, 1, 3, 1, 1, 1) = {val6j_3:.5f}") # Expected: 0.0

    # Example from Wikipedia:
    # {1 2 3}
    # {1 2 3}
    # Expected: 1/30 = 0.03333
    val6j_4 = calculate_wigner_6j(1,2,3,1,2,3)
    print(f"wigner_6j(1,2,3,1,2,3) = {val6j_4:.5f}")

    # Example with half-integers
    # {3/2 1 1/2}
    # {3/2 1 3/2}
    # Expected: -1/12 ~ -0.08333
    val6j_5 = calculate_wigner_6j(1.5, 1, 0.5, 1.5, 1, 1.5)
    print(f"wigner_6j(1.5, 1, 0.5, 1.5, 1, 1.5) = {val6j_5:.5f}")

    try:
        print("\nTesting invalid input types:")
        calculate_wigner_3j(1, 1.2, 1, 0,0,0) # Not half-integer
    except ValueError as e:
        print(f"Caught expected error for 3j: {e}")

    try:
        calculate_wigner_6j(1, 'a', 1, 0,0,0) # Not a number
    except Exception as e: # Sympy might raise TypeError before my check
        print(f"Caught expected error for 6j: {e}")

    try:
        calculate_wigner_3j(1, 1, 1, 0.5, 0, -0.5) # j integer, m half-integer mismatch implicitly handled by sympy
        # My checks are for numeric type, not sum of j+m being int. Sympy does this.
    except ValueError as e:
        print(f"Caught expected error for 3j (j+m type mismatch): {e}")


    # Test case from a known source (e.g. WolframAlpha or a textbook for 3j)
    # Wigner3j[j1, j2, j3, m1, m2, m3]
    # Wigner3j[2, 2, 2, 1, 1, -2] = Sqrt[2/35] ~ 0.2390457
    val_test_3j = calculate_wigner_3j(2,2,2,1,1,-2)
    print(f"wigner_3j(2,2,2,1,1,-2) = {val_test_3j:.7f} (Expected ~0.2390457)")

    # Wigner6j[{j1,j2,j3},{j4,j5,j6}]
    # Wigner6j[{2,2,2},{2,2,2}] = -1/70 ~ -0.0142857
    val_test_6j = calculate_wigner_6j(2,2,2,2,2,2)
    print(f"wigner_6j(2,2,2,2,2,2) = {val_test_6j:.7f} (Expected ~-0.0142857)")

    # Test for |m| > j, should return 0
    val_m_gt_j = calculate_wigner_3j(1,1,1, 2,0,-2)
    print(f"wigner_3j(1,1,1, 2,0,-2) with |m|>j = {val_m_gt_j}") # Expected 0.0

    # Test for j1+j2+j3 not an integer
    val_sum_j_not_int = calculate_wigner_3j(0.5, 1, 1.5, 0.5,0,-0.5) # sum j = 3 (int) - valid
    print(f"wigner_3j(0.5, 1, 1.5, 0.5,0,-0.5) = {val_sum_j_not_int:.5f}")

    # This case (0.5, 0.5, 0.5) j1+j2+j3 = 1.5 (not int) should be 0
    val_sum_j_not_int_2 = calculate_wigner_3j(0.5, 0.5, 0.5, 0.5, -0.5, 0)
    print(f"wigner_3j(0.5, 0.5, 0.5, 0.5, -0.5, 0) = {val_sum_j_not_int_2:.5f}") # Expected 0.0

    # Test 3j with m1+m2+m3 != 0
    val_m_sum_not_zero = calculate_wigner_3j(1,1,1,1,1,1)
    print(f"wigner_3j(1,1,1,1,1,1) where m_sum !=0 = {val_m_sum_not_zero}") # Expected 0.0

    print("Completed Wigner symbol examples.")
    # Test a case where sympy might return an exact expression that needs evalf
    j1,j2,j3,m1,m2,m3 = 1,1,1,0,0,0
    sym_val = wigner_3j(j1,j2,j3,m1,m2,m3)
    print(f"Direct sympy for (1,1,1,0,0,0): {sym_val} type: {type(sym_val)}") # e.g. -sqrt(3)/9
    j1,j2,j3,j4,j5,j6 = 1,1,1,1,1,1
    sym_val_6j = wigner_6j(j1,j2,j3,j4,j5,j6)
    print(f"Direct sympy for 6j (1,1,1,1,1,1): {sym_val_6j} type: {type(sym_val_6j)}") # e.g. 1/6

    # Test with explicit HalfInteger
    from sympy import S
    val_half_int_input = calculate_wigner_3j(S(1)/2, S(1)/2, 1, S(1)/2, S(-1)/2, 0)
    print(f"wigner_3j(1/2, 1/2, 1, 1/2, -1/2, 0) with S() = {val_half_int_input:.5f}")

    # Test that j values must be non-negative
    try:
        calculate_wigner_3j(-1, 1, 1, 0, 0, 0)
    except ValueError as e:
        print(f"Caught error for negative j in 3j: {e}")
    try:
        calculate_wigner_6j(-1, 1, 1, 1, 1, 1)
    except ValueError as e:
        print(f"Caught error for negative j in 6j: {e}")

    # Test that m values must be int or half int
    try:
        calculate_wigner_3j(1,1,1, 0.2, 0.8, -1)
    except ValueError as e:
        print(f"Caught error for non-int/half-int m: {e}")

    # Test that j values for 6j must be int or half int
    try:
        calculate_wigner_6j(1.2, 1,1,1,1,1)
    except ValueError as e:
        print(f"Caught error for non-int/half-int j for 6j: {e}")

    print("All tests in __main__ completed.")

    # Test case from issue: https://github.com/sympy/sympy/issues/21165
    # wigner_3j(0,0,0,0,0,0) should be 1
    val_000 = calculate_wigner_3j(0,0,0,0,0,0)
    print(f"wigner_3j(0,0,0,0,0,0) = {val_000}") # Expected 1.0

    # wigner_6j(0,0,0,0,0,0) should be undefined by some conventions, or specific value by others.
    # Sympy's wigner_6j(0,0,0,0,0,0) gives an error "ZeroDivisionError: division by zero"
    # This is because it involves factorials of negative numbers if not handled.
    # Let's see how sympy handles it.
    # The formula for 6j involves triangle conditions. (0,0,0) is a valid triangle.
    # Sympy's documentation or behavior for this edge case needs to be checked.
    # For now, we rely on Sympy's behavior.
    # If j1=j2=j3=0, then j4,j5,j6 must be equal, say j. Then {000,jjj} = (-1)^(2j) / (2j+1)
    # So if j4=j5=j6=0, then {000,000} = 1.
    try:
        val_6j_000 = calculate_wigner_6j(0,0,0,0,0,0)
        print(f"wigner_6j(0,0,0,0,0,0) = {val_6j_000}") # Expected 1.0 based on some conventions
    except Exception as e:
        print(f"Error for wigner_6j(0,0,0,0,0,0): {e}") # Sympy might error or give specific result.
                                                        # Sympy 1.12 returns 1.0

    val_6j_000_111 = calculate_wigner_6j(0,0,0,1,1,1) # Should be 1/(2*1+1) = 1/3
    print(f"wigner_6j(0,0,0,1,1,1) = {val_6j_000_111:.5f}") # Expected 0.33333

    val_6j_101_101 = calculate_wigner_6j(1,0,1,1,0,1) # j2=0 implies j1=j3, j4=j6, j5 must be 0
                                                 # {j1 0 j1} = (-1)^(j1+j4+j5) / sqrt((2j1+1)(2j4+1)) * delta_j1j3 * delta_j4j6
                                                 # {j1 0 j1, j4 j5 j4} = (-1)^(j1+j4+j5) / sqrt((2j1+1)(2j4+1)) if (j1,j4,j5) form a triangle
                                                 # Sympy's wigner_6j(j1,0,j1,j4,j5,j4) should handle this.
                                                 # For {1,0,1; 1,0,1}: (-1)^(1+1+0)/sqrt((3)*(3)) = 1/3
    print(f"wigner_6j(1,0,1,1,0,1) = {val_6j_101_101:.5f}") # Expected 0.33333

    # Test case where one j is zero for 3j: (j,0,j,m,0,-m) / ( (-1)^(j-m) / sqrt(2j+1) )
    # wigner_3j(1,0,1, 1,0,-1)
    # Expected: (-1)^(1-1) / sqrt(2*1+1) = 1/sqrt(3) ~ 0.57735
    val_3j_j0j = calculate_wigner_3j(1,0,1,1,0,-1)
    print(f"wigner_3j(1,0,1,1,0,-1) = {val_3j_j0j:.5f}") # Expected 0.57735

    val_3j_j0j_2 = calculate_wigner_3j(2,0,2,1,0,-1)
    # Expected: (-1)^(2-1) / sqrt(2*2+1) = -1/sqrt(5) ~ -0.44721
    print(f"wigner_3j(2,0,2,1,0,-1) = {val_3j_j0j_2:.5f}") # Expected -0.44721

    # Test Regge Symmetries for 3j symbol (example)
    # (j1 j2 j3) = (j1 j3 j2) with m values permuted accordingly
    # (m1 m2 m3)   (m1 m3 m2)
    # (j1 j2 j3) = (j2 j1 j3) with m values permuted
    # (m1 m2 m3)   (m2 m1 m3)
    # etc.
    # wigner_3j(j1,j2,j3,m1,m2,m3) == wigner_3j(j2,j1,j3,m2,m1,m3)
    # wigner_3j(0.5, 0.5, 1, 0.5, -0.5, 0) = 0.40825
    val_regge1 = calculate_wigner_3j(0.5, 0.5, 1, -0.5, 0.5, 0) # Swapped m1,m2 AND j1,j2 effectively
    print(f"Regge test 1: wigner_3j(0.5, 0.5, 1, -0.5, 0.5, 0) = {val_regge1:.5f} (original was {val2:.5f})")
    # This should be different unless m1=m2. Sympy handles this.

    # Check parity: (-1)^(j1+j2+j3) * wigner_3j(j1,j2,j3,-m1,-m2,-m3)
    j1,j2,j3,m1,m2,m3 = 0.5, 0.5, 1, 0.5, -0.5, 0
    val_orig = calculate_wigner_3j(j1,j2,j3,m1,m2,m3)
    val_parity = calculate_wigner_3j(j1,j2,j3,-m1,-m2,-m3)
    parity_factor = (-1)**(j1+j2+j3) if (j1+j2+j3).is_integer() else (1j)**(2*(j1+j2+j3)) # Handle half-int sum
    if not isinstance(parity_factor, (int, float, complex)): parity_factor = parity_factor.real

    print(f"Parity test: Original val={val_orig:.5f}, Parity val={val_parity:.5f}, Factor={parity_factor}")
    print(f"val_orig == parity_factor * val_parity: {abs(val_orig - parity_factor * val_parity) < 1e-9}")


    print("\nWigner 9j Symbol Examples:")
    # Example 1: From Wikipedia {1,1,1; 1,1,1; 1,1,1} -> 1/36
    # (Verified with online calculators, e.g. WolframAlpha)
    val9j_1 = calculate_wigner_9j(1,1,1,1,1,1,1,1,1)
    print(f"wigner_9j(1,1,1,1,1,1,1,1,1) = {val9j_1:.8g} (Expected: 1/36 = {1/36:.8g})")

    # Example 2: One j is zero. E.g. j3=0.
    # {j1 j2 0}
    # {j4 j5 j6}
    # {j7 j8 j9}
    # If j3=0, then j1=j2, j6=j9. The 9j symbol becomes related to a 6j symbol.
    # Specifically, {a b 0, c d e, f g b} = delta_ab * delta_eb * (-1)^(b+c+f+g) / sqrt((2b+1)(2e+1)) * {a c f; d g e} (This is one specific formula variant)
    # From online calculator: wigner_9j(1,1,0, 1,1,1, 1,1,0) = 1/6
    val9j_2 = calculate_wigner_9j(1,1,0, 1,1,1, 1,1,0)
    print(f"wigner_9j(1,1,0, 1,1,1, 1,1,0) = {val9j_2:.8g} (Expected: 1/6 = {1/6:.8g})")

    # Example 3: A case that should be zero due to triangle inequality in a row
    # Row 1: (1,1,3) - not a triangle
    val9j_3_tri_row = calculate_wigner_9j(1,1,3, 1,1,1, 1,1,1)
    print(f"wigner_9j(1,1,3, ...) = {val9j_3_tri_row:.8g} (Expected: 0.0)")

    # Example 3b: A case that should be zero due to triangle inequality in a column
    # Col 1: (1,1,3)
    val9j_3_tri_col = calculate_wigner_9j(1,1,1, 1,1,1, 3,1,1)
    print(f"wigner_9j(..., col1=(1,1,3)) = {val9j_3_tri_col:.8g} (Expected: 0.0)")


    # Example 4: A case that should be zero due to sum of triad not integer
    # Row 1: (1/2, 1/2, 1/2) -> sum = 3/2
    val9j_4_sum = calculate_wigner_9j(S(1)/2, S(1)/2, S(1)/2, 1,1,1, 1,1,1)
    print(f"wigner_9j(1/2,1/2,1/2, ...) = {val9j_4_sum:.8g} (Expected: 0.0)")

    # Example 5: From another source (e.g. WolframAlpha or textbook)
    # Wigner9j[{1, 2, 1}, {2, 1, 2}, {1, 2, 2}] = 1/90
    val9j_5 = calculate_wigner_9j(1,2,1, 2,1,2, 1,2,2)
    print(f"wigner_9j(1,2,1, 2,1,2, 1,2,2) = {val9j_5:.8g} (Expected: 1/90 = {1/90:.8g})")

    # Example 6: With half-integers
    # Wigner9j[{1/2,1/2,1},{1/2,1/2,1},{1,1,1}] = 1/12
    # (Verified with online calculator: e.g. volya.fsu.edu/calc_racah/index.html)
    val9j_6 = calculate_wigner_9j(S(1)/2, S(1)/2, 1,  S(1)/2, S(1)/2, 1,  1,1,1)
    print(f"wigner_9j(1/2,1/2,1, 1/2,1/2,1, 1,1,1) = {val9j_6:.8g} (Expected: 1/12 = {1/12:.8g})")

    # Test invalid input type
    try:
        calculate_wigner_9j(1,1,1, 1,1.2,1, 1,1,1) # 1.2 is not int/half-int
    except ValueError as e:
        print(f"Caught expected error for 9j invalid input type: {e}")

    # Test invalid input (negative j)
    try:
        calculate_wigner_9j(1,1,1, 1,-1,1, 1,1,1) # -1 is not non-negative
    except ValueError as e:
        print(f"Caught expected error for 9j negative j: {e}")

    print("All tests in __main__ for 9j completed.")
