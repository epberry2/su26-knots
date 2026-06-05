from helpers.findloops import findloops
from helpers.findloops import findpath
from helpers.findloops import findpathwithloops
from sympy import symbols, Poly, gcd
from jones.JonesPolyAlg import AlgBracket, AlgJonesKnot
from jones.GeoJonesKnotHelpers import *
from helpers.tanglestate import get_state
from jones.GeoJonesKnot import geoJonesKnot

q = symbols('q')
a = symbols('a')


def geoHomflyKnot(num, denom):
    # Strategy: Walk along loop with right hand to wall, keeping track of how many times you have looped
    # around each point. At the end, handle the figure-8 loop, and then return.
    # If a portion of the figure-8 loop ends with a half C-arc, then the curve
    # obtained by following the figure-8 half-loop with RH to wall for that portion
    # will wind an extra time around the middle point due to the half-loop iff the half-loop is CCW
    # when RH is to wall; likewise, for a half R-arc, iff the half-loop is CW.
    # Therefore, we include the arcs in the monomials at the midpoints depending on their orientation.
    
    assert num >= denom, "only rationals >= 1 considered at the moment"
    assert num % 2 == 1, "only odd numerators considered at the moment"
    
    orientation, points_ = get_state(num, denom)
    assert orientation == "UP" or orientation == "OP", "tangle has RI orientation"
    
    path, loops = findpathwithloops(num, denom)
    curr_state = LoopWalkState("H", points_) # stores winding number around each vertex, as well as lowest power to normalize
    monomial_array = [(0, 0)] * num
    tracePath(num, denom, path, loops, curr_state, monomial_array)
    # handle endpoint
    if denom % 2 == 0: 
        idx = 0
        parity = 1
    else:
        curr_state.increment_r(1) # loop around endpoint first
        idx = num - 1
        parity = -1
    store_monomial(monomial_array, curr_state, idx, parity) # parity = 1: RH is to wall
    if denom % 2 == 0: curr_state.increment_c(1) 
    # walk backwards
    tracePath(num, denom, list(reversed(path)), list(reversed(loops)), curr_state, monomial_array)
    # done!
    homfly_polynomial = Poly(0, q, a)
    for coeff, m_exp, n_exp in monomial_array:
        homfly_polynomial += Poly(coeff * q ** (n_exp - curr_state.lowest_n_power)
                                  * a ** (m_exp - curr_state.lowest_m_power), q, a)
    return homfly_polynomial

def homflyToJones(homfly_polynomial):
    jones_polynomial = Poly(homfly_polynomial.as_expr().subs(a, q**2), q)
    lowest_jones_power = min(exp[0] for exp in jones_polynomial.monoms())
    jones_polynomial = Poly(jones_polynomial.as_expr() * q ** (-lowest_jones_power), q)
    return jones_polynomial


homfly = geoHomflyKnot(13, 3)
jones = homflyToJones(homfly)
jonestest = geoJonesKnot(13, 3)
print(homfly.as_expr())
print(jones.as_expr())
print(jonestest.as_expr())