from helpers.findloops import findloops
from helpers.findloops import findpath
from helpers.findloops import findpathwithloops
from sympy import symbols, Poly, gcd
from jones.JonesPolyAlg import AlgBracket, AlgJonesKnot
from jones.GeoJonesKnotHelpers import *
from helpers.tanglestate import get_state

q = symbols('q')


def geoJonesKnot(num, denom):
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
    curr_state = LoopWalkState("J", points_) # stores winding number around each vertex, as well as lowest power to normalize
    monomial_array = [(0, 0)] * num
    tracePath(num, denom, path, loops, curr_state, monomial_array)
    # handle endpoint
    if denom % 2 == 0: 
        curr_state.increment_c(-1) # loop around endpoint first
        idx = 0
        parity = -1
    else:
        idx = num - 1
        parity = 1
    store_monomial(monomial_array, curr_state, idx, parity) # parity = 1: RH is to wall
    if denom % 2 == 1: curr_state.increment_r(-1) 
    # walk backwards
    tracePath(num, denom, list(reversed(path)), list(reversed(loops)), curr_state, monomial_array)
    # done!
    jones_polynomial = Poly(0, q)
    for coeff, _, exp in monomial_array:
        jones_polynomial += Poly(coeff * q ** (exp - curr_state.lowest_n_power), q)
    return jones_polynomial
