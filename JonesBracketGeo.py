from findloops import findloops
from findloops import findpath
from findloops import findpathwithloops
from sympy import symbols, Poly, gcd
from JonesPolyAlg import AlgBracket, AlgJonesKnot
from JonesBracketGeoHelpers import *
from tanglestate import get_state

q = symbols('q')


def JonesBracketGeo(num, denom):
    path, loops = findpathwithloops(num, denom)    
    startOnLeft = (path[0] <= num)
    state = tangleRotateState("J", startOnLeft)
    poly_left = Poly(0, q)
    poly_right = Poly(0, q)
    left_monomials = []
    right_monomials = []
    
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        parity, _, n_power = state.powers()
        if state.on_left:
            left_monomials.append((parity, n_power))
        else:
            right_monomials.append((parity, n_power))
        state.traverse(rotateCCW(num, denom, curr_point, next_point, path_type), path_type)
        
    parity, _, n_power = state.powers() # didn't hit last point
    if state.on_left:
        left_monomials.append((parity, n_power))
    else:
        right_monomials.append((parity, n_power))
    
    min_power = state.min_n_power
    for parity, power in left_monomials:
        poly_left += Poly(parity * q ** (power - min_power), q)
    for parity, power in right_monomials:
        poly_right += Poly(parity * q ** (power - min_power), q)
    
    return (poly_left, poly_right)

for poly in JonesBracketGeo(5, 2):
    print(poly.as_expr())