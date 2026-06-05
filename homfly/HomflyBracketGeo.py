import helpers.findloops as fl
import sympy as sp
from helpers.tanglestate import get_state
from jones.JonesBracketGeoHelpers import *

q = sp.symbols("q")
a = sp.symbols("a")

def HomflyBracketGeo(num, denom):
    path, loops = fl.findpathwithloops(num, denom)    
    startOnLeft = (path[0] <= num)
    state = tangleRotateState("H", startOnLeft, get_state(num, denom)[1])
    poly_left = sp.Poly(0, q)
    poly_right = sp.Poly(0, q)
    left_monomials = []
    right_monomials = []
    
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        parity, m_power, n_power = state.powers()
        if state.on_left:
            left_monomials.append((parity, m_power, n_power))
        else:
            right_monomials.append((parity, m_power, n_power))
        state.traverse(rotateCCW(num, denom, curr_point, next_point, path_type), path_type)
        
    parity, m_power, n_power = state.powers() # didn't hit last point
    if state.on_left:
        left_monomials.append((parity, m_power, n_power))
    else:
        right_monomials.append((parity, m_power, n_power))
    
    min_m_power = state.min_m_power
    min_n_power = state.min_n_power
    for parity, m_power, n_power in left_monomials:
        poly_left += sp.Poly(parity * q ** (n_power - min_n_power) *  a ** (m_power - min_m_power), q, a)
    for parity, m_power, n_power in right_monomials:
        poly_right += sp.Poly(parity * q ** (n_power - min_n_power) * a ** (m_power - min_m_power), q, a)
    
    return (poly_left, poly_right)


for poly in HomflyBracketGeo(5, 1):
    print(poly.as_expr())
    
for poly in HomflyBracketGeo(5, 2):
    print(poly.as_expr())