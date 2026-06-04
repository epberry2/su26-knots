from findloops import findloops
from findloops import findpath
from findloops import findpathwithloops
from sympy import symbols, Poly, gcd
from JonesPolyAlg import AlgBracket, AlgJonesKnot
from JonesPolyGeoHelpers import *
from tanglestate import get_state

q = symbols('q')


def JonesPolyGeo(num, denom):
    path = findpathwithloops(num, denom)[0]
    loops = findpathwithloops(num, denom)[1]
    if num >= denom:
        startOnLeft = (num % 2 == 1)
    else:
        startOnLeft = (denom % 2 == 0)
    # This matches whether or not the path and loops array start on the left-hand side
    leftLine = startOnLeft
    if startOnLeft:
        PolyLeft = Poly(1, q)
        PolyRight = Poly(0, q)
    else:
        PolyLeft = Poly(0, q)
        PolyRight = Poly(1, q)
    currPow = 0
    currPar = 1
    for i in range(len(loops)):
        currPoint = path[i]
        pathType = loops[i]
        aboveZero = intersectionAboveAxis(num, denom, currPoint)
        if pathType == 'T':
            if leftLine:
                currPow -= 1
            else:
                currPow += 1
            leftLine = not leftLine
        else:
            nextPoint = path[i + 1]
            if rotateCCW(num, denom, currPoint, nextPoint, pathType):
                currPow += 2
            else:
                currPow -= 2
        currPar = coeff(pathType)
        if currPow < 0:
            PolyLeft *= Poly(q**(-currPow), q)
            PolyRight *= Poly(q**(-currPow), q)
            currPow = 0
        if leftLine:
            PolyLeft += Poly(currPar * (q ** currPow), q)
        else:
            PolyRight += Poly(currPar * (q ** currPow), q)
    return (PolyLeft, PolyRight)


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
    curr_state = CurrState(points = points_) # stores winding number around each vertex, as well as lowest power to normalize
    monomial_array = [(0, 0)] * num
    tracePath(num, denom, path, loops, curr_state, monomial_array)
    # handle endpoint
    if denom % 2 == 0: 
        curr_state.increment_r(-1) # loop around endpoint first
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
    for coeff, exp in monomial_array:
        jones_polynomial += Poly(coeff * q ** (exp - curr_state.lowest_power), q)
    return jones_polynomial
    
    
print(geoJonesKnot(17, 4).as_expr())
print(geoJonesKnot(5,2).as_expr())   # works!
# print(AlgJonesKnot(3,1).as_expr())
# print(geoJonesKnot(19,11).as_expr()) # works!
# print(AlgJonesKnot(19,11).as_expr())



# for i in range(50):
#     for j in range(50):
#         if i < j: continue
#         if i % 2 == 0 or j % 2 == 0: continue
#         if gcd(i, j) > 1: continue
#         skip = False
#         for k in range(j):
#             if (k * j) % i == 1 or (k * j) % i == i - 1:
#                 if k % 2 == 1:
#                     skip = True
#         if skip: continue
#         if(geoJonesKnot(i, j) == AlgJonesKnot(i, j)):
#             print(f"True {{{i}/{j}}}:")
#         else: print(f"False {{{i}/{j}}}:")
