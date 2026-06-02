from findloops import findloops
from findloops import findpath
from findloops import findpathwithloops
from sympy import symbols, Poly
import sympy as sp
from JonesPolyAlg import AlgBracket

q = symbols('q')

def intersectionAboveAxis(num, denom, point):
    if point <= num // 2:
        return False
    elif point <= num and point > (num + 1) // 2:
        return True
    elif point < num:
        return (num >= denom % (2 * num))
    elif point - num <= denom // 2:
        return False
    return True

def rotateCCW(num, denom, point, nextPoint, pathType):
    if pathType == "L":
        return intersectionAboveAxis(num, denom, point)
    if pathType == "R":
        return not intersectionAboveAxis(num, denom, point)
    if num >= denom:
        return point < nextPoint
    return point > nextPoint

def coeff(pathType):
    if pathType == "L" or pathType == "R":
        return 1
    return -1

def JonesPolyGeo(num, denom):
    path = findpathwithloops(num, denom)[0]
    loops = findpathwithloops(num, denom)[1]
    if num < denom and denom % 2 == 1:
        startOnLeft = False
    if num < denom and denom % 2 == 0:
        startOnLeft = True
    if num >= denom and num % 2 == 1:
        startOnLeft = True
    if num >= denom and num % 2 == 0:
        startOnLeft = False
    # This matches whether or not the path and loops array start on the left-hand side when drawing.
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

# Test case

# for p in range(1, 20):
#     for r in range(1, 20):
#         if sp.gcd(p, r) > 1:
#             continue
#         jonesalg = AlgBracket(p, r)
#         jonesgeo = JonesPolyGeo(p, r)
#         a0 = jonesalg[0, 0].as_poly(q)
#         a1 = jonesalg[1, 0].as_poly(q)
#         j0 = jonesgeo[0]
#         j1 = jonesgeo[1]

#         if  (a0 == j0 and a1 == j1) or (-a0 == j0 and -a1 == j1):
#             continue
#         else:
#             print("broke at (" + str(p) + ", " + str(r) + ")")
#             break
