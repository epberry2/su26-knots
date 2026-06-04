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
