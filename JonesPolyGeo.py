from findloops import findloops
from findloops import findpath
from findloops import findpathwithloops
from sympy import symbols, Poly
import sympy as sp

q = symbols('q')

def JonesPolyGeo(num, denom):
    path = findpathwithloops(num, denom)
    startOnLeft = num >= denom
    leftLine = startOnLeft
    # PolyLeft, PolyRight
    if startOnLeft:
        PolyLeft = Poly(1, q)
        PolyRight = Poly(0, q)
    else:
        PolyLeft = Poly(0, q)
        PolyRight = Poly(1, q)
    CurrPow = 0
    CurrPar = 1
    lowestPow = 0
    for i in range(1, len(path) - 1):
        # loops_i = findloops(num, denom, path[i], path[i+1])
        pathType = path[i]
        if pathType == 'C' or pathType == 'L' or pathType == 'R' or pathType == 'T':
            if pathType == 'C':
                currPoint = (int(path[i - 1]))
                nextPoint = (int(path[i + 1]))
                if startOnLeft:
                    if aboveZero:
                        if (nextPoint < currPoint):
                            CurrPow = CurrPow - 2
                            CurrPar = -1
                        else: 
                            CurrPow = CurrPow + 2
                            CurrPar = -1
                    else:
                        if (nextPoint > currPoint):
                            CurrPow = CurrPow + 2
                            CurrPar = -1
                        else: 
                            CurrPow = CurrPow - 2
                            CurrPar = -1
                else:
                    if aboveZero:
                        if (nextPoint < currPoint):
                            CurrPow = CurrPow + 2
                            CurrPar = -1
                        else:
                            CurrPow = CurrPow - 2
                            CurrPar = -1
                    else:
                        if (nextPoint > currPoint):
                            CurrPow = CurrPow - 2
                            CurrPar = -1
                        else:
                            CurrPow = CurrPow + 2
                            CurrPar = -1
            elif pathType == 'L':
                if aboveZero:
                    CurrPow = CurrPow + 2
                    CurrPar = 1
                else:
                    CurrPow = CurrPow - 2
                    CurrPar = 1
            elif pathType == 'R':
                if aboveZero:
                    CurrPow = CurrPow - 2
                    CurrPar = 1
                else:
                    CurrPow = CurrPow + 2
                    CurrPar = 1
            elif pathType == 'T':
                if leftLine:
                    CurrPow -= 1
                    CurrPar = -1
                else:
                    CurrPow += 1
                    CurrPar = -1
                leftLine = not leftLine
            if CurrPow < 0:
                PolyLeft *= Poly(q**(-CurrPow), q)
                PolyRight *= Poly(q**(-CurrPow), q)
                CurrPow = 0
            if leftLine:
                PolyLeft += Poly(CurrPar * (q ** CurrPow), q)
                print("adding " + str( CurrPar * (q ** CurrPow)))
            else:
                PolyRight += Poly(CurrPar * (q ** CurrPow), q)
                print("adding " + str( CurrPar * (q ** CurrPow)))
        else: 
            point = int(path[i])
            if point <= num // 2:
                aboveZero = False
                print(str(point) + " below zero")
            elif point <= num and point > (num + 1) // 2:
                aboveZero = True
                print(str(point) + " above zero")
            elif point < num:
                aboveZero = num >= denom % (2 * num)
            elif point - num <= denom // 2:
                aboveZero = False
                print(str(point) + " below zero")
            else:
                aboveZero = True
                print(str(point) + " above zero")
        if CurrPow < lowestPow:
            lowestPow = CurrPow
    return (PolyLeft, PolyRight)



print(tuple(p.as_expr() for p in JonesPolyGeo(19, 17)))
#print(tuple(p.as_expr() for p in JonesPolyGeo(2, 5)))