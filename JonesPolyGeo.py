from findloops import findloops
from findloops import findpath
from findloops import findpathwithloops
from sympy import symbols, Poly
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

# Test case

# for p in range(1, 20):
#     for r in range(1, 20):
#         if sp.gcd(p, r) > 1:
#             continue
#         jonesalg = AlgBracket(p, r)
#         jonesgeo = JonesPolyGeo(p, r)
#         a0 = jonesalg[0].as_poly(q)
#         a1 = jonesalg[1].as_poly(q)
#         j0 = jonesgeo[0]
#         j1 = jonesgeo[1]

#         if  (a0 == j0 and a1 == j1) or (-a0 == j0 and -a1 == j1):
#             continue
#         else:
#             print("broke at (" + str(p) + ", " + str(r) + ")")
#             break
# print("done")

# def printDebug(intersectionNum, adjust, windL, windC, windR, parity):
#     # print(f"hit intersection {intersectionNum}")
#     print(f"monomialArr[{intersectionNum}] = {parity * Poly(q**(2 * (adjust + windL + windC + windR)), q).as_expr()}")
#     print(f"windL = {windL}, windC = {windC} windR = {windR}")
    
# def adjustArr(arr, adjustBy):
#     for pol in arr:
#         pol *= Poly(q**adjustBy, q)

def geoJonesKnot(num, denom):
    # Strategy: Go through the path. You hit a point on the right arc when you 
    # do a C arc or an R arc. You swap left and right when you do a C arc;
    # if it is CCW you swap immediately, otherwise, you swap afterwards.
    # When you do a CW L arc, your winding number around
    # the left point increases; a CW C arc increases it around the center point;
    # and a CW R arc increases it around the rightmost point. When you reach the
    # end, swap left and right and trace backwards to hit the remaining points.
    assert num >= denom, "only rationals >= 1 considered at the moment"
    assert num % 2 == 1, "only odd numerators considered at the moment"
    
    path = findpathwithloops(num, denom)[0]
    loops = findpathwithloops(num, denom)[1]
    
    if num >= denom:
        startOnLeft = (num % 2 == 1)
    else:
        startOnLeft = (denom % 2 == 0)
        
    if not startOnLeft:
        path.reverse()
        loops.reverse()
    
    whichSide = True   # False = Left, True = Right, start on
    lowestPow = 0
    windL = 0
    windC = 0
    windR = 0
    monomialArr = [(0, 0)] * (num)
    
    for i in range(len(loops)):   # len(loops) + 1 = len(path)
        pathType = loops[i]
        if pathType == "T": continue
        point = path[i]
        nextPoint = path[i + 1]
        isCCW = rotateCCW(num, denom, point, nextPoint, pathType)
        isCW = not isCCW
        parity = -1 if (point > nextPoint) ^ (pathType == "C") else 1
        if pathType == "L":
            if isCCW:
                windL -= 1
            else: windL += 1
        if pathType == "C":
            intersectionNum = (min(point, nextPoint) - 1) * 2
            if isCCW:
                whichSide = not whichSide
                if whichSide: intersectionNum += 1
                windC -= 1 # increment beforehand
                power = 2 * (windL + windC + windR)
                if lowestPow > power: lowestPow = power
                monomialArr[intersectionNum] = (parity, power)
                # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
            else:
                if whichSide: intersectionNum += 1
                power = 2 * (windL + windC + windR)
                if lowestPow > power: lowestPow = power
                monomialArr[intersectionNum] = (parity, power)
                # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
                windC += 1 # increment aftwerwards
                whichSide = not whichSide
        if pathType == "R":
            intersectionNum = 2 * (min(point, nextPoint) - 1 - num) + (num - denom)
            if whichSide: intersectionNum += 1
            if isCW: # increment beforehand
                windR += 1
            power = 2 * (windL + windC + windR)
            if lowestPow > power: lowestPow = power
            monomialArr[intersectionNum] = (parity, power)
            # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
            if isCCW: # increment afterwards
                windR -= 1
    # end loop
    

    if intersectionAboveAxis(num, denom, path[-1]):
        if whichSide:
            windR += 1
            parity = 1
    elif not whichSide: # ?????? denom > num
        windR -= 1
        parity = -1
    intersectionNum = num - 1
    power = 2 * (windL + windC + windR)
    if lowestPow > power: lowestPow = power
    monomialArr[intersectionNum] = (parity, power)
    # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
    whichSide = not whichSide
    if intersectionAboveAxis(num, denom, path[-1]):
        if whichSide:
            windR -= 1
            parity = -1
    elif not whichSide: # ????? denom > num
            windR += 1
            parity = 1

    path.reverse()
    loops.reverse()

    for i in range(len(loops)):   # len(loops) + 1 = len(path)
        pathType = loops[i]
        if pathType == "T": continue
        point = path[i]
        nextPoint = path[i + 1]
        isCCW = rotateCCW(num, denom, point, nextPoint, pathType)
        isCW = not isCCW
        parity = -1 if (point > nextPoint) ^ (pathType == "C") else 1
        if pathType == "L":
            if isCCW:
                windL -= 1
            else: windL += 1
        if pathType == "C":
            intersectionNum = (min(point, nextPoint) - 1) * 2
            if isCCW:
                whichSide = not whichSide
                if whichSide: intersectionNum += 1
                windC -= 1 # increment beforehand
                power = 2 * (windL + windC + windR)
                if lowestPow > power: lowestPow = power
                monomialArr[intersectionNum] = (parity, power)
                # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
            else:
                if whichSide: intersectionNum += 1
                power = 2 * (windL + windC + windR)
                if lowestPow > power: lowestPow = power
                monomialArr[intersectionNum] = (parity, power)
                # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
                windC += 1 # increment aftwerwards
                whichSide = not whichSide
        if pathType == "R":
            intersectionNum = 2 * (min(point, nextPoint) - 1 - num) + (num - denom)
            if whichSide: intersectionNum += 1
            if isCW: # increment beforehand
                windR += 1
            power = 2 * (windL + windC + windR)
            if lowestPow > power: lowestPow = power
            monomialArr[intersectionNum] = (parity, power)
            # printDebug(intersectionNum, adjust, windL, windC, windR, parity)
            if isCCW: # increment afterwards
                windR -= 1
    jonesPolynomial = Poly(0, q)
    for (coeff, exp) in monomialArr:
        jonesPolynomial += Poly(coeff * q ** (exp - lowestPow), q)
    return jonesPolynomial
    
print(geoJonesKnot(3,1).as_expr())
print(geoJonesKnot(5,3).as_expr())
