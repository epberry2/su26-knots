import findloops as fl
import sympy as sp
from collections import defaultdict

q = sp.symbols("q")

def GeoJonesKnot(num,denom):
    u = num - denom # take inverse top twist, so u/v = (u-v)/v
    v = denom
    path, loops = fl.findpathwithloops(u, v)
    #print(path)
    jones = 1
    currPow = 0
    
    sign = 0 # we change sign every L or R loop
    qsign = 1 # the sign of the coefficients


    
    for i in range(len(loops)):
        
        if loops[i] == 'L':
            if path[i+1] > path[i]:
                currPow -= 2
            else:
                currPow += 2
            sign = 1 - sign

        elif loops[i] == 'R':
            if path[i+1] < path[i]:
                currPow -= 2
            else:
                currPow += 2
            sign = 1 - sign
        
        elif loops[i] == 'C':
            if path[i+1] < path[i]:
                currPow -= 2
                if path[i] > u // 2 and path[i+1] < u // 2:
                    sign = 1 - sign
            else:
                currPow += 2
                if path[i+1] > u // 2 and path[i] < u // 2:
                    sign = 1 - sign
        else:
            # T arc
            if path[i+1] > path[i]:
                currPow -= 2            
                if (path[i+1] - u) < (v + 1) / 2 and sign == 1 or (path[i+1] - u) > (v + 1) / 2 and sign == 0:
                    sign = 1 - sign
            else:
                currPow += 2
                if path[i+1] <= (u + 1) / 2 and sign == 1 or path[i+1] > (u + 1) / 2 and sign == 0:
                    sign = 1 - sign


        if sign == 0:
            jones += (-1)**qsign * q**currPow

        
        qsign = 1 - qsign
        #print(currPow, sign)
        
            
    # now we have traversed the outside blue path, so traverse the inside segment in reverse
    
    if u < v:
        currPow += 2
    else:
        currPow -= 2
        sign = 1 - sign
    
    if sign == 0:
        jones += (-1)**qsign * q**currPow

    qsign = 1 - qsign
    loops.reverse()
    path.reverse()

    #print(currPow, sign)

    for i in range(len(loops)):
        
        if loops[i] == 'L':
            if path[i+1] > path[i]:
                currPow -= 2
            else:
                currPow += 2
            sign = 1 - sign

        elif loops[i] == 'R':
            if path[i+1] < path[i]:
                currPow -= 2
            else:
                currPow += 2
            sign = 1 - sign

        elif loops[i] == 'C':
            if path[i+1] < path[i]:
                currPow -= 2
                if path[i] > u // 2 and path[i+1] < u // 2:
                    sign = 1 - sign
            else:
                currPow += 2
                if path[i+1] > u // 2 and path[i] < u // 2:
                    sign = 1 - sign

        else:
            # T arc
            if path[i+1] > path[i]:
                currPow -= 2            
                if (path[i+1] - u) < (v + 1) / 2 and sign == 1 or (path[i+1] - u) > (v + 1) / 2 and sign == 0:
                    sign = 1 - sign
            else:
                currPow += 2
                if path[i+1] <= (u + 1) / 2 and sign == 1 or path[i+1] > (u + 1) / 2 and sign == 0:
                    sign = 1 - sign
            
        if sign == 0:

            jones += (-1)**qsign * q**currPow

        qsign = 1 - qsign
        #print(currPow, sign)

    powers = []
    for term in jones.as_coefficients_dict().keys():
        if term == 1:  # It's a constant term (x^0)
            powers.append(0)
        elif term.is_Pow:
            powers.append(term.exp)
        elif term == q:
            powers.append(1)

    min_power = min(powers)
    jones = (q**(-min_power) * jones).expand()

    return sp.Poly(jones, q)

print(GeoJonesKnot(5,3))