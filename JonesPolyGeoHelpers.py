from dataclasses import dataclass

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

@dataclass
class CurrState:
    lowest_power: int = 0
    wind_l: int = 0
    wind_c: int = 0
    wind_r: int = 0



def store_monomial(monomials, state, index, parity):
    power = 2 * (state.wind_l + state.wind_c + state.wind_r)
    state.lowest_power = min(state.lowest_power, power)
    monomials[index] = (parity, power)
    
def intersection_index(num, denom, point, next_point, path_type, parity):
    if path_type == "C":
        index = 2 * (min(point, next_point) - 1)
    elif path_type == "R":
        index = 2 * (min(point, next_point) - 1 - num) + (num - denom)
    else:
        raise ValueError(f"Only intersect at C or R arcs")
    return index + int(parity == 1) # We walk with RH to wall, so on right iff parity = 1

def tracePath(num, denom, path, loops, state, monomial_array):
    for point, next_point, path_type in zip(path, path[1:], loops):   # fancy for loop
        if path_type == "T": continue
        is_ccw = rotateCCW(num, denom, point, next_point, path_type)
        parity = -1 if (point > next_point) ^ (path_type == "C") else 1
        if path_type == "L":
            state.wind_l += (-1 if is_ccw else 1)
        if path_type == "C":
            if is_ccw:
                state.wind_c -= 1 # increment beforehand
                intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                store_monomial(monomial_array, state, intersection_num, parity)
            else:
                intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                store_monomial(monomial_array, state, intersection_num, parity)
                state.wind_c += 1 # increment aftwerwards
        if path_type == "R":
            intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
            if not is_ccw: # increment beforehand
                state.wind_r += 1
            store_monomial(monomial_array, state, intersection_num, parity)
            if is_ccw: # increment afterwards
                state.wind_r -= 1

