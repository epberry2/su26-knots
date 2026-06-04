from dataclasses import dataclass, field

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
    points: tuple[str, str, str]
    lowest_power: int = 0
    winding_nums: list[int] = field(default_factory=lambda: [0, 0, 0])
    arc_right: bool = field(init = False)

    def __post_init__(self):
        self.arc_right = (self.points[1] == "X+")
        
    def increment_l(self, amount = 1):
        self.winding_nums[0] += amount
        
    def increment_c(self, amount = 1):
        self.winding_nums[1] += amount
        
    def increment_r(self, amount = 1):
        self.winding_nums[2] += amount
    
    def power(self):
        pow = -1 * 2 * (self.winding_nums[0] + self.winding_nums[1] + self.winding_nums[2])
        self.lowest_power = min(self.lowest_power, pow)
        return pow
        

def store_monomial(monomials, state, index, parity):
    if state.arc_right: parity *= -1
    monomials[index] = (parity, state.power())
    
def intersection_index(num, denom, point, next_point, path_type, parity):
    if path_type == "C":
        index = 2 * (min(point, next_point) - 1) + int(denom % 2 == 0)
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
            state.increment_l(-1 if is_ccw else 1)
        if path_type == "C":
            if is_ccw:
                state.increment_c(-1) # increment beforehand
                intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                store_monomial(monomial_array, state, intersection_num, parity)
            else:
                intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                store_monomial(monomial_array, state, intersection_num, parity)
                state.increment_c(1) # increment aftwerwards
        if path_type == "R":
            intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
            if not is_ccw: # increment beforehand
                state.increment_r(1)
            store_monomial(monomial_array, state, intersection_num, parity)
            if is_ccw: # increment afterwards
                state.increment_r(-1)

