from dataclasses import dataclass, field
from jones.JonesBracketGeoHelpers import rotateCCW

@dataclass
class LoopWalkState:
    poly_type: str
    points: tuple[str, str, str]
    winding_nums: list[int] = field(default_factory=lambda: [0, 0, 0])
    lowest_m_power: int = field(init = False, default = 0)
    lowest_n_power: int = field(init = False, default = 0)
    arc_right: bool = field(init = False)
    x_plus_idx: int = field(init = False, default = 0)
    x_minus_idx: int = field(init = False, default = 0)
    y_idx: int = field(init = False, default = 0)

    def __post_init__(self):
        self.x_plus_idx = self.points.index("X+")
        self.x_minus_idx = self.points.index("X-") # should be 0
        self.y_idx = self.points.index("Y")
        self.arc_right = (self.x_plus_idx == 1)
        
    def increment_l(self, amount = 1):
        self.winding_nums[0] += amount
        
    def increment_c(self, amount = 1):
        self.winding_nums[1] += amount
        
    def increment_r(self, amount = 1):
        self.winding_nums[2] += amount
    
    def powers(self):
        match self.poly_type:
            case "J":
                pow = -1 * 2 * (self.winding_nums[0] + self.winding_nums[1] + self.winding_nums[2])
                self.lowest_n_power = min(self.lowest_n_power, pow)
                return 0, pow
            case "H":
                n_pow = -2 * (self.winding_nums[self.x_minus_idx] + self.winding_nums[self.y_idx]
                              - self.winding_nums[self.x_plus_idx])
                m_pow = -2 * (self.winding_nums[self.x_plus_idx])
                self.lowest_m_power = min(self.lowest_m_power, m_pow)
                self.lowest_n_power = min(self.lowest_n_power, n_pow)
                return m_pow, n_pow
    
    def quiver_diag_term(self):
        return (self.winding_nums[self.x_minus_idx] + self.winding_nums[self.y_idx]
                    - 3 * self.winding_nums[self.x_plus_idx])

def store_monomial(monomials, state, index, parity):
    if state.poly_type == "H": parity *= -1
    if state.arc_right: parity *= -1
    monomials[index] = (parity, *state.powers())
    
def intersection_index(num, denom, point, next_point, path_type):
    """
    Inputs an arc on the loop. Assuming we walk with RH to tangle.
    
    If the arc is of C or R type, it intersects the beta arc.
    Returns the 0-indexed point of intersection of the inputted arc.
    
    If the arc is of T type, it can be helpful to consider the extension
    of the arc on the right until it hits beta. \\ 
    In this case, the arc may loop around the right point, or be a "shortened
    T arc" which crosses beta immediately to the right. \\
    Returns num (to the right of all intersections) if not shortened,
    or the index of the immediate rightward intersection if shortened.

    If the arc is of L type, returns -1 (to the left of all intersections)
    """
    parity = -1 if (point > next_point) ^ (path_type == "C") else 1
    if path_type == "C":
        index = 2 * (min(point, next_point) - 1) + int(denom % 2 == 0)
    elif path_type == "R":
        index = 2 * (min(point, next_point) - 1 - num) + (num - denom)
    elif path_type == "T": #return right intersection
        higher_point = point if point > next_point else next_point
        if higher_point > num + (denom // 2):
            return num # goes over all intersections
        index = 2 * (higher_point - 1 - num) + (num - denom)
        parity = 1 if point > next_point else -1
    else:
        return -1 # left arc to left of all intersections
    return index + int(parity == 1) # We walk with RH to wall, so on right iff parity = 1

def tracePath(num, denom, path, loops, state, monomial_array = [], store_array = True):
    for point, next_point, path_type in zip(path, path[1:], loops):   # fancy for loop
        is_ccw = rotateCCW(num, denom, point, next_point, path_type)
        parity = -1 if (point > next_point) ^ (path_type == "C") else 1
        match path_type:
            case "L":
                state.increment_l(-1 if is_ccw else 1)
            case "C":
                if is_ccw:
                    state.increment_c(-1) # increment beforehand
                    intersection_num = intersection_index(num, denom, point, next_point, path_type)
                    if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                else:
                    intersection_num = intersection_index(num, denom, point, next_point, path_type)
                    if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                    state.increment_c(1) # increment aftwerwards
            case "R":
                intersection_num = intersection_index(num, denom, point, next_point, path_type)
                if not is_ccw: # increment beforehand
                    state.increment_r(1)
                if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                if is_ccw: # increment afterwards
                    state.increment_r(-1)
