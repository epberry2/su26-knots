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
    
def intersection_index(num, denom, point, next_point, path_type, parity):
    if path_type == "C":
        index = 2 * (min(point, next_point) - 1) + int(denom % 2 == 0)
    elif path_type == "R":
        index = 2 * (min(point, next_point) - 1 - num) + (num - denom)
    else:
        raise ValueError(f"Only intersect at C or R arcs")
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
                    intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                    if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                else:
                    intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                    if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                    state.increment_c(1) # increment aftwerwards
            case "R":
                intersection_num = intersection_index(num, denom, point, next_point, path_type, parity)
                if not is_ccw: # increment beforehand
                    state.increment_r(1)
                if store_array: store_monomial(monomial_array, state, intersection_num, parity)
                if is_ccw: # increment afterwards
                    state.increment_r(-1)
