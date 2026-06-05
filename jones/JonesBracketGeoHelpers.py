from dataclasses import dataclass, field, InitVar

@dataclass
class tangleRotateState:
        # can pass into constructor
    poly_type: str
    on_left: bool
    points_state: InitVar[tuple[str, str, str] | None] = None
        # can't pass into constructor
    min_m_power: int = field(init = False, default = 0)
    min_n_power: int = field(init = False, default = 0)
    quarter_loops: list[int] = field(init = False, default_factory = lambda: [0, 0, 0])
    total_loops: int = field(init = False, default = 0)
    x_plus_idx: int = field(init = False, default = 0)
    x_minus_idx: int = field(init = False, default = 0)
    y_idx: int = field(init = False, default = 0)
    
    def __post_init__(self, points):
        if points is not None:
            self.x_plus_idx = points.index("X+")
            self.x_minus_idx = points.index("X-")
            self.y_idx = points.index("Y")
        else: 
            assert self.poly_type != "H", "For HOMFLY bracked must initialize location of points"

    def traverse(self, is_ccw, path_type):
        match path_type:
            case "T":
                self.quarter_loops[1] -= 1 if self.on_left else -1
                self.total_loops += 1
                self.on_left = not self.on_left
            case _:
                self.quarter_loops[path_type_index(path_type)] += 2 if is_ccw else -2
                self.total_loops += 1       
                
    def powers(self):
        match self.poly_type:
            case "J":
                parity = -1 if self.total_loops % 2 == 0 else 1
                m_power = 0
                n_power = self.quarter_loops[0] + self.quarter_loops[1] + self.quarter_loops[2]
                self.min_n_power = min(self.min_n_power, n_power)
                return parity, m_power, n_power
            case "H":
                parity = -1 if self.total_loops % 2 == 0 else 1
                m_power = self.quarter_loops[self.x_plus_idx]
                n_power = (self.quarter_loops[self.x_minus_idx] + self.quarter_loops[self.y_idx]
                                - self.quarter_loops[self.x_plus_idx])
                self.min_m_power = min(self.min_m_power, m_power)
                self.min_n_power = min(self.min_n_power, n_power)
                return parity, m_power, n_power
    
    

def path_type_index(path_type):
    match path_type:
        case "L":
            return 0
        case "C":
            return 1
        case "R":
            return 2
        case "T":
            return 3


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


