from dataclasses import dataclass, field

@dataclass
class tangleRotateState:
    poly_type: str
    on_left: bool
    quarter_loops: list[int] = field(default_factory=lambda: [0, 0, 0])
    total_loops: int = 0
    min_power: int = 0

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
                self.min_power = min(self.min_power, n_power)
                return parity, m_power, n_power
            case _: # e.g. H for HOMFLY
                return 0, 0, 0
                
    
    

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

def coeff(pathType):
    if pathType == "L" or pathType == "R":
        return 1
    return -1


