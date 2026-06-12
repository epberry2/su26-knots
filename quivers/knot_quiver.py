from knot_vectors import knot_vectors
from helpers.findloops import findpathwithloops
from jones.JonesBracketGeoHelpers import rotateCCW
from jones.GeoJonesKnotHelpers import intersection_index
from helpers.tanglestate import get_state
from dataclasses import dataclass
from typing import List, Set




def winds(u, v):
    path, loops = findpathwithloops(u, v)
    
    states = [[0] * u for _ in range(u)]
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        "Unfinished"
        
        
@dataclass
class winding_tracker:
    num: int
    denom: int
    diag_winds: List[List[int]]
    x_plus_states: List[int]
    orientation: str
    has_hit: Set[int] = set()
    curr_x_plus: int = 0
    def __post_init__(self):
        self.diag_winds = [[0] * self.num for _ in range(self.num)]
        self.x_plus_states = [0] * self.num
        self.orientation = get_state(self.num, self.denom)[0]
        if self.orientation == "RI":
            raise ValueError("RI orientation")
        
    def step(self, curr_point, next_point, path_type):
        is_ccw = rotateCCW(self.num, self.denom, curr_point, next_point, path_type)
        match path_type:
            case "R":
                intersect = intersection_index(self.num, self.denom, curr_point, next_point, path_type)
                if self.orientation == "UP":
                    self.curr_x_plus -= 0 if is_ccw else -2
                self.x_plus_states[intersect] = self.curr_x_plus
                self.has_hit.add(intersect)
                if self.orientation == "UP":
                    self.curr_x_plus -= 2 if is_ccw else 0
            case "C":
                intersect = intersection_index(self.num, self.denom, curr_point, next_point, path_type)
                self.curr_x_plus -= 1 if is_ccw else -1
                for j in self.has_hit:
                    for i in range(intersect):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] -= 2 if is_ccw else -2
                self.x_plus_states[intersect] = self.curr_x_plus
                self.has_hit.add(intersect)
                self.curr_x_plus -= 1 if is_ccw else -1
            case "T":
                if self.orientation == "OP":
                    self.curr_x_plus -= 1 if curr_point > self.num else -1
                for j in self.has_hit:
                    for i in range(self.num):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] -= 1 if curr_point > self.num else -1
    def wind_diag(self, j, i):
        if j > i:
            raise ValueError("j must be less than i")
        return self.diag_winds[j][i]
    def wind_x_plus(self, j, i):
        if j > i:
            raise ValueError("j must be less than i")
        return (self.x_plus_states[i] - self.x_plus_states[j]) // 2