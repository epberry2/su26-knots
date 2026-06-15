from quivers.knot_vectors import knot_vectors
from helpers.findloops import findpathwithloops, findpath
from jones.JonesBracketGeoHelpers import rotateCCW
from jones.GeoJonesKnotHelpers import intersection_index
from helpers.tanglestate import get_state
from dataclasses import dataclass, field
from typing import List, Set
import sympy as sp


def knot_quiver(u, v):
    winder = winding_tracker(u, v)
    winder.trace_path()
    Q = sp.Matrix.zeros(u, u)
    q_diag = knot_vectors(u, v)[2]
    path = winder.intersection_path
    for idx_i in range(u):
        i = path[idx_i]
        Q[i, i] = q_diag[i]
        for idx_j in range(idx_i):
            j = path[idx_j]
            Q[i, j] = Q[i, i] + winder.wind_diag(j, i) - 2 * winder.wind_x_plus(j, i)
            Q[j, i] = Q[i, j]
    return Q
        
@dataclass
class winding_tracker:
    num: int
    denom: int
    diag_winds: List[List[int]] = field(init = False)
    x_plus_states: List[int] = field(init = False)
    orientation: str = field(init = False)
    intersection_path: List[int] = field(default_factory=list)
    has_hit: Set[int] = field(default_factory=set)
    curr_x_plus: int = 0
    
    def __post_init__(self):
        self.diag_winds = [[0] * self.num for _ in range(self.num)]
        self.x_plus_states = [0] * self.num
        self.orientation = get_state(self.num, self.denom)[0]
        if self.orientation == "RI":
            raise ValueError("RI orientation")
        
    def step(self, curr_point, next_point, path_type):
        is_cw = not rotateCCW(self.num, self.denom, curr_point, next_point, path_type)
        intersect = intersection_index(self.num, self.denom, curr_point, next_point, path_type)
        print(curr_point, next_point, path_type, is_cw, intersect)
        match path_type:
            case "R":
                if self.orientation == "UP" and is_cw: # CW loop goes around right point before hitting beta
                    self.curr_x_plus += 2 

                # R loop passes under all the points to the right of intersect
                for j in self.has_hit:
                    for i in range(intersect + 1, self.num):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] += 1 if is_cw else -1
                self.x_plus_states[intersect] = self.curr_x_plus
                self.has_hit.add(intersect)
                self.intersection_path.append(intersect)
                if not is_cw: # R loop passes underneath rightward points after hitting intersection
                    for i in range(intersect + 1, self.num):
                        if i not in self.has_hit:
                            self.diag_winds[intersect][i] -= 1
                if self.orientation == "UP" and not is_cw: # CCW loop goes around right point after hitting beta 
                    self.curr_x_plus -= 2
            case "C":
                if self.orientation == "OP":
                    self.curr_x_plus += 1 if is_cw else -1
                for j in self.has_hit:
                    for i in range(intersect):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] += 2 if is_cw else -2
                    if not is_cw: # If γ_{j,i} ends with a CCW center loop, the writhe ends up decreasing by 1 due to the square.
                        self.diag_winds[j][intersect] -= 1
                self.x_plus_states[intersect] = self.curr_x_plus
                self.has_hit.add(intersect)
                self.intersection_path.append(intersect)
                for i in range(intersect):
                    if i not in self.has_hit:
                        self.diag_winds[intersect][i] += 1 if is_cw else -1 # loops starting at intersection point only get half the writhe
                if self.orientation == "OP":
                    self.curr_x_plus += 1 if is_cw else -1
            case "T":
                if self.orientation == "OP":
                    self.curr_x_plus -= 1 if curr_point > self.num else -1
                for j in self.has_hit:
                    for i in range(intersect): # If the T arc is shortened, it only affects writhe where i is less than the intersection point to the right
                        if i not in self.has_hit:
                            self.diag_winds[j][i] -= 1 if curr_point > self.num else -1
                                
    def turn_around_CW(self): # Will always turn CW because we walk with RH to wall
        match self.orientation:
            case "UP": # Path ends on rightmost point
                self.curr_x_plus += 2
                self.x_plus_states[self.num - 1] = self.curr_x_plus
                self.has_hit.add(self.num - 1)
                self.intersection_path.append(self.num - 1)
            case "OP": # Path ends on central point Z
                self.curr_x_plus += 1
                self.x_plus_states[0] = self.curr_x_plus
                self.has_hit.add(0)
                self.intersection_path.append(0)
                self.curr_x_plus += 1
                
    def trace_path(self):
        path, loops = findpathwithloops(self.num, self.denom)
        for curr_point, next_point, path_type in zip(path, path[1:], loops):
            self.step(curr_point, next_point, path_type)
        self.turn_around_CW()
        for curr_point, next_point, path_type in zip(reversed(path), reversed(path[:-1]), reversed(loops)):
            self.step(curr_point, next_point, path_type)
    
    def wind_diag(self, j, i):
        "j must come before i in intersection_path"
        return self.diag_winds[j][i]
    def wind_x_plus(self, j, i):
        return (self.x_plus_states[i] - self.x_plus_states[j]) // 2


# knot_quiver(11, 8)
# print(knot_vectors(11, 8))