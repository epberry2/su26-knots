import quivers.mu_funcs as mu_func
from quivers.knot_vectors import knot_vectors
from helpers.findloops import findpathwithloops
from jones.JonesBracketGeoHelpers import rotateCCW
from jones.GeoJonesKnotHelpers import intersection_index
from helpers.tanglestate import get_state
from dataclasses import dataclass, field
from typing import List, Set, Tuple
import sympy as sp


def intersection_path(u, v):
    winder = winding_tracker(u, v)
    return [idx + 1 for idx in winder.intersection_path]

def diag_wind_matrix(u, v):
    winder = winding_tracker(u, v)
    _, _, q_diag = winder.homfly_vectors()
    path = winder.intersection_path
    W = sp.Matrix.zeros(u, u)
    for idx_i in range(u):
        i = path[idx_i]
        for idx_j in range(idx_i):
            j = path[idx_j]
            W[j, i] = winder.wind_diag(j, i)
            W[i, j] = q_diag[i] - q_diag[j] + winder.wind_diag(j, i) - 4 * winder.wind_x_plus(j, i)
    return W

def x_plus_wind_matrix(u, v):
    winder = winding_tracker(u, v)
    _, _, q_diag = winder.homfly_vectors()
    path = winder.intersection_path
    X = sp.Matrix.zeros(u, u)
    for idx_i in range(u):
        i = path[idx_i]
        for idx_j in range(idx_i):
            j = path[idx_j]
            X[j, i] = winder.wind_x_plus(j, i)
            X[i, j] = -X[j, i]
    return X

def colored_homfly_vectors_and_quiver(
    u: int, v: int) -> Tuple[List[int], List[int], sp.Matrix]:
    """Calculates the Homfly quiver data for the knot K_u/v
    
    Input: u, v
    
    Output: (S, A, Q) of colored Homfly polynomial for the knot K_{u/v}
    """
    if u < v:
        raise ValueError("U")
    winder = winding_tracker(u, v)
    Q = sp.Matrix.zeros(u, u)
    s_vec, a_vec, q_diag = winder.homfly_vectors()
    path = winder.intersection_path
    for idx_i in range(u):
        i = path[idx_i]
        Q[i, i] = q_diag[i]
        for idx_j in range(idx_i):
            j = path[idx_j]
            Q[i, j] = Q[i, i] + winder.wind_diag(j, i) - 2 * winder.wind_x_plus(j, i)
            Q[j, i] = Q[i, j] # = Q[j, j] + winddiag(ij) - 2winderwindxplus(ji)
    return s_vec, a_vec, Q

def colored_jones_vector_and_quiver(
    u, v) -> Tuple[List[int], sp.Matrix]:
    """Calculates the Jones quiver data for the knot K_u/v

    Input: u, v
    
    Output: (H, Q') of colored Jones polynomial for the knot K_{u/v}
    """
    winder = winding_tracker(u, v)
    Q = sp.Matrix.zeros(u, u)
    h_vec, q_diag = winder.colored_jones_vectors()
    path = winder.intersection_path
    for idx_i in range(u):
        i = path[idx_i]
        Q[i, i] = q_diag[i]
        for idx_j in range(idx_i):
            j = path[idx_j]
            Q[i, j] = Q[i, i] + 2 * winder.wind_x_plus(j, i) - winder.wind_diag(j, i)
            Q[j, i] = Q[i, j]
    return h_vec, Q
        
        
@dataclass
class winding_tracker:
    """Tracks the winding numbers for loops on alpha_u/v bar"""
    num: int
    denom: int
    diag_winds: List[List[int]] = field(init = False)
    writhe_states: List[Tuple[int]] = field(init = False)
    curr_writhes: List[int] = field(init = False)
    intersection_path: List[int] = field(default_factory=list)
    has_hit: Set[int] = field(default_factory=set)
    permute: Tuple[int] = field(init = False)
    writhes_k_w: Tuple[int] = (0, 0, 0)
    
    def __post_init__(self):
        if self.num % 2 == 0:
            raise ValueError("num must be odd")
        self.diag_winds = [[0] * self.num for _ in range(self.num)]
        self.writhe_states = [[0, 0, 0] for _ in range(self.num)]
        self.curr_writhes = [0, 0, 0]
        orientation = get_state(self.num, self.denom)[0]
        if orientation == "RI":
            raise ValueError("RI orientation")
        self.permute = (1, 2 if orientation == "UP" else 0, 0 if orientation == "UP" else 2)
        self._trace_path()
        
    def _step(self, curr_point, next_point, path_type):
        is_cw = not rotateCCW(self.num, self.denom, curr_point, next_point, path_type)
        intersect = intersection_index(self.num, self.denom, curr_point, next_point, path_type)
        match path_type:
            case "R":
                if is_cw: # CW loop goes around right point before hitting beta
                    self.curr_writhes[self.permute[2]] += 2
                # R loop passes under all the points to the right of intersect
                for j in self.has_hit:
                    for i in range(intersect + 1, self.num):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] += 1 if is_cw else -1

                self.writhe_states[intersect] = tuple(self.curr_writhes)
                self.has_hit.add(intersect)
                self.intersection_path.append(intersect)
                if not is_cw: # R loop passes underneath rightward points after hitting intersection
                    for i in range(intersect + 1, self.num):
                        if i not in self.has_hit:
                            self.diag_winds[intersect][i] -= 1
                if not is_cw:  # CCW loop goes around right point after hitting beta 
                    self.curr_writhes[self.permute[2]] -= 2
                    
            case "C":
                self.curr_writhes[self.permute[1]] += 1 if is_cw else -1
                for j in self.has_hit:
                    for i in range(intersect):
                        if i not in self.has_hit:
                            self.diag_winds[j][i] += 2 if is_cw else -2
                    if not is_cw: # If γ_{j,i} ends with a CCW center loop, the writhe ends up decreasing by 1 due to the square.
                        self.diag_winds[j][intersect] -= 1
                self.writhe_states[intersect] = tuple(self.curr_writhes)
                self.has_hit.add(intersect)
                self.intersection_path.append(intersect)
                for i in range(intersect):
                    if i not in self.has_hit:
                        self.diag_winds[intersect][i] += 1 if is_cw else -1 # loops starting at intersection point only get half the writhe
                self.curr_writhes[self.permute[1]] += 1 if is_cw else -1

            case "T":
                self.curr_writhes[self.permute[1]] += 1 if curr_point <= self.num else -1
                for j in self.has_hit:
                    for i in range(intersect): # If the T arc is shortened, it only affects writhe where i is less than the intersection point to the right
                        if i not in self.has_hit:
                            self.diag_winds[j][i] -= 1 if curr_point > self.num else -1
            case "L":
                self.curr_writhes[self.permute[0]] += 2 if is_cw else -2
                
                                
    def _turn_around_CW(self): # Will always turn CW because we walk with RH to wall
        if self.permute[2] == 0: # Path end (X+) is rightmost point
            self.curr_writhes[self.permute[2]] += 2
            self.writhes_k_w = self.writhe_states[self.num - 1] = tuple(self.curr_writhes)
            self.has_hit.add(self.num - 1)
            self.intersection_path.append(self.num - 1)
        else: # Path end (X+) is central point Z
            self.curr_writhes[self.permute[1]] += 1
            self.writhes_k_w = self.writhe_states[0] = tuple(self.curr_writhes)
            self.has_hit.add(0)
            self.intersection_path.append(0)
            self.curr_writhes[self.permute[1]] += 1

       
    def _trace_path(self):
        path, loops = findpathwithloops(self.num, self.denom)
        for curr_point, next_point, path_type in zip(path, path[1:], loops):
            self._step(curr_point, next_point, path_type)
        self._turn_around_CW()
        for curr_point, next_point, path_type in zip(reversed(path), reversed(path[:-1]), reversed(loops)):
            self._step(curr_point, next_point, path_type)
    
    def wind_diag(self, j, i):
        #j MUST come before i in intersection_path
        return self.diag_winds[j][i]
    def wind_x_plus(self, j, i):
        "j MUST come before i in intersection_path"
        return (self.writhe_states[i][0] - self.writhe_states[j][0]) // 2
    def homfly_vectors(self):
        "return s_vec, a_vec, q_diag"
        mu1, mu2, mu3 = mu_func.mus(self.num, self.denom)
        s_vec = [0] * self.num
        a_vec = [0] * self.num
        q_diag = [0] * self.num
        for i in range(self.num): # u = len(writhes) = len(s_vec) = len(a_vec) = len(q_diag)
            s_vec[i] = mu1 + (
                    self.writhes_k_w[0] - self.writhe_states[i][0]
                    + self.writhes_k_w[1] - self.writhe_states[i][1]
                    + self.writhes_k_w[2] - self.writhe_states[i][2]
                ) // 2
            a_vec[i] = mu2 + self.writhes_k_w[0] - self.writhe_states[i][0]
            q_diag[i] = mu3 + (
                    -3 * (self.writhes_k_w[0] - self.writhe_states[i][0])
                    + self.writhes_k_w[1] - self.writhe_states[i][1]
                    + self.writhes_k_w[2] - self.writhe_states[i][2]
                ) // 2
        return s_vec, a_vec, q_diag
    
    def colored_jones_vectors(self):
        mu1, mu2, mu3 = mu_func.mus(self.num, self.denom)
        h_vec = [0] * self.num
        q_diag = [0] * self.num
        for i in range(self.num):
            wind_term = (
                3 * (self.writhes_k_w[0] - self.writhe_states[i][0])
                - self.writhes_k_w[1] + self.writhe_states[i][1]
                - self.writhes_k_w[2] + self.writhe_states[i][2]
            ) // 2
            h_vec[i] = wind_term + 2 * mu2 - mu1
            q_diag[i] = wind_term - mu3
        return h_vec, q_diag

