from helpers.findloops import findpath, findpathwithloops_sector
from quivers.tangle_vectors import tangle_vectors, q_delta_terms
from helpers.tanglestate import get_state
from jones.GeoJonesKnotHelpers import *
from jones.JonesBracketGeoHelpers import rotateCCW
import sympy as sp

def quiver(num, denom):
    Q = sp.Matrix.zeros(num + denom, num + denom)
    points_ = get_state(num, denom)[1]
    path = findpath(num, denom)
    x_w = path[-1]
    x_plus_idx = points_.index("X+")
    for a in range(1, num + denom + 1):
        path_a, loops_a = findpathwithloops_sector(num, denom, a, x_w)
        state = LoopWalkState(poly_type = "", points = points_)
        tracePath(num, denom, path_a, loops_a, state, store_array = False)
        Q[a-1, a-1] = state.quiver_diag_term() + -2 * q_delta_terms(num, points_, a, x_w)
        for b in range(1, a):
            entry = Q[a-1, a-1] + writhe_of_path(num, denom, b, a)
            a_b_path, a_b_loops = findpathwithloops_sector(num, denom, b, a)
            wind_x_plus = 0
            for curr_point, next_point, path_type in zip(a_b_path, a_b_path[1:], a_b_loops):
                if (x_plus_idx, path_type) in ((0, "L"), (1, "C"), (2, "R")):
                    wind_x_plus += 1 if rotateCCW(num, denom, curr_point, next_point, path_type) else -1
            entry -= 2 * wind_x_plus
            entry += q_delta_terms(num, points_, a, x_w)
            Q[a-1, b-1] = entry
            Q[b-1, a-1] = entry
    return Q
                


def writhe_of_path(num, denom, j, i):
    "assume i < j and num >= denom"
    if i == j: return 0
    path, loops = findpathwithloops_sector(num, denom, j, i)
    i_left = i <= num
    j_left = j <= num
    writhe = 0
    on_left = j_left
    

    if i_left and j_left:
        if loops[0] == "L":
            writhe += -1 if j > i else 0
        elif loops[0] == "C":
            writhe += 0 if j > i else -1
        elif loops[0] == "T":
            writhe += 0 if j > i else -1
    elif (not i_left) and (not j_left):
        if loops[0] == "R":
            writhe += -1 if i > j else 0
        elif loops[0] == "T":
            writhe += 0 if i > j else -1
            
    for next_point, path_type in zip(path[1:-1], loops[0:-1]):
        if i_left:
            match path_type:
                case "L":
                    writhe += 1 if next_point > i else -1
                case "C":
                    writhe += 1 if next_point < i else -1
                case "T":
                    if not on_left:
                        writhe += 1 if next_point < i else -1
                    on_left = not on_left
        else:
            match path_type:
                case "R":
                    writhe += 1 if next_point < i else -1
                case "T":
                    if on_left:
                        writhe += -1 if next_point < i else 1
                    on_left = not on_left
    return writhe

# print(quiver(3, 1))