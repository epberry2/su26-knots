from helpers.findloops import findpathwithloops
from collections import defaultdict
from jones.JonesBracketGeoHelpers import rotateCCW
from jones.GeoJonesKnotHelpers import intersection_index
from helpers.tanglestate import get_state
from quivers.mu_funcs import mus

def knot_vectors(u, v):
    writhes = writhe_list(u, v)
    mu1, mu2, mu3 = mus(u, v)
    
    s_vec = [0] * u
    a_vec = [0] * u
    q_diag = [0] * u
    for i in range(u): # u = len(writhes) = len(s_vec) = len(a_vec) = len(q_diag)
        s_vec[i] = mu1 + (writhes[i][0] + writhes[i][1] + writhes[i][2]) // 2
        a_vec[i] = mu2 + writhes[i][0]
        q_diag[i] = mu3 + (-3 * writhes[i][0] + writhes[i][1] + writhes[i][2]) // 2
    return s_vec, a_vec, q_diag

def writhe_list(num, denom):
    "writhe_list[i] = [w_{X+}(gamma_i), w_{X-}(gamma_i), w_Y(gamma_i)]"
    orientation, _ = get_state(num, denom)
    path, loops = findpathwithloops(num, denom)
    writhe_list = [[0, 0, 0] for _ in range(num)]
    writhe_counter = [0, 0, 0]
    permute = (1, 2 if orientation == "UP" else 0, 0 if orientation == "UP" else 2)
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter, permute)
    match orientation:
        case "UP": # X+ on right
            writhe_counter[permute[2]] += 2
            writhe_list[num - 1] = writhe_counter.copy()
            offset = writhe_counter.copy()
        case "OP": # X+ in middle
            writhe_counter[permute[1]] += 1
            writhe_list[0] = writhe_counter.copy()
            offset = writhe_counter.copy()
            writhe_counter[permute[1]] += 1
        case "RI":
            raise ValueError("RI orientation")
    for curr_point, next_point, path_type in zip(reversed(path), reversed(path[:-1]), reversed(loops)):
        add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter, permute)
    normalized_list = [[off - writhe for writhe, off in zip(intersection_point, offset)] for intersection_point in writhe_list]
    return normalized_list
    
def add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter, permute):
    match path_type:
        case "T":
            writhe_counter[permute[1]] -= 1 if curr_point > num else -1
        case "L":
            writhe_counter[permute[0]] -= 2 if rotateCCW(num, denom, curr_point, next_point, path_type) else -2
        case "C":
            intersection_point = intersection_index(num, denom, curr_point, next_point, path_type)
            writhe_counter[permute[1]] -= 1 if rotateCCW(num, denom, curr_point, next_point, path_type) else -1
            writhe_list[intersection_point] = writhe_counter.copy()
            writhe_counter[permute[1]] -= 1 if rotateCCW(num, denom, curr_point, next_point, path_type) else -1
        case "R":
            intersection_point = intersection_index(num, denom, curr_point, next_point, path_type)
            writhe_counter[permute[2]] -= 0 if rotateCCW(num, denom, curr_point, next_point, path_type) else -2
            writhe_list[intersection_point] = writhe_counter.copy()
            writhe_counter[permute[2]] -= 2 if rotateCCW(num, denom, curr_point, next_point, path_type) else 0

# print(knot_vectors(7, 5))