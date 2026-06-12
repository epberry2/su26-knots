from helpers.findloops import findpathwithloops
from jones.JonesBracketGeoHelpers import rotateCCW
from helpers.tanglestate import get_state


def tangle_vectors(u, v):
    "Returns (s_vec, a_vec)"
    path, loops = findpathwithloops(u, v)
    curr_wind_total = 0
    curr_wind_x_plus = 0
    s_vec = [0] * (u + v)
    a_vec = [0] * (u + v)
    x_plus_index = get_state(u, v)[1].index("X+")
    # idx = u + v - 1
    x_omega = path[-1]
    for curr_point, next_point, path_type in reversed(list(zip(path, path[1:], loops))):
        if path_type != "T":
            curr_wind_total -= 1 if rotateCCW(u, v, curr_point, next_point, path_type) else -1
            if (path_type, x_plus_index) in (("L", 0), ("C", 1), ("R", 2)):
                curr_wind_x_plus += 1
        s_vec[curr_point - 1] = curr_wind_total + s_delta_terms(u, curr_point, x_omega)
        a_vec[curr_point - 1] = 2 * curr_wind_x_plus + (s_delta_terms(u, curr_point, x_omega) if x_plus_index == 1 else 0)
        # idx -= 1
        # s_vec[idx] = curr_wind + delta_terms(u, curr_point, x_omega)
        # for his notation
    return s_vec, a_vec

def s_delta_terms(num, i, j):
    return int(i <= num and j > num) - int(i > num and j <= num)

def q_delta_terms(num, points, i, j):
    if points[1] != "X+": return 0
    return s_delta_terms(num, i, j)

# print(tangle_vectors(5, 2))