from helpers.findloops import findpathwithloops
from collections import defaultdict
from jones.JonesBracketGeoHelpers import rotateCCW
from jones.GeoJonesKnotHelpers import intersection_index
from helpers.tanglestate import get_state


def writhe_list(num, denom):
    path, loops = findpathwithloops(num, denom)
    orientation, _ = get_state(num, denom)
    writhe_list = [[0, 0, 0] for _ in range(num)]
    writhe_counter = [0, 0, 0]
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter)
    match orientation:
        case "UP": # X+ on right
            writhe_counter[2] += 2
            writhe_list[num - 1] = writhe_counter.copy()
            offset = writhe_counter.copy()
        case "OP": # X+ in middle
            writhe_counter[1] += 1
            writhe_list[0] = writhe_counter.copy()
            offset = writhe_counter.copy()
            writhe_counter[1] += 1
        case "RI":
            raise ValueError("RI orientation")
    for curr_point, next_point, path_type in zip(reversed(path), reversed(path[:-1]), reversed(loops)):
        add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter)
    normalized_list = [[off - writhe for writhe, off in zip(intersection_point, offset)] for intersection_point in writhe_list]
    return normalized_list
    
def add_arc(writhe_list, num, denom, curr_point, next_point, path_type, writhe_counter):
    match path_type:
        case "L":
            writhe_counter[0] -= 2 if rotateCCW(num, denom, curr_point, next_point, path_type) else -2
        case "C":
            intersection_point = intersection_index(num, denom, curr_point, next_point, path_type)
            writhe_counter[1] -= 1 if rotateCCW(num, denom, curr_point, next_point, path_type) else -1
            writhe_list[intersection_point] = writhe_counter.copy()
            writhe_counter[1] -= 1 if rotateCCW(num, denom, curr_point, next_point, path_type) else -1
        case "R":
            intersection_point = intersection_index(num, denom, curr_point, next_point, path_type)
            writhe_counter[2] -= 0 if rotateCCW(num, denom, curr_point, next_point, path_type) else -2
            writhe_list[intersection_point] = writhe_counter.copy()
            writhe_counter[2] -= 2 if rotateCCW(num, denom, curr_point, next_point, path_type) else 0

print(writhe_list(5, 2))