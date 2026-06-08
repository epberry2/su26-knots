from helpers.findloops import findpathwithloops
from helpers.tanglestate import get_state
from functools import cache
from collections import defaultdict
from sympy import symbols, Poly

q = symbols('q')
a = symbols('a')


def two_colored_geo(num, denom):
    path, loops = findpathwithloops(num, denom)
    points_config = get_state(num, denom)[1]
    m_graph = defaultdict(list)
    n_graph = defaultdict(list)
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        add_arc(m_graph, n_graph, num, denom, curr_point, next_point, path_type, points_config)
    coeff_2_0 = Poly(0, a, q)
        
    

def add_arc(m_graph, n_graph, num, denom, p1, p2, path_type, points_config, j = 2): # j = 2
    match path_type:
        case "L":
            m_weight = 2 if points_config[0] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[0] == "X+" else 2
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(m_graph, (p1, i), (p2, i), m_weight)
                    add_relation(n_graph, (p1, i), (p2, i), n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(m_graph, (p1, i), (p2, i), -m_weight)
                    add_relation(n_graph, (p1, i), (p2, i), -n_weight)
        case "C":
            m_weight = 2 if points_config[1] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[1] == "X+" else 2
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(m_graph, (i, p1), (i, p2), -m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), -n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(m_graph, (i, p1), (i, p2), m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), n_weight)
        case "R":
            m_weight = 2 if points_config[2] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[2] == "X+" else 2
            if p1 < p2:
                for i in range(num + 1, denom + 1):
                    add_relation(m_graph, (i, p1), (i, p2), -m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), -n_weight)
            else:
                for i in range(1, denom + 1):
                    add_relation(m_graph, (i, p1), (i, p2), m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), n_weight)
        case "T":
            m_weight = 1 if points_config[1] == "X+" else 0
            n_weight = 1 - (2 * j) if points_config[1] == "X+" else 1
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(m_graph, (i, p1), (i, p2), m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), n_weight)
                for i in range(num + 1, denom + 1):
                    add_relation(m_graph, (p1, i), (p2, i), m_weight)
                    add_relation(m_graph, (p1, i), (p2, i), n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(m_graph, (i, p1), (i, p2), -m_weight)
                    add_relation(n_graph, (i, p1), (i, p2), -n_weight)
                for i in range(num + 1, denom + 1):
                    add_relation(m_graph, (p1, i), (p2, i), -m_weight)
                    add_relation(n_graph, (p1, i), (p2, i), -n_weight)

def add_relation(graph, u, v, weight, j = 2): # u, v = tuple, j = 2. #(u) - #(v) = weight
    u_sorted, u_offset = sorted_and_offset(u)
    v_sorted, v_offset = sorted_and_offset(v)
    new_weight = v_offset + weight - u_offset
    graph[v_sorted].append((u_sorted, new_weight))
    graph[u_sorted].append((v_sorted, -new_weight))
    

@cache
def sorted_and_offset(t: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    inv = 0
    for i in range(len(t)):
        for k in range(i + 1, len(t)):
            if t[i] > t[k]:
                inv += 1
    return tuple(sorted(t)), 2 * inv