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
    graph = defaultdict(list)
    parities = {}
    parity = 1
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        add_arc(graph, num, denom, curr_point, next_point, path_type, points_config)
        parities[curr_point] = parity
        parity *= -1
    parities[path[num + denom - 1]] = parity
    
    values = compute_pair_values(graph)
        
    min_m = min(m for m, n in values.values())
    min_n = min(n for m, n in values.values())
    coeff_2_2 = coeff_2_1 = coeff_2_0 = 0
    for key, (m, n) in values.items():
        a_pow, q_pow = m - min_m, n - min_n
        parity = parities[key[0]] * parities[key[1]]
        if key[1] <= num:
            coeff_2_2 += parity * (a**a_pow) * (q**q_pow)
            if (key[0] != key[1]):
                coeff_2_2 += parity * (a**a_pow) * (q**(q_pow + 2))
        elif key[0] > num:
            coeff_2_0 += parity * (a**a_pow) * (q**q_pow)
            if (key[0] != key[1]):
                coeff_2_0 += parity * (a**a_pow) * (q**(q_pow + 2))
        else:
            coeff_2_1 += parity * (a**a_pow) * (q**q_pow)
    
    return coeff_2_2, coeff_2_1, coeff_2_0
    
    
    

def add_arc(graph, num, denom, p1, p2, path_type, points_config, j = 2): # j = 2
    match path_type:
        case "L":
            m_weight = 2 if points_config[0] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[0] == "X+" else 2
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(graph, (p1, i), (p2, i), m_weight, n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(graph, (p1, i), (p2, i), -m_weight, -n_weight)
        case "C":
            m_weight = 2 if points_config[1] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[1] == "X+" else 2
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(graph, (i, p1), (i, p2), -m_weight, n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(graph, (i, p1), (i, p2), m_weight, n_weight)
        case "R":
            m_weight = 2 if points_config[2] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[2] == "X+" else 2
            if p1 < p2:
                for i in range(num + 1, num + denom + 1):
                    add_relation(graph, (i, p1), (i, p2), -m_weight, -n_weight)
            else:
                for i in range(1, denom + 1):
                    add_relation(graph, (i, p1), (i, p2), m_weight, n_weight)
        case "T":
            m_weight = 1 if points_config[1] == "X+" else 0
            n_weight = 1 - (2 * j) if points_config[1] == "X+" else 1
            if p1 < p2:
                for i in range(1, num + 1):
                    add_relation(graph, (i, p1), (i, p2), m_weight, n_weight)
                for i in range(num + 1, num + denom + 1):
                    add_relation(graph, (p1, i), (p2, i), m_weight, n_weight)
            else:
                for i in range(1, num + 1):
                    add_relation(graph, (i, p1), (i, p2), -m_weight, -n_weight)
                for i in range(num + 1, num + denom + 1):
                    add_relation(graph, (p1, i), (p2, i), -m_weight, -n_weight)

def add_relation(graph, u, v, m_weight, n_weight, j = 2): # u, v = tuple, j = 2. #(u) - #(v) = weight
    u_sorted, u_offset = sorted_and_offset(u)
    v_sorted, v_offset = sorted_and_offset(v)
    new_n_weight = v_offset + n_weight - u_offset
    graph[v_sorted].append((u_sorted, m_weight, new_n_weight))
    graph[u_sorted].append((v_sorted, -m_weight, -new_n_weight))
    

def compute_pair_values(graph):
    """
    graph[cur] contains (nxt, m_diff, n_diff), meaning:
        M(nxt) = M(cur) + m_diff
        N(nxt) = N(cur) + n_diff
    """

    if not graph:
        return {}
    start = next(iter(graph))
    values = {start: (0, 0)}
    stack = [start]
    while stack: # whatever-first search
        cur = stack.pop()
        cur_m, cur_n = values[cur]

        for nxt, m_diff, n_diff in graph[cur]:
            proposed = (cur_m + m_diff, cur_n + n_diff)

            if nxt not in values:
                values[nxt] = proposed
                stack.append(nxt)
            elif values[nxt] != proposed:
                raise ValueError("Inconsistent relations")

    if len(values) != len(graph):
        raise ValueError("Graph is not connected")

    return values


@cache
def sorted_and_offset(t: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    inv = 0
    for i in range(len(t)):
        for k in range(i + 1, len(t)):
            if t[i] > t[k]:
                inv += 1
    return tuple(sorted(t)), 2 * inv

print(two_colored_geo(5, 2))