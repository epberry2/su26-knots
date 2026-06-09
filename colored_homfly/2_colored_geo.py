from helpers.findloops import findpathwithloops
from helpers.tanglestate import get_state
from functools import cache
from collections import defaultdict
from sympy import symbols, Poly
from itertools import combinations_with_replacement

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

def colored_homfly_geo(num, denom, j):
    path, loops = findpathwithloops(num, denom)
    points_config = get_state(num, denom)[1]
    graph = defaultdict(list)
    parities = {}
    parity = 1
    for curr_point, next_point, path_type in zip(path, path[1:], loops):
        add_arc(graph, num, denom, curr_point, next_point, path_type, points_config, j)
        parities[curr_point] = parity
        parity *= -1
    parities[path[num + denom - 1]] = parity
    values = compute_pair_values(graph)
    min_m = min(m for m, n in values.values())
    min_n = min(n for m, n in values.values())
    
    # TODO: figure out coefficients from values dict
    
    
    

def add_arc(graph, num, denom, p1, p2, path_type, points_config, j = 2): # defaults to 2-colored case
    match path_type:
        case "L":
            m_weight = 2 if points_config[0] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[0] == "X+" else 2
            if p1 < p2:
                for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                    add_relation(graph, (p1,) + indices, (p2,) + indices, m_weight, n_weight) # (p1,) is a 1-tuple
            else:
                for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                    add_relation(graph, (p1,) + indices, (p2,) + indices, -m_weight, -n_weight)
        case "C":
            m_weight = 2 if points_config[1] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[1] == "X+" else 2
            if p1 < p2:
                for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                    add_relation(graph, indices + (p1,), indices + (p2,), -m_weight, n_weight)
            else:
                for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                    add_relation(graph, indices + (p1,), indices + (p2,), m_weight, n_weight)
        case "R":
            m_weight = 2 if points_config[2] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[2] == "X+" else 2
            if p1 < p2:
                for indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - 1):
                    add_relation(graph, indices + (p1,), indices + (p2,), -m_weight, -n_weight)
            else:
                for indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - 1):
                    add_relation(graph, indices + (p1,), indices + (p2,), m_weight, n_weight)
        case "T":
            m_weight = 1 if points_config[1] == "X+" else 0
            n_weight = 1 - (2 * j) if points_config[1] == "X+" else 1
            if p1 < p2:
                for i in range(j):
                    for left_indices in combinations_with_replacement(range(1, num + 1), i):
                        for right_indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - i - 1):
                            add_relation(graph, left_indices + (p1,) + right_indices,
                                            left_indices + (p2,) + right_indices, m_weight, n_weight)
            else:
                for i in range(j):
                    for left_indices in combinations_with_replacement(range(1, num + 1), i):
                        for right_indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - i - 1):
                            add_relation(graph, left_indices + (p1,) + right_indices,
                                            left_indices + (p2,) + right_indices, -m_weight, -n_weight)


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
                inv += 1                    # calculate how many inversions
    return tuple(sorted(t)), 2 * inv        # with inv inversions, n(t) - n(sorted(t)) = 2 * inv

print(two_colored_geo(5, 2))
# should print (a**6*q**2 - a**4*q**8 + 2*a**4*q**4 - a**4 + a**2*q**12 - 2*a**2*q**10 - a**2*q**8 + 4*a**2*q**6 - a**2*q**4 - 2*a**2*q**2 + a**2, -a**3*q**5 + a**3*q**3 + a*q**11 - 2*a*q**9 + 2*a*q**5 - a*q**3, q**12 - q**10 - q**8 + q**6)
# (unless it shouldn't...)