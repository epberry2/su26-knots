from helpers.findloops import findpathwithloops
from helpers.tanglestate import get_state
from functools import cache
from collections import defaultdict, deque, Counter
from sympy import symbols, Poly
from itertools import combinations_with_replacement
import math

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
    bases_dicts = [defaultdict(int) for _ in range(j+1)] # creates a list of dictionaries mapping powers to coefficients (a,q) -> k
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

    ab = 0
    for key, (m, n) in values.items():
        total_sign = math.prod(parities.get(x, 1) for x in key) # compute sign of tuple
        
        # split tuple into right and left (i_1,...,i_k) -> (i_1,...,i_r), (i_r+1,...,i_k)
        split = find_split(key, num)
        left = key[:split]
        right = key[split:]

        compute_permutation_powers(left, right, (m - min_m,n - min_n), bases_dicts[split], total_sign)  
    homfly = []
    for i in range(j + 1):
        homfly.append(Poly(dict(bases_dicts[i]), (a,q)))
    
    homfly.reverse()
    return homfly
    
def find_split(point, num):
    split = 0
    for i in range(len(point)):
        if point[i] > num:
            break 
        split += 1
    return split

def add_arc(graph, num, denom, p1, p2, path_type, points_config, j = 2): # defaults to 2-colored case
    match path_type:
        case "L": # L loop flips around point on left. m_weight, n_weight are the resp. weights of that point.
            m_weight = 2 if points_config[0] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[0] == "X+" else 2
            for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                add_relation(
                    graph,
                    (p1,) + indices,    # (p1,) is a 1-tuple
                    (p2,) + indices,  
                    m_weight if p1 < p2 else -m_weight, 
                    n_weight if p1 < p2 else -n_weight,
                ) 
        case "C":
            m_weight = 2 if points_config[1] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[1] == "X+" else 2
            for indices in combinations_with_replacement(range(1, num + 1), j - 1):
                add_relation(
                    graph,
                    indices + (p1,),
                    indices + (p2,), 
                    -m_weight if (p1 < p2) ^ (num < denom) else m_weight, 
                    -n_weight if (p1 < p2) ^ (num < denom) else n_weight,
                )
        case "R":
            m_weight = 2 if points_config[2] == "X+" else 0
            n_weight = 2 - (4 * j) if points_config[2] == "X+" else 2
            for indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - 1):
                add_relation(
                    graph,
                    indices + (p1,),
                    indices + (p2,), 
                    -m_weight if p1 < p2 else m_weight, 
                    -n_weight if p1 < p2 else n_weight,
                )
        case "T":
            m_weight = 1 if points_config[1] == "X+" else 0
            n_weight = 1 - (2 * j) if points_config[1] == "X+" else 1
            for i in range(j):
                for left_indices in combinations_with_replacement(range(1, num + 1), i):
                    for right_indices in combinations_with_replacement(range(num + 1, num + denom + 1), j - i - 1):
                        add_relation(
                            graph,
                            left_indices + (p1,) + right_indices,
                            left_indices + (p2,) + right_indices, 
                            m_weight if p1 < p2 else -m_weight,
                            n_weight if p1 < p2 else -n_weight,
                        )

def add_relation(graph, u, v, m_weight, n_weight, j = 2): # u, v = tuple, #(u) - #(v) = weight
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

def sign(n):
    return (n > 0) - (n < 0)

def compute_permutation_powers(t1, t2, start_value, p_counter, sgn):
    # updates the p_counter for all permutations of initial_tuple

    # tuple is of the form (i_1,...,i_k) in sorted order
    # the sign stays the same for every distinct permutation 
 
    k1 = len(t1)
    k2 = len(t2)
    
    # p_counter represents our polynomial as a dictionary {(a_exp, q_exp): coefficient}
    # initial tuple contributes 1 to the start powers
    p_counter[start_value] += sgn 
    
    # Queue for BFS: stores (current_tuple, current_power)
    queue = deque([(t1, t2, start_value[1])])
    
    # Set to keep track of unique permutations we've already visited
    visited = {(t1, t2)}
    
    while queue:
        curr_t1, curr_t2, curr_pow = queue.popleft()
        
        # Try all possible elementary permutations (adjacent swaps) in tuple 1
        for i in range(k1 - 1):
            # Element at i and element at i+1
            val_i, val_j = curr_t1[i], curr_t1[i+1]
            
            # If they are the same, swapping doesn't create a new permutation
            if val_i == val_j:
                continue
                
            # Create the new tuple by swapping them
            next_t1 = list(curr_t1)
            next_t1[i], next_t1[i+1] = next_t1[i+1], next_t1[i]
            next_t1 = tuple(next_t1)
            
            # If we haven't seen this permutation yet, we move to it
            if (next_t1, curr_t2) not in visited:
                visited.add((next_t1, curr_t2))
                
                # Calculate the new power by swapping 2 points at a time
                new_pow = curr_pow - sign(val_i - val_j) * 2
                
                # Add q^new_pow to our polynomial p
                p_counter[(start_value[0], new_pow)] += sgn
                
                queue.append((next_t1, curr_t2, new_pow))

        for i in range(k2 - 1):
            # Element at i and element at i+1
            val_i, val_j = curr_t2[i], curr_t2[i+1]
            
            # If they are the same, swapping doesn't create a new permutation
            if val_i == val_j:
                continue
                
            # Create the new tuple by swapping them
            next_t2 = list(curr_t2)
            next_t2[i], next_t2[i+1] = next_t2[i+1], next_t2[i]
            next_t2 = tuple(next_t2)
            
            # If we haven't seen this permutation yet, we move to it
            if (curr_t1, next_t2) not in visited:
                visited.add((curr_t1, next_t2))
                
                # Calculate the new power by swapping 2 points at a time
                new_pow = curr_pow - sign(val_i - val_j) * 2
                
                # Add q^new_pow to our polynomial p
                p_counter[(start_value[0], new_pow)] += sgn
                
                queue.append((curr_t1, next_t2, new_pow))
                
    return True


@cache
def sorted_and_offset(t: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    inv = 0
    for i in range(len(t)):
        for k in range(i + 1, len(t)):
            if t[i] > t[k]:
                inv += 1                    # calculate how many inversions
    return tuple(sorted(t)), 2 * inv        # with inv inversions, n(t) - n(sorted(t)) = 2 * inv

#print(two_colored_geo(5, 2))
# should print (a**6*q**2 - a**4*q**8 + 2*a**4*q**4 - a**4 + a**2*q**12 - 2*a**2*q**10 - a**2*q**8 + 4*a**2*q**6 - a**2*q**4 - 2*a**2*q**2 + a**2, -a**3*q**5 + a**3*q**3 + a*q**11 - 2*a*q**9 + 2*a*q**5 - a*q**3, q**12 - q**10 - q**8 + q**6)
# (unless it shouldn't...)
print(colored_homfly_geo(5, 2, 7))