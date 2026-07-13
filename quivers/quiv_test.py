from quivers.knot_quiver import colored_jones_vector_and_quiver, colored_homfly_vectors_and_quiver, intersection_path, winding_tracker, diag_wind_matrix
from helpers.quantum_algebra import QuantumCombinatorics
import numpy as np
import sympy as sp
from bisect import insort
from quivers.evaluate_quiver import get_tuples
from helpers.arc_directions import arc_directions

def print_min_vecs(u, v, n, m=1):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v) # returns Tuple[List[int], List[int], Matrix]
    S = np.array(S)
    A = np.array(A)
    for j in range(1, n + 1):
        combos = get_tuples(u, j) # returns list[tuple] of vectors of length u which partition j
        combo_list = list(combos)
        min_d = []
        count = 0
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(S, d)        
            p2 = np.dot(A, d)
            p3 = np.einsum('i,ij,j', d, Q, d)
            term = p3
            if not flag or p1 - 2 * p2 + p3 < min_d[m - 1][1]:
                insort(min_d, (d.tolist(), term, p2, int((-1) ** (p1 % 2))), key=lambda x: x[1])
                if flag:
                    min_d.pop()
                else:
                    count += 1
                    if count >= m:
                        flag = True
        print([(vec, val, int(val2), sgn) for vec, val, val2, sgn in min_d])
        print("-------------------------")
    print(S)
    print(A)
    sp.pprint(Q)

def print_max_vecs(u, v, n, m=1):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    print(S)
    sp.pprint(Q)
    S = np.array(S)
    A = np.array(A)
    for j in range(1, n + 1):
        qc = QuantumCombinatorics()
        combos = get_tuples(u, j)
        combo_list = list(combos)
        max_d = [(np.zeros(u, dtype=int), 0, 0) for _ in range(m)]
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(S, d)        
            p2 = np.dot(A, d)
            p3 = np.einsum('i,ij,j', d, Q, d)
            multinom = qc.get_multinomial(d)
            for i, _ in enumerate(multinom):
                term = 2 * i + p3
                # if not flag or p1 + p3 + 2 * i > max_d[m - 1][1]:
                #     insort(max_d, (d, p1 + p3 + 2 * i, int((-1) ** (p1 % 2))), key=lambda x: -(x[1]))
                #     max_d.pop()
                if not flag or term > max_d[m - 1][1]:
                    insort(max_d, (d.tolist(), term, int((-1) ** (p1 % 2))), key=lambda x: -(x[1]))
                    max_d.pop()
                flag = True
        print([(vec, val, sgn) for vec, val, sgn in max_d])
        print("-------------------------")

def print_max_jones_vecs(u, v, n, m=1):
    H, Q = colored_jones_vector_and_quiver(u, v)
    H = np.array(H)
    for j in range(1, n + 1):
        qc = QuantumCombinatorics()
        combos = get_tuples(u, j)
        combo_list = list(combos)
        max_d = []
        count = 0
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(H, d)        
            p3 = np.einsum('i,ij,j', d, Q, d)
            multinom = qc.get_multinomial(d)
            for i, coeff in enumerate(multinom):
                term = p1 + p3 - 2 * i
                if (not flag) or term > max_d[len(max_d) - 1][1]:
                    insort(max_d, (d.tolist(), term, coeff * int((-1) ** (p1 % 2))), key=lambda x: -(x[1]))
                    if flag:
                        max_d.pop()
                    else:
                        count += 1
                        if count >= m:
                            flag = True
        print([(vec, val, sgn) for vec, val, sgn in max_d])
        print("-------------------------")
    print(H)
    sp.pprint(Q)
    
def print_intersection_paths(max_num):
    for u in range(1, max_num):
        for v in range(1, u):
            try:
                path = intersection_path(u, v)
                print(f"path {u}/{v}:")
                print(path)
            except:
                pass  
            
# print(arc_directions(17, 5))
# sp.pprint(colored_homfly_vectors_and_quiver(7, 3))
# print(intersection_path(11, 4))
# sp.pprint(diag_wind_matrix(11, 4))
