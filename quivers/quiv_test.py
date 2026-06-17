from quivers.knot_quiver import colored_jones_vector_and_quiver, colored_homfly_vectors_and_quiver
from helpers.quantum_algebra import QuantumCombinatorics
import numpy as np
import sympy as sp
from bisect import insort
from quivers.evaluate_quiver import get_tuples

def print_min_vecs(u, v, n, m=1):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v) # returns Tuple[List[int], List[int], Matrix]
    print(S)
    sp.pprint(Q)
    S = np.array(S)
    A = np.array(A)
    for j in range(1, n + 1):
        combos = get_tuples(u, j) # returns list[tuple] of vectors of length u which partition j
        combo_list = list(combos)
        min_d = [(np.zeros(u, dtype=int), 0) for _ in range(m)]
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(S, d)        
            # p2 = np.dot(A, d)
            p3 = np.einsum('i,ij,j', d, Q, d)
            if not flag or p1 + p3 < min_d[m - 1][1]:
                insort(min_d, (d, p1 + p3), key=lambda x: x[1])
                min_d.pop()
            flag = True
        print([(vec.tolist(), val) for vec, val in min_d])
        print("-------------------------")
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
        max_pow = 0
        max_d = (0) * u
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(S, d)        
            # p2 = np.dot(A, d)
            p3 = np.einsum('i,ij,j', d, Q, d)
            multinom = qc.get_multinomial(d)
            for i, _ in enumerate(multinom):
                if not flag or p1 + p3 + 2 * i > max_pow:
                    max_pow = p1 + p3 + 2 * i
                    max_d = d
                flag = True
        print(max_d)

u, v = 13, 8
printnum = 10

S, A, Q = colored_homfly_vectors_and_quiver(u, v)
sp.pprint(Q)



# print_min_vecs(u, v, 30, 10)