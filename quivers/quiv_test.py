from quivers.knot_quiver import colored_jones_vector_and_quiver, colored_homfly_vectors_and_quiver
from helpers.quantum_algebra import QuantumCombinatorics
import numpy as np
import sympy as sp
from quivers.evaluate_quiver import get_tuples

def print_min_vecs(u, v, n):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    print(S)
    sp.pprint(Q)
    S = np.array(S)
    A = np.array(A)
    for j in range(1, n + 1):
        combos = get_tuples(u, j)
        combo_list = list(combos)
        min_pow = 0
        min_d = (0) * u
        flag = False
        for d in combo_list:
            d = np.array(d)
            p1 = np.dot(S, d)        
            # p2 = np.dot(A, d)
            p3 = np.einsum('i,ij,j', d, Q, d)
            if not flag or p1 + p3 < min_pow:
                min_pow = p1 + p3
                min_d = d
            flag = True
        print(min_d)      
def print_max_vecs(u, v, n):
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

u, v = 3, 1
printnum = 10

print_min_vecs(u, v, printnum)