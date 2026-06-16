from quivers.knot_quiver import colored_jones_vector_and_quiver, colored_homfly_vectors_and_quiver
from helpers.quantum_algebra import QuantumCombinatorics
import numpy as np
import sympy as sp
from quivers.evaluate_quiver import get_tuples


u, v = 55, 2

S, A, Q = colored_homfly_vectors_and_quiver(u, v)
print(S)
sp.pprint(Q)
S = np.array(S)
A = np.array(A)
for j in range(10):
    combos = get_tuples(u, j)
    combo_list = list(combos)
    min_pow = 0
    min_d = (0) * u
    flag = False
    for d in combo_list:
        d = np.array(d)
        p1 = np.dot(S, d)        
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        if not flag or p1 + p3 < min_pow:
            min_pow = p1 + p3
            min_d = d
        flag = True
    print(min_d, np.einsum('i,ij,j', min_d, Q, min_d))