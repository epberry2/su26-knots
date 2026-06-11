import sympy as sp
import numpy as np
from itertools import combinations
from collections import defaultdict
from helpers.quantum_nums import qmultinom

q = sp.symbols('q')
a = sp.symbols('a')

def get_tuples(k, n):
    # uses stars and bars method to get all k-tuples with values adding up to n
    # If k is 1, there's only one valid tuple: (n,)
    if k == 1:
        return [(n,)]
        
    result = []
    # We choose k-1 divider positions out of a total of n + k - 1 slots
    for dividers in combinations(range(n + k - 1), k - 1):
        # Add boundaries at the start (-1) and end (n + k - 1)
        full_dividers = (-1,) + dividers + (n + k - 1,)
        
        # Calculate the distance between dividers to get the tuple values
        tup = tuple(full_dividers[i+1] - full_dividers[i] - 1 for i in range(k))
        result.append(tup)
        
    return result


def evaluate_quiver(Q, S, A, u, v, j):
    # Evaluates the j colored homfly polynomial from the quiver
    S = np.array(S)
    A = np.array(A)
    bases_polys = [0 for _ in range(j+1)]
    combos = get_tuples(u+v, j)
    for d in combos:        
        d = np.array(d)
        weight = d[:u].sum()
        p1 = np.dot(S, d)
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        multi_nom = qmultinom(weight, list(d[:u])) * qmultinom(j - weight, list(d[u:]))
        bases_polys[weight] += a**p2 * (-q)**p1 * q**p3 * multi_nom
    for i in range(len(bases_polys)):
        bases_polys[i] = sp.expand(sp.simplify(bases_polys[i]))
    return bases_polys

Q = np.array([[1, 1, 0, 0],
              [1, 2, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0]])
S = [2, 3, 1, 0]
A = [0, 0, 0, 0]
p = evaluate_quiver(Q, S, A, 3, 1, 2)
print(p)