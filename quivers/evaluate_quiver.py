import sympy as sp
import numpy as np
from itertools import combinations
from collections import defaultdict
from helpers.quantum_nums import qmultinom
from quivers.tangle_quiver import quiver
from quivers.tangle_vectors import tangle_vectors
from helpers.normalize_laurent import normalize_laurent_2var, normalize_laurents_2var

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
    # Evaluates the j colored homfly polynomial from the quiver tangle
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
    return normalize_laurents_2var(bases_polys, a, q)

def evaluate_quiver_knot(Q, S, A, u, v, j):
    # Evaluates the j colored homfly polynomial from the quiver knot
    S = np.array(S)
    A = np.array(A)
    homfly = 0
    combos = get_tuples(u, j)
    for d in combos:        
        d = np.array(d)
        p1 = np.dot(S, d)
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        multi_nom = qmultinom(d.sum(), list(d))
        homfly += a**p2 * (-q)**p1 * q**p3 * multi_nom

    return normalize_laurent_2var(homfly, a, q)

knot = (5, 2)

Q = quiver(knot[0], knot[1])

S, A = tangle_vectors(knot[0], knot[1])

Q_numpy = np.array(Q.tolist(), dtype=int)

quiv = evaluate_quiver(Q_numpy, S, A, knot[0], knot[1], 2)



# print(p)