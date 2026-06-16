import sympy as sp
import numpy as np
from tqdm import tqdm
import time
from itertools import combinations, chain
from collections import defaultdict
from quivers.tangle_quiver import quiver
from quivers.tangle_vectors import tangle_vectors
from helpers.normalize_laurent import normalize_laurent_2var, normalize_laurents_2var
from colored_homfly.knot_closure import colored_homfly_knot
from colored_homfly.colored_geo import colored_homfly_geo
from helpers.quantum_algebra import QuantumCombinatorics
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
    # Evaluates the j colored homfly polynomial from the quiver tangle
    qc = QuantumCombinatorics() # creates cache for quantum multinomials
    S = np.array(S)
    A = np.array(A)
    bases_polys = [0 for _ in range(j+1)]
    bases_dicts = [defaultdict(int) for _ in range(j+1)] # creates a list of dictionaries mapping powers to coefficients (a,q) -> k
    combos = get_tuples(u+v, j)
    combo_list = list(combos)

    for d in tqdm(combo_list, desc="Processing tuples", unit="tuple"):      
        d = np.array(d)
        weight = d[:u].sum()
        p1 = np.dot(S, d)
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        #multi_nom = qmultinom(weight, list(d[:u])) * qmultinom(j - weight, list(d[u:]))
        multinom = qc._multiply_polynomials(qc.get_multinomial(d[:u]), qc.get_multinomial(d[u:]))
        sign = (-1) ** (p1 % 2)
        
        for i in range(len(multinom)):
            bases_dicts[weight][p2, p1+p3+(2*i)] += sign * multinom[i]
    keys = chain.from_iterable(bases_dicts)
    min_a, min_q = map(min, zip(*keys))
    normalized = [
        {
            (k[0] - min_a, k[1] - min_q): val 
            for k, val in d.items()
        }
        for d in bases_dicts
    ]
    for i in range(len(bases_dicts)):
        bases_polys[i] = sp.Poly(dict(normalized[i]), (a,q))
    
    return bases_polys

def evaluate_quiver_knot(Q, S, A, u, v, j):
    # Evaluates the j colored homfly polynomial given the K_(u/v) quiver
    qc = QuantumCombinatorics() # creates cache for quantum multinomials
    poly = defaultdict(int)
    S = np.array(S)
    A = np.array(A)
    combos = get_tuples(u, j) # get all u-tuples with values adding up to j
    combo_list = list(combos)

    for d in tqdm(combo_list, desc="Processing tuples", unit="tuple"):        
        d = np.array(d)
        p1 = np.dot(S, d)        
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        multinom = qc.get_multinomial(d)
        #multinom = qc.q_multinomial(d.sum(), list(d))
        sign = (-1) ** (p1 % 2)
        for i, coeff in enumerate(multinom):
            if coeff != 0:
                poly[(p2, p1+p3+(2*i))] += sign * coeff

    # normalize polynomial so lowest powers are 0
    min_a, min_q = map(min, zip(*poly))
    normalized = {
        (k[0] - min_a, k[1] - min_q): val 
        for k, val in poly.items()
    }

    return sp.Poly(dict(normalized), (a, q))

def evaluate_quiver_knot_dict(Q, S, A, u, v, j):
    # Evaluates the j colored homfly polynomial given the K_(u/v) quiver
    qc = QuantumCombinatorics() # creates cache for quantum multinomials
    poly = defaultdict(int)
    S = np.array(S)
    A = np.array(A)
    combos = get_tuples(u, j) # get all u-tuples with values adding up to j
    combo_list = list(combos)

    for d in tqdm(combo_list, desc="Processing tuples", unit="tuple"):        
        d = np.array(d)
        p1 = np.dot(S, d)        
        p2 = np.dot(A, d)
        p3 = np.einsum('i,ij,j', d, Q, d)
        multinom = qc.get_multinomial(d)
        #multinom = qc.q_multinomial(d.sum(), list(d))
        sign = (-1) ** (p1 % 2)
        for i, coeff in enumerate(multinom):
            if coeff != 0:
                poly[(p2, p1+p3+(2*i))] += sign * coeff
    
    return (dict(poly))

def evaluate_quiver_jones(Q, H, u, v, j):
    # Evaluates the j colored jones polynomial given the K_(u/v) quiver
    qc = QuantumCombinatorics() # creates cache for quantum multinomials
    poly = defaultdict(int)
    H = np.array(H)
    combos = get_tuples(u, j) # get all u-tuples with values adding up to j
    combo_list = list(combos)

    for d in tqdm(combo_list, desc="Processing tuples", unit="tuple"):        
        d = np.array(d)
        p1 = np.dot(H, d)   
        p2 = np.einsum('i,ij,j', d, Q, d)
        multinom = qc.get_multinomial(d)
        sign = (-1) ** (p1 % 2)
        for i, coeff in enumerate(multinom):
            if coeff != 0:
                poly[int(p1+p2-(2*i))] += int(sign * coeff)

    # normalize polynomial so lowest powers are 0

    min_value_key = min(poly)
    normalized = {
        k - min_value_key: val 
        for k, val in poly.items()
    }

    return sp.Poly(dict(normalized), q)

def evaluate_quiver_jones_from_homfly(Q, S, A, u, v, j):
    # Evaluates the j colored jones polynomial given the K_(u/v) quiver (homfly quiver)
    qc = QuantumCombinatorics() # creates cache for quantum multinomials
    poly = defaultdict(int)
    S = np.array(S)
    A = np.array(A)
    combos = get_tuples(u, j) # get all u-tuples with values adding up to j
    combo_list = list(combos)

    for d in tqdm(combo_list, desc="Processing tuples", unit="tuple"):        
        d = np.array(d)
        p1 = np.dot(-S, d)   
        p2 = np.dot(2 * A, d)
        p3 = np.einsum('i,ij,j', d, -Q, d)
        multinom = qc.get_multinomial(d)
        sign = (-1) ** (p1 % 2)
        for i, coeff in enumerate(multinom):
            if coeff != 0:
                poly[int(p1+p2+p3-(2*i))] += int(sign * coeff)

    # normalize polynomial so lowest powers are 0

    min_value_key = min(poly)
    normalized = {
        k - min_value_key: val 
        for k, val in poly.items()
    }

    return sp.Poly(dict(normalized), q)

def evaluate_quiver_knot_old(Q, S, A, u, v, j):
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

