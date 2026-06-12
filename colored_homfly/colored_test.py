import sympy as sp
import time
import numpy as np

from colored_homfly.colored_alg import closure_polynomial_matrix, evaluate_tangle_vector
from colored_homfly.colored_geo import colored_homfly_geo
from colored_homfly.knot_closure import colored_homfly_knot
from helpers.normalize_laurent import *
from helpers.tanglestate import get_state
from quivers.tangle_quiver import quiver
from quivers.tangle_vectors import tangle_vectors
from quivers.evaluate_quiver import evaluate_quiver

a = sp.symbols('a')
q = sp.symbols('q')

def test_jcolor_polynomial(u, v, j, times=True):
    print("")
    print(f"Testing {(u, v, j)} knot")
    t1 = time.perf_counter()
    geo = colored_homfly_knot(u, v, j)
    t2 = time.perf_counter()
    if times:
        print(f"Geo computation took {t2 - t1:.6f} seconds")
    alg = normalize_laurent_2var(closure_polynomial_matrix(u, v, j), a, q)
    t3 = time.perf_counter()
    if times:
        print(f"Alg computation took {t3 - t2:.6f} seconds")
    return alg == geo

def test_multiple(knots, n):
    begin = time.perf_counter()
    correct = 0
    total_knots = len(knots) * n
    for knot in knots:
        state = get_state(knot[0], knot[1])
        if state[0] == "RI":
            print("")
            print(f"Cannot close tangle, skipping {(knot[0], knot[1])}")
            total_knots -= n
            continue
        for i in range(1,n+1):
            x = test_jcolor_polynomial(knot[0], knot[1], i)
            print(x)
            if x:
                correct += 1
    end = time.perf_counter()
    print("")
    print(f"Testing completed in {end - begin:.2f} seconds")
    print(f"{correct} knots match out of {total_knots}")


def test_quiver(u, v, j, times=True):
    print("")
    print(f"Testing {(u, v, j)} knot")

    t1 = time.perf_counter()
    geo = colored_homfly_geo(u, v, j)
    t2 = time.perf_counter()

    if times:
        print(f"Geo computation took {t2 - t1:.6f} seconds")
    
    Q = quiver(u,v)
    S, A = tangle_vectors(u,v)
    tq = time.perf_counter()
    if times:
        print(f"Computing quiver took {tq - t2:.6f} seconds")
    Q_numpy = np.array(Q.tolist(), dtype=int)    
    quiv = normalize_laurents_2var(evaluate_quiver(Q_numpy, S, A, u, v, j), a, q)
    t3 = time.perf_counter()
    if times:
        print(f"Evaluating quiver took {t3 - tq:.6f} seconds")
    for i in range(j):
        if geo[i].as_expr() != quiv[i]:
            print(geo[i].as_expr())
            print(quiv[i])
            return False
    return True

print(test_jcolor_polynomial(7, 5, 3))

# knots = [(3, 1), (7, 3)]
# test_multiple(knots, 3)
