import sympy as sp
import time
import numpy as np

from colored_homfly.colored_alg import closure_polynomial_matrix, evaluate_tangle_vector
from colored_homfly.colored_geo import colored_homfly_geo
from colored_homfly.knot_closure import colored_homfly_knot
from helpers.normalize_laurent import normalize_laurent_2var
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
        print(f"Alg computation took {t3 - t1:.6f} seconds")
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


knots = [(3, 1), (7, 3)]
#test_multiple(knots, 3)

knot = (5, 2)

Q = quiver(knot[0], knot[1])
S, A = tangle_vectors(knot[0], knot[1])
Q_numpy = np.array(Q.tolist(), dtype=int)

quiv = evaluate_quiver(Q_numpy, S, A, knot[0], knot[1], 2)
geo = colored_homfly_geo(knot[0], knot[1], 2)
print(quiv)
print(geo)

