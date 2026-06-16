import sympy as sp
import numpy as np
import time

from helpers.visualize import visualize_lattice
from quivers.knot_quiver import colored_homfly_vectors_and_quiver
from quivers.knot_vectors import knot_vectors
from quivers.evaluate_quiver import evaluate_quiver_knot


knot = (23, 21)
j = 6
t1 = time.perf_counter()
S, A, Q = colored_homfly_vectors_and_quiver(knot[0],knot[1])
Q_numpy = np.array(Q.tolist(), dtype=int)    
tm = time.perf_counter()
quiv = evaluate_quiver_knot(Q_numpy, S, A, knot[0], knot[1], j)
t2 = time.perf_counter()
print(f"Computed quiver in {tm - t1:.6f} seconds")
print(f"Computed polynomial in {t2 - tm:.6f} seconds")

visualize_lattice(quiv, path="knot", show_plot=True, log_scale=True)