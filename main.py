import sympy as sp
import numpy as np

from helpers.visualize import visualize_lattice
from colored_homfly.knot_closure import colored_homfly_knot
from quivers.knot_quiver import knot_quiver
from quivers.knot_vectors import knot_vectors
from quivers.evaluate_quiver import evaluate_quiver_knot


knot = (19,17)
Q = knot_quiver(knot[0],knot[1])
Q_numpy = np.array(Q.tolist(), dtype=int)    
S, A, _ = knot_vectors(knot[0], knot[1])

quiv = evaluate_quiver_knot(Q_numpy, S, A, knot[0], knot[1], 2)

geo = sp.Poly(colored_homfly_knot(knot[0], knot[1], 2))
# print(Q_numpy)
# print(quiv.as_expr())
# print(geo.as_expr())
# print(geo.as_expr() == quiv.as_expr())
# print(geo.as_expr() - quiv.as_expr())

#visualize_lattice(geo, path="knot", show_plot=False)