import sympy as sp

from helpers.visualize import visualize_lattice
from colored_homfly.knot_closure import colored_homfly_knot

poly = sp.Poly(colored_homfly_knot(11, 8, 3))
visualize_lattice(poly, path="9-7-3", show_plot=False)