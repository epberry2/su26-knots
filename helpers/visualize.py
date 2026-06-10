import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial):
    p_dict = polynomial.as_dict()
    keys = p_dict.keys()
    x, y = zip(*keys)
    plt.scatter(x, y)
    plt.xlabel('q powers')
    plt.ylabel('a powers')
    plt.title('Lattice points of Homfly Polynomial')
    plt.show()
    

#polynomial = sp.Poly(colored_homfly_knot(3, 1, 2), q, a)
#visualize_lattice(polynomial)

