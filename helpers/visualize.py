import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial, path=None, show_plot=True):
    p_dict = polynomial.as_dict()
    keys = p_dict.keys()
    x, y = zip(*keys)
    z = p_dict.values()
    scatter = plt.scatter(x, y, c=z, cmap='viridis', s=50, edgecolor='black')
    plt.xlabel('q powers')
    plt.ylabel('a powers')
    colorbar = plt.colorbar(scatter)
    colorbar.set_label('coefficient of monomial')
    plt.title('Lattice points of Homfly Polynomial')
    if path is not None:
        plt.savefig(f"{path}.png")
    if show_plot:
        plt.show()
    
    

#polynomial = sp.Poly(colored_homfly_knot(3, 1, 2), q, a)
#visualize_lattice(polynomial)

