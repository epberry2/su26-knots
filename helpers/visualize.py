import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sympy as sp


q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial, path=None, show_plot=True, log_scale=True, edgecolor=None):
    p_dict = polynomial.as_dict()
    keys = p_dict.keys()
    x, y = zip(*keys)
    z = p_dict.values()
    norm = colors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=min(z), vmax=max(z))
    if log_scale:
        scatter = plt.scatter(y, x, c=z, cmap='viridis', norm=norm, s=50, edgecolor=edgecolor)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label('coefficient of monomial (log scale)')
    else:
        scatter = plt.scatter(y, x, c=z, cmap='viridis', s=50, edgecolor=edgecolor)
        colorbar = plt.colorbar(scatter)
        colorbar.set_label('coefficient of monomial')
    plt.xlabel('q powers')
    plt.ylabel('a powers')
    
    plt.title('Lattice points of Homfly Polynomial')
    if path is not None:
        plt.savefig(f"{path}.png")
    if show_plot:
        plt.show()
    
    

#polynomial = sp.Poly(colored_homfly_knot(3, 1, 2), q, a)
#visualize_lattice(polynomial)

