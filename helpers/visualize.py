import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import sympy as sp

from quivers.knot_quiver import colored_homfly_vectors_and_quiver


q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial, path=None, show_plot=True, log_scale=True, edgecolor=None, thick_dots = False, cols = 'viridis'):
    p_dict = polynomial.as_dict()
    keys = p_dict.keys()
    x, y = zip(*keys)
    z = p_dict.values()
    sizes = [float(np.clip(abs(c), 5, 35)) for c in z] if thick_dots else 10
    norm = colors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=min(z), vmax=max(z))
    if log_scale:
        if cols in ["black", "red", "blue", "green"]:
            scatter = plt.scatter(y, x, s = sizes, color=cols, norm=norm, edgecolor=edgecolor)
        else:
            scatter = plt.scatter(y, x, c=z, s = sizes, cmap=cols, norm=norm, edgecolor=edgecolor)
            colorbar = plt.colorbar(scatter)
            colorbar.set_label('coefficient of monomial (log scale)')
    else:
        if cols in ["black", "red", "blue", "green"]:
            scatter = plt.scatter(y, x, s = sizes, color=cols, edgecolor=edgecolor)
        else:
            scatter = plt.scatter(y, x, c=z, s = sizes, cmap=cols, edgecolor=edgecolor)
            colorbar = plt.colorbar(scatter)
            colorbar.set_label('coefficient of monomial')
    plt.xlabel('q powers')
    plt.ylabel('a powers')
    
    plt.title('Lattice points of Homfly Polynomial')
    if path is not None:
        plt.savefig(f"{path}.png")
    if show_plot:
        plt.show()

def visualize_lattice_from_dict(p_dict, path=None, show_plot=True, log_scale=True, edgecolor=None, thick_dots = False, cols = 'viridis'):
    keys = p_dict.keys()
    x, y = zip(*keys)
    z = p_dict.values()
    sizes = [float(min(abs(c), 50)) for c in z] if thick_dots else 10
    norm = colors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=min(z), vmax=max(z))
    if log_scale:
        if cols in ["black", "red", "blue", "green"]:
            scatter = plt.scatter(y, x, s = sizes, color=cols, norm=norm, edgecolor=edgecolor)
        else:
            scatter = plt.scatter(y, x, c=z, s = sizes, cmap=cols, norm=norm, edgecolor=edgecolor)
            colorbar = plt.colorbar(scatter)
            colorbar.set_label('coefficient of monomial (log scale)')
    else:
        if cols in ["black", "red", "blue", "green"]:
            scatter = plt.scatter(y, x, s = sizes, color=cols, edgecolor=edgecolor)
        else:
            scatter = plt.scatter(y, x, c=z, s = sizes, cmap=cols, edgecolor=edgecolor)
            colorbar = plt.colorbar(scatter)
            colorbar.set_label('coefficient of monomial')
    plt.xlabel('q powers')
    plt.ylabel('a powers')
    
    plt.title('Lattice points of Homfly Polynomial')
    if path is not None:
        plt.savefig(f"{path}.png")
    if show_plot:
        plt.show()
        
def diff_quiv(i, j):
    _, _, Q = colored_homfly_vectors_and_quiver(i, j)
    D = Q[1:, 1:] - Q[1:, :-1] - Q[:-1, 1:] + Q[:-1, :-1]
    return D

def quiv(i, j):
    return colored_homfly_vectors_and_quiver(i, j)[2]

    
def create_quiver_table(n, path="table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    Q = diff_quiv(i,j)
                    Q_numpy = np.array(Q.tolist(), dtype=int) 
                    f.write(f"Quiver for K_{i}/{j}\n")
                    np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    f.write("\n")
                except:
                    pass
                
# create_quiver_table(30)


