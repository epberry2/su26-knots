import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import math
import sympy as sp

from quivers.knot_quiver import colored_homfly_vectors_and_quiver
from helpers.tanglestate import get_state
from visualizers.heatmap import HybridNormalize

colors = [
    (0.0, '#084594'),   # Far negative: Deep Blue
    (0.4, '#80b1d3'),   # Approaches -1: Smooth Light Blue
    (0.45, '#b3de69'),  # EXACTLY -1: Solid Sage Green
    (0.5, '#e0e0e0'),   # EXACTLY 0: Solid Neutral Gray
    (0.55, '#fdb462'),  # EXACTLY 1: Solid Pastel Orange
    (0.6, '#fb8072'),   # Leaves 1: Smooth Light Red
    (1.0, '#b10026')    # Far positive: Deep Red
]
custom_cmap = mcolors.LinearSegmentedColormap.from_list('custom_hybrid', colors)

q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial, path=None, show_plot=True, log_scale=True, edgecolor=None, thick_dots = False, cols = 'viridis'):
    """Visualizes the colored Homfly polynomial of some rational knot"""
    p_dict = polynomial.as_dict()
    keys = p_dict.keys()
    x, y = zip(*keys)
    z = p_dict.values()
    sizes = [float(np.clip(abs(c), 5, 35)) for c in z] if thick_dots else 10
    norm = mcolors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=min(z), vmax=max(z))
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
        plt.savefig(f"{path}.svg", format="svg")
        print(f"plot saved to {path}.svg")
    if show_plot:
        plt.show()

def visualize_lattice_colormap(polynomial, path=None, show_plot=True):
    """Visualizes the colored Homfly polynomial of some rational knot"""
    coeff_dict = polynomial.as_dict()
    keys = coeff_dict.keys()
    x, y = zip(*keys)
    z = coeff_dict.values()
    sizes = 10
    #norm = colors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=min(z), vmax=max(z))
    norm = HybridNormalize(vmin=min(coeff_dict.values()), vmax=max(coeff_dict.values()))

    scatter = plt.scatter(y, x, c=z, s = sizes, cmap=custom_cmap, norm=norm)
    colorbar = plt.colorbar(scatter)
    colorbar.set_label('coefficient of monomial')
    

    plt.xlabel('q powers')
    plt.ylabel('a powers')
    
    plt.title('Lattice points of Homfly Polynomial')
    if path is not None:
        plt.savefig(f"{path}.svg", format="svg")
        print(f"plot saved to {path}.svg")
    if show_plot:
        plt.show()


def visualize_lattice_from_dict(p_dict, path=None, show_plot=True, log_scale=True, edgecolor=None, thick_dots = False, cols = 'viridis'):
    plt.figure()
    
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
        


  




