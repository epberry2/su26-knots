import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
import sympy as sp

from quivers.knot_quiver import colored_homfly_vectors_and_quiver, colored_jones_vector_and_quiver, winding_tracker
from helpers.tanglestate import get_state
from quivers.tail_quivs import colored_homfly_head_vectors_and_quiver, colored_homfly_tail_vectors_and_quiver


q = sp.symbols('q')
a = sp.symbols('a')

def visualize_lattice(polynomial, save_path=None, show_plot=True, log_scale=True, edgecolor=None, thick_dots = False, cols = 'viridis'):
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

    if save_path is not None:
        plt.savefig(f"plots/{save_path}.png")
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

def jones_quiv(i, j):
    return colored_jones_vector_and_quiver(i, j)[1]

    
def create_quiver_table(n, path="table"):
    # writes to text file all quivers of knots up to numerator n
    with open(f"{path}.txt", 'w') as f:
        for i in range(1,n+1,2):
            for j in range(1,i):
                try:                   
                    # _, Q = colored_jones_vector_and_quiver(i,j)
                    # u = i - j
                    # v = j
                    # if v <= u // 2 :
                    #     raise ValueError(" ")
                    # orientation = get_state(u, v)[0]
                    # if not orientation == "UP":
                    #     raise ValueError(" ")
                    # S, A, Q = colored_homfly_vectors_and_quiver(i, j)
                    
                    # diag_terms = [(Q[ii, ii], ii) for ii in range(i)]
                    # diag_terms.sort(key=lambda x: -x[0])
                    # if diag_terms[0][0] > diag_terms[1][0] and diag_terms[1][0] > diag_terms[2][0]:
                    #     # raise ValueError(" ")
                    #     # f.write(f"c_1 - c_2 = {H[diag_terms[0][1]] - H[diag_terms[1][1]]} \n")
                    #     "do nothiong"
                    # else:
                    #     f.write("WEIRD KNOT...\n")
                    
                
                    # S, A, Q = colored_homfly_tail_vectors_and_quiver(i, j)
                    S, A, Q = colored_homfly_head_vectors_and_quiver(i, j)
                    
                    # np.savetxt(
                    #     f,
                    #     Q,
                    #     fmt="%d",
                    #     delimiter=" & ",
                    #     newline=r" \\" + "\n"
                    # )
                    # f.write("\n")
                    
                    f.write(f"S, A, Q for K_{i}/{j}:\n")
                    f.write(f"{S}\n")
                    f.write(f"{A}\n")
                    # f.write(f"{[diag_terms[ii][0] for ii in range(i)]}\n")
                    # w = winding_tracker(i, j)
                    # _, _, vec = w.homfly_vectors()
                    # sorted_vec = [vec[idx] for idx in w.intersection_path]
                    # f.write(f"K_{i}/{j} diag els in order:\n")
                    # f.write(f"{sorted_vec}\n")

                    # np.savetxt(f, Q_numpy, fmt="%3d", delimiter=" ")
                    
                    
                    np.savetxt(
                        f,
                        Q,
                        fmt="%3d",
                        delimiter=" ",
                        # newline=r" \\" + "\n"
                    )
                    f.write("\n")
                except:
                    pass
                
create_quiver_table(50, path = "tail_minors_latex")


