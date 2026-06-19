import sympy as sp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import time

from helpers.visualize import visualize_lattice_from_dict, create_quiver_table
from quivers.knot_quiver import colored_homfly_vectors_and_quiver, colored_jones_vector_and_quiver
from quivers.knot_vectors import knot_vectors
from quivers.evaluate_quiver import evaluate_quiver_knot_dict, evaluate_quiver_jones_from_homfly, evaluate_quiver_knot

q = sp.symbols('q')
a = sp.symbols('a')


def homfly_tail(u, v, j, trunc=10):
    S, A, Q = colored_homfly_vectors_and_quiver(u,v)
    Q_numpy = np.array(Q.tolist(), dtype=int)       
    polies = []
    xs, ys, zs = [], [], []

    for color in range(1,j+1):
        quiv_dict = evaluate_quiver_knot_dict(Q_numpy, S, A, u, v, color)
        min_qp = min(quiv_dict, key=lambda k: k[1]) # get key with smallest q power
        sign = 1
        if quiv_dict[min_qp] < 0:
            sign = -1
        normalized = {
            (int(k[0] - min_qp[0]), int(k[1] - min_qp[1])): int(sign * val)
            for k, val in quiv_dict.items()
        }
        sorted_items = sorted(normalized.items(), key=lambda item: (item[0][1], item[0][0]))
        trunced = dict(sorted_items[:trunc])
        polies.append(trunced)
        x, y = zip(*trunced)
        xs.append(x)
        ys.append(y)
        zs.append(list(trunced.values()))

    return xs, ys, zs

def homfly_tail_animation(u, v, j, trunc):
    ys, xs, zs = homfly_tail(u, v, j, trunc)
    num_frames = len(xs)
    all_z = [val for frame in zs for val in frame]
    global_min_z = min(all_z)
    global_max_z = max(all_z)

    all_x = [val for frame in xs for val in frame]
    all_y = [val for frame in ys for val in frame]

    shared_norm = mcolors.SymLogNorm(linthresh=1.0, linscale=1.0, vmin=global_min_z, vmax=global_max_z)

    # --- 3. MATPLOTLIB CANVAS SETUP ---
    fig, ax = plt.subplots(figsize=(7, 6))

    # Lock limits to global min/max so the frame window doesn't vibrate/shake
    ax.set_xlim(min(xs[-1]) - 0.5, max(xs[-1]) + 0.5)
    ax.set_ylim(min(ys[-1]) - 0.5, max(ys[-1]) + 0.5)

    ax.set_xlabel('q powers')
    ax.set_ylabel('a powers')
    title = ax.set_title('j=1')

    scatter = ax.scatter([], [], c=[], cmap='viridis', norm=shared_norm, s=80, edgecolor=None)
    fig.colorbar(scatter, ax=ax, label='Coefficient Intensity')

    def update(frame_idx):
        # Get the x, y, and z lists for the current frame
        x = xs[frame_idx]
        y = ys[frame_idx]
        z = zs[frame_idx]
        
        # 1. Update the coordinates (requires an Nx2 numpy array shape)
        positions = np.column_stack((x, y))
        scatter.set_offsets(positions)
        
        # 2. Update the color mapping values
        z_numeric = np.array(z, dtype=float)
        scatter.set_array(z)
        
        # 3. Update the title text dynamically
        title.set_text(f'j={frame_idx+1}')
        
        return scatter, title
    ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=500, blit=False)

    output_path = f"lattice_evolution_{u}_{v}.gif"
    ani.save(output_path, writer='pillow')

    print(f"Successfully generated and saved {output_path}")
    plt.close()
    return True

u, v = 7, 3

homfly_tail_animation(u, v, 15, trunc=100)