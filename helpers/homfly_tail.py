import sympy as sp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import time

from quivers.knot_quiver import colored_homfly_vectors_and_quiver, colored_jones_vector_and_quiver
from quivers.evaluate_quiver import evaluate_quiver_knot_dict, evaluate_quiver_jones_from_homfly, evaluate_quiver_knot

q = sp.symbols('q')
a = sp.symbols('a')


def homfly_tail(u, v, j, trunc=10, dir="se"):
    if dir not in {"se", "sw", "nw", "ne"}:
        raise ValueError("direction is se, sw, nw, ne")
    S, A, Q = colored_homfly_vectors_and_quiver(u,v)
    Q_numpy = np.array(Q.tolist(), dtype=int)       
    polies = []
    xs, ys, zs = [], [], []

    for color in range(1,j+1):
        quiv_dict = evaluate_quiver_knot_dict(Q_numpy, S, A, u, v, color)
        if dir == "se":
            min_qp = max(quiv_dict, key=lambda k: k[1]) # get key with highest q power
            min_ap = min(quiv_dict, key=lambda k: k[0]) # get key with smallest a power
        elif dir == "nw":
            min_qp = min(quiv_dict, key=lambda k: k[1]) # get key with smallest q power
            min_ap = max(quiv_dict, key=lambda k: k[0]) # get key with highest a power
        elif dir == "sw":
            min_qp = min(quiv_dict, key=lambda k: k[1]) # get key with smallest q power
            min_ap = min(quiv_dict, key=lambda k: k[0]) # get key with highest a power
        else:
            min_qp = max(quiv_dict, key=lambda k: k[1]) # get key with smallest q power
            min_ap = max(quiv_dict, key=lambda k: k[0]) # get key with highest a power
        sign = 1
        if quiv_dict[min_qp] < 0:
            sign = -1
        normalized = {
            (int(k[0] - min_ap[0]), int(k[1] - min_qp[1])): int(sign * val)
            for k, val in quiv_dict.items()
        }
        if dir == "se":
            sorted_items = sorted(normalized.items(), key=lambda item: (-item[0][1], item[0][0]))
        elif dir == "nw":
            sorted_items = sorted(normalized.items(), key=lambda item: (item[0][1], item[0][0]))
        elif dir == "sw":
            sorted_items = sorted(normalized.items(), key=lambda item: (item[0][1], -item[0][0]))
        else:
            sorted_items = sorted(normalized.items(), key=lambda item: (-item[0][1], -item[0][0]))
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

    output_path = f"homfly_tail_lattices/lattice_evolution_{u}_{v}.mp4"
    ani.save(output_path, writer='ffmpeg')

    print(f"Successfully generated and saved {output_path}")
    plt.close()
    return True

def homfly_tail_heatmap(u, v, j, trunc, dir="se"):
    xs, ys, zs = homfly_tail(u, v, j, trunc, dir)

    num_frames = j

    all_z = [val for frame in zs for val in frame]
    global_min_z = min(all_z)
    global_max_z = max(all_z)

    all_unique_x = sorted(list(set([val for frame in xs for val in frame])))
    all_unique_y = sorted(list(set([val for frame in ys for val in frame])))

    fig, ax = plt.subplots(figsize=(10, 5))
    
    def update(frame_idx):
        # Clear the previous frame's heatmap and colorbar axis
        ax.clear()

        x = xs[frame_idx]
        y = ys[frame_idx]
        z = zs[frame_idx]

        df = pd.DataFrame({'x': x, 'y': y, 'z': z})
        matrix = df.pivot(index='x', columns='y', values='z')
        matrix = matrix.reindex(index=all_unique_x, columns=all_unique_y).fillna(0).astype(int)
        matrix = matrix.iloc[::-1]
        
        # Get the matrix for the current frame
        sns.heatmap(
            matrix, 
            ax=ax, 
            cmap='seismic', 
            vmin=global_min_z, 
            vmax=global_max_z, 
            cbar=False,
            center=0,
            annot=True, # Optional: Set to True if you want to see the 0s written out
            fmt="d"   # Text formatting for annotations
        )
        
        ax.set_title(f"Lattice Density (j = {frame_idx + 1})")
        ax.set_xlabel("q powers")
        ax.set_ylabel("a powers")

    mappable = plt.cm.ScalarMappable(
        norm=plt.Normalize(vmin=global_min_z, vmax=global_max_z), 
        cmap='seismic'
    )
    fig.colorbar(mappable, ax=ax, label='Coefficient Intensity')

    ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=200, blit=False)

    output_path = f"homfly_tail_lattices/heatmap_evolution_{dir}_{u}_{v}.mp4"
    ani.save(output_path, writer='ffmpeg')

    print(f"Heatmap video successfully saved to {output_path}")
    plt.close()

homfly_tail_heatmap(13, 3, 6, trunc=100, dir="se")