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
from helpers.heatmap import HybridNormalize

q = sp.symbols('q')
a = sp.symbols('a')

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


def homfly_tail(u, v, j, trunc=10, dir="se"):
    if dir not in {"se", "sw", "nw", "ne"}:
        raise ValueError("direction is se, sw, nw, ne")
    S, A, Q = colored_homfly_vectors_and_quiver(u,v)
    Q_numpy = np.array(Q.tolist(), dtype=int)       
    polies = []
    xs, ys, zs = [], [], []

    for color in range(1,j+1):
        quiv_dict = evaluate_quiver_knot_dict(Q_numpy, S, A, u, v, color)
        sign = 1
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
        if quiv_dict[min_ap] < 0:
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
    min_x = min(all_unique_y)
    max_a = max(all_unique_x)
    q_crop = min(ys[-1])
    cropped_q = [val for val in all_unique_y if q_crop <= val]

    custom_norm = HybridNormalize(vmin=global_min_z, vmax=global_max_z)

    fig, ax = plt.subplots(figsize=(15, 5))
    
    def update(frame_idx):
        # Clear the previous frame's heatmap and colorbar axis
        ax.clear()
        fdx = frame_idx + 1
        x = xs[frame_idx]
        y = ys[frame_idx]
        z = zs[frame_idx]

        df = pd.DataFrame({'x': x, 'y': y, 'z': z})
        matrix = df.pivot(index='x', columns='y', values='z')
        matrix = matrix.reindex(index=all_unique_x, columns=cropped_q).fillna(0).astype(int)
        matrix = matrix.iloc[::-1]
        
        # Get the matrix for the current frame
        sns.heatmap(
            matrix, 
            ax=ax, 
            cmap=custom_cmap, 
            norm=custom_norm,
            vmin=global_min_z, 
            vmax=global_max_z, 
            cbar=False,
            annot=True, # Optional: Set to True if you want to see the 0s written out
            fmt="d"   # Text formatting for annotations
        )

        if dir == "se":
            # plot stable zone
            line_color = 'black'
            _, A, _ = colored_homfly_vectors_and_quiver(u, v)
            if v == 1:
                shift = q_crop // 2 - 1
                yshift = max_a // 2 + 1
                ax.plot([-fdx - shift - 1,-shift], [yshift,yshift], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift - 1,-fdx - shift - 1], [yshift,yshift - 1], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift - 1,-fdx - shift], [yshift - 1,yshift - 1], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift,-fdx - shift], [yshift - 1,yshift - 2], color=line_color, linewidth=2.5)
                for i in range(2, yshift):
                    xl = -(fdx * i) - ((i)**2 + (i)) // 2 + 1 - shift + (i-1)
                    xr = -(fdx * (i-1)) - ((i-1)**2 + (i-1)) // 2 + 1 - shift + (i - 2)
                    ax.plot([xl,xr], [yshift - i,yshift - i], color=line_color, linewidth=2.5)
                    ax.plot([xl,xl], [yshift - i,yshift - i-1], color=line_color, linewidth=2.5)         
                    pass
            elif v % 2 == 0:
                shift = q_crop // 2 - 1
                yshift = max_a // 2 + 1
                if A.count(min(A)) == 1:
                    ax.plot([0,- shift], [yshift,yshift], color=line_color, linewidth=2.5)
                    #ax.plot([0,0], [yshift, yshift - 1], color=line_color, linewidth=2.5)
                    ax.plot([0,- fdx - shift], [yshift - 1,yshift - 1], color=line_color, linewidth=2.5)
                    #ax.plot([-fdx - shift,-fdx - shift], [yshift - 1,yshift - 2], color=line_color, linewidth=2.5)
                else:
                    ax.plot([-fdx - 2 - shift,- shift], [yshift,yshift], color=line_color, linewidth=2.5)
                    ax.plot([-fdx-2-shift,-fdx-2-shift], [yshift, yshift - 1], color=line_color, linewidth=2.5)
                    ax.plot([-fdx-2-shift,- fdx - shift], [yshift - 1,yshift - 1], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift,-fdx - shift], [yshift - 1,yshift - 2], color=line_color, linewidth=2.5)                
                ax.plot([-fdx- shift,-fdx+1 - shift], [yshift - 2,yshift - 2], color=line_color, linewidth=2.5)
                ax.plot([-fdx+1 - shift,-fdx+1- shift], [yshift - 2,yshift - 4], color=line_color, linewidth=2.5)
                for i in range(4, yshift):
                    xl = -fdx - ((i-3)**2 + (i-3)) // 2 + 1 - shift
                    xr = -fdx - ((i-4)**2 + (i-4)) // 2 + 1 - shift
                    ax.plot([xl,xr], [yshift - i,yshift - i], color=line_color, linewidth=2.5)
                    ax.plot([xl,xl], [yshift - i,yshift - i-1], color=line_color, linewidth=2.5)
                #ax.plot([-fdx - ((yshift - 4)**2 + (yshift - 4)) // 2 + 1 - shift, - shift],[0,0], color=line_color, linewidth=2.5)
            elif v % 2 == 1:
                shift = q_crop // 2 - 1
                yshift = max_a // 2 + 1
                ax.plot([-fdx - shift - 1,-shift], [yshift,yshift], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift - 1,-fdx - shift - 1], [yshift,yshift - 1], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift - 1,-fdx - shift], [yshift - 1,yshift - 1], color=line_color, linewidth=2.5)
                ax.plot([-fdx - shift,-fdx - shift], [yshift - 1,yshift - 3], color=line_color, linewidth=2.5)
                for i in range(3, yshift):
                    xl = -fdx - ((i-2)**2 + (i-2) + 1) + 1 - shift
                    xr = -fdx - ((i-3)**2 + (i-3) + 1) + 1 - shift
                    ax.plot([xl,xr], [yshift - i,yshift - i], color=line_color, linewidth=2.5)
                    ax.plot([xl,xl], [yshift - i,yshift - i-1], color=line_color, linewidth=2.5)
        ax.set_title(f"Lattice Density (j = {frame_idx + 1})")
        ax.set_xlabel("q powers")
        ax.set_ylabel("a powers")

    mappable = plt.cm.ScalarMappable(
        norm=custom_norm,
        cmap=custom_cmap
    )
    fig.colorbar(mappable, ax=ax, label='Coefficient Intensity')

    ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=200, blit=False)

    output_path = f"homfly_tail_lattices/heatmap_evolution_{dir}_{u}_{v}.mp4"
    ani.save(output_path, writer='ffmpeg')

    print(f"Heatmap video successfully saved to {output_path}")
    plt.close()

def homfly_plots(u, v, trunc=100, dir="se"):
    xs, ys, zs = homfly_tail(u, v, 10, trunc, dir)
    
    fig, ax = plt.subplots(2, 5, figsize=(20,12), sharex=True, sharey=True)

    num_plots = 10 

    all_z = [val for frame in zs for val in frame]
    global_min_z = min(all_z)
    global_max_z = max(all_z)

    all_unique_x = sorted(list(set([val for frame in xs for val in frame])))
    all_unique_y = sorted(list(set([val for frame in ys for val in frame])))

    for j in range(2):
        for frame_idx in range(5):
            x = xs[frame_idx + j * 5]
            y = ys[frame_idx + j * 5]
            z = zs[frame_idx + j * 5]

            df = pd.DataFrame({'x': x, 'y': y, 'z': z})
            matrix = df.pivot(index='x', columns='y', values='z')
            matrix = matrix.reindex(index=all_unique_x, columns=all_unique_y).fillna(0).astype(int)
            matrix = matrix.iloc[::-1]
            
            # Get the matrix for the current frame
            sns.heatmap(
                matrix, 
                ax=ax[j, frame_idx], 
                cmap='seismic', 
                vmin=global_min_z, 
                vmax=global_max_z, 
                cbar=False,
                center=0,
                annot=True, # Optional: Set to True if you want to see the 0s written out
                fmt="d"   # Text formatting for annotations
            )

            ax[j, frame_idx].set_title(f"j={j * 5 + frame_idx + 1}", fontweight="bold")
            ax[j, frame_idx].set_xlabel("")
            ax[j, frame_idx].set_ylabel("")

    mappable = plt.cm.ScalarMappable(
        norm=mcolors.Normalize(vmin=global_min_z, vmax=global_max_z), 
        cmap='seismic'
    )

        # --- ADD OUTLINES TO EACH SUBPLOT ---
    # If your axes array is 2D (e.g. 2x2 grid), use axes.flat to loop through all of them
    for axi in ax.flat:
        # 1. Force the spine borders to be visible (Seaborn heatmaps sometimes turn them off)
        for spine in axi.spines.values():
            spine.set_visible(True)
            spine.set_color('black')       # Set your outline color
            spine.set_linewidth(3.0)       # Set thickness of the outline
            
        # Optional: If you want a small gap between the heatmap and the border line
        # ax.set_frame_on(True)
    fig.supxlabel("q powers", fontweight="bold")
    fig.supylabel("a powers", fontweight="bold")
    fig.suptitle(f"Heatmaps for K_{u}/{v}", fontsize=16, fontweight='bold')
    # Passing ax=axes.tolist() spans the colorbar nicely across both subplots
    cbar = fig.colorbar(mappable, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
    cbar.set_label('Coefficient Value')


    output_path = f"homfly_tail_lattices/heatmap_subplots_{dir}_{u}_{v}.png"
    plt.savefig(output_path)

    print(f"Heatmap plots successfully saved to {output_path}")
    plt.close()

#homfly_plots(5, 1, trunc=25, dir="se")

homfly_tail_heatmap(9, 4, j=10, trunc=300, dir="nw")