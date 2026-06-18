import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import pandas as pd
import sympy as sp
import seaborn as sns
from quivers.evaluate_quiver import evaluate_quiver_knot, evaluate_quiver_qsub
from quivers.knot_quiver import colored_homfly_vectors_and_quiver

# 1. Define the anchor colors for your gradient scale
# Positions scale from 0.0 (most negative) to 1.0 (most positive)
# 0.5 is exactly zero.
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

# 2. Define a custom Normalization class to create flat plateaus around [-1, 0, 1]
class HybridNormalize(mcolors.Normalize):
    def __init__(self, vmin, vmax, clip=False):
        super().__init__(vmin, vmax, clip)
        
    def __call__(self, value, clip=None):
        # Convert inputs to numpy array safely
        val = np.atleast_1d(value).astype(float)
        result = np.zeros_like(val)
        
        # Determine the dynamic max extension to keep the colorbar perfectly symmetric
        max_abs = max(abs(self.vmin), abs(self.vmax), 2.0)
        
        # Walk through the zones and map them manually to the cmap segments
        # Zone A: Highly negative numbers down to -1
        mask_neg = val < -1
        if np.any(mask_neg):
            result[mask_neg] = 0.0 + 0.4 * (val[mask_neg] - (-max_abs)) / (-1 - (-max_abs))
            
        # Zone B: Fixed plateaus for structural anchors
        result[val == -1] = 0.45
        result[val == 0] = 0.5
        result[val == 1] = 0.55
        
        # Zone C: Slopes between anchors to prevent harsh edges
        mask_sl1 = (val > -1) & (val < 0)
        result[mask_sl1] = 0.45 + 0.05 * (val[mask_sl1] - (-1))
        
        mask_sl2 = (val > 0) & (val < 1)
        result[mask_sl2] = 0.50 + 0.05 * (val[mask_sl2] - 0)
        
        # Zone D: Highly positive numbers starting from 1
        mask_pos = val > 1
        if np.any(mask_pos):
            result[mask_pos] = 0.6 + 0.4 * (val[mask_pos] - 1) / (max_abs - 1)
            
        # Bound between 0 and 1 for matplotlib standard mapping
        return np.clip(result, 0.0, 1.0)



q = sp.symbols('q')
a = sp.symbols('a')

def substitute(u, v, j, p):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    Q_numpy = np.array(Q.tolist(), dtype=int) 
    poly = evaluate_quiver_qsub(Q_numpy, S, A, u, v, j, p)
    
    if poly.eval(1) < 0:
        return -poly
    return poly

def truncate(p, max_power):
    filtered_terms = {
        deg: coeff 
        for deg, coeff in p.as_dict().items() 
        if deg[0] <= max_power
    }
    truncated_poly = sp.Poly.from_dict(filtered_terms, p.gens)
    return truncated_poly

def create_heatmap(u, v, j=10, qsub=2, trunc=50):
    # creates a heatmap to visualize the tail of the colored jones polynomial
    polies = []
    for i in range(1, j+1):
        sub = substitute(u, v, i, qsub)
        polies.append(truncate(sub, trunc))

    data = []
    for idx, p in enumerate(polies):
        # p.to_dict() returns {(degree,): coeff}
        for (deg,), coeff in p.as_dict().items():
            data.append({'Poly_Index': idx + 1, 'Degree': deg, 'Coeff': int(coeff)})


    df = pd.DataFrame(data)

    # Step 2: Pivot so Degrees are columns and Poly_Index are rows
    pivot_df = df.pivot(index='Poly_Index', columns='Degree', values='Coeff').fillna(0)

    # 1. Find the true integer range of your degrees
    min_deg = pivot_df.columns.min()
    max_deg = pivot_df.columns.max()

    # 2. Create a complete, unbroken sequence of integers
    complete_range = range(min_deg, max_deg + 1)

    # 3. Force Pandas to include all columns, filling the missing ones with 0

    pivot_df_fixed = pivot_df.reindex(columns=complete_range, fill_value=0)
    even_columns = [c for c in pivot_df_fixed.columns if c % 2 == 0]
    pivot_df_even = pivot_df_fixed[even_columns]

    norm = HybridNormalize(vmin=pivot_df.values.min(), vmax=pivot_df.values.max())

    # Step 3: Plot the heatmap
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(pivot_df_even, cmap=custom_cmap, norm=norm, annot=False, linewidths=0.5 ,fmt='.0f', cbar=True)
    ax.plot([0,1], [0,0], color='black', linewidth=2.5)
    for row_idx in range(len(pivot_df_even)):
        ax.plot([row_idx+1, row_idx+2 ], [row_idx, row_idx], color='black', linewidth=2.5)
        ax.plot([row_idx+2, row_idx+2], [row_idx, row_idx+1], color='black', linewidth=2.5)
    plt.title("Polynomial Coefficients Heatmap")
    plt.xlabel("Degree (q^d)")
    plt.ylabel("jth Colored Jones Polynomial")
    dir_path = Path("tails")
    dir_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(f"tails/K_{u}_{v}_{qsub}")

create_heatmap(11, 8, 8, 3, 50)