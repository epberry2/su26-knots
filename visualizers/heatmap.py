import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import pandas as pd
import sympy as sp
import seaborn as sns
from quivers.evaluate_quiver import evaluate_quiver_knot, evaluate_quiver_qsub, evaluate_quiver_fraction
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
        # Convert inputs to a clean floating-point numpy array
        val = np.atleast_1d(value).astype(float)
        
        # Establish symmetric max bounds
        max_abs = max(abs(self.vmin), abs(self.vmax), 2.0)
        
        # 1. Define distinct, mutually exclusive conditions
        conditions = [
            (val < -1),                                # Zone A: Extreme negatives
            (val == -1),                               # Anchor: Exactly -1
            (val > -1) & (val < 0),                    # Zone C1: Slope between -1 and 0
            (val == 0),                                # Anchor: Exactly 0
            (val > 0) & (val < 1),                     # Zone C2: Slope between 0 and 1
            (val == 1),                                # Anchor: Exactly 1
            (val > 1)                                  # Zone D: Extreme positives
        ]
        
        # 2. Define corresponding mathematical formulas for each condition
        choices = [
            0.0 + 0.4 * (val - (-max_abs)) / (-1 - (-max_abs)), # Zone A Formula
            0.45,                                               # Anchor -1
            0.45 + 0.05 * (val - (-1)),                         # Slope C1 Formula
            0.50,                                               # Anchor 0 (Guaranteed Neutral Gray)
            0.50 + 0.05 * (val - 0),                            # Slope C2 Formula
            0.55,                                               # Anchor 1
            0.60 + 0.4 * (val - 1) / (max_abs - 1)              # Zone D Formula
        ]
        
        # np.select matches conditions to choices cleanly; default handles fallback
        result = np.select(conditions, choices, default=0.5)
        
        return np.clip(result, 0.0, 1.0)



q = sp.symbols('q')
a = sp.symbols('a')

def substitute(u, v, j, p):
    # computes the jth colored homfly polynomial and then substitutes q -> q^-1 and q -> q^p
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    Q_numpy = np.array(Q.tolist(), dtype=int) 
    poly = evaluate_quiver_qsub(Q_numpy, S, A, u, v, j, p)
    
    # normalize so that the first term is positive
    if poly(0) < 0:
        return -poly
    return poly

def sub_frac(u, v, j, b, c):
    # computes the jth colored homfly polynomial and then substitutes q -> q^-c and q -> q^b
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    Q_numpy = np.array(Q.tolist(), dtype=int) 
    poly = evaluate_quiver_fraction(Q_numpy, S, A, u, v, j, b, c)
    
    # normalize so that the first term is positive
    if poly(0) < 0:
        return -poly
    return poly

def invert_and_renormalize(expr, variable):
    # Step 1: Invert all powers by substituting variable with 1/variable
    inverted_expr = expr.subs(variable, 1/variable)
    
    # Expand to simplify expressions like (1/q)**2 into q**-2
    inverted_expr = sp.expand(inverted_expr)
    
    # Get the components of the inverted polynomial
    coeff_dict = inverted_expr.as_coefficients_dict()
    
    # Helper to extract the numerical exponent from any term
    def get_exponent(term):
        if term == 1:
            return 0
        if term == variable:
            return 1
        if isinstance(term, sp.Pow) and term.base == variable:
            return term.exp
        return 0
    
    # Find the lowest exponent present in the inverted dictionary
    exponents = [get_exponent(term) for term in coeff_dict.keys()]
    lowest_power = min(exponents)
    
    # Step 2: Renormalize by multiplying by variable**(-lowest_power)
    # This shifts the lowest power perfectly to 0
    renormalization_factor = variable ** (-lowest_power)
    normalized_expr = sp.expand(inverted_expr * renormalization_factor)
    
    return sp.Poly(normalized_expr, q)

def truncate(p, max_power):
    # truncates polynomial up to max_power
    filtered_terms = {
        deg: coeff 
        for deg, coeff in p.as_dict().items() 
        if deg[0] <= max_power
    }
    truncated_poly = sp.Poly.from_dict(filtered_terms, p.gens)
    return truncated_poly

def create_heatmap(u, v, j=10, qsub=2, trunc=50, tail=True):
    """Creates heatmap for the tail of the colored jones polynomial of a rational knot

    Args:
        u (int): Numerator.
        v (int): Denominator.
        j (int): Highest color to calculate.
        qsub (int): q specialization, 2 gives the Jones polynomial.
        trunc (int): Highest degree for each polynomial.
        tail (bool): True to compute tail, False to compute head.

    Returns:
        None (NoneType): Saves image to path.
    """
    
    polies = [] # create list of j colored jones polynomials
    for i in range(1, j+1):
        sub = substitute(u, v, i, qsub)
        if not tail:
            sub = sub.as_expr()
            sub = invert_and_renormalize(sub, q)
            if sub(0) < 0:
                sub = -sub
        polies.append(truncate(sub, trunc))

    data = []
    for idx, p in enumerate(polies):
        # p.as_dict() returns {(degree,): coeff}
        for (deg,), coeff in p.as_dict().items():
            data.append({'Poly_Index': idx + 1, 'Degree': deg, 'Coeff': int(coeff)})


    df = pd.DataFrame(data) # saves data in a data frame

    # Pivot so Degrees are columns and Poly_Index are rows
    pivot_df = df.pivot(index='Poly_Index', columns='Degree', values='Coeff').fillna(0)

    # Find the integer range of degrees for heatmap
    min_deg = pivot_df.columns.min()
    max_deg = pivot_df.columns.max()

    complete_range = range(min_deg, max_deg + 1)

    # Force Pandas to include all even columns, filling the missing ones with 0

    pivot_df_fixed = pivot_df.reindex(columns=complete_range, fill_value=0)
    even_columns = [c for c in pivot_df_fixed.columns if c % 2 == 0]
    pivot_df_even = pivot_df_fixed[even_columns]

    norm = HybridNormalize(vmin=pivot_df.values.min(), vmax=pivot_df.values.max())

    # Plot the heatmap
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(pivot_df_even, cmap=custom_cmap, norm=norm, annot=False, linewidths=0.5 ,fmt='.0f', cbar=True)
    
    # calculate shift in staircase
    shift = 0 
    _, A, _ = colored_homfly_vectors_and_quiver(u, v)
    if A.count(min(A)) == 1:
        shift = qsub - 2

    # plot staircase
    ax.plot([0,1+shift], [0,0], color='black', linewidth=2.5)
    for row_idx in range(len(pivot_df_even)):
        ax.plot([row_idx+1+shift, row_idx+2+shift], [row_idx, row_idx], color='black', linewidth=2.5)
        ax.plot([row_idx+2+shift, row_idx+2+shift], [row_idx, row_idx+1], color='black', linewidth=2.5)

    plt.title("Polynomial Coefficients Heatmap")
    plt.xlabel("Degree (q^d)")
    plt.ylabel("jth Colored Jones Polynomial")
    dir_path = Path("tails")
    dir_path.mkdir(parents=True, exist_ok=True)
    if tail:
        plt.savefig(f"tails/K_{u}_{v}_{qsub}")
    else:
        plt.savefig(f"tails/H_{u}_{v}_{qsub}")

def create_heatmap_frac(u, v, j=10, qnum=2, qdenom=1, trunc=50, tail=True):
    """Creates heatmap for the tail of the colored HOMFLY-PT polynomial of a rational knot specialized to qnum/qdenom

    Args:
        u (int): Numerator.
        v (int): Denominator.
        j (int): Highest color to calculate.
        qnum (int): specializes a - > q^qnum
        qdenom (int): specializes q -> q^-qdenom
        trunc (int): Highest degree for each polynomial.
        tail (bool): True to compute tail, False to compute head.

    Returns:
        None (NoneType): Saves image to path.
    """
    
    polies = [] # create list of j colored jones polynomials
    for i in range(1, j+1):
        sub = sub_frac(u, v, i, qnum, qdenom)
        if not tail:
            sub = sub.as_expr()
            sub = invert_and_renormalize(sub, q)
            if sub(0) < 0:
                sub = -sub
        polies.append(truncate(sub, trunc))

    data = []
    for idx, p in enumerate(polies):
        # p.as_dict() returns {(degree,): coeff}
        for (deg,), coeff in p.as_dict().items():
            data.append({'Poly_Index': idx + 1, 'Degree': deg, 'Coeff': int(coeff)})


    df = pd.DataFrame(data) # saves data in a data frame

    # Pivot so Degrees are columns and Poly_Index are rows
    pivot_df = df.pivot(index='Poly_Index', columns='Degree', values='Coeff').fillna(0)

    # Find the integer range of degrees for heatmap
    min_deg = pivot_df.columns.min()
    max_deg = pivot_df.columns.max()

    complete_range = range(min_deg, max_deg + 1)

    # Force Pandas to include all even columns, filling the missing ones with 0

    pivot_df_fixed = pivot_df.reindex(columns=complete_range, fill_value=0)
    even_columns = [c for c in pivot_df_fixed.columns if c % 2 == 0]
    pivot_df_even = pivot_df_fixed[even_columns]

    norm = HybridNormalize(vmin=pivot_df.values.min(), vmax=pivot_df.values.max())

    # Plot the heatmap
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(pivot_df_even, cmap=custom_cmap, norm=norm, annot=True, linewidths=0.5 ,fmt='.0f', cbar=True)
    
    # calculate shift in staircase
    shift = 0 
    _, A, _ = colored_homfly_vectors_and_quiver(u, v)
    #if A.count(min(A)) == 1:
    #    shift = qsub - 2

    # plot staircase
    ax.plot([0,1+shift], [0,0], color='black', linewidth=2.5)
    for row_idx in range(len(pivot_df_even)):
        ax.plot([row_idx+1+shift, row_idx+2+shift], [row_idx, row_idx], color='black', linewidth=2.5)
        ax.plot([row_idx+2+shift, row_idx+2+shift], [row_idx, row_idx+1], color='black', linewidth=2.5)

    plt.title("Polynomial Coefficients Heatmap")
    plt.xlabel("Degree (q^d)")
    plt.ylabel("jth Colored Jones Polynomial")
    dir_path = Path("tails")
    dir_path.mkdir(parents=True, exist_ok=True)
    if tail:
        plt.savefig(f"tails/K_{u}_{v}_{qnum}_{qdenom}")
    else:
        plt.savefig(f"tails/H_{u}_{v}_{qnum}_{qdenom}")

#create_heatmap_frac(5,2,j=10,qnum=1,qdenom=5,trunc=50,tail=True)