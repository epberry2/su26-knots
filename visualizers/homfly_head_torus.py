import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from visualizers.heatmap import HybridNormalize
import matplotlib.colors as mcolors

q, a = sp.symbols('q a')

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

def truncate_q(expr, N):
    p = sp.Poly(sp.expand(expr), q)
    return sum(c*q**e[0] for e, c in p.terms() if 0 <= e[0] <= N)

def theta_sum_terms(n_val, N):
    n_val = sp.Integer(n_val)
    poly = sp.Integer(0)
    j = 0
    while True:
        found_any = False
        for jj in ({j, -j} if j != 0 else {0}):
            exp = sp.Rational(jj * ((2*n_val+1)*jj - (2*n_val-1)), 2)
            if exp.is_integer and 0 <= exp <= N:
                poly += (-1)**jj * q**int(exp)
                found_any = True
        if not found_any and j > 0:
            break
        j += 1
    return sp.expand(poly)

def poch_trunc(x, step, N):
    result = sp.Integer(1)
    k = 0
    while step*k <= N:
        result = sp.expand(result * (1 - x * q**(step*k)))
        result = truncate_q(result, N)
        k += 1
    return result

def full_expression_dict(n_val, N):
    n_val = sp.Integer(n_val)
    theta = theta_sum_terms(n_val, N)
    a_inv_poch = poch_trunc(1/a, 1, N)
    qq_poch = poch_trunc(q, 1, N)
    inv_qq_poch = sp.series(1/qq_poch, q, 0, N+1).removeO()
    inv_qq_poch = truncate_q(inv_qq_poch, N)

    combined = sp.expand(theta * a_inv_poch * inv_qq_poch)
    combined = truncate_q(combined, N + 1)
    combined = sp.expand(combined * (1 - 1/(a*q)))
    combined = sp.expand(combined)

    a_shift = N + 2
    q_shift = 1
    cleared = sp.expand(combined * a**a_shift * q**q_shift)
    cleared = sp.Poly(cleared, a, q)

    result = {}
    for (ae, qe), coeff in cleared.terms():
        a_pow = ae - a_shift
        q_pow = qe - q_shift
        if q_pow <= N:
            result[(a_pow, q_pow)] = result.get((a_pow, q_pow), 0) + coeff

    # sp.series can leave float-ish coeffs; snap to exact integers/rationals
    result = {k: sp.nsimplify(v) for k, v in result.items() if v != 0}
    return result


def plot_heatmap(coeff_dict, title="Coefficient heatmap", save_path=None):
    """
    coeff_dict: {(a_power, q_power): coefficient}
    Rows = a_power, Columns = q_power.
    """
    max_q = max(k[1] for k in coeff_dict.keys())
    coeff_dict = {k: v for k,v in coeff_dict.items() if k[1] != max_q}
    a_powers = sorted(set(k[0] for k in coeff_dict))
    q_powers = sorted(set(k[1] for k in coeff_dict))

    a_index = {v: i for i, v in enumerate(a_powers)}
    q_index = {v: i for i, v in enumerate(q_powers)}

    grid = np.full((len(a_powers), len(q_powers)), np.nan)
    for (ap, qp), coeff in coeff_dict.items():
        grid[a_index[ap], q_index[qp]] = float(coeff)

    fig, ax = plt.subplots(figsize=(max(8, len(q_powers)*0.5),
                                     max(6, len(a_powers)*0.4)))

    norm = HybridNormalize(vmin=min(coeff_dict.values()), vmax=max(coeff_dict.values()))

    vmax = np.nanmax(np.abs(grid))
    im = ax.imshow(grid, cmap=custom_cmap, norm=norm, aspect="auto", origin="upper")

    ax.set_xticks(range(len(q_powers)))
    ax.set_xticklabels(q_powers)
    ax.set_yticks(range(len(a_powers)))
    ax.set_yticklabels(a_powers)
    ax.set_xlabel("power of q")
    ax.set_ylabel("power of a")
    ax.invert_yaxis()
    ax.set_title(title)

    for i in range(len(a_powers)):
        for j in range(len(q_powers)):
            val = grid[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{int(val)}", ha="center", va="center",
                         fontsize=7,
                         color="white" if abs(val) > vmax*0.6 else "black")

    fig.colorbar(im, ax=ax, label="coefficient")
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
    return fig


if __name__ == "__main__":
    n_val = 3
    N = 30
    coeff_dict = full_expression_dict(n_val, N)
    print(f"Number of nonzero terms: {len(coeff_dict)}")
    print("a_pow range:", min(k[0] for k in coeff_dict), "to", max(k[0] for k in coeff_dict))
    print("q_pow range:", min(k[1] for k in coeff_dict), "to", max(k[1] for k in coeff_dict))

    fig = plot_heatmap(
        coeff_dict,
        title=f"Homfly head for K_{2*n_val+1}/{1}",
        save_path=f"homfly_heads/head_{2*n_val+1}_1.png"
    )
