import sympy as sp
from typing import Dict, List, Tuple, Optional

q, a = sp.symbols("q a")

def polynomial_lattice_points(poly: sp.Expr) -> List[Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]]:
    """
    Convert an expanded Laurent polynomial in q and a into lattice points.

    Each monomial coeff*q^m*a^n becomes one point (m, n).
    The returned list contains (q_exponent, a_exponent, coeff, monomial).
    """
    expanded = sp.expand(poly)
    terms = sp.Add.make_args(expanded)
    points: List[Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]] = []

    for term in terms:
        powers = term.as_powers_dict()
        q_exp = powers.get(q, sp.Integer(0))
        a_exp = powers.get(a, sp.Integer(0))
        coeff = sp.simplify(term / (q ** q_exp * a ** a_exp))

        if coeff.has(q) or coeff.has(a):
            raise ValueError(
                "The expression contains a term that is not a Laurent monomial "
                f"in q and a: {term}"
            )

        if not (q_exp.is_integer and a_exp.is_integer):
            raise ValueError(
                f"Non-integer exponent found in term {term}: q^{q_exp}, a^{a_exp}"
            )

        points.append((int(q_exp), int(a_exp), sp.factor(coeff), term))

    points.sort(key=lambda item: (item[0], item[1], str(item[2])))
    return points

def plot_polynomial_lattice_points(
    poly: sp.Expr, *,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show: bool = True,
    label_coefficients: bool = True,
) -> List[Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]]:
    """
    Plot the monomial support of a Laurent polynomial in q and a.

    Horizontal axis: exponent of q.
    Vertical axis: exponent of a.
    Every monomial coeff*q^m*a^n is plotted as one point (m,n).
    """
    import matplotlib.pyplot as plt

    points = polynomial_lattice_points(poly)

    if not points:
        raise ValueError("Cannot plot the zero polynomial.")

    xs = [q_exp for q_exp, _, _, _ in points]
    ys = [a_exp for _, a_exp, _, _ in points]

    fig, ax = plt.subplots()
    ax.scatter(xs, ys)

    if label_coefficients:
        for q_exp, a_exp, coeff, _ in points:
            ax.annotate(str(coeff), (q_exp, a_exp), textcoords="offset points", xytext=(5, 5))

    ax.set_xlabel("q exponent")
    ax.set_ylabel("a exponent")
    ax.set_xticks(sorted(set(xs)))
    ax.set_yticks(sorted(set(ys)))
    ax.grid(True)

    if title is not None:
        ax.set_title(title)

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"\nSaved polynomial lattice plot to: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return points
