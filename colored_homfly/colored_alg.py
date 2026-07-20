from __future__ import annotations
from math import gcd
from typing import Dict, List, Tuple, Optional
from helpers.quantum_nums import *
from helpers.continuedfrac import *
from visualizers.plot_polynomial import *
import helpers.closure_formulas as cf
from helpers.twist_matrix import twist_matrix, basis_index
from helpers.normalize_laurent import normalize_laurent_2var
import helpers.tanglestate as tanglestate

import sympy as sp

q, a = sp.symbols("q a")

def cf_to_chronological_word(cf: List[int]) -> str:
    """
    e.g. [1, 2, 2] -> TTRRT
    """
    blocks: List[str] = []
    
    for i, power in enumerate(cf):
        if power < 0:
            raise ValueError("Only nonnegative continued fraction entries are supported.")
        op = "T" if i % 2 == 0 else "R"
        blocks.append(op * power)

    return "".join(reversed(blocks))

def word_from_fraction(num: int, denom: int) -> str:
    cf = continuedfrac(num, denom)
    word = cf_to_chronological_word(cf)
    return word

def initial_vector(j: int) -> sp.Matrix:
    size = 3 * (j + 1)
    v = sp.zeros(size, 1)
    v[0, 0] = sp.Integer(1)
    return v

def evaluate_tangle_vector(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = True,
):
    """
    Returns vector of size (j + 1)
    """
    if j < 0:
        raise ValueError("j must be nonnegative")

    word = word_from_fraction(num, denom)
    Tmat = twist_matrix("T", j)
    Rmat = twist_matrix("R", j)
    v = initial_vector(j)

    for op in word:
        if op == "T":
            v = Tmat * v
        elif op == "R":
            v = Rmat * v
        else:
            raise ValueError("word must contain only T and R")

        if simplify_each_step:
            v = v.applyfunc(lambda x: sp.factor(sp.simplify(x)))
    return v

def closure_polynomial_matrix(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = True,
) -> sp.Expr:
    if denom == 0:
        raise ZeroDivisionError("denom cannot be zero")
    if j < 0:
        raise ValueError("j must be nonnegative")
    if gcd(num, denom) != 1:
        raise ValueError("Please input a reduced fraction: gcd(num, denom) must be 1")
    if num < denom:
        raise ValueError(
            "This numerator-closure algorithm is for u/v > 1. "
            "For u/v < 1 the natural closure in the paper is denominator closure, "
            "which is not implemented here."
        )

    vec = evaluate_tangle_vector(num, denom, j, simplify_each_step=simplify_each_step)
    final_orientation = tanglestate.get_state(num, denom)[0]
    value = sp.Integer(0)
    for k in range(j + 1):
        coeff = vec[basis_index(final_orientation, j, k), 0]
        value += coeff * cf.cl_num(final_orientation, j, k)
    value = sp.expand(sp.factor(sp.simplify(value)))
    return value

def plot_closure_lattice_points(
    expr: sp.Expr,
    *,
    save_path: Optional[str] = None,
    show: bool = True,
    label_coefficients: bool = True,
) -> List[Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]]:
    """Plot the exponent points of a ClosureEvaluation output."""
    title = f"Monomial support for {out.num}/{out.denom}, j={out.j}"

    return plot_polynomial_lattice_points(
        expr,
        title=title,
        save_path=save_path,
        show=show,
        label_coefficients=label_coefficients,
    )
    
if __name__ == "__main__":
    num = 3
    denom = 1
    j = 2

    # Tangle evaluation of tau_{num/denom}.
    ev = evaluate_tangle_vector(num, denom, j)
    print(ev)

    # Numerator closure value of Cl(tau_{num/denom}).
    out = closure_polynomial_matrix(num, denom, j)
    
    print(normalize_laurent_2var(out, a, q))

    # plot_closure_lattice_points(
    #     out,
    #     save_path=f"homfly_points_{num}_{denom}_j{j}.png",
    #     show=True,
    #     label_coefficients=True,
    # )
