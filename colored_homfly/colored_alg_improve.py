from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd
from typing import Dict, List, Tuple, Optional

from helpers.continuedfrac import *
from helpers.plot_polynomial import *
import helpers.closure_formulas as cf
from helpers.twist_matrix import basis_index
from helpers.normalize_laurent import normalize_laurent_2var
import helpers.tanglestate as tanglestate

import sympy as sp

q, a = sp.symbols("q a")

ORIENTATIONS = ("UP", "OP", "RI")
BasisKey = Tuple[str, int]

# Sparse Laurent polynomial in q,a:
#     {(q_exp, a_exp): integer coefficient}
SparsePoly = Dict[Tuple[int, int], int]
SparseState = Dict[BasisKey, SparsePoly]


# ============================================================
# Continued fraction / word
# ============================================================

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
    cf_list = continuedfrac(num, denom)
    return cf_to_chronological_word(cf_list)


# ============================================================
# Sparse Laurent polynomial arithmetic
# ============================================================

def _clean(p: SparsePoly) -> SparsePoly:
    """Remove zero coefficients."""
    return {m: c for m, c in p.items() if c != 0}


def poly_zero() -> SparsePoly:
    return {}


def poly_one() -> SparsePoly:
    return {(0, 0): 1}


def poly_monomial(q_exp: int = 0, a_exp: int = 0, coeff: int = 1) -> SparsePoly:
    if coeff == 0:
        return {}
    return {(int(q_exp), int(a_exp)): int(coeff)}


def poly_add_into(target: SparsePoly, source: SparsePoly) -> None:
    """In-place target += source."""
    for mon, coeff in source.items():
        new_coeff = target.get(mon, 0) + coeff
        if new_coeff:
            target[mon] = new_coeff
        elif mon in target:
            del target[mon]


def poly_shift_mul(
    p: SparsePoly,
    *,
    q_shift: int = 0,
    a_shift: int = 0,
    coeff: int = 1,
) -> SparsePoly:
    """Return coeff*q^q_shift*a^a_shift*p."""
    if coeff == 0 or not p:
        return {}

    out: SparsePoly = {}
    for (qe, ae), c in p.items():
        nc = coeff * c
        if nc:
            out[(qe + q_shift, ae + a_shift)] = out.get((qe + q_shift, ae + a_shift), 0) + nc
    return _clean(out)


def poly_mul(p: SparsePoly, r: SparsePoly) -> SparsePoly:
    """Sparse multiplication. Mostly used for optional checks."""
    if not p or not r:
        return {}

    out: SparsePoly = {}
    for (q1, a1), c1 in p.items():
        for (q2, a2), c2 in r.items():
            mon = (q1 + q2, a1 + a2)
            out[mon] = out.get(mon, 0) + c1 * c2
    return _clean(out)


def sparse_to_sympy(p: SparsePoly) -> sp.Expr:
    """Convert sparse Laurent polynomial to a SymPy expression."""
    if not p:
        return sp.Integer(0)

    expr = sp.Integer(0)
    for (qe, ae), coeff in p.items():
        expr += sp.Integer(coeff) * q ** qe * a ** ae
    return expr


def sympy_to_sparse(expr: sp.Expr) -> SparsePoly:
    """
    Convert a Laurent polynomial in q,a into sparse form.

    This is not used in the main recursion. It is useful for debugging or for
    comparing with old SymPy-based code.
    """
    expr = sp.expand(expr)
    out: SparsePoly = {}

    for term in sp.Add.make_args(expr):
        powers = term.as_powers_dict()
        qe = powers.get(q, sp.Integer(0))
        ae = powers.get(a, sp.Integer(0))

        if not (qe.is_integer and ae.is_integer):
            raise ValueError(f"Non-integer exponent in term: {term}")

        coeff = sp.simplify(term / (q ** qe * a ** ae))
        if coeff.has(q) or coeff.has(a):
            raise ValueError(f"Not a Laurent monomial in q,a: {term}")

        mon = (int(qe), int(ae))
        out[mon] = out.get(mon, 0) + int(coeff)

    return _clean(out)


# ============================================================
# Cached quantum binomials as sparse q-polynomials
# ============================================================

@lru_cache(maxsize=None)
def qbinom_sparse(n: int, k: int, base_exp: int) -> Tuple[Tuple[int, int], ...]:
    """
    Sparse q-polynomial for Gaussian binomial [n choose k]_{q^base_exp}.

    Return format is immutable tuple:
        ((q_exp, coeff), ...)

    This avoids SymPy inside the twist recursion.

    Recurrence:
        [n k]_Q = [n-1 k]_Q + Q^(n-k) [n-1 k-1]_Q.
    """
    if k < 0 or k > n:
        return tuple()
    if k == 0 or k == n:
        return ((0, 1),)

    left = dict(qbinom_sparse(n - 1, k, base_exp))
    right = dict(qbinom_sparse(n - 1, k - 1, base_exp))

    shift = base_exp * (n - k)
    out: Dict[int, int] = dict(left)

    for qe, coeff in right.items():
        q_new = qe + shift
        out[q_new] = out.get(q_new, 0) + coeff

    out = {qe: c for qe, c in out.items() if c != 0}
    return tuple(sorted(out.items()))


def qbinom_as_poly(n: int, k: int, base_exp: int) -> SparsePoly:
    """
    Convert cached one-variable q-binomial into two-variable sparse polynomial.
    """
    return {(qe, 0): coeff for qe, coeff in qbinom_sparse(n, k, base_exp)}


# ============================================================
# Twist rules without matrices
# ============================================================

def twist_terms_sparse(j: int, X: str, k: int, op: str) -> List[Tuple[str, int, SparsePoly]]:
    """
    The six twist rules, but coefficients are sparse Laurent polynomials.

    Returns:
        [(target_orientation, h, coeff_poly), ...]
    """
    if not (0 <= k <= j):
        raise ValueError(f"k must satisfy 0 <= k <= j, got k={k}, j={j}")

    out: List[Tuple[str, int, SparsePoly]] = []

    if op == "T":
        if X == "UP":
            # TUP[j,k] = sum_{h=k}^j (-q)^h q^{k^2} [h k]_{q^2} UP[j,h]
            for h in range(k, j + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(h, k, 2)
                coeff = poly_shift_mul(coeff, q_shift=h + k * k, coeff=sign)
                out.append(("UP", h, coeff))

        elif X == "OP":
            # TOP[j,k] = sum_{h=k}^j (-q)^h a^k q^{k(k-2j)} [h k]_{q^2} RI[j,h]
            for h in range(k, j + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(h, k, 2)
                coeff = poly_shift_mul(
                    coeff,
                    q_shift=h + k * (k - 2 * j),
                    a_shift=k,
                    coeff=sign,
                )
                out.append(("RI", h, coeff))

        elif X == "RI":
            # TRI[j,k] = sum_{h=k}^j (-q)^h a^h q^{h(h-2j)} [h k]_{q^2} OP[j,h]
            for h in range(k, j + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(h, k, 2)
                coeff = poly_shift_mul(
                    coeff,
                    q_shift=h + h * (h - 2 * j),
                    a_shift=h,
                    coeff=sign,
                )
                out.append(("OP", h, coeff))

        else:
            raise ValueError(f"Unknown orientation: {X}")

    elif op == "R":
        if X == "UP":
            # RUP[j,k] = sum_{h=0}^k (-q)^h a^h q^{k(2j-k)-2hj}
            #            [j-h, k-h]_{q^-2} OP[j,h]
            for h in range(0, k + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(j - h, k - h, -2)
                coeff = poly_shift_mul(
                    coeff,
                    q_shift=h + k * (2 * j - k) - 2 * h * j,
                    a_shift=h,
                    coeff=sign,
                )
                out.append(("OP", h, coeff))

        elif X == "OP":
            # ROP[j,k] = sum_{h=0}^k (-q)^h a^k q^{-k^2}
            #            [j-h, k-h]_{q^-2} UP[j,h]
            for h in range(0, k + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(j - h, k - h, -2)
                coeff = poly_shift_mul(
                    coeff,
                    q_shift=h - k * k,
                    a_shift=k,
                    coeff=sign,
                )
                out.append(("UP", h, coeff))

        elif X == "RI":
            # RRI[j,k] = sum_{h=0}^k (-q)^h q^{-k(k-2j)}
            #            [j-h, k-h]_{q^-2} RI[j,h]
            for h in range(0, k + 1):
                sign = -1 if h % 2 else 1
                coeff = qbinom_as_poly(j - h, k - h, -2)
                coeff = poly_shift_mul(
                    coeff,
                    q_shift=h - k * (k - 2 * j),
                    coeff=sign,
                )
                out.append(("RI", h, coeff))

        else:
            raise ValueError(f"Unknown orientation: {X}")

    else:
        raise ValueError("op must be T or R")

    return out


def initial_state(j: int) -> SparseState:
    if j < 0:
        raise ValueError("j must be nonnegative")
    return {("UP", 0): poly_one()}


def apply_twist_state(state: SparseState, j: int, op: str) -> SparseState:
    """
    Apply one T/R twist directly to the sparse state.

    This replaces:
        v = Tmat * v
        v = Rmat * v
    and avoids full matrix multiplication.
    """
    new_state: SparseState = {}

    for (X, k), coeff_poly in state.items():
        if not coeff_poly:
            continue

        for Y, h, rule_poly in twist_terms_sparse(j, X, k, op):
            contribution = poly_mul(coeff_poly, rule_poly)
            if not contribution:
                continue

            key = (Y, h)
            if key not in new_state:
                new_state[key] = {}
            poly_add_into(new_state[key], contribution)

    return {key: poly for key, poly in new_state.items() if poly}


def evaluate_tangle_state(num: int, denom: int, j: int) -> SparseState:
    """
    Fast internal evaluator.

    Returns sparse state:
        {("UP", k): sparse_poly, ("OP", k): ..., ("RI", k): ...}

    No SymPy matrices and no simplify/factor calls are used here.
    """
    if j < 0:
        raise ValueError("j must be nonnegative")

    word = word_from_fraction(num, denom)
    state = initial_state(j)

    for op in word:
        state = apply_twist_state(state, j, op)

    return state


def state_to_vector(state: SparseState, j: int) -> sp.Matrix:
    """
    Convert sparse state to the old SymPy vector format.

    This preserves the old interface of evaluate_tangle_vector().
    """
    size = 3 * (j + 1)
    v = sp.zeros(size, 1)

    for (X, k), poly in state.items():
        v[basis_index(X, j, k), 0] = sparse_to_sympy(poly)

    return v


def initial_vector(j: int) -> sp.Matrix:
    """
    Kept for compatibility with the old code.
    """
    size = 3 * (j + 1)
    v = sp.zeros(size, 1)
    v[0, 0] = sp.Integer(1)
    return v


def evaluate_tangle_vector(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = False,
) -> sp.Matrix:
    """
    Fast replacement of the old matrix evaluator.

    The argument simplify_each_step is kept for compatibility, but is ignored
    in the fast path. The whole point is to avoid SymPy simplification during
    the recursion.

    Returns the same type as before:
        SymPy column vector of size 3*(j+1).
    """
    state = evaluate_tangle_state(num, denom, j)
    return state_to_vector(state, j)


# ============================================================
# Closure
# ============================================================

def closure_polynomial_matrix(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = False,
    final_simplify: bool = True,
) -> sp.Expr:
    """
    Same mathematical output as the old closure_polynomial_matrix(), but faster.

    Main optimization:
        - evaluate_tangle_state() avoids SymPy matrices and simplification.
        - SymPy is only used at the final closure step, because cf.cl_num()
          is already a SymPy expression.
    """
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

    state = evaluate_tangle_state(num, denom, j)
    final_orientation = tanglestate.get_state(num, denom)[0]

    value = sp.Integer(0)
    for k in range(j + 1):
        coeff_sparse = state.get((final_orientation, k), {})
        if not coeff_sparse:
            continue
        coeff = sparse_to_sympy(coeff_sparse)
        value += coeff * cf.cl_num(final_orientation, j, k)

    if final_simplify:
        return sp.expand(sp.factor(sp.simplify(value)))

    return sp.expand(value)


# ============================================================
# Plot wrapper
# ============================================================

def plot_closure_lattice_points(
    expr: sp.Expr,
    *,
    num: Optional[int] = None,
    denom: Optional[int] = None,
    j: Optional[int] = None,
    save_path: Optional[str] = None,
    show: bool = True,
    label_coefficients: bool = True,
) -> List[Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]]:
    """
    Plot the exponent points of a closure polynomial.

    The previous version referenced an undefined variable `out`.
    This keeps the same plotting behavior but fixes the title bug.
    """
    if num is not None and denom is not None and j is not None:
        title = f"Monomial support for {num}/{denom}, j={j}"
    else:
        title = "Monomial support"

    return plot_polynomial_lattice_points(
        expr,
        title=title,
        save_path=save_path,
        show=show,
        label_coefficients=label_coefficients,
    )


# ============================================================
# Optional comparison helper
# ============================================================

def compare_with_old_vector(old_vector: sp.Matrix, new_vector: sp.Matrix) -> bool:
    """
    Return True if two vectors are algebraically equal.

    Useful when testing against the old SymPy-matrix implementation.
    """
    if old_vector.shape != new_vector.shape:
        return False

    for i in range(old_vector.rows):
        if sp.simplify(old_vector[i, 0] - new_vector[i, 0]) != 0:
            return False

    return True


if __name__ == "__main__":
    num = 3
    denom = 1
    j = 2

    ev = evaluate_tangle_vector(num, denom, j)
    print(ev)

    out = closure_polynomial_matrix(num, denom, j)
    print(normalize_laurent_2var(out, a, q))

    # plot_closure_lattice_points(
    #     out,
    #     num=num,
    #     denom=denom,
    #     j=j,
    #     save_path=f"homfly_points_{num}_{denom}_j{j}.png",
    #     show=True,
    #     label_coefficients=True,
    # )
