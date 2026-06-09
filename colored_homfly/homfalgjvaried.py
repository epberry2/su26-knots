"""
homflag_matrix_general_j.py

General-j matrix implementation for antisymmetric-colored HOMFLY-PT rational
tangle closures, following the conventions used in Jonathan A. Higgins,
"A Geometric Approach to the Links-Quivers Correspondence II: Rational Links".

What this file computes:
    1. The j-colored skein module evaluation <tau_{u/v}>_j of a positive
       rational tangle tau_{u/v}.
    2. The numerator closure value

           P^{V_j}_{u/v}(q,a) = Cl(< T tau_{(u-v)/v} >_j)

       for u/v > 1, using the paper's Lemma 2.7 closure formulas.

Important convention:
    The closure formulas used here are the paper's formulas, i.e. up to
    framing shift. If your group uses a different normalization, an overall
    monomial factor may be needed.

The matrix size depends on j:
    There are three orientation sectors UP, OP, RI, each with k = 0,...,j.
    Hence the full twist matrices T and R have size 3*(j+1) by 3*(j+1).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Dict, List, Tuple, Optional

import sympy as sp


q, a = sp.symbols("q a")

ORIENTATIONS = ("UP", "OP", "RI")
BasisKey = Tuple[str, int]
CollectedState = Dict[BasisKey, sp.Expr]


@dataclass
class MatrixTangleEvaluation:
    j: int
    num: int
    denom: int
    frac: Fraction
    cf: List[int]
    word: str
    final_orientation: str
    matrix_size: int
    vector: sp.Matrix
    collected: CollectedState


@dataclass
class ClosureEvaluation:
    j: int
    num: int
    denom: int
    frac: Fraction
    closure_type: str
    inner_num: int
    inner_denom: int
    inner: MatrixTangleEvaluation
    closure_basis: str
    value: sp.Expr
    expanded_value: sp.Expr


# ============================================================
# Continued fraction part
# ============================================================

def basic_continued_fraction(frac: Fraction) -> List[int]:
    """
    Ordinary continued fraction [a0, a1, ..., an] for a positive rational.
    """
    if frac <= 0:
        raise ValueError("Only positive rational numbers are supported.")

    u = frac.numerator
    v = frac.denominator
    cf: List[int] = []

    while v != 0:
        a0 = u // v
        cf.append(a0)
        u, v = v, u - a0 * v

    return cf


def cf_value(cf: List[int]) -> Fraction:
    """
    Evaluate a finite continued fraction.
    """
    if not cf:
        raise ValueError("Continued fraction cannot be empty.")

    x = Fraction(cf[-1], 1)

    for a0 in reversed(cf[:-1]):
        if x == 0:
            raise ZeroDivisionError("Invalid continued fraction with zero tail.")
        x = a0 + Fraction(1, x)

    return x


def fix_parity(cf: List[int], want_odd: bool) -> List[int]:
    """
    Change parity of continued fraction length without changing its value.

    Uses:
        [..., m] = [..., m - 1, 1]  if m > 1
        [..., b, 1] = [..., b + 1]
    """
    if not cf:
        raise ValueError("Continued fraction cannot be empty.")

    cf = cf[:]
    is_odd = (len(cf) % 2 == 1)

    if is_odd == want_odd:
        return cf

    if cf[-1] > 1:
        cf[-1] -= 1
        cf.append(1)
    else:
        if len(cf) < 2:
            raise ValueError("Cannot change parity for this continued fraction.")
        cf[-2] += 1
        cf.pop()

    return cf


def continued_fraction_for_tangle(num: int, denom: int) -> List[int]:
    """
    Paper convention:
        u/v >= 1  -> odd length continued fraction
        0<u/v<1  -> even length continued fraction
    """
    if denom == 0:
        raise ZeroDivisionError("denom cannot be zero.")

    frac = Fraction(num, denom)

    if frac <= 0:
        raise ValueError("Only positive rational numbers are supported.")

    cf = basic_continued_fraction(frac)
    cf = fix_parity(cf, want_odd=(frac >= 1))

    if cf_value(cf) != frac:
        raise ValueError(
            f"Continued fraction error: {cf} gives {cf_value(cf)}, not {frac}"
        )

    return cf


def cf_to_chronological_word(cf: List[int]) -> str:
    """
    If formally

        tau = T^a0 R^a1 T^a2 ... tau_0,

    then program execution starts at tau_0, so the chronological word is
    the reversed block word.
    """
    blocks: List[str] = []

    for i, power in enumerate(cf):
        if power < 0:
            raise ValueError("Only nonnegative continued fraction entries are supported.")

        op = "T" if i % 2 == 0 else "R"
        blocks.append(op * power)

    return "".join(reversed(blocks))


def word_from_fraction(num: int, denom: int) -> Tuple[Fraction, List[int], str]:
    frac = Fraction(num, denom)
    cf = continued_fraction_for_tangle(num, denom)
    word = cf_to_chronological_word(cf)
    return frac, cf, word


# ============================================================
# Basis indexing and orientation tracking
# ============================================================

def basis_index(j: int, X: str, k: int) -> int:
    """
    Convert basis key X[j,k] into a matrix/vector index.
    """
    if X not in ORIENTATIONS:
        raise ValueError(f"Unknown orientation {X}")

    if not (0 <= k <= j):
        raise ValueError(f"k must satisfy 0 <= k <= j, got k={k}, j={j}")

    return ORIENTATIONS.index(X) * (j + 1) + k


def index_basis(j: int, idx: int) -> BasisKey:
    """
    Convert matrix/vector index back to basis key X[j,k].
    """
    n = j + 1

    if not (0 <= idx < 3 * n):
        raise ValueError(f"index out of range: {idx}")

    X = ORIENTATIONS[idx // n]
    k = idx % n

    return X, k


def apply_orientation_word(word: str, start: str = "UP") -> str:
    """
    Orientation action read from the six twist rules:

        T: UP -> UP, OP -> RI, RI -> OP
        R: UP -> OP, OP -> UP, RI -> RI
    """
    orient = start

    for op in word:
        if op == "T":
            orient = {
                "UP": "UP",
                "OP": "RI",
                "RI": "OP",
            }[orient]

        elif op == "R":
            orient = {
                "UP": "OP",
                "OP": "UP",
                "RI": "RI",
            }[orient]

        else:
            raise ValueError("word must contain only T and R")

    return orient


# ============================================================
# Quantum algebra helpers
# ============================================================

def qbinom(n: int, k: int, base: sp.Expr) -> sp.Expr:
    """
    Gaussian binomial coefficient [n choose k]_base.
    """
    if k < 0 or k > n:
        return sp.Integer(0)

    num = sp.Integer(1)
    den = sp.Integer(1)

    for i in range(1, k + 1):
        num *= 1 - base ** (n - k + i)
        den *= 1 - base ** i

    return sp.factor(sp.simplify(num / den))


def qpochhammer(x: sp.Expr, base: sp.Expr, n: int) -> sp.Expr:
    """
    (x; base)_n = product_{i=0}^{n-1} (1 - x base^i).
    """
    if n < 0:
        raise ValueError("Pochhammer length must be nonnegative")

    out = sp.Integer(1)

    for i in range(n):
        out *= 1 - x * base**i

    return sp.factor(sp.simplify(out))


# ============================================================
# Twist coefficients and matrices
# ============================================================

def twist_terms_unexpanded(
    j: int,
    X: str,
    k: int,
    op: str,
) -> List[Tuple[str, int, sp.Expr]]:
    """
    The six twist rules.

    Return list of:

        (Y, h, coeff)

    meaning:

        op applied to X[j,k] contains coeff * Y[j,h].

    This version keeps the quantum binomial unexpanded, which is better
    for larger j.
    """
    if not (0 <= k <= j):
        raise ValueError(f"k must satisfy 0 <= k <= j, got k={k}, j={j}")

    out: List[Tuple[str, int, sp.Expr]] = []

    if op == "T":

        if X == "UP":
            # TUP[j,k]
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** h
                    * q ** (k * k)
                    * qbinom(h, k, q**2)
                )
                out.append(("UP", h, sp.factor(sp.simplify(coeff))))

        elif X == "OP":
            # TOP[j,k]
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** h
                    * a ** k
                    * q ** (k * (k - 2 * j))
                    * qbinom(h, k, q**2)
                )
                out.append(("RI", h, sp.factor(sp.simplify(coeff))))

        elif X == "RI":
            # TRI[j,k]
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** h
                    * a ** h
                    * q ** (h * (h - 2 * j))
                    * qbinom(h, k, q**2)
                )
                out.append(("OP", h, sp.factor(sp.simplify(coeff))))

        else:
            raise ValueError(f"Unknown basis type: {X}")

    elif op == "R":

        if X == "UP":
            # RUP[j,k]
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** h
                    * a ** h
                    * q ** (k * (2 * j - k) - 2 * h * j)
                    * qbinom(j - h, k - h, q**-2)
                )
                out.append(("OP", h, sp.factor(sp.simplify(coeff))))

        elif X == "OP":
            # ROP[j,k]
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** h
                    * a ** k
                    * q ** (-k * k)
                    * qbinom(j - h, k - h, q**-2)
                )
                out.append(("UP", h, sp.factor(sp.simplify(coeff))))

        elif X == "RI":
            # RRI[j,k]
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** h
                    * q ** (-k * (k - 2 * j))
                    * qbinom(j - h, k - h, q**-2)
                )
                out.append(("RI", h, sp.factor(sp.simplify(coeff))))

        else:
            raise ValueError(f"Unknown basis type: {X}")

    else:
        raise ValueError("op must be T or R")

    return out


def twist_matrix(j: int, op: str) -> sp.Matrix:
    """
    Full 3(j+1) by 3(j+1) matrix for T or R.

    Columns are source basis vectors X[j,k].
    Rows are target basis vectors Y[j,h].
    """
    if j < 0:
        raise ValueError("j must be nonnegative")

    if op not in {"T", "R"}:
        raise ValueError("op must be T or R")

    size = 3 * (j + 1)
    M = sp.zeros(size, size)

    for X in ORIENTATIONS:
        for k in range(j + 1):
            col = basis_index(j, X, k)

            for Y, h, coeff in twist_terms_unexpanded(j, X, k, op):
                row = basis_index(j, Y, h)
                M[row, col] += coeff

    return M


def initial_vector(j: int) -> sp.Matrix:
    """
    Vector for <tau_{0/1}>_j = UP[j,0].
    """
    size = 3 * (j + 1)
    v = sp.zeros(size, 1)
    v[basis_index(j, "UP", 0), 0] = sp.Integer(1)
    return v


def vector_to_collected(
    j: int,
    v: sp.Matrix,
    *,
    simplify: bool = True,
) -> CollectedState:
    """
    Convert a matrix vector to {(X,k): coefficient}.
    """
    out: CollectedState = {}

    for idx in range(3 * (j + 1)):
        coeff = v[idx, 0]

        if simplify:
            coeff = sp.factor(sp.simplify(coeff))

        if coeff != 0:
            out[index_basis(j, idx)] = coeff

    return out


def collect_orientation(collected: CollectedState) -> Optional[str]:
    """
    Check which orientation sector contains nonzero terms.
    """
    actual = {X for (X, k), coeff in collected.items() if coeff != 0}

    if len(actual) == 0:
        return None

    if len(actual) == 1:
        return next(iter(actual))

    raise ValueError(f"More than one orientation present: {sorted(actual)}")


def evaluate_tangle_matrix(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = True,
    verbose: bool = False,
) -> MatrixTangleEvaluation:
    """
    Matrix version of <tau_{u/v}>_j.

    This is the main general-j tangle evaluator. It uses matrices of size
    3*(j+1), so the matrix size changes automatically with j.
    """
    if j < 0:
        raise ValueError("j must be nonnegative")

    frac, cf, word = word_from_fraction(num, denom)
    expected_orientation = apply_orientation_word(word)

    Tmat = twist_matrix(j, "T")
    Rmat = twist_matrix(j, "R")
    v = initial_vector(j)

    if verbose:
        print("start:")
        print_vector_summary(j, v)

    for op in word:
        if op == "T":
            v = Tmat * v
        elif op == "R":
            v = Rmat * v
        else:
            raise ValueError("word must contain only T and R")

        if simplify_each_step:
            v = v.applyfunc(lambda x: sp.factor(sp.simplify(x)))

        if verbose:
            print(f"\nafter {op}:")
            print_vector_summary(j, v)

    collected = vector_to_collected(j, v, simplify=True)
    actual_orientation = collect_orientation(collected)

    if actual_orientation != expected_orientation:
        raise ValueError(
            f"Orientation mismatch: expected {expected_orientation}, got {actual_orientation}"
        )

    return MatrixTangleEvaluation(
        j=j,
        num=num,
        denom=denom,
        frac=frac,
        cf=cf,
        word=word,
        final_orientation=expected_orientation,
        matrix_size=3 * (j + 1),
        vector=v,
        collected=collected,
    )


# ============================================================
# Closure formulas, arbitrary j
# ============================================================

def closure_TUP(j: int, k: int) -> sp.Expr:
    """
    Closure formula for Cl(T UP[j,k]), up to framing shift.
    """
    if k < 0 or k > j:
        return sp.Integer(0)

    expr = (
        (-q) ** k
        * q ** (2 * k * k)
        * qbinom(j, k, q**2)
        * qpochhammer(a**2 * q ** (2 - 2 * j - 2 * k), q**2, k)
        / qpochhammer(q**2, q**2, k)
    )

    return sp.factor(sp.simplify(expr))


def closure_TRI(j: int, k: int) -> sp.Expr:
    """
    Closure formula for Cl(T RI[j,k]), up to framing shift.
    """
    if k < 0 or k > j:
        return sp.Integer(0)

    r = j - k

    expr = (
        (-q) ** (k - j)
        * a ** (2 * (k - j))
        * q ** (2 * (k - j) ** 2)
        * qbinom(j, k, q**2)
        * qpochhammer(a**2 * q ** (2 - 2 * j - 2 * r), q**2, r)
        / qpochhammer(q**2, q**2, r)
    )

    return sp.factor(sp.simplify(expr))


def closure_polynomial_matrix(
    num: int,
    denom: int,
    j: int,
    *,
    simplify_each_step: bool = True,
    verbose: bool = False,
    show_inner: bool = True,
) -> ClosureEvaluation:
    """
    Compute the general-j numerator closure value

        P^{V_j}_{u/v}(q,a) = Cl(< T tau_{(u-v)/v} >_j)

    for u/v > 1, using matrix tangle evaluation and closure formulas.

    The output is up to framing shift.
    """
    if denom == 0:
        raise ZeroDivisionError("denom cannot be zero")

    if j < 0:
        raise ValueError("j must be nonnegative")

    if gcd(num, denom) != 1:
        raise ValueError("Please input a reduced fraction: gcd(num, denom) must be 1")

    frac = Fraction(num, denom)

    if frac <= 1:
        raise ValueError(
            "This numerator-closure algorithm is for u/v > 1. "
            "For u/v < 1 the natural closure in the paper is denominator closure, "
            "which is not implemented here."
        )

    inner_num = num - denom
    inner_denom = denom

    inner = evaluate_tangle_matrix(
        inner_num,
        inner_denom,
        j,
        simplify_each_step=simplify_each_step,
        verbose=verbose,
    )

    if inner.final_orientation == "UP":
        basis = "UP"
        closure_func = closure_TUP

    elif inner.final_orientation == "RI":
        basis = "RI"
        closure_func = closure_TRI

    else:
        raise ValueError(
            "Unexpected inner orientation for numerator closure. "
            f"Expected UP or RI after removing final T, got {inner.final_orientation}."
        )

    value = sp.Integer(0)

    for k in range(j + 1):
        coeff = inner.collected.get((basis, k), sp.Integer(0))

        if coeff != 0:
            value += coeff * closure_func(j, k)

    value = sp.factor(sp.simplify(value))
    expanded_value = sp.expand(value)

    closure_type = "knot" if num % 2 == 1 else "2-component link"

    out = ClosureEvaluation(
        j=j,
        num=num,
        denom=denom,
        frac=frac,
        closure_type=closure_type,
        inner_num=inner_num,
        inner_denom=inner_denom,
        inner=inner,
        closure_basis=f"T{basis}[{j},k]",
        value=value,
        expanded_value=expanded_value,
    )

    if show_inner:
        print_closure_summary(out)

    return out


# ============================================================
# Printing helpers
# ============================================================

def print_vector_summary(j: int, v: sp.Matrix) -> None:
    """
    Print which k-values are nonzero in each orientation sector.
    """
    collected = vector_to_collected(j, v, simplify=False)

    for X in ORIENTATIONS:
        nonzero = []

        for k in range(j + 1):
            coeff = collected.get((X, k), sp.Integer(0))

            if coeff != 0:
                nonzero.append(k)

        if nonzero:
            print(f"{X}[{j},k] nonzero k = {nonzero}")


def print_collected_coefficients(collected: CollectedState, j: int) -> None:
    """
    Print nonzero collected coefficients.
    """
    print("\nCollected coefficients:")

    for X in ORIENTATIONS:
        for k in range(j + 1):
            coeff = collected.get((X, k), sp.Integer(0))

            if coeff != 0:
                print(f"{X}[{j},{k}] = {sp.factor(coeff)}")
                print(f"expanded = {sp.expand(coeff)}")
                print()


def print_tangle_summary(ev: MatrixTangleEvaluation) -> None:
    """
    Print a tangle evaluation summary.
    """
    print(f"j = {ev.j}")
    print(f"fraction = {ev.frac}")
    print(f"continued fraction = {ev.cf}")
    print(f"chronological T/R word = {ev.word}")
    print(f"matrix size = {ev.matrix_size} x {ev.matrix_size}")
    print(f"final orientation = {ev.final_orientation}")
    print_collected_coefficients(ev.collected, ev.j)


def print_closure_summary(out: ClosureEvaluation) -> None:
    """
    Print a closure computation summary.
    """
    print("\nClosure computation")
    print(f"j = {out.j}")
    print(f"input fraction = {out.frac}")
    print(f"closure type = {out.closure_type}")
    print(f"compute as Cl(< T tau_{{({out.num}-{out.denom})/{out.denom}}} >_{out.j})")
    print(f"inner fraction = {out.inner.frac}")
    print(f"inner continued fraction = {out.inner.cf}")
    print(f"inner chronological T/R word = {out.inner.word}")
    print(f"inner matrix size = {out.inner.matrix_size} x {out.inner.matrix_size}")
    print(f"inner orientation = {out.inner.final_orientation}")
    print(f"closure basis = {out.closure_basis}")

    print("\nInner collected coefficients:")
    print_collected_coefficients(out.inner.collected, out.j)

    print("\nClosure value, up to framing shift:")
    print(out.value)

    print("\nExpanded:")
    print(out.expanded_value)


# ============================================================
# Tests and examples
# ============================================================

def smoke_test() -> None:
    """
    Basic tests for continued fraction, word, orientation, matrix dimensions,
    and representative closure computations.
    """
    examples = [
        (3, 1, "TTT", "UP"),
        (5, 2, "TRTT", "OP"),
        (4, 3, "TRRT", "UP"),
        (8, 3, "TTRTT", "OP"),
        (3, 2, "TRT", "RI"),
    ]

    for j in [1, 2, 3, 4]:
        for num, denom, expected_word, expected_orientation in examples:
            ev = evaluate_tangle_matrix(num, denom, j)

            assert ev.word == expected_word, (
                j,
                num,
                denom,
                ev.word,
                expected_word,
            )

            assert ev.final_orientation == expected_orientation, (
                j,
                num,
                denom,
                ev.final_orientation,
                expected_orientation,
            )

            assert ev.matrix_size == 3 * (j + 1)

    for j in [1, 2, 3]:
        for num, denom in [(3, 1), (5, 2), (7, 3), (4, 3), (8, 3)]:
            _ = closure_polynomial_matrix(num, denom, j, show_inner=False)

    print("smoke_test passed.")


if __name__ == "__main__":
    # Change these three numbers to test.
    num = 3
    denom = 1
    j = 3

    # Tangle evaluation of tau_{num/denom}.
    ev = evaluate_tangle_matrix(num, denom, j, verbose=False)
    print_tangle_summary(ev)

    # Numerator closure value of Cl(tau_{num/denom}).
    closure_polynomial_matrix(num, denom, j, verbose=False, show_inner=True)

    # Uncomment for sanity checks.
    # smoke_test()