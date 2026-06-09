from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Tuple

import sympy as sp


q, a = sp.symbols("q a")


# ============================================================
# Data containers
# ============================================================

BasisTerm = Tuple[str, int, sp.Expr]       # (X, k, coefficient)
CollectedState = Dict[Tuple[str, int], sp.Expr]


@dataclass
class TangleEvaluation:
    """
    Output package for the rational tangle skein evaluator.

    frac:
        The rational number u/v.

    cf:
        Continued fraction with the paper's parity convention.

    word:
        Chronological T/R word applied to tau_{0/1}.

    final_orientation:
        One of "UP", "OP", "RI".

    raw_state:
        List of raw summands (X, k, coefficient), with no collection.

    collected:
        Dictionary mapping (X, k) to the collected coefficient of X[j,k].
    """
    j: int
    num: int
    denom: int
    frac: Fraction
    cf: List[int]
    word: str
    final_orientation: str
    raw_state: List[BasisTerm]
    collected: CollectedState


# ============================================================
# Continued fraction part
# ============================================================

def basic_continued_fraction(frac: Fraction) -> List[int]:
    """
    Return ordinary continued fraction [a0, a1, ..., an]
    for a positive rational number.
    """
    if frac <= 0:
        raise ValueError("Only positive rational numbers are supported.")

    u = frac.numerator
    v = frac.denominator

    cf = []

    while v != 0:
        a0 = u // v
        cf.append(a0)
        u, v = v, u - a0 * v

    return cf


def cf_value(cf: List[int]) -> Fraction:
    """
    Compute the value of a continued fraction.
    Used as a safety check.
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
    Change the parity of the continued fraction length
    without changing its value.

    If want_odd=True, return odd length.
    If want_odd=False, return even length.
    """
    if not cf:
        raise ValueError("Continued fraction cannot be empty.")

    cf = cf[:]

    is_odd = (len(cf) % 2 == 1)

    if is_odd == want_odd:
        return cf

    # Use [..., m] = [..., m - 1, 1] when m > 1.
    if cf[-1] > 1:
        cf[-1] -= 1
        cf.append(1)

    # Use [..., b, 1] = [..., b + 1].
    else:
        if len(cf) >= 2:
            cf[-2] += 1
            cf.pop()
        else:
            raise ValueError("Cannot change parity for this continued fraction.")

    return cf


def continued_fraction_for_tangle(num: int, denom: int) -> List[int]:
    """
    Paper convention:

    If num/denom >= 1, use odd length continued fraction.
    If num/denom < 1, use even length continued fraction.
    """
    frac = Fraction(num, denom)

    if frac <= 0:
        raise ValueError("Only positive rational numbers are supported.")

    cf = basic_continued_fraction(frac)

    if frac >= 1:
        cf = fix_parity(cf, want_odd=True)
    else:
        cf = fix_parity(cf, want_odd=False)

    # Safety check.
    if cf_value(cf) != frac:
        raise ValueError(
            f"Continued fraction error: {cf} gives {cf_value(cf)}, not {frac}"
        )

    return cf


def cf_to_chronological_word(cf: List[int]) -> str:
    """
    Convert continued fraction to chronological T/R word.

    If cf = [a0, a1, a2, ...], then formally:

        tau = T^a0 R^a1 T^a2 R^a3 ... tau_0.

    But our program applies operations to tau_0 in chronological order,
    so the rightmost block acts first.

    Example:
        5/2 = [2, 1, 1]

        Formal expression:
            tau_{5/2} = T^2 R T tau_0

        Chronological word:
            T R T T

        So the program uses:
            "TRTT"
    """
    blocks = []

    for i, power in enumerate(cf):
        if power < 0:
            raise ValueError("Only nonnegative continued fraction entries are supported.")

        op = "T" if i % 2 == 0 else "R"
        blocks.append(op * power)

    word = "".join(reversed(blocks))

    return word


def word_from_fraction(num: int, denom: int) -> Tuple[Fraction, List[int], str]:
    """
    Return:
        frac: Fraction object
        cf: continued fraction with correct parity
        word: chronological T/R word
    """
    frac = Fraction(num, denom)
    cf = continued_fraction_for_tangle(num, denom)
    word = cf_to_chronological_word(cf)

    return frac, cf, word


# ============================================================
# Orientation / state part
# ============================================================

def apply_orientation_word(word: str, start: str = "UP") -> str:
    """
    Track only the basis-orientation type under T/R.

    From the twist rules:
        T: UP -> UP, OP -> RI, RI -> OP
        R: UP -> OP, OP -> UP, RI -> RI

    This is useful because after applying a rational tangle word,
    all collected terms should land in the same orientation type.
    """
    orient = start

    for op in word:
        if op == "T":
            if orient == "UP":
                orient = "UP"
            elif orient == "OP":
                orient = "RI"
            elif orient == "RI":
                orient = "OP"
            else:
                raise ValueError(f"Unknown orientation: {orient}")

        elif op == "R":
            if orient == "UP":
                orient = "OP"
            elif orient == "OP":
                orient = "UP"
            elif orient == "RI":
                orient = "RI"
            else:
                raise ValueError(f"Unknown orientation: {orient}")

        else:
            raise ValueError("word must contain only T and R.")

    return orient


# ============================================================
# q-binomial part
# ============================================================

def qbinom(n: int, k: int, base: sp.Expr) -> sp.Expr:
    """
    Gaussian binomial coefficient:

        [n choose k]_base.
    """
    if k < 0 or k > n:
        return sp.Integer(0)

    num = sp.Integer(1)
    den = sp.Integer(1)

    for i in range(1, k + 1):
        num *= 1 - base ** (n - k + i)
        den *= 1 - base ** i

    return sp.simplify(num / den)


def expanded_terms(poly: sp.Expr) -> List[sp.Expr]:
    """
    Return additive terms of an expanded polynomial.

    This preserves raw summand counts before later collection.
    """
    poly = sp.expand(poly)
    return list(sp.Add.make_args(poly))


def qbinom_plus_terms(n: int, k: int) -> List[sp.Expr]:
    """
    [n choose k]_+ means base q^2.
    """
    return expanded_terms(qbinom(n, k, q**2))


def qbinom_minus_terms(n: int, k: int) -> List[sp.Expr]:
    """
    [n choose k]_- means base q^(-2).
    """
    return expanded_terms(qbinom(n, k, q**-2))


# ============================================================
# Twist rules
# ============================================================

def twist_raw_terms(j: int, X: str, k: int, op: str) -> List[BasisTerm]:
    """
    Return list of raw summands:

        (new_X, new_k, coefficient)

    We expand q-binomials into separate additive terms
    and do not combine same target basis terms.

    Basis type X is one of:
        "UP", "OP", "RI".

    k is between 0 and j.
    """
    if not (0 <= k <= j):
        raise ValueError(f"k must satisfy 0 <= k <= j, got k={k}, j={j}.")

    out = []

    if op == "T":

        if X == "UP":
            # TUP[j,k]
            for h in range(k, j + 1):
                base = (-q)**h * q**(k * k)
                for term in qbinom_plus_terms(h, k):
                    out.append(("UP", h, sp.expand(base * term)))

        elif X == "OP":
            # TOP[j,k]
            for h in range(k, j + 1):
                base = (-q)**h * a**k * q**(k * (k - 2 * j))
                for term in qbinom_plus_terms(h, k):
                    out.append(("RI", h, sp.expand(base * term)))

        elif X == "RI":
            # TRI[j,k]
            for h in range(k, j + 1):
                base = (-q)**h * a**h * q**(h * (h - 2 * j))
                for term in qbinom_plus_terms(h, k):
                    out.append(("OP", h, sp.expand(base * term)))

        else:
            raise ValueError(f"Unknown basis type: {X}")

    elif op == "R":

        if X == "UP":
            # RUP[j,k]
            for h in range(0, k + 1):
                base = (
                    (-q)**h
                    * a**h
                    * q**(k * (2 * j - k) - 2 * h * j)
                )
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("OP", h, sp.expand(base * term)))

        elif X == "OP":
            # ROP[j,k]
            for h in range(0, k + 1):
                base = (-q)**h * a**k * q**(-k * k)
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("UP", h, sp.expand(base * term)))

        elif X == "RI":
            # RRI[j,k]
            for h in range(0, k + 1):
                base = (-q)**h * q**(-k * (k - 2 * j))
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("RI", h, sp.expand(base * term)))

        else:
            raise ValueError(f"Unknown basis type: {X}")

    else:
        raise ValueError("op must be T or R")

    return out


# ============================================================
# Apply word and collect result
# ============================================================

def raw_counts(state: List[BasisTerm]) -> Dict[Tuple[str, int], int]:
    """
    Return how many raw summands land in each basis web X[j,k].
    """
    counts = defaultdict(int)

    for X, k, coeff in state:
        counts[(X, k)] += 1

    return dict(counts)


def print_raw_counts(state: List[BasisTerm], j: int) -> None:
    """
    Print how many raw summands land in each basis web X[j,k].
    """
    counts = raw_counts(state)

    for X in ["UP", "OP", "RI"]:
        for k in range(j + 1):
            c = counts.get((X, k), 0)
            if c:
                print(f"{X}[{j},{k}] : {c} raw terms")


def apply_word_raw(word: str, j: int = 2, verbose: bool = True) -> List[BasisTerm]:
    """
    Start from:

        <tau_{0/1}>_j = UP[j,0].

    State is a list of raw summands:

        (X, k, coefficient).

    No collecting during the process.
    """
    if j < 0:
        raise ValueError("j must be nonnegative.")

    state = [("UP", 0, sp.Integer(1))]

    if verbose:
        print("start:")
        print_raw_counts(state, j)

    for op in word:
        new_state = []

        for X, k, coeff in state:
            for Y, h, twist_coeff in twist_raw_terms(j, X, k, op):
                new_state.append((Y, h, sp.expand(coeff * twist_coeff)))

        state = new_state

        if verbose:
            print(f"\nafter {op}:")
            print_raw_counts(state, j)

    return state


def collect_state(state: List[BasisTerm]) -> CollectedState:
    """
    Collect raw summands into actual Laurent polynomial coefficients.
    """
    collected = defaultdict(lambda: sp.Integer(0))

    for X, k, coeff in state:
        collected[(X, k)] += coeff

    simplified = {}

    for key, coeff in collected.items():
        coeff = sp.factor(sp.simplify(coeff))
        if coeff != 0:
            simplified[key] = coeff

    return simplified


def print_collected_coefficients(collected: CollectedState, j: int) -> None:
    """
    Print nonzero collected coefficients.
    """
    print("\nCollected coefficients:")

    for X in ["UP", "OP", "RI"]:
        for k in range(j + 1):
            coeff = collected.get((X, k), sp.Integer(0))

            if coeff != 0:
                print(f"{X}[{j},{k}] = {sp.factor(coeff)}")
                print(f"expanded = {sp.expand(coeff)}")
                print()


def assert_single_orientation(collected: CollectedState, expected_orientation: str) -> None:
    """
    Safety check: after a full rational tangle word, all nonzero basis terms
    should have the final orientation predicted by the T/R orientation action.
    """
    actual_orientations = {X for (X, k), coeff in collected.items() if coeff != 0}

    if actual_orientations != {expected_orientation}:
        raise ValueError(
            "Orientation mismatch. "
            f"Expected only {expected_orientation}, got {sorted(actual_orientations)}."
        )


def evaluate_tangle(
    num: int,
    denom: int,
    j: int = 2,
    *,
    verbose: bool = False,
    check_orientation: bool = True,
) -> TangleEvaluation:
    """
    Main evaluator for Task A.

    Input:
        num, denom:
            Positive rational number num/denom.

        j:
            Color.

        verbose:
            If True, print raw term counts after each twist.

        check_orientation:
            If True, verify that collected terms land only in the predicted
            final basis orientation.

    Output:
        TangleEvaluation object containing:
            frac, cf, chronological word, final_orientation,
            raw_state, collected coefficients.
    """
    if denom == 0:
        raise ZeroDivisionError("denom cannot be zero.")
    if j < 0:
        raise ValueError("j must be nonnegative.")

    frac, cf, word = word_from_fraction(num, denom)
    final_orientation = apply_orientation_word(word)

    raw_state = apply_word_raw(word, j=j, verbose=verbose)
    collected = collect_state(raw_state)

    if check_orientation:
        assert_single_orientation(collected, final_orientation)

    return TangleEvaluation(
        j=j,
        num=num,
        denom=denom,
        frac=frac,
        cf=cf,
        word=word,
        final_orientation=final_orientation,
        raw_state=raw_state,
        collected=collected,
    )


def print_evaluation_summary(ev: TangleEvaluation) -> None:
    """
    Human-readable summary for one tangle evaluation.
    """
    print(f"j = {ev.j}")
    print(f"fraction = {ev.frac}")
    print(f"continued fraction = {ev.cf}")
    print(f"chronological T/R word = {ev.word}")
    print(f"final orientation = {ev.final_orientation}")
    print(f"number of raw terms = {len(ev.raw_state)}")
    print_collected_coefficients(ev.collected, ev.j)




# ============================================================
# j = 2 closure formulas
# ============================================================

def pochhammer_q2(x: sp.Expr, n: int) -> sp.Expr:
    """
    Pochhammer symbol (x; q^2)_n.

    (x; q^2)_n = product_{i=0}^{n-1} (1 - x q^{2i}).
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    out = sp.Integer(1)

    for i in range(n):
        out *= 1 - x * q**(2 * i)

    return sp.expand(out)


def q2_pochhammer_den(n: int) -> sp.Expr:
    """
    Denominator (q^2; q^2)_n.
    """
    return pochhammer_q2(q**2, n)


def closure_TUP_j2(k: int) -> sp.Expr:
    """
    Closure formula for Cl(T UP[2,k]), up to the framing shift
    used in Lemma 2.7 of the paper.

    General formula:

        Cl(T UP[j,k]) ~
            (-q)^k q^{2k^2} [j choose k]_+
            (a^2 q^{2-2j-2k}; q^2)_k / (q^2; q^2)_k.

    Here j = 2.
    """
    j = 2

    if k < 0 or k > j:
        return sp.Integer(0)

    coeff = (
        (-q)**k
        * q**(2 * k * k)
        * qbinom(j, k, q**2)
        * pochhammer_q2(a**2 * q**(2 - 2 * j - 2 * k), k)
        / q2_pochhammer_den(k)
    )

    return sp.factor(sp.simplify(coeff))


def closure_TRI_j2(k: int) -> sp.Expr:
    """
    Closure formula for Cl(T RI[2,k]), up to the framing shift
    used in Lemma 2.7 of the paper.

    General formula:

        Cl(T RI[j,k]) ~
            (-q)^{k-j} a^{2(k-j)} q^{2(k-j)^2} [j choose k]_+
            (a^2 q^{2-2j-2(j-k)}; q^2)_{j-k}
            / (q^2; q^2)_{j-k}.

    Here j = 2.
    """
    j = 2

    if k < 0 or k > j:
        return sp.Integer(0)

    coeff = (
        (-q)**(k - j)
        * a**(2 * (k - j))
        * q**(2 * (k - j) * (k - j))
        * qbinom(j, k, q**2)
        * pochhammer_q2(a**2 * q**(2 - 2 * j - 2 * (j - k)), j - k)
        / q2_pochhammer_den(j - k)
    )

    return sp.factor(sp.simplify(coeff))


def closure_polynomial_j2(
    num: int,
    denom: int,
    *,
    verbose: bool = False,
    show_inner: bool = True,
) -> sp.Expr:
    """
    Compute the j = 2 numerator closure polynomial for the rational
    knot/link represented as Cl(tau_{num/denom}), using the paper's convention:

        P^{V_2}_{num/denom}(q,a) = Cl(< T tau_{(num-denom)/denom} >_2).

    Important:
        This uses the closure formulas from Lemma 2.7, which are stated
        up to framing shift in the paper.

    Input requirements:
        num/denom > 1.

    Output:
        A SymPy expression in q and a.

    If num is odd, the closure is a rational knot K_{num/denom}.
    If num is even, the closure is a 2-component rational link L_{num/denom}.
    """
    frac = Fraction(num, denom)

    if frac <= 1:
        raise ValueError(
            "This numerator-closure routine expects num/denom > 1. "
            "For fractions < 1, the paper uses denominator closure instead."
        )

    # Undo the final top twist T.
    inner_num = num - denom
    inner_denom = denom

    inner = evaluate_tangle(
        inner_num,
        inner_denom,
        j=2,
        verbose=verbose,
        check_orientation=True,
    )

    if inner.final_orientation == "UP":
        basis_type = "UP"
        closure_function = closure_TUP_j2
    elif inner.final_orientation == "RI":
        basis_type = "RI"
        closure_function = closure_TRI_j2
    else:
        raise ValueError(
            "After undoing the final top twist, the inner tangle should have "
            f"orientation UP or RI, but got {inner.final_orientation}."
        )

    result = sp.Integer(0)

    for k in range(3):
        coeff = inner.collected.get((basis_type, k), sp.Integer(0))
        result += coeff * closure_function(k)

    result = sp.factor(sp.simplify(result))

    if show_inner:
        link_type = "knot" if num % 2 == 1 else "2-component link"
        print("\nClosure computation, j = 2")
        print(f"input fraction = {frac}")
        print(f"closure type = {link_type}")
        print(f"compute as Cl(< T tau_{{({num}-{denom})/{denom}}} >_2)")
        print(f"inner fraction = {inner.frac}")
        print(f"inner continued fraction = {inner.cf}")
        print(f"inner chronological T/R word = {inner.word}")
        print(f"inner orientation = {inner.final_orientation}")
        print(f"closure basis = T{basis_type}[2,k]")
        print("\nInner collected coefficients:")
        print_collected_coefficients(inner.collected, 2)
        print("Closure value, up to framing shift:")
        print(result)
        print("expanded =")
        print(sp.expand(result))

    return result


# ============================================================
# Small test suite
# ============================================================

def smoke_test() -> None:
    """
    Quick internal tests.

    These are not proof of mathematical correctness, but they catch:
        - bad continued fraction parity,
        - bad chronological word construction,
        - orientation inconsistency,
        - syntax/runtime errors.
    """
    examples = [
        (3, 1, "TTT", "UP"),
        (5, 2, "TRTT", "OP"),
        (4, 3, "TRRT", "UP"),
        (8, 3, "TTRTT", "OP"),
        (3, 2, "TRT", "RI"),
    ]

    for num, denom, expected_word, expected_orientation in examples:
        ev = evaluate_tangle(num, denom, j=2, verbose=False)
        assert ev.word == expected_word, (num, denom, ev.word, expected_word)
        assert ev.final_orientation == expected_orientation, (
            num,
            denom,
            ev.final_orientation,
            expected_orientation,
        )

    print("smoke_test passed.")


# ============================================================
# Main test
# ============================================================

if __name__ == "__main__":
    # Change these two numbers to test different numerator closures.
    # For example:
    #   5/2 gives the figure-8 knot K_{5/2}.
    #   4/3 gives a 2-component rational link L_{4/3}.
    num = 5
    denom = 2

    # First: ordinary tangle evaluation for tau_{num/denom}.
    ev = evaluate_tangle(num, denom, j=2, verbose=True)
    print()
    print_evaluation_summary(ev)

    # Second: j = 2 numerator closure value for Cl(tau_{num/denom}).
    closure_polynomial_j2(num, denom, verbose=False, show_inner=True)

    # Uncomment this if you want quick sanity tests.
    # smoke_test()
