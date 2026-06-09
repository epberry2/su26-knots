from collections import defaultdict
from fractions import Fraction

import sympy as sp


q, a = sp.symbols("q a")


# ============================================================
# Continued fraction part
# ============================================================

def basic_continued_fraction(frac):
    """
    Return ordinary continued fraction [a0, a1, ..., an]
    for a positive rational number.
    """
    u = frac.numerator
    v = frac.denominator

    cf = []

    while v != 0:
        a0 = u // v
        cf.append(a0)
        u, v = v, u - a0 * v

    return cf


def cf_value(cf):
    """
    Compute the value of a continued fraction.
    Used as a safety check.
    """
    x = Fraction(cf[-1], 1)

    for a0 in reversed(cf[:-1]):
        x = a0 + Fraction(1, x)

    return x


def fix_parity(cf, want_odd):
    """
    Change the parity of the continued fraction length
    without changing its value.

    If want_odd=True, return odd length.
    If want_odd=False, return even length.
    """
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


def continued_fraction_for_tangle(num, denom):
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


def cf_to_chronological_word(cf):
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
        op = "T" if i % 2 == 0 else "R"
        blocks.append(op * power)

    word = "".join(reversed(blocks))

    return word


def word_from_fraction(num, denom):
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
# q-binomial part
# ============================================================

def qbinom(n, k, base):
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


def expanded_terms(poly):
    """
    Return additive terms of an expanded polynomial.

    This preserves raw summand counts before later collection.
    """
    poly = sp.expand(poly)
    return list(sp.Add.make_args(poly))


def qbinom_plus_terms(n, k):
    """
    [n choose k]_+ means base q^2.
    """
    return expanded_terms(qbinom(n, k, q**2))


def qbinom_minus_terms(n, k):
    """
    [n choose k]_- means base q^(-2).
    """
    return expanded_terms(qbinom(n, k, q**-2))


# ============================================================
# Twist rules
# ============================================================

def twist_raw_terms(j, X, k, op):
    """
    Return list of raw summands:

        (new_X, new_k, coefficient)

    We expand q-binomials into separate additive terms
    and do not combine same target basis terms.

    Basis type X is one of:
        "UP", "OP", "RI".

    k is between 0 and j.
    """

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

def print_raw_counts(state, j):
    """
    Print how many raw summands land in each basis web X[j,k].
    """
    counts = defaultdict(int)

    for X, k, coeff in state:
        counts[(X, k)] += 1

    for X in ["UP", "OP", "RI"]:
        for k in range(j + 1):
            c = counts[(X, k)]
            if c:
                print(f"{X}[{j},{k}] : {c} raw terms")


def apply_word_raw(word, j=2):
    """
    Start from:

        <tau_{0/1}>_j = UP[j,0].

    State is a list of raw summands:

        (X, k, coefficient).

    No collecting during the process.
    """
    state = [("UP", 0, sp.Integer(1))]

    print("start:")
    print_raw_counts(state, j)

    for op in word:
        new_state = []

        for X, k, coeff in state:
            for Y, h, twist_coeff in twist_raw_terms(j, X, k, op):
                new_state.append((Y, h, sp.expand(coeff * twist_coeff)))

        state = new_state

        print(f"\nafter {op}:")
        print_raw_counts(state, j)

    return state


def collect_state(state):
    """
    Collect raw summands into actual Laurent polynomial coefficients.
    """
    collected = defaultdict(lambda: sp.Integer(0))

    for X, k, coeff in state:
        collected[(X, k)] += coeff

    for key in collected:
        collected[key] = sp.factor(sp.simplify(collected[key]))

    return collected


def print_collected_coefficients(collected, j):
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


# ============================================================
# Main test
# ============================================================

if __name__ == "__main__":
    # Change these three numbers to test different examples.
    j = 2
    num = 5
    denom = 2

    frac, cf, word = word_from_fraction(num, denom)

    print(f"j = {j}")
    print(f"fraction = {frac}")
    print(f"continued fraction = {cf}")
    print(f"chronological T/R word = {word}")
    print()

    state = apply_word_raw(word, j=j)
    collected = collect_state(state)

    print_collected_coefficients(collected, j)