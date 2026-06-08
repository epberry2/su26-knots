from collections import defaultdict
import sympy as sp

q, a = sp.symbols("q a")


def qbinom(n, k, base):
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
    Return additive terms of expanded polynomial.
    This preserves the number of summands before later collection.
    """
    poly = sp.expand(poly)
    return list(sp.Add.make_args(poly))


def qbinom_plus_terms(n, k):
    return expanded_terms(qbinom(n, k, q**2))


def qbinom_minus_terms(n, k):
    return expanded_terms(qbinom(n, k, q**-2))


def twist_raw_terms(j, X, k, op):
    """
    Return list of raw summands:
        (new_X, new_k, coefficient)

    Important:
    we expand q-binomials into separate terms
    and DO NOT combine same target basis terms.
    """

    out = []

    if op == "T":

        if X == "UP":
            for h in range(k, j + 1):
                base = (-q)**h * q**(k*k)
                for term in qbinom_plus_terms(h, k):
                    out.append(("UP", h, sp.expand(base * term)))

        elif X == "OP":
            for h in range(k, j + 1):
                base = (-q)**h * a**k * q**(k * (k - 2*j))
                for term in qbinom_plus_terms(h, k):
                    out.append(("RI", h, sp.expand(base * term)))

        elif X == "RI":
            for h in range(k, j + 1):
                base = (-q)**h * a**h * q**(h * (h - 2*j))
                for term in qbinom_plus_terms(h, k):
                    out.append(("OP", h, sp.expand(base * term)))

    elif op == "R":

        if X == "UP":
            for h in range(0, k + 1):
                base = (-q)**h * a**h * q**(k * (2*j - k) - 2*h*j)
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("OP", h, sp.expand(base * term)))

        elif X == "OP":
            for h in range(0, k + 1):
                base = (-q)**h * a**k * q**(-k*k)
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("UP", h, sp.expand(base * term)))

        elif X == "RI":
            for h in range(0, k + 1):
                base = (-q)**h * q**(-k * (k - 2*j))
                for term in qbinom_minus_terms(j - h, k - h):
                    out.append(("RI", h, sp.expand(base * term)))

    else:
        raise ValueError("op must be T or R")

    return out


def apply_word_raw(word, j=2):
    """
    Start from <tau_{0/1}>_j = UP[j,0].

    State is a list of raw summands:
        (X, k, coefficient)

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


def print_raw_counts(state, j):
    counts = defaultdict(int)

    for X, k, coeff in state:
        counts[(X, k)] += 1

    for X in ["UP", "OP", "RI"]:
        for k in range(j + 1):
            c = counts[(X, k)]
            if c:
                print(f"{X}[{j},{k}] : {c} raw terms")


def collect_state(state, j):
    """
    Collect raw summands into actual Laurent polynomial coefficients.
    """
    collected = defaultdict(lambda: sp.Integer(0))

    for X, k, coeff in state:
        collected[(X, k)] += coeff

    for key in collected:
        collected[key] = sp.factor(sp.simplify(collected[key]))

    return collected


if __name__ == "__main__":
    j = 2

    # tau_{5/2} = T^2 R T tau_{0/1}
    # chronological word: T, R, T, T
    state = apply_word_raw("TTT", j=j)

    collected = collect_state(state, j)

    print("\nCollected coefficients:")
    for X in ["UP", "OP", "RI"]:
        for k in range(j + 1):
            coeff = collected.get((X, k), 0)
            if coeff != 0:
                print(f"{X}[{j},{k}] = {sp.factor(coeff)}")
                print(f"expanded = {sp.expand(coeff)}")
                print()