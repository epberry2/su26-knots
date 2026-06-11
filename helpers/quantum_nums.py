import sympy as sp

a = sp.symbols('a')
q = sp.symbols('q')


def squiggle(n):
    return (a * q ** (n) - a ** (-1) * q ** (-n)) / (q - q ** (-1))

def quantum(N):
    return (q ** N - q ** (-N)) / (q - q **(-1))

def qbinom(n: int, k: int) -> sp.Expr:
    """
    return shifted binomial
    """
    num = qpochhammer(q**2, q**2, n)
    denom = qpochhammer(q**2, q**2, k) * qpochhammer(q**2, q**2, n-k)
    return sp.factor(sp.simplify(num / denom))

def qmultinom(n: int, k: tuple[int, ...]) -> sp.Expr:
    """ 
    returns shifted multinomial
    """
    num = qpochhammer(q**2, q**2, n)
    denom = 1
    for i in k:
        denom *= qpochhammer(q**2, q**2, i)
    return sp.factor(sp.simplify(num / denom))


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

