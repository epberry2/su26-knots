import numpy
from fractions import Fraction

def continuedfrac(a, b): 
    X = []
    num = a
    denom = b
    if num >= b:
        c = num // denom
        X.append(c)
        num -= (c * denom)
    parity = 0
    while (num != 0):
        num, denom = denom, num
        if parity == 0:
            num -= 1
        c = num // denom
        X.append(c)
        if parity == 0:
            num += 1
        num -= (c * denom)
        parity = 1 - parity
    return X

def cf_value(cf):
    if not cf:
        raise ValueError("Continued fraction cannot be empty.")
    x = Fraction(cf[-1], 1)
    for a0 in reversed(cf[:-1]):
        if x == 0:
            raise ZeroDivisionError("Continued fraction cannot contain zero.")
        x = a0 + Fraction(1, x)
    return x