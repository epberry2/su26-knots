import helpers.quantum_nums as qnums
import sympy as sp

q = sp.symbols('q')
a = sp.symbols('a')

def cl_num(orientation, j, k):
    "Numerator closure"
    match orientation:
        case "UP":
            out_num = q ** (j**2 + k**2)
            out_denom = a ** j
            out_num *= qnums.qpochhammer(a**2 * q**(2 - 2 * j - 2 * k), q**2, j)
            out_denom *= qnums.qpochhammer(q**2, q**2, j)
            return sp.factor(sp.simplify(out_num * qnums.qbinom(j, k) / out_denom))
        case "OP":
            out_num = q ** ((j - k) ** 2)
            out_denom = a ** (j - k)
            out_num *= qnums.qpochhammer(a**2 * q**(2 - 2 * j), q**2, j - k)
            out_denom *= qnums.qpochhammer(q**2, q**2, j - k)
            return sp.factor(sp.simplify(out_num * qnums.qbinom(j, k) / out_denom))
        case _:
            raise TypeError("Orientation not UP or OP")
            

def cl_denom(orientation, j, k):
    "Warning: Denominator closure"
    return

