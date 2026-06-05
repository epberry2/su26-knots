from sympy import Matrix, symbols
import sympy as sp
from continuedfrac import continuedfrac

q = symbols('q')

T = Matrix([
    [-q**2, -q],
    [0, 1]
])

R = Matrix([
    [1, 0],
    [-q**-1,-q**-2]
])

def AlgBracket(num,denom):
    cfrac = continuedfrac(num,denom)
    cfrac.reverse()
    v = Matrix([[0],[1]])
    par = 0

    for k in cfrac:
        if par == 0:
            for _ in range(k):
                v = T * v

        else:
            for _ in range(k):
                v = R * v

        par = 1 - par
    v = v.applyfunc(sp.expand)
    powers = []
    for term in v[0].as_coefficients_dict().keys():
        if term == 1:  # It's a constant term (x^0)
            powers.append(0)
        elif term.is_Pow:
            powers.append(term.exp)
        elif term == q:
            powers.append(1)

    for term in v[1].as_coefficients_dict().keys():
        if term == 1:  # It's a constant term (x^0)
            powers.append(0)
        elif term.is_Pow:
            powers.append(term.exp)
        elif term == q:
            powers.append(1)        

    min_power = min(powers)
    #min_power=0
    v = (q**(1-min_power) * v).applyfunc(sp.expand)
    
    return (sp.Poly(v[0], q), sp.Poly(v[1], q))

def AlgJonesKnot(u,v):
    x, y = AlgBracket(u,v)
    y = y.as_expr() * (q + q**-1)
    jones_poly = sp.Poly((x.as_expr() + y).expand(),q)
    if (jones_poly(1) == -1):
        jones_poly *= sp.Poly(-1, q)
    return jones_poly

print(AlgBracket(5,1))