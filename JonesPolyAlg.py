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

def AlgBracket(u,v):
    cfrac = continuedfrac(u,v)
    cfrac.reverse()
    v = Matrix([[0],[1]])
    par = 0

    for k in cfrac:
        if par == 0:
            v = (T ** k) * v

        else:
            v = (R ** k) * v

        par = 1 - par
    return v.applyfunc(sp.expand)

print(AlgBracket(2,5))
