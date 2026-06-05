from sympy import Matrix, symbols
import sympy as sp
from continuedfrac import continuedfrac

q = symbols('q')
a = symbols('a')


# Order the basis as (UP(1,1), UP(1,0), OP[1,1], OP[1,0], RI[1,1], RI[1,0])

T = Matrix([
    [-q**2, -q, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, -a, -a*q**-1],
    [0, 0, 0, 0, 0, 1],
    [0, 0, -a, -q, 0, 0],
    [0, 0, 0, 1, 0, 0]
])


R = Matrix([
    [0, 0, -a, 0, 0, 0],
    [0, 0, a*q**-1, 1, 0, 0],
    [-a, 0, 0, 0, 0, 0],
    [q, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, -q**2, 0],
    [0, 0, 0, 0, 1, 1]
])


def AlgBracketHom(num,denom):
    cfrac = continuedfrac(num,denom)
    cfrac.reverse()
    v = Matrix([[0],[1],[0],[0],[0],[0]])
    par = 0

    for k in cfrac:
        if par == 0:
            for _ in range(k):
                v = T * v

        else:
            for _ in range(k):
                v = R * v

        par = 1 - par
    v = v
    v = v.applyfunc(sp.expand)

    return v

def AlgHomKnot(num,denom):
    v = AlgBracketHom(num,denom)
    x = Matrix([[((a * q**-1) - (a**-1 * q)) / (q - q**-1)],[(a - a**-1)/(q-q**-1)],[1],[(a - a**-1)/(q-q**-1)],[0],[0]])
    homfly_poly = v.dot(x)
    return homfly_poly.simplify().expand()

print(AlgHomKnot(5,2))