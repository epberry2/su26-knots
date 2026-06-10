from colored_homfly.colored_geo import colored_homfly_geo
import helpers.closure_formulas as cf
import helpers.tanglestate as ts
from helpers.normalize_laurent import normalize_laurent_2var
import sympy as sp

a = sp.symbols('a')
q = sp.symbols('q')

def colored_homfly_knot(num, denom, j):
    orientation = ts.get_state(num, denom)[0]
    if orientation == "RI":
        raise ValueError("tangle has RI orientation")
    tangle_homfly = colored_homfly_geo(num, denom, j)
    homfly = 0
    for k in range(len(tangle_homfly)):
        homfly += tangle_homfly[k].as_expr() * cf.cl_num(orientation, j, k)
    return normalize_laurent_2var(homfly, a, q)
