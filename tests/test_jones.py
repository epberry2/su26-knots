import pandas as pd
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
from quivers.evaluate_quiver import colored_jones_polynomial


two_bridge = {
    (3,1) : (3,1),
    (4,1) : (5,2),
    (5,1) : (5,1),
    (5,2) : (7,3),
    (6,1) : (9,4),
    (6,2) : (11,4), # (11, 3)
    (6,3) : (13,8), # (13, 5)
    (7,1) : (7,1),
    (7,2) : (11,5),
    (7,3) : (13,3),
    (7,4) : (15,11),
    (7,5) : (17,5),
    (7,6) : (19,7),
    (7,7) : (21,13) # (21, 6)
}

q = sp.symbols('q')

def normalize_mathematica(expr):
    expr = sp.denom(expr) * expr

    poly = sp.Poly(expr, q)

    lowest_deg = min(poly.as_dict())[0]
    poly = poly.as_expr() * q ** -lowest_deg
    poly = poly.subs(q, q**2)

    return sp.Poly(poly, q)

def invert_poly(poly):
    highest_deg = max(poly.as_dict())[0]
    poly = poly.as_expr().subs(q,q**-1) * q**highest_deg
    return sp.Poly(poly, q)


def test_colored_jones():
    dfs = [0, 0, 0]
    dfs[0] = pd.read_csv("tests/jones_polynomials.csv")
    dfs[0]["sympy"] = dfs[0]["jones_polynomial"].apply(parse_mathematica)
    dfs[1] = pd.read_csv("tests/jones2_polynomials.csv")
    dfs[1]["sympy"] = dfs[1]["jones2_polynomial"].apply(parse_mathematica)
    dfs[2] = pd.read_csv("tests/jones3_polynomials.csv")
    dfs[2]["sympy"] = dfs[2]["jones3_polynomial"].apply(parse_mathematica)

    correct = 0
    total = (len(dfs[0]) - 1) * 3
    for i in range(1,len(dfs[0])):
        for j in range(3):
            row = dfs[j].loc[i]
            u, v = int(row["knot"][5]), int(row["knot"][8])
            x, y = two_bridge[(u, v)]
            mathematica_poly = normalize_mathematica(row["sympy"])
            quiver_poly = colored_jones_polynomial(x, y, j + 1)
            eq = quiver_poly == mathematica_poly or -quiver_poly == mathematica_poly or invert_poly(quiver_poly) == mathematica_poly or -invert_poly(quiver_poly) == mathematica_poly
            print(f"{(x,y,j+1)}")
            print(f"{eq}")
            if not eq:
                print(quiver_poly)
                print(mathematica_poly)
            else:
                correct += 1
            print()
    print(f"{correct} polynomials correct out of {total} ")

test_colored_jones()
