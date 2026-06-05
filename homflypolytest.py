from HomflyPolyAlg import AlgHomKnot
from GeoHomflyKnot import geoHomflyKnot
from sympy import expand, simplify, symbols, together, cancel, Poly
from numpy import gcd
from tanglestate import get_state

q = symbols('q')
a = symbols('a')

def normalize(expr):
    expr = expand(cancel(together(expr)))

    if expr == 0:
        return expr

    terms = expr.as_ordered_terms()

    min_q = min(term.as_coeff_exponent(q)[1] for term in terms)
    min_a = min(term.as_coeff_exponent(a)[1] for term in terms)

    normalized = Poly(expand(cancel(expr / (q**min_q * a**min_a))), q, a)
    if (homflyToAlexander(normalized)(1) == -1):
        normalized *= Poly(-1, q)
    return normalized

def test(u, v):
    print("")
    print(f"------------------------------ Testing {u}/{v} ------------------------------")
    print("")
    gcd_ = gcd(u, v)
    if gcd_ > 1:
        new_u = u // gcd_
        new_v = v // gcd_
        print(f"{u}/{v} is not simplified. Simplifying to {new_u}/{new_v}")
        u = new_u
        v = new_v
    if get_state(u, v)[0] == "RI":
        print(f"Error: {u}/{v} has RI orientation")
    else:
        alg = normalize(AlgHomKnot(u, v)).as_expr()
        geo = geoHomflyKnot(u, v).as_expr()
        too_large_to_print = (u + v > 50)
        if not too_large_to_print:
            print(alg)
            print("")
            print(geo)
            print("")
        print(alg == geo)
    print("")
    
def homflyToJones(homfly_polynomial):
    jones_polynomial = Poly(homfly_polynomial.as_expr().subs(a, q**2), q)
    lowest_jones_power = min(exp[0] for exp in jones_polynomial.monoms())
    jones_polynomial = Poly(jones_polynomial.as_expr() * q ** (-lowest_jones_power), q)
    return jones_polynomial

def homflyToAlexander(homfly_polynomial):
    alexander_polynomial = Poly(homfly_polynomial.as_expr().subs(a, 1), q)
    lowest_alexander_power = min(exp[0] for exp in alexander_polynomial.monoms())
    alexander_polynomial = Poly(alexander_polynomial.as_expr() * q ** (-lowest_alexander_power), q)
    return alexander_polynomial


# u, v = 201, 47

# test(u, v)