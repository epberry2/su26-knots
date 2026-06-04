import findloops as fl
import sympy as sp
from tanglestate import get_state

q = sp.symbols("q")
a = sp.symbols("a")

def HomflyPolyGeo(num, denom):
    path, loops = fl.findpathwithloops(num, denom)
    state = get_state(num, denom)

    return