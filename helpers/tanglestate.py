from fractions import Fraction

def continuedfrac(a, b): 
    """Computes the continued fraction expansion of a/b."""
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

class TangleState():
    """Keeps track of the orientation of a tangle.
    
    Orientation is one of UP, OP, RI.
    The points are a permutation of Y, X-, X+.
    """
    def __init__(self):
        self.orient = "UP"
        self.points = ("Y", "X-", "X+")

    def set_points(self, points):
        self.points = points
        if self.points[0] == "X+":
            self.orient = "RI"
        if self.points[1] == "X+":
            self.orient = "OP"
        if self.points[2] == "X+":
            self.orient = "UP"
    
    def t_twist(self):
        """Updates the state after a top twist """
        self.points = (self.points[1], self.points[0], self.points[2])
        if self.points[0] == "X+":
            self.orient = "RI"
        if self.points[1] == "X+":
            self.orient = "OP"
        if self.points[2] == "X+":
            self.orient = "UP"
    
    def r_twist(self):
        """Updates the state after a right twist"""
        self.points = (self.points[0], self.points[2], self.points[1])
        if self.points[0] == "X+":
            self.orient = "RI"
        if self.points[1] == "X+":
            self.orient = "OP"
        if self.points[2] == "X+":
            self.orient = "UP"

# Returns the permutation of [Y, X-, X+] given tau_(u/v)
def get_state(u,v):
    """Computes the state of K_u/v.

    Returns (orientation, points).
    """
    x = continuedfrac(u,v)
    tangle = TangleState()
    x.reverse()
 
    par = 0
    for k in x:
        if par == 0:
            if k % 2 != 0:
                tangle.t_twist()

        else:
            if k % 2 != 0:
                tangle.r_twist()
        par = 1 - par
    return (tangle.orient, tangle.points)
