from helpers.continuedfrac import continuedfrac

class TangleState():
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
        self.points = (self.points[1], self.points[0], self.points[2])
        if self.points[0] == "X+":
            self.orient = "RI"
        if self.points[1] == "X+":
            self.orient = "OP"
        if self.points[2] == "X+":
            self.orient = "UP"
    
    def r_twist(self):
        self.points = (self.points[0], self.points[2], self.points[1])
        if self.points[0] == "X+":
            self.orient = "RI"
        if self.points[1] == "X+":
            self.orient = "OP"
        if self.points[2] == "X+":
            self.orient = "UP"

# Returns the permutation of [Y, X-, X+] given tau_(u/v)
def get_state(u,v):
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

# print(get_state(5,2))