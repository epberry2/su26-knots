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

tangle = TangleState()
print(tangle.orient, tangle.points)
tangle.t_twist()
print(tangle.orient, tangle.points)
tangle.r_twist()
print(tangle.orient, tangle.points)