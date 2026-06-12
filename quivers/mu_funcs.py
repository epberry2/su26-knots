from helpers.tanglestate import TangleState
from colored_homfly.colored_alg import word_from_fraction

def mus(u, v):
    word = word_from_fraction(u - v, v)
    pos = True
    tangle = TangleState()
    mu1 = mu2 = mu3 = 0
    for op in word:
        if op == "T":
            match tangle.orient:
                case "UP":
                    mu1 -= 1
                case "OP":
                    if pos:
                        mu1 += 1
                        mu2 += 1
                        mu3 -= 1
                case "RI":
                    if pos:
                        mu1 += 1
                        mu2 += 1
                        mu3 -= 1
                      
            tangle.t_twist()
            pos = True
        else:
            match tangle.orient:
                case "UP":
                    mu1 -= 1
                    mu2 -= 1
                    mu3 += 1
                    
                case "OP":
                    if pos:
                        mu1 -= 1
                        mu2 -= 1
                        mu3 += 1

                case "RI":
                    if pos:
                        mu3 -= 1

            tangle.r_twist()
            pos = False

    return mu1, mu2, mu3