import numpy

def continuedfrac(a, b): 
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

# something
print("hello")
print("world")