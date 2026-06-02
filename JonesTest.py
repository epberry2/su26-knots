import sympy as sp
import math
from itertools import permutations

import JonesPolyAlg as jpa
import JonesPolyGeo as jpg

def coprime_pairs(n):
    numbers = range(1, n + 1)
    all_pairs = permutations(numbers, 2)
    coprime_pairs = [pair for pair in all_pairs if math.gcd(pair[0], pair[1]) == 1]
    
    return coprime_pairs

num_correct = 0
pairs = coprime_pairs(20)
for i, j in pairs:
    va = jpa.AlgBracket(i, j)
    vb = jpg.JonesPolyGeo(i, j)

    correct = va[0] == vb[0] or -va[0] == vb[0] and va[1] == vb[1] or -va[1] == vb[1]
    print(f"({i}, {j}): {correct}")
    if correct:
        num_correct += 1
    else:
        print(va)
        print(vb)

print(f"Num correct: {num_correct} / {len(pairs)}, {(num_correct / len(pairs)):.2f} ")
    