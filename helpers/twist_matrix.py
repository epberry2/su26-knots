import sympy as sp
from helpers.quantum_nums import qbinom

q, a = sp.symbols("q a")

def twist_matrix(twist: str, j: int) -> sp.Matrix:
    size = 3 * (j + 1)
    M = sp.zeros(size, size)
    for orientation in ("UP", "OP", "RI"):
        for k in range(j + 1):
            idx = basis_index(orientation, j, k)
            M[:, idx] = twist_map(twist, orientation, j, k)
    return M

def twist_map(twist: str, orientation: str, j: int, k: int):
    """
    e.g. "T", "UP", j, k -> TUP[j:k] as a vector
    """
    if not (0 <= k <= j):
        raise ValueError(f"k must satisfy 0 <= k <= j, got k={k}, j={j}")
    out = sp.zeros(3 * (j + 1), 1)
    match twist, orientation:
        case "T", "UP":
            new_orientation = "UP"
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** (h - j)
                    * q ** (k * k)
                    * qbinom(h, k)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case "T", "OP":
            new_orientation = "RI"
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** h
                    * a ** k
                    * q ** (k * (k - 2 * j))
                    * qbinom(h, k)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case "T", "RI":
            new_orientation = "OP"
            for h in range(k, j + 1):
                coeff = (
                    (-q) ** h
                    * a ** h
                    * q ** (k ** 2 - 2 * j * h)
                    * qbinom(h, k)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case "R", "UP":
            new_orientation = "OP"
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** (h - j)
                    * a ** (h - j)
                    * q ** (k ** 2 + j ** 2 - 2 * k * h)
                    * qbinom(j - h, k - h)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case "R", "OP":
            new_orientation = "UP"
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** (h - j)
                    * a ** (k - j)
                    * q ** (2 * h * (j - k) + (k - j) ** 2)
                    * qbinom(j - h, k - h)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case "R", "RI":
            new_orientation = "RI"
            for h in range(0, k + 1):
                coeff = (
                    (-q) ** h
                    * q ** (h * (2 * j - 2 * k) + k ** 2 - j ** 2)
                    * qbinom(j - h, k - h)
                )
                out[basis_index(new_orientation, j, h)] = coeff
        case _:
            raise ValueError("""Twist must be "T" or "R", orientation must be "UP", "OP", or "RI".""")
    return out

def basis_index(X: str, j: int, k: int) -> int:
    """
    Convert X, j, k into index of basis element X[j, k].
    """
    orientations = ("UP", "OP", "RI")
    return orientations.index(X) * (j + 1) + k

