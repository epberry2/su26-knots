import sympy as sp

def normalize_laurent_2var(
    expr: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[tuple[int, int], sp.Expr]:
    expr = sp.cancel(sp.together(expr))
    expr = sp.expand(expr)
    if expr == 0:
        return (0, 0), sp.Integer(0)
    terms = sp.Add.make_args(expr)
    x_powers = []
    y_powers = []
    for term in terms:
        powers = term.as_powers_dict()
        x_pow = powers.get(x, sp.Integer(0))
        y_pow = powers.get(y, sp.Integer(0))
        if not x_pow.is_integer or not y_pow.is_integer:
            raise ValueError("Expression is not a Laurent polynomial with integer powers.")
        x_powers.append(int(x_pow))
        y_powers.append(int(y_pow))
    x_shift = min(x_powers)
    y_shift = min(y_powers)
    normalized = sp.expand(expr * x**(-x_shift) * y**(-y_shift))
    return normalized