import sympy as sp
from typing import List

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

def normalize_laurents_2var(
    exprs: List[sp.Expr],
    x: sp.Symbol,
    y: sp.Symbol,
) -> List[sp.Expr]:
    flag = False
    for term in exprs:
        expr = sp.cancel(sp.together(term))
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
        term_x_shift = min(x_powers)
        term_y_shift = min(y_powers)
        if flag:
            x_shift = min(x_shift, term_x_shift)
            y_shift = min(y_shift, term_y_shift)
        else:
            x_shift = term_x_shift
            y_shift = term_y_shift
        flag = True
    out = []
    for expr in exprs:
        normalized = sp.expand(expr * x**(-x_shift) * y**(-y_shift))
        out.append(normalized)
    return out