r"""
Two integrals that physics needs and Sage's integrate() does not do by itself.

    delta_integrate(f, g, x, a, b)      Int_a^b f(x) delta(g(x)) dx = sum over the roots x_i of g in (a, b)
                                        of f(x_i)/|g'(x_i)|
    principal_value(f, x, a, b, c)      PV Int_a^b f(x) dx with a simple pole at c: the limit of
                                        Int_a^(c-h) + Int_(c+h)^b as h -> 0 (the window must be symmetric)
"""
from sage.all import SR, solve, integrate, limit


def delta_integrate(f, g, x, a, b, assumptions=()):
    """Int_a^b f(x) delta(g(x)) dx, summed over the simple roots of g between a and b."""
    f, g = SR(f), SR(g)
    tot = SR(0)
    for sol in solve(g == 0, x):
        r = sol.rhs()
        inside = True
        for lo_hi in ((a, r), (r, b)):
            try:
                if bool(SR(lo_hi[0]) >= SR(lo_hi[1])):      # a root on an end point is left out
                    inside = False
            except TypeError:
                pass
        if not inside:
            continue
        gp = g.diff(x).subs({x: r})
        tot += f.subs({x: r}) / gp.abs()
    return tot.simplify_full()


def principal_value(f, x, a, b, c):
    """PV Int_a^b f(x) dx for a simple pole of f at x = c (a < c < b), from the antiderivative F:
    the limit of F(c - h) - F(a) + F(b) - F(c + h) as h -> 0+ (an i pi from log of a negative number
    cancels inside each of the two pieces)."""
    h = SR.var('pv_h')
    F = integrate(SR(f), x)
    expr = (F.subs({x: c - h}) - F.subs({x: a})) + (F.subs({x: b}) - F.subs({x: c + h}))
    out = limit(expr, pv_h=0, dir='+')
    return out.real_part().simplify_full()        # f is real on the line, so is its principal value


def integrate_termwise(expr, x, *bounds):
    """Integrate a sum term by term (Maxima is much more reliable on the single terms of an expanded sum
    than on the whole, where it may join fractions and lose track of logarithm branches)."""
    e = SR(expr).expand()
    op = e.operator()
    terms = e.operands() if op is not None and 'add' in getattr(op, '__name__', str(op)) else [e]
    return sum((integrate(t, x, *bounds) for t in terms), SR(0))
