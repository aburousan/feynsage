r"""
Integrals that physics needs and Sage's integrate() does not do by itself.

    delta_integrate(f, g, x, a, b)      Int_a^b f(x) delta(g(x)) dx = sum over the roots x_i of g in (a, b)
                                        of f(x_i)/|g'(x_i)|
    principal_value(f, x, a, b, c)      PV Int_a^b f(x) dx with a simple pole at c: the limit of
                                        Int_a^(c-h) + Int_(c+h)^b as h -> 0 (the window must be symmetric)
    integrate_termwise(expr, x, a, b)   a sum integrated term by term

The integrals themselves are done by feynsage.integrators.integrate_checked: several backends tried in
turn, every answer checked numerically.
"""
from sage.all import SR, CC, solve, limit, log


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


def principal_value(f, x, a, b, c, check=True):
    """PV Int_a^b f(x) dx for a simple pole of f at x = c (a < c < b), from the antiderivative F:
    the limit of F(c - h) - F(a) + F(b) - F(c + h) as h -> 0+ (an i pi from log of a negative number
    cancels inside each of the two pieces).  F comes from integrate_checked (checked by differentiating
    it), and the result is checked against PV = Int_a^b [f - R/(x - c)] dx + R log((b - c)/(c - a)),
    R the residue, done numerically at random values of the other symbols."""
    from .integrators import integrate_checked, IntegrationError, _points, _numeric_definite, _close
    import random
    from sage.all import assumptions
    f, x = SR(f), SR(x)
    h = SR.var('pv_h')
    F = integrate_checked(f, x)
    expr = (F.subs({x: c - h}) - F.subs({x: a})) + (F.subs({x: b}) - F.subs({x: c + h}))
    out = limit(expr, pv_h=0, dir='+').real_part().simplify_full()   # f is real on the line, so is its PV
    if check:
        syms = set(f.variables()) | set(SR(a).variables()) | set(SR(b).variables()) | set(SR(c).variables())
        syms = sorted(syms - {x}, key=str)
        pts = _points(syms, list(assumptions()), 2, random.Random(1)) if syms else [{}]
        for p in pts:
            fa, fb, fc = SR(a).subs(p), SR(b).subs(p), SR(c).subs(p)
            fp = f.subs(p)
            R = limit((x - fc)*fp, x=fc)
            g = (fp - R/(x - fc)).simplify_full()            # the pole cancelled exactly, not by rounding
            num = _numeric_definite(g, x, fa, fb) + complex(CC((R*log((fb - fc)/(fc - fa))).n(80)))
            if not _close(complex(CC(out.subs(p).n(80))), num, 1e-8):
                raise IntegrationError("principal_value: the limit of the antiderivative does not agree with "
                                       "the numerical principal value at %s" % p)
    return out


def integrate_termwise(expr, x, *bounds):
    """Integrate a sum term by term, each term with integrate_checked (several backends, each answer checked
    numerically).  Single terms are much easier than the whole sum, where a backend may join fractions and
    lose track of logarithm branches."""
    from .integrators import integrate_checked
    e = SR(expr).expand()
    op = e.operator()
    terms = e.operands() if op is not None and 'add' in getattr(op, '__name__', str(op)) else [e]
    return sum((integrate_checked(t, x, *bounds) for t in terms), SR(0))


def tidy(expr, radicals=True):
    """
    Make a result from integrate() readable: the argument of every log is factored and the logs are
    expanded (so log((a + b)^2) - 2 log(a + b) cancels), radicals are simplified with the assumptions
    in force (assume(x > 0) first) and the rest is simplified.  The value is unchanged.
    """
    from sage.all import SR, log
    e = SR(expr)
    w0 = SR.wild(0)
    for L in set(e.find(log(w0))):
        arg = L.operands()[0]
        try:
            e = e.subs({L: log(arg.factor())})
        except Exception:
            pass
    e = e.expand_log(algorithm='powers')
    if radicals:
        try:
            e = e.canonicalize_radical()
        except Exception:
            pass
    try:
        return e.simplify_full()
    except Exception:
        return e
