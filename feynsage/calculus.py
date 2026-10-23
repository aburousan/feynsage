r"""
Integrals that physics needs and Sage's integrate() does not do by itself.

    delta_integrate(f, g, x, a, b)      Int_a^b f(x) delta(g(x)) dx = sum over the roots x_i of g in (a, b)
                                        of f(x_i)/|g'(x_i)|
    principal_value(f, x, a, b, c)      PV Int_a^b f(x) dx with simple poles at c (one pole or a list): the
                                        limit of Int_a^(c-h) + Int_(c+h)^b as h -> 0 (symmetric windows)
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


def _residue(f, x, c):
    """The residue of f at a simple pole x = c (f without free symbols other than x): (x - c) f with the
    pole cancelled exactly, at x = c.  limit() is the fallback, since some Maxima versions (5.49) leave
    lim (x - c)/(c^2 - x^2) unevaluated."""
    g = (x - c)*SR(f)
    for simp in (lambda e: e.simplify_rational(), lambda e: e.factor(), lambda e: e.simplify_full()):
        try:
            R = simp(g).subs({x: c})
            v = complex(CC(R.n(80)))
            if v == v and abs(v) != float('inf'):     # not NaN, not infinite
                return R

        except (ValueError, TypeError, ZeroDivisionError, ArithmeticError, RuntimeError):
            pass
    return limit(g, **{str(x): c})


def principal_value(f, x, a, b, c, check=True):
    """PV Int_a^b f(x) dx for simple poles of f at x = c (a < c < b); c may be a list of poles, and every
    pole inside (a, b) must be in it.  From the antiderivative F: with the poles in order c_1 < ... < c_n,
    the limit h -> 0+ of F(c_1 - h) - F(a) + F(c_2 - h) - F(c_1 + h) + ... + F(b) - F(c_n + h) (an i pi
    from log of a negative number cancels inside each piece).  F comes from integrate_checked (checked by
    differentiating it), and the result is checked against
        PV = Int_a^b [f - sum_k R_k/(x - c_k)] dx + sum_k R_k log((b - c_k)/(c_k - a)),
    R_k the residues, done numerically at random values of the other symbols.  The check also looks for the
    real roots of the denominator of a rational f, so a pole left out of c is an error, not a wrong number."""
    from .integrators import integrate_checked, IntegrationError, _points, _numeric_definite, _close
    import random
    from sage.all import assumptions
    f, x = SR(f), SR(x)
    cs = [SR(e) for e in (c if isinstance(c, (list, tuple)) else [c])]
    syms = set(f.variables()) | set(SR(a).variables()) | set(SR(b).variables())
    for e in cs:
        syms |= set(e.variables())
    syms = sorted(syms - {x}, key=str)
    pts = _points(syms, list(assumptions()), 2, random.Random(1)) if syms else [{}]
    if len(cs) > 1:                                      # order the poles at a point that obeys the assumptions
        p0 = pts[0] if pts else {}
        cs = sorted(cs, key=lambda e: float(CC(e.subs(p0).n(60)).real()))
    h = SR.var('pv_h')
    F = integrate_checked(f, x)
    ends = [SR(a)] + cs + [SR(b)]
    expr = SR(0)
    for k in range(len(ends) - 1):                       # the piece from ends[k] (+ h) to ends[k + 1] (- h)
        lo = ends[k] if k == 0 else ends[k] + h
        hi = ends[k + 1] if k == len(ends) - 2 else ends[k + 1] - h
        expr += F.subs({x: hi}) - F.subs({x: lo})
    out = limit(expr, pv_h=0, dir='+').real_part().simplify_full()   # f is real on the line, so is its PV
    if check:
        for p in pts:
            fa, fb = SR(a).subs(p), SR(b).subs(p)
            fcs = [e.subs(p) for e in cs]
            fp = f.subs(p)
            _all_poles_listed(fp, x, fa, fb, fcs)
            Rs = [_residue(fp, x, fc) for fc in fcs]
            g = (fp - sum(R/(x - fc) for R, fc in zip(Rs, fcs))).simplify_full()   # poles cancelled exactly
            num = _numeric_definite(g, x, fa, fb) + sum(
                complex(CC((R*log((fb - fc)/(fc - fa))).n(80))) for R, fc in zip(Rs, fcs))
            if not _close(complex(CC(out.subs(p).n(80))), num, 1e-8):
                raise IntegrationError("principal_value: the limit of the antiderivative does not agree with "
                                       "the numerical principal value at %s" % p)
    return out


def _all_poles_listed(f, x, a, b, cs):
    """For a rational f (numbers only besides x): every real root of its denominator inside (a, b) must be
    one of the poles cs, otherwise the principal value is not what the caller asked for."""
    from .integrators import IntegrationError
    from sage.all import QQbar
    try:
        den = SR(f).denominator().polynomial(QQbar)
    except (TypeError, ValueError, NotImplementedError):
        return                                           # not a rational function: nothing to look for
    fa, fb = CC(SR(a).n(60)).real(), CC(SR(b).n(60)).real()
    listed = [CC(SR(e).n(60)) for e in cs]
    for r in den.roots(multiplicities=False):
        if r.imag() != 0 or not (fa < r.real() < fb):
            continue
        if not any(abs(CC(r) - e) < 1e-12*(1 + abs(e)) for e in listed):
            raise IntegrationError("principal_value: f has a pole at x = %s inside (%s, %s) that is not in c; "
                                   "pass every pole, e.g. c = [%s, ...]" % (r, a, b, r))


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
