r"""
Laurent expansions that work at the poles of the Gamma function.

    laurent(expr, eps, 0)                  expr = sum_k c_k eps^k up to eps^0, poles included
    laurent(expr, [(delta, 0), (eps, 0)])  first in delta, then every coefficient in eps

Sage's series() stops at Gamma(-k + c eps) with k = 0, 1, 2, ...  Here every such Gamma is first
rewritten with Gamma(z + 1) = z Gamma(z) as Gamma(1 + c eps) / [(c eps - k)(c eps - k + 1) ... (c eps)],
so what is left is regular and series() can do the rest.  No factor e^(eps gamma_E) is put in (that is
what oneloop.expand_eps does for loop integrals).
"""
from sage.all import SR, ZZ, gamma, prod, psi, factorial


def _at_zero(z, x):
    """z at x = 0, also when plain substitution meets a removable singularity such as (e^x - 1)/x."""
    try:
        return z.subs({x: 0})
    except (ValueError, ZeroDivisionError, RuntimeError):
        from sage.all import limit
        return limit(z, **{str(x): 0})


def _regularise(e, x):
    op = e.operator()
    if op is None:
        return e
    args = [_regularise(a, x) for a in e.operands()]
    if op == gamma:
        z = args[0]
        z0 = _at_zero(z, x)
        if z0 in ZZ and z0 <= 0:
            k = -ZZ(z0)
            return gamma(z + k + 1) / prod(z + j for j in range(k + 1))
        return gamma(z)
    name = getattr(op, '__name__', str(op))
    if name in ('psi', 'polygamma') or str(op) in ('psi', 'polygamma'):
        # psi(m, z) = psi(m, z + k + 1) - sum_j (-1)^m m!/(z + j)^(m + 1), regular at the poles z = -k
        if len(args) == 1:
            m, z = 0, args[0]
        else:
            m, z = args
        z0 = _at_zero(z, x)
        if z0 in ZZ and z0 <= 0 and m in ZZ:
            k = -ZZ(z0)
            shifted = psi(z + k + 1) if m == 0 else psi(m, z + k + 1)
            return shifted - sum((-1)**m * factorial(m) / (z + j)**(m + 1) for j in range(k + 1))
    try:
        return op(*args)
    except TypeError:
        return e


def laurent(expr, x, order=0, simplify=True):
    r"""
    The Laurent expansion of expr in x around 0 up to and including x^order (poles included).
    x can also be a list [(x1, order1), (x2, order2), ...]: expand in x1 first, then each coefficient
    in x2, and so on (for an analytic regulator delta and then eps).
    """
    if isinstance(x, (list, tuple)):
        out = SR(expr)
        for v, o in x:
            out = laurent(out, v, o, simplify=False)
        return out.simplify_full() if simplify else out
    e = _regularise(SR(expr), x)
    try:
        s = e.series(x, order + 1).truncate()
    except (ValueError, ZeroDivisionError):
        # GiNaC's series stops at a removable singularity such as (e^x - 1)/x; Maxima's taylor does not
        s = e.taylor(x, 0, order)
    if simplify:
        s = s.collect(x)
        s = sum((c.simplify_full() * x**k for c, k in s.coefficients(x)), SR(0))
    return s


def coefficient(expr, x, k):
    """The coefficient of x^k in the Laurent expansion of expr (poles of Gamma included)."""
    return laurent(expr, x, k, simplify=False).expand().coefficient(x, k)


def residue(expr, z, z0, order=None):
    """The residue of expr at z = z0, also at poles of Gamma and psi functions (z0 may be an
    expression such as -1 - eps).  It is the coefficient of 1/h in the Laurent series in h = z - z0."""
    h = SR.var('res_h')
    e = SR(expr).subs({SR(z): SR(z0) + h})
    L = laurent(e, h, -1, simplify=False)
    r = L.expand().coefficient(h, -1)
    try:
        return r.simplify_rational()           # simplify_full would turn Gamma(-eps) into factorials
    except Exception:
        return r
