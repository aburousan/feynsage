r"""
One-loop tensor integrals and Passarino-Veltman functions, in the conventions of
Package-X and LoopTools.

Conventions
-----------
Minkowski metric (+,-,-,-).  A propagator (q, m) is 1/(q^2 - m^2 + i0).  The measure is

    mu^(2 eps) e^(eps gamma_E) (4 pi)^(-eps) * (1/(i pi^2)) * Int d^d l ,   d = 4 - 2 eps

so that   A0(m) = m^2 (1/eps + log(mu^2/m^2) + 1)   and   B0 = 1/eps + (finite).
These are exactly the normalisations of Package-X (LoopRefine with Eps, Mu) and of
LoopTools.  The +i0 is kept by two special functions:

    LogM(z)          = log(z - i0)   (z is always a "mass^2 minus momenta" combination)
    DiscB(s, m1, m2) = sqrt(lam)/s * log((m1^2 + m2^2 - s + sqrt(lam))/(2 m1 m2)),  s -> s + i0

C0 and D0: IR-finite ones are symbolic functions whose .n() uses the closed forms of
feynsage.scalar (dilogarithms, > 30 digits) and explicit() writes them out; IR-divergent ones
come out as Laurent series in eps with the poles explicit (feynsage.ir, Ellis-Zanderighi).

Quick start
-----------
    from feynsage.pv import *
    loop("1", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})             # B0(s, m, m)
    loop("l^mu l^nu", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})   # g^{mu nu} B00 + p^mu p^nu B11
    B0(s, m, m).explicit()   # closed form, like Package-X's LoopRefine
"""
import itertools
import re
from sage.all import (SR, RR, var, function, log, sqrt, pi, I, gamma, matrix, factorial, PolynomialRing, lcm,
                      exp, euler_gamma, beta, integrate,
                      Integer, QQ, prod, binomial)

eps, mu = var('eps mu')
Dim = var('d')                      # the space-time dimension while reducing, d = 4 - 2 eps at the end

# ---------------------------------------------------------------------------- special functions
def _mp():
    import mpmath
    return mpmath


def _cplx(z):
    return complex(z)


def _logm_evalf(self, z, parent=None, algorithm=None):
    mp = _mp()
    z = mp.mpc(_cplx(z))
    if z.imag == 0 and z.real < 0:
        val = mp.log(-z.real) - 1j * mp.pi
    else:
        val = mp.log(z)
    return parent(complex(val)) if parent is not None else complex(val)


LogM = function('LogM', nargs=1, evalf_func=_logm_evalf,
                print_latex_func=lambda self, z: r"\ln\left(%s - i0\right)" % z._latex_())


def _discb_num(s, m1, m2):
    mp = _mp()
    with mp.workdps(40):
        s, m1, m2 = mp.mpf(float(s)), mp.mpf(float(m1)), mp.mpf(float(m2))
        z = s + 1j * mp.mpf('1e-25') * max(1, abs(s), m1 ** 2, m2 ** 2)
        lam = z ** 2 + m1 ** 4 + m2 ** 4 - 2 * z * m1 ** 2 - 2 * z * m2 ** 2 - 2 * m1 ** 2 * m2 ** 2
        r = mp.sqrt(lam)
        return complex(r / z * mp.log((m1 ** 2 + m2 ** 2 - z + r) / (2 * m1 * m2)))


def _discb_evalf(self, s, m1, m2, parent=None, algorithm=None):
    v = _discb_num(s, m1, m2)
    return parent(v) if parent is not None else v


def _discb_deriv(self, s, m1, m2, diff_param=None):
    if diff_param != 0:
        raise NotImplementedError("DiscB is differentiated only in s")
    lam = s ** 2 + m1 ** 4 + m2 ** 4 - 2 * s * m1 ** 2 - 2 * s * m2 ** 2 - 2 * m1 ** 2 * m2 ** 2
    dlam = 2 * s - 2 * m1 ** 2 - 2 * m2 ** 2
    return (s * dlam - 2 * lam) / (2 * s * lam) * DiscB(s, m1, m2) - 1 / s


DiscB = function('DiscB', nargs=3, evalf_func=_discb_evalf, derivative_func=_discb_deriv,
                 print_latex_func=lambda self, s, a, b: r"\Lambda(%s; %s, %s)" % (s._latex_(), a._latex_(), b._latex_()))


def _exact_args(args):
    """Numbers -> exact rationals (floats by their simplest rational), exact algebraic numbers kept."""
    out = []
    for a in args:
        a = SR(a)
        if a.is_numeric():
            try:
                out.append(QQ(a))
            except (TypeError, ValueError):
                out.append(RR(a).simplest_rational())
        else:
            out.append(a)
    return out


def _to_parent(v, parent):
    """An mpmath value into Sage's requested field, keeping all its digits."""
    from sage.all import RealField, ComplexField
    if parent is None:
        return complex(v)
    import mpmath
    prec = getattr(parent, 'prec', lambda: 53)()
    R = RealField(prec)
    z = ComplexField(prec)(R(mpmath.nstr(mpmath.re(v), 75)), R(mpmath.nstr(mpmath.im(v), 75)))
    try:
        return parent(z)
    except (TypeError, ValueError):
        return z


def _c0_evalf(self, *args, parent=None, algorithm=None):
    """Closed form (scalar.c0_value, about 40 digits); contour integration only where the closed
    form does not apply (lambda(s1, s12, s2) = 0 with a nonzero invariant)."""
    from .scalar import c0_value
    try:
        v = c0_value(*_exact_args(args), full=True)
    except ZeroDivisionError:
        v = c0_numeric(*[float(a) for a in args])
    return _to_parent(v, parent)


def _d0_evalf(self, *args, parent=None, algorithm=None):
    """Closed form (scalar.d0_value: Denner, Denner-Dittmaier, about 40 digits); contour integration
    only at exceptional points where no labelling gives a regular closed form."""
    from .scalar import d0_value
    try:
        v = d0_value(*_exact_args(args), full=True)
    except ZeroDivisionError:
        v = d0_numeric(*[float(a) for a in args])
    return _to_parent(v, parent)


_C0f = function('C0', nargs=6, evalf_func=_c0_evalf)
_D0f = function('D0', nargs=10, evalf_func=_d0_evalf)


def C0(s1, s12, s2, m0, m1, m2):
    """Scalar triangle C0(s1, s12, s2; m0, m1, m2), propagators (l, m0), (l + p1, m1), (l + p2, m2),
    s1 = p1^2, s2 = p2^2, s12 = (p1 - p2)^2 (Package-X order).
    Output: IR finite -> the symbol C0(...) (.n() gives > 30 digits, explicit() the dilogarithms);
    IR divergent (soft or collinear) -> c_-2/eps^2 + c_-1/eps + c_0 with mu, poles explicit."""
    from .ir import ir_coefficients, with_mu
    c = ir_coefficients('C', (s1, s12, s2, m0, m1, m2))
    if c is not None:
        return with_mu(c)
    return _C0f(s1, s12, s2, m0, m1, m2)


def D0(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3):
    """Scalar box D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3) (Package-X / LoopTools order),
    propagators (l, m0), (l + p1, m1), (l + p1 + p2, m2), (l + p1 + p2 + p3, m3).
    Output: IR finite -> the symbol D0(...) (.n() gives > 30 digits, explicit() the dilogarithms);
    IR divergent -> c_-2/eps^2 + c_-1/eps + c_0 with mu, poles explicit."""
    from .ir import ir_coefficients, with_mu
    c = ir_coefficients('D', (s1, s2, s3, s4, s12, s23, m0, m1, m2, m3))
    if c is not None:
        return with_mu(c)
    return _D0f(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3)


# ---------------------------------------------------------------------------- scalar functions
def _is_zero(x):
    """Structural zero test (no Maxima call): numbers are compared exactly, symbols are not zero."""
    x = SR(x)
    return x.is_trivial_zero() or (x.is_numeric() and x == 0)


def A0(m):
    """Tadpole A0(m) = m^2 (1/eps + log(mu^2/m^2) + 1) in the Package-X normalisation; A0(0) = 0.
    Output: a symbolic expression in eps (the 1/eps UV pole) and mu (the renormalisation scale)."""
    m = SR(m)
    if _is_zero(m):
        return SR(0)
    return m ** 2 * (1 / eps + log(mu ** 2 / m ** 2) + 1)


def B0(s, m1, m2):
    """Two-point function B0(p^2 = s; m1, m2) in the Package-X / LoopTools normalisation.
    Output: 1/eps + finite part, the finite part in closed form (as Package-X's LoopRefine):
    logs, LogM(z) = log(z - i0) and DiscB(s, m1, m2) = sqrt(lam)/s * log((m1^2+m2^2-s+sqrt(lam))/(2 m1 m2))
    with s -> s + i0.  Numbers in, a number with .n(); symbols in, a formula (masses assumed nonzero)."""
    s, m1, m2 = SR(s), SR(m1), SR(m2)
    z1, z2, zs = _is_zero(m1), _is_zero(m2), _is_zero(s)
    if z1 and not z2:
        m1, m2, z1, z2 = m2, m1, z2, z1          # put the massless line second
    if zs:
        if z1 and z2:
            return SR(0)                          # scaleless
        if z2:                                    # B0(0; m, 0)
            return 1 / eps + 1 + log(mu ** 2 / m1 ** 2)
        if (m1 - m2).is_trivial_zero():
            return 1 / eps + log(mu ** 2 / m1 ** 2)
        return 1 / eps + 1 + (m1 ** 2 * log(mu ** 2 / m1 ** 2) - m2 ** 2 * log(mu ** 2 / m2 ** 2)) / (m1 ** 2 - m2 ** 2)
    if z1 and z2:
        return 1 / eps + 2 + log(mu ** 2) - LogM(-s)
    if z2:                                        # B0(s; m, 0)
        return 1 / eps + 2 + log(mu ** 2 / m1 ** 2) + (m1 ** 2 - s) / s * (LogM(m1 ** 2 - s) - log(m1 ** 2))
    return (1 / eps + 2 + DiscB(s, m1, m2) - (m1 ** 2 - m2 ** 2 + s) / (2 * s) * log(m1 ** 2 / m2 ** 2)
            + log(mu ** 2 / m2 ** 2))


def uv_part(expr):
    """The coefficient of 1/eps (the UV pole if the integral is IR finite; UV and IR together otherwise)."""
    return SR(expr).expand().coefficient(eps, -1)


def pole_parts(expr):
    """(coefficient of 1/eps^2, coefficient of 1/eps)."""
    e = SR(expr).expand()
    return e.coefficient(eps, -2), e.coefficient(eps, -1)


def finite_part(expr, mu_value=None):
    """The eps^0 part; with mu_value the scale mu is set to that number."""
    e = SR(expr).expand().coefficient(eps, 0)
    return e.subs(mu=mu_value) if mu_value is not None else e


# ---------------------------------------------------------------------------- numerical C0 and D0
def _feynman_numeric(npoint, masses, sij, n):
    r"""
    (-1)^N Gamma(N-2) Int_simplex d^(N-1)x  Delta^(2-N),   Delta = sum x_i m_i^2 - sum_{i<j} x_i x_j s_ij - i0.
    The simplex is mapped to the unit cube (x1 = t1, x2 = (1-t1) t2, ...) and the contour
    is deformed, t -> t - i lam t(1-t) dDelta/dt, so that Im Delta < 0 everywhere.
    """
    import numpy as np
    N = npoint
    k = N - 1
    msq = np.array([float(m) ** 2 for m in masses])
    S = np.zeros((N, N))
    for (i, j), v in sij.items():
        S[i, j] = S[j, i] = float(v)
    xg, wg = np.polynomial.legendre.leggauss(n)
    xg, wg = (xg + 1) / 2, wg / 2
    grids = np.meshgrid(*([xg] * k), indexing='ij')
    T = np.stack([g.ravel() for g in grids])                     # k x P
    W = np.prod(np.stack(np.meshgrid(*([wg] * k), indexing='ij')), axis=0).ravel()

    def xs_of(t):
        x, rest = [], 1
        for a in range(k):
            x.append(rest * t[a])
            rest = rest * (1 - t[a])
        return [rest] + x, None                                  # x0 = what is left

    def jac_simplex(t):
        j = 1
        for a in range(k):
            j = j * (1 - t[a]) ** (k - 1 - a)
        return j

    def delta(t):
        x, _ = xs_of(t)
        val = sum(x[i] * msq[i] for i in range(N))
        for i in range(N):
            for j in range(i + 1, N):
                if S[i, j]:
                    val = val - x[i] * x[j] * S[i, j]
        return val

    D0 = delta(T)
    if np.all(D0 > 0):
        Z, J = T, np.ones(T.shape[1])
    else:
        h = 1e-6
        grad = np.zeros_like(T)
        for a in range(k):
            e = np.zeros((k, 1)); e[a] = h
            grad[a] = (delta(T + e) - delta(T - e)) / (2 * h)
        scale = max(1.0, float(np.max(np.abs(grad))))
        lam = 0.5 / scale
        h2 = 1e-4
        hess = np.zeros((k, k, T.shape[1]))
        for a in range(k):
            for b in range(k):
                ea = np.zeros((k, 1)); ea[a] = h2
                eb = np.zeros((k, 1)); eb[b] = h2
                hess[a, b] = (delta(T + ea + eb) - delta(T + ea - eb) - delta(T - ea + eb) + delta(T - ea - eb)) / (4 * h2 * h2)
        f = T * (1 - T)
        Z = T - 1j * lam * f * grad
        Jm = np.zeros((k, k, T.shape[1]), dtype=complex)
        for a in range(k):
            for b in range(k):
                Jm[a, b] = (1.0 if a == b else 0.0) - 1j * lam * ((1 - 2 * T[a]) * grad[a] * (a == b) + f[a] * hess[a, b])
        J = np.linalg.det(np.moveaxis(Jm, -1, 0))
    Dz = delta(Z)
    pref = (-1) ** N * float(gamma(N - 2))
    return complex(pref * np.sum(W * J * jac_simplex(Z) * Dz ** (2 - N)))


_checked = {}


def _twice(npoint, masses, sij, n, tol, name):
    """Evaluate on two grids; warn if they differ by more than tol (relative)."""
    import warnings
    a = _feynman_numeric(npoint, masses, sij, n)
    b = _feynman_numeric(npoint, masses, sij, int(n * 1.3))
    err = abs(a - b) / max(1.0, abs(b))
    _checked[(name, tuple(masses), tuple(sorted(sij.items())))] = err
    if err > tol:
        warnings.warn("%s: two grids differ by %.1e; raise n or check for an IR divergence" % (name, err))
    return b


def c0_numeric(s1, s12, s2, m0, m1, m2, n=260, tol=1e-6):
    """C0(s1, s12, s2; m0, m1, m2) numerically (Package-X/LoopTools argument order), IR-finite cases.
    Computed on two Gauss-Legendre grids (n and 1.3 n points per direction) as an error check."""
    return _twice(3, (m0, m1, m2), {(0, 1): s1, (0, 2): s2, (1, 2): s12}, n, tol, "C0")


def d0_numeric(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, n=90, tol=1e-5):
    """D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3) numerically, IR-finite cases (two grids)."""
    sij = {(0, 1): s1, (1, 2): s2, (2, 3): s3, (0, 3): s4, (0, 2): s12, (1, 3): s23}
    return _twice(4, (m0, m1, m2, m3), sij, n, tol, "D0")


# ---------------------------------------------------------------------------- exact kinematics
def _parse_key(key):
    key = key.replace(' ', '')
    if key.endswith('^2'):
        a = key[:-2]
        return a, a
    a, b = key.split('.')
    return a, b


class Kin:
    r"""
    External momenta and their scalar products.  All arithmetic of the reduction is done
    exactly in the field K = Q(d, invariants, masses), which is much faster than general
    symbolic expressions.  kin = {"p^2": "s", "p1.p2": "(s1 + s2 - s12)/2", ...};
    a missing product p.q becomes a new symbol pq (p^2 becomes p2).
    """

    def __init__(self, vectors, kin=None, symbols=(), euclidean=False, fixed=None):
        self.vectors = list(vectors)
        self.euclidean = euclidean
        raw = {}
        for key, val in (kin or {}).items():
            a, b = _parse_key(key)
            raw[tuple(sorted((a, b)))] = SR(val)
        for key, name in (fixed or {}).items():
            raw[tuple(sorted(key))] = SR.var(name)
        names = {'d'} | {str(x) for x in symbols}
        for v in raw.values():
            names |= {str(x) for x in v.variables()}
        for i, a in enumerate(self.vectors):
            for b in self.vectors[i:]:
                k = tuple(sorted((a, b)))
                if k not in raw:
                    raw[k] = SR.var(a + '2' if a == b else a + b)
                    names.add(str(raw[k]))
        self.R = PolynomialRing(QQ, sorted(names))
        self.K = self.R.fraction_field()
        self.d = self.K('d')
        self.sp = {k: self.to_K(v) for k, v in raw.items()}

    def to_K(self, x):
        return self.K(str(SR(x)))

    def dot(self, u, v):
        K = self.K
        return sum((K(cu) * K(cv) * self.sp[tuple(sorted((a, b)))] for a, cu in u.items() for b, cv in v.items()), K(0))


def _vec(s, vectors):
    """'l + 2 p1 - p2' -> ({'p1': 2, 'p2': -1}, coefficient of l)."""
    e = SR(s)
    out, cl = {}, 0
    for v in list(vectors) + ['l']:
        c = e.coefficient(SR.var(v))
        if c != 0:
            if v == 'l':
                cl = c
            else:
                out[v] = QQ(c)
    return out, cl


def _key(v):
    return tuple(sorted((a, c) for a, c in v.items() if c != 0))


def _add(u, v, cv=1):
    w = dict(u)
    for a, c in v.items():
        w[a] = w.get(a, 0) + cv * c
    return {a: c for a, c in w.items() if c != 0}


def _acc(total, part, c):
    """total += c * part for {master: coefficient} dictionaries."""
    if c == 0:
        return total
    for k, v in part.items():
        w = total.get(k, 0) + c * v
        if w == 0:
            total.pop(k, None)
        else:
            total[k] = w
    return total


# ---------------------------------------------------------------------------- the reduction engine
_XREG = []          # values of the Feynman-parameter masters ('X', n): (full, (1/eps^2, 1/eps))


def _divided_difference(nodes, deriv):
    """f[z_0, ..., z_n] with repeated nodes allowed; deriv(j, z) = f^(j)(z)."""
    nodes = sorted(nodes, key=str)
    memo = {}

    def dd(t):
        if t in memo:
            return memo[t]
        if all((t[0] - z).is_trivial_zero() for z in t):
            v = deriv(len(t) - 1, t[0]) / factorial(len(t) - 1)
        else:
            i1 = max(i for i in range(len(t)) if not (t[i] - t[0]).is_trivial_zero())
            a = t[:i1] + t[i1 + 1:]          # without z_i1
            b = t[1:]                        # without z_0
            v = (dd(b) - dd(a)) / (t[i1] - t[0])
        memo[t] = v
        return v
    return dd(tuple(SR(z) for z in nodes))


def _manifestly_nonneg(e):
    """Every term of the expanded e has a positive coefficient and even powers of all symbols."""
    e = SR(e).expand()
    if e.is_numeric():
        return bool(e >= 0)
    terms = e.operands() if e.operator() is not None and e.operator().__name__ == 'add' else [e]
    for t in terms:
        c = t
        for v in t.variables():
            dgr = t.degree(v)
            if int(dgr) != dgr or int(dgr) % 2:
                return False
            c = c.coefficient(v, dgr)
        if not (c.is_numeric() and bool(c > 0)):
            return False
    return True


def _bernstein_nonneg(Delta, y):
    """Quadratic Delta(y) >= 0 on [0, 1] if its Bernstein coefficients are manifestly >= 0."""
    D0, D1 = Delta.subs({y: 0}), Delta.subs({y: 1})
    b1 = D0 + Delta.diff(y).subs({y: 0}) / 2
    return all(_manifestly_nonneg(b) for b in (D0, b1, D1))


def _fast(expr, y):
    return lambda t: float(expr.subs({y: t}))


def _simplex_integral(groups, n, Msq, sfun, N, k):
    r"""Int over the (N-1)-simplex of prod_i x_i^n_i Delta^s, s = d/2 + k - N = 2 - eps + k - N, as an
    exact expression in eps.  sfun(i, j) is the coefficient of x_i x_j in Delta (i in one group, j in another)."""
    s = 2 - eps + k - N
    if len(groups) == 1:
        Kd = N - 1 + sum(n)

        def deriv(j, z):
            if z.is_trivial_zero():
                return SR(0)                 # 0^(c - eps) = 0 in dimensional regularisation
            return z ** (s + Kd - j) / prod([s + i for i in range(1, Kd - j + 1)], SR(1))
        nodes = [Msq[i] for i in range(N) for _ in range(n[i] + 1)]
        return prod([factorial(x) for x in n], SR(1)) * _divided_difference(nodes, deriv)
    for G in groups:
        if any(not (Msq[i] - Msq[G[0]]).is_trivial_zero() for i in G):
            raise ZeroDivisionError("vanishing Gram determinant with unequal masses on kinematically identical "
                                    "lines: not covered; move the kinematics slightly")
    if len(groups) != 2:
        raise ZeroDivisionError("vanishing Gram determinant with more than two independent line groups: "
                                "not covered; move the kinematics slightly")
    weight, alpha = SR(1), []
    for G in groups:
        nG = sum(n[i] for i in G)
        weight *= prod([factorial(n[i]) for i in G], SR(1)) / factorial(nG + len(G) - 1)
        alpha.append(nG + len(G) - 1)
    y = SR.var('fs_y')
    G0, G1 = groups
    M0, M1 = Msq[G0[0]], Msq[G1[0]]
    c01 = sfun(G0[0], G1[0])
    Delta = (c01 * y * (1 - y) + M0 * (1 - y) + M1 * y).expand()      # y = group-1 variable
    a, b = alpha[1], alpha[0]                                          # y^a (1-y)^b
    # Delta = c y^p (1-y)^q: Beta functions, exact in eps
    for p in (0, 1, 2):
        for q in (0, 1, 2):
            c = (Delta / (y ** p * (1 - y) ** q)).simplify_rational()
            if not c.has(y):
                u, v = a + p * s + 1, b + q * s + 1
                if c.is_numeric() and bool(c < 0):          # (c - i0)^s: the +i0 of the propagators
                    cs = (-c) ** s * exp(-I * pi * s)
                else:
                    cs = c ** s
                return weight * cs * gamma(u) * gamma(v) / gamma(u + v)
    # general quadratic: expand in eps (the prefactor has at most a simple pole) and integrate.
    # Only where Delta > 0 on [0, 1] (below threshold), so that log(Delta - i0) = log(Delta).
    if not _bernstein_nonneg(Delta, y) and set(Delta.variables()) - {y}:
        raise ZeroDivisionError("vanishing Gram determinant with symbolic kinematics where Delta > 0 on [0, 1] "
                                "is not manifest (it may be above threshold): give numbers")
    if not Delta.variables() or set(Delta.variables()) == {y}:
        import mpmath as _m
        f = _fast(Delta, y)
        vals = [f(t / 64) for t in range(65)]
        roots = [r for r in Delta.roots(y, ring=RR, multiplicities=False) if 0 <= r <= 1] if Delta.degree(y) > 0 else []
        if min(vals) <= 0 or roots:
            raise ZeroDivisionError("vanishing Gram determinant above threshold (Delta changes sign): "
                                    "not covered; move the kinematics slightly")
    s0 = 2 + k - N
    L = log(Delta)
    out = SR(0)
    jmax = 1 if N - k - 2 <= 0 else 0              # Gamma(N - k - d/2) has a pole only then
    for j in range(jmax + 1):
        integrand = (y ** a * (1 - y) ** b * Delta ** s0 * (-L) ** j / factorial(j)).expand()
        try:
            val = integrate(integrand, y, 0, 1)
        except ValueError:                           # Maxima asks for the sign of a symbolic combination
            raise ZeroDivisionError("vanishing Gram determinant: the Feynman-parameter integral needs the "
                                    "sign of a symbolic combination; give numbers")
        out += val * eps ** j
    return weight * out


class _Engine:
    r"""
    Reduces one-loop integrals with numerators built from l^2 and l.v to the masters
    A0, B0, C0, D0.  Results are {master key: coefficient in K}.

    Minkowski (default): D = (l + q)^2 - m^2.   Euclidean: D = (l + q)^2 + m^2.
    Masters are stored with the invariants of the chosen metric; they are converted to
    the Package-X (Minkowski) functions only at the end, using I_M = (-1)^N I_E, p_M^2 = -p_E^2.
    """

    def __init__(self, kin):
        self.kin = kin
        self.K = kin.K
        self.d = kin.d
        self.e = -1 if kin.euclidean else 1          # D = q^2 - e m^2
        self.cache = {}

    def dot(self, u, v):
        return self.kin.dot(u, v)

    def msq(self, m):
        return self.K(m) ** 2

    # masters ----------------------------------------------------------------
    def master(self, props):
        N = len(props)
        q = [p[0] for p in props]
        m = [p[1] for p in props]
        sq = lambda a, b: self.dot(_add(q[a], q[b], -1), _add(q[a], q[b], -1))
        if N == 1:
            return {('A', m[0]): self.K(1)}
        if N == 2:
            a, b = sorted([m[0], m[1]], key=str)
            return {('B', sq(0, 1), a, b): self.K(1)}
        if N == 3:
            best = None
            for pm in itertools.permutations(range(3)):
                args = (sq(pm[0], pm[1]), sq(pm[1], pm[2]), sq(pm[0], pm[2]), m[pm[0]], m[pm[1]], m[pm[2]])
                if best is None or str(args) < str(best):
                    best = args
            return {('C',) + best: self.K(1)}
        if N == 4:
            return {('D', sq(0, 1), sq(1, 2), sq(2, 3), sq(0, 3), sq(0, 2), sq(1, 3), m[0], m[1], m[2], m[3]): self.K(1)}
        raise NotImplementedError("more than four propagators: reduce to boxes first")

    # span of the momenta ---------------------------------------------------------
    def span_coeffs(self, props, v):
        Q = [p[0] for p in props[1:]]
        if not v:
            return [0] * len(Q)
        if not Q:
            return None
        names = sorted({a for q in Q for a in q} | set(v))
        A = matrix(QQ, [[q.get(a, 0) for q in Q] for a in names])
        b = matrix(QQ, [[v.get(a, 0)] for a in names])
        try:
            x = A.solve_right(b)
        except ValueError:
            return None
        return [x[i, 0] for i in range(len(Q))]

    # scalar integrals ----------------------------------------------------------
    def scalar(self, props, lv, l2pow=0):
        """Int (l^2)^l2pow prod_j (l.v_j) / prod_i D_i  (props[0] has q = 0)."""
        if not props:
            return {}
        key = ('S', tuple((_key(q), m) for q, m in props), tuple(sorted(_key(v) for v in lv)), l2pow)
        if key not in self.cache:
            self.cache[key] = self._scalar(props, list(lv), l2pow)
        return self.cache[key]

    def _cancel(self, props, i, lv, l2pow):
        rest = props[:i] + props[i + 1:]
        if not rest:
            return {}
        qr = rest[0][0]
        new = tuple((_add(q, qr, -1), m) for q, m in rest)
        return self.poly_numerator(new, lv, l2pow, qr)

    def poly_numerator(self, props, lv, l2pow, qr):
        """Integrand prod (l.v - qr.v) (l^2 - 2 l.qr + qr^2)^l2pow over props (the shift l -> l - qr)."""
        total = {}
        terms = [([], self.K(1))]
        for v in lv:
            nt = []
            qv = self.dot(qr, v)
            for fac, c in terms:
                nt.append((fac + [v], c))
                if qv != 0:
                    nt.append((fac, -c * qv))
            terms = nt
        qq = self.dot(qr, qr)
        for fac, c in terms:
            for a in range(l2pow + 1):
                rest = l2pow - a
                for b in range(rest + 1):
                    cb = binomial(l2pow, a) * binomial(rest, b) * (-2) ** b * qq ** (rest - b)
                    if cb != 0:
                        _acc(total, self.scalar(props, fac + [qr] * b, a), c * cb)
        return total

    def _scalar(self, props, lv, l2pow):
        e = self.e
        m0sq = self.msq(props[0][1])
        if l2pow > 0:                                    # l^2 = D_0 + e m0^2
            out = dict(self._cancel(props, 0, lv, l2pow - 1))
            return _acc(out, self.scalar(props, lv, l2pow - 1), e * m0sq)
        if not lv:
            return self.master(props)
        for j, v in enumerate(lv):                        # l.q_i = (D_i - D_0 - q_i^2 + e m_i^2 - e m0^2)/2
            c = self.span_coeffs(props, v)
            if c is None:
                continue
            others = lv[:j] + lv[j + 1:]
            out = {}
            for i, ci in enumerate(c, start=1):
                if ci == 0:
                    continue
                qi, mi = props[i]
                f = self.dot(qi, qi) - e * self.msq(mi) + e * m0sq
                _acc(out, self._cancel(props, i, others, 0), QQ(ci) / 2)
                _acc(out, self._cancel(props, 0, others, 0), -QQ(ci) / 2)
                _acc(out, self.scalar(props, others, 0), -QQ(ci) / 2 * f)
            return out
        T = self.tensor(props, len(lv))                   # no factor in the span: tensor reduction
        return self.contract_tensor(T, lv, props)

    # tensor integrals by projection ----------------------------------------------
    def independent(self, props):
        """Labels 1..N-1 of a linearly independent subset of the momenta (zero momenta dropped)."""
        names = sorted({a for q, _ in props for a in q})
        labs, rows = [], []
        for i in range(1, len(props)):
            row = [props[i][0].get(a, 0) for a in names]
            if any(row) and matrix(QQ, rows + [row]).rank() == len(rows) + 1:
                labs.append(i)
                rows.append(row)
        return labs

    def structures(self, props, r):
        labs = self.independent(props)
        return [(k, ms) for k in range(r // 2 + 1)
                for ms in itertools.combinations_with_replacement(labs, r - 2 * k)]

    def tensor2_lightlike(self, props, r):
        r"""Two-point tensor with q^2 = 0, from Feynman parameters (no Gram matrix):
        Delta(x) = (1-x) m0^2 + x m1^2 is linear in x, so every coefficient is an exact
        combination of the tadpoles A(m0), A(m1).  J(n, M) = Int (k^2)^n/(k^2 - e M)^2
        = e^(n-1) (d/2 - 1 + n) M^(n-1) A(M)."""
        K, d, e = self.K, self.d, self.e
        m0, m1 = props[0][1], props[1][1]
        M0, M1 = K(m0) ** 2, K(m1) ** 2
        res = {}
        for k in range(r // 2 + 1):
            a = r - 2 * k                              # power of x from (-x q)^a
            norm = prod([d + 2 * i for i in range(k)], K(1))
            pref = K((-1) ** a) * K(e) ** (k - 1) * (d / 2 - 1 + k) / norm
            p = d / 2 + k - 2                          # Delta^p
            out = {}
            if M1 == M0:
                if m0 != 0:
                    _acc(out, {('A', m0): K(1)}, pref * M0 ** (k - 1) / (a + 1))
            else:
                delta = M1 - M0
                for b in range(a + 1):
                    c = pref * binomial(a, b) * (-M0) ** (a - b) / delta ** (a + 1) / (p + b + 1)
                    # [y^(p+b+1)]_{M0}^{M1},  y^(p+b+1) A-normalised: A(m) m^(2(k+b))
                    if m1 != 0:
                        _acc(out, {('A', m1): K(1)}, c * M1 ** (k + b))
                    if m0 != 0:
                        _acc(out, {('A', m0): K(1)}, -c * M0 ** (k + b))
            res[(k, (1,) * a)] = out
        return res

    def terms_of(self, struct, r):
        key = ('terms', struct, r)
        if key in self.cache:
            return self.cache[key]
        k, ms = struct
        seen, out = set(), []
        for perm in itertools.permutations(range(r)):
            fac = [('g',) + tuple(sorted((perm[2 * a], perm[2 * a + 1]))) for a in range(k)]
            fac += [('q', lab, perm[2 * k + b]) for b, lab in enumerate(ms)]
            sig = (tuple(sorted(f for f in fac if f[0] == 'g')), tuple(sorted((f[2], f[1]) for f in fac if f[0] == 'q')))
            if sig not in seen:
                seen.add(sig)
                out.append(fac)
        self.cache[key] = out
        return out

    def contract_terms(self, t1, t2, props):
        slots = {}
        for side, t in ((0, t1), (1, t2)):
            for fi, f in enumerate(t):
                for i in ((f[1], f[2]) if f[0] == 'g' else (f[2],)):
                    slots.setdefault(i, []).append((side, fi))
        tt = (t1, t2)

        def other(node, i):
            f = tt[node[0]][node[1]]
            return (f[2] if f[1] == i else f[1]) if f[0] == 'g' else None
        visited, val = set(), self.K(1)
        for side0 in (0, 1):
            for fi0, f0 in enumerate(tt[side0]):
                if (side0, fi0) in visited or f0[0] != 'q':
                    continue
                visited.add((side0, fi0))
                cur, i = (side0, fi0), f0[2]
                while True:
                    nxt = [s for s in slots[i] if s != cur][0]
                    visited.add(nxt)
                    f = tt[nxt[0]][nxt[1]]
                    if f[0] == 'q':
                        val *= self.dot(props[f0[1]][0], props[f[1]][0])
                        break
                    i, cur = other(nxt, i), nxt
        for side0 in (0, 1):
            for fi0, f0 in enumerate(tt[side0]):
                if (side0, fi0) in visited:
                    continue
                visited.add((side0, fi0))
                cur, i = (side0, fi0), f0[2]
                while True:
                    nxt = [s for s in slots[i] if s != cur][0]
                    if nxt == (side0, fi0):
                        break
                    visited.add(nxt)
                    i, cur = other(nxt, i), nxt
                val *= self.d
        return val

    def tensor(self, props, r):
        key = ('T', tuple((_key(q), m) for q, m in props), r)
        if key in self.cache:
            return self.cache[key]
        if len(props) == 2 and props[1][0] and self.dot(props[1][0], props[1][0]) == 0:
            res = self.tensor2_lightlike(props, r)
            self.cache[key] = res
            return res
        basis = self.structures(props, r)
        terms = [self.terms_of(b, r) for b in basis]
        n = len(basis)
        G = matrix(self.K, n, n, lambda a, b: sum((self.contract_terms(t1, t2, props)
                                                   for t1 in terms[a] for t2 in terms[b]), self.K(0)))
        # invert through the polynomial ring: det and adjugate there are fast (Singular)
        R = self.kin.R
        den = lcm([x.denominator() for x in G.list()]) if n else R(1)
        GR = matrix(R, n, n, [R(x * den) for x in G.list()])
        det = GR.det()
        if det == 0:
            res = self.tensor_feynman(props, r, basis)
            self.cache[key] = res
            return res
        adj = GR.adjugate()
        Ginv = matrix(self.K, n, n, [self.K(a) * self.K(den) / self.K(det) for a in adj.list()])
        rhs = [{k: len(terms[b]) * v for k, v in self.scalar(props, [props[lab][0] for lab in basis[b][1]], basis[b][0]).items()}
               for b in range(n)]
        res = {}
        for a in range(n):
            out = {}
            for b in range(n):
                _acc(out, rhs[b], Ginv[a, b])
            res[basis[a]] = out
        self.cache[key] = res
        return res

    def tensor_feynman(self, props, r, basis):
        r"""
        Tensor coefficients straight from Feynman parameters, for a vanishing Gram determinant:
            T(k, ms) = (-1)^|ms| / prod_{i<k}(d + 2i) * [(-1)^(N+k)] Gamma(k + d/2) Gamma(N - k - d/2)/Gamma(d/2)
                       * e^(eps gamma_E) mu^(2 eps) * Int_simplex prod_{l in ms} x_l  Delta(x)^(d/2 + k - N)
        (the bracket only in Minkowski), Delta = -e sum_{i<j} x_i x_j s_ij + sum_i x_i m_i^2.
        Lines with identical kinematics (s_ij = 0 and s_ik = s_jk for all k) are grouped and their
        parameters integrated out exactly (Dirichlet).  One group left (e.g. all invariants zero):
        divided differences of z^s, exact for any masses.  Two groups: a 1D integral, done exactly
        (Beta functions if Delta is c y^a (1-y)^b, otherwise expanded in eps and integrated).
        Each coefficient becomes a master ('X', n) whose value is stored in _XREG.
        """
        N = len(props)
        K, e = self.K, self.e
        sij = lambda i, j: self.dot(_add(props[i][0], props[j][0], -1), _add(props[i][0], props[j][0], -1))
        groups = []
        for i in range(N):
            for G in groups:
                j = G[0]
                if sij(i, j) == 0 and all(sij(i, kk) == sij(j, kk) for kk in range(N) if kk not in (i, j)):
                    G.append(i)
                    break
            else:
                groups.append([i])
        Msq = [_SR(K(m) ** 2) for _, m in props]
        # the shift P = sum_i x_i q_i written in the basis momenta: q_i = sum_b c_ib q_b
        labs = self.independent(props)
        names = sorted({a for q, _ in props for a in q})
        A = matrix(QQ, [[props[b][0].get(a, 0) for b in labs] for a in names]) if labs else None
        X = PolynomialRing(QQ, ['fsx%d' % i for i in range(N)])
        xs = X.gens()
        lam = {b: X(0) for b in labs}
        for i in range(1, N):
            if not props[i][0]:
                continue
            c = A.solve_right(matrix(QQ, [[props[i][0].get(a, 0)] for a in names]))
            for jb, b in enumerate(labs):
                lam[b] += c[jb, 0] * xs[i]
        res = {}
        for struct in basis:
            k, ms = struct
            poly = prod([lam[lab] for lab in ms], X(1))
            F = SR(0)
            for coef, mon in zip(poly.coefficients(), poly.monomials()):
                n = list(mon.exponents()[0])
                F += coef * _simplex_integral(groups, n, Msq, lambda i, j: -e * _SR(sij(i, j)), N, k)
            pref = (-1) ** len(ms) / prod([(4 - 2 * eps) + 2 * i for i in range(k)], SR(1))
            dh = 2 - eps
            pref *= gamma(k + dh) * gamma(N - k - dh) / gamma(dh)
            if e == 1:
                pref *= (-1) ** (N + k)
            val = pref * F * exp(eps * euler_gamma) * mu ** (2 * eps)
            ser = val.series(eps, 1).truncate().expand()
            c2, c1, c0 = (ser.coefficient(eps, -2), ser.coefficient(eps, -1), ser.coefficient(eps, 0))
            full = c2 / eps ** 2 + c1 / eps + c0
            _XREG.append((full, (c2, c1)))
            res[struct] = {('X', len(_XREG) - 1): K(1)}
        return res

    def contract_tensor(self, T, lv, props):
        r = len(lv)
        total = {}
        for struct, coef in T.items():
            if not coef:
                continue
            for fac in self.terms_of(struct, r):
                val = self.K(1)
                for f in fac:
                    val *= self.dot(lv[f[1]], lv[f[2]]) if f[0] == 'g' else self.dot(props[f[1]][0], lv[f[2]])
                _acc(total, coef, val)
        return total


# ---------------------------------------------------------------------------- masters -> formulas
def _SR(x):
    return SR(str(x))


def _master_formula(key, euclidean):
    """(full value with its poles, (1/eps^2 coefficient, 1/eps coefficient)) in the Package-X
    normalisation.  IR-divergent C0 and D0 come from feynsage.ir."""
    from .ir import ir_coefficients, with_mu
    kind = key[0]
    sign = -1 if euclidean else 1                    # Minkowski invariant = sign * Euclidean one
    zero = SR(0)
    if kind == 'A':
        m = _SR(key[1])
        full, pole = A0(m), (zero, m ** 2)
        n = 1
    elif kind == 'B':
        s, m1, m2 = sign * _SR(key[1]), _SR(key[2]), _SR(key[3])
        full, pole = B0(s, m1, m2), (zero, SR(1))
        if _is_zero(s) and _is_zero(m1) and _is_zero(m2):
            pole = (zero, zero)
        n = 2
    else:
        if kind == 'X':                              # Feynman-parameter value (vanishing Gram), already final
            return _XREG[key[1]]
        if kind == 'C':
            args = [sign * _SR(x) for x in key[1:4]] + [_SR(x) for x in key[4:]]
            n, f = 3, _C0f
        else:
            args = [sign * _SR(x) for x in key[1:7]] + [_SR(x) for x in key[7:]]
            n, f = 4, _D0f
        c = ir_coefficients(kind, args)
        if c is None:
            full, pole = f(*args), (zero, zero)
        else:
            full, pole = with_mu(c), (c[0], c[1] + log(mu ** 2) * c[0])
    if euclidean and n % 2:
        full, pole = -full, (-pole[0], -pole[1])
    return full, pole


def _at4(c, K):
    """c(d = 4), dc/dd(d = 4) and d^2c/dd^2(d = 4) for c in K."""
    R = K.ring()
    dgen = R('d')
    num, den = c.numerator(), c.denominator()
    n0, d0 = num.subs({dgen: 4}), den.subs({dgen: 4})
    if d0 == 0:
        raise ZeroDivisionError("a coefficient has a pole at d = 4; the eps expansion needs more terms")
    n1, d1 = num.derivative(dgen).subs({dgen: 4}), den.derivative(dgen).subs({dgen: 4})
    n2, d2 = num.derivative(dgen, 2).subs({dgen: 4}), den.derivative(dgen, 2).subs({dgen: 4})
    n0, d0, n1, d1, n2, d2 = [K(x) for x in (n0, d0, n1, d1, n2, d2)]
    f0 = n0 / d0
    f1 = (n1 * d0 - n0 * d1) / d0 ** 2
    f2 = (n2 - 2 * f1 * d1 - f0 * d2) / d0             # from (f d)'' = n''
    return f0, f1, f2


def _finalize(dct, K, euclidean):
    """sum_M c_M(d) M  ->  expression in eps (up to eps^0), with the rational terms from
    d = 4 - 2 eps:  c(d) = c(4) - 2 eps c'(4) + 2 eps^2 c''(4)."""
    K = K.K if isinstance(K, Kin) else K
    total = SR(0)
    for key, c in dct.items():
        full, (p2, p1) = _master_formula(key, euclidean)
        c4, dc4, d2c4 = _at4(c, K)
        total += _SR(c4) * full
        if dc4 != 0 and not (p2.is_trivial_zero() and p1.is_trivial_zero()):
            total += -2 * _SR(dc4) * (p2 / eps + p1)
        if d2c4 != 0 and not p2.is_trivial_zero():
            total += 2 * _SR(d2c4) * p2
    return total


def _masters_view(dct, euclidean):
    """sum_M c_M(d) M with the masters as plain symbols (like Package-X's LoopIntegrate)."""
    total = SR(0)
    for key, c in dct.items():
        k = key[0]
        if k == 'X':
            sym = SR.symbol('Xfeyn%d' % key[1], latex_name=r'X_{%d}' % key[1])
            total += _SR(c) * _XREG[key[1]][0]
            continue
        name = {'A': 'A0', 'B': 'B0', 'C': 'C0', 'D': 'D0'}[k]
        sym = SR.symbol('%s(%s)' % (name, ', '.join(str(x) for x in key[1:])))
        total += _SR(c) * sym
    return total


# ---------------------------------------------------------------------------- results
class LoopResult:
    """What loop() returns.  Print it, or use .coefficients(), .masters(), .info()."""

    def __init__(self, parts, master_parts, euclidean, description):
        self.parts = parts                      # [(structure, expression in eps)]
        self.master_parts = master_parts        # [(structure, sum c(d) * master)]
        self.euclidean = euclidean
        self.description = description

    def coefficients(self):
        return dict(self.parts)

    def __getitem__(self, structure):
        return dict(self.parts)[structure]

    def __iter__(self):
        return iter(self.parts)

    def masters(self):
        """The same result before the eps expansion: exact coefficients in d times the masters."""
        return dict(self.master_parts)

    def info(self):
        print(self.description)

    def __repr__(self):
        if not self.parts:
            return "0"
        return "\n".join("%s :  %s" % (st, c) for st, c in self.parts)


_INFO = r"""loop(): a one-loop integral reduced to scalar functions.
  Normalisation (Package-X / LoopTools): mu^(2 eps) e^(eps gamma_E) (4 pi)^(-eps) Int d^d l/(i pi^2),
  d = 4 - 2 eps, propagators 1/(q^2 - m^2 + i0) [Minkowski] or, with euclidean=True,
  Int d^d l/pi^(d/2) with 1/(q^2 + m^2) [Euclidean; the functions are then the Minkowski
  ones at p_M^2 = -p_E^2 times (-1)^N].
  Output: one line per tensor structure (g^{mu nu}, p^mu, ... or 1 for a scalar numerator),
  each with its coefficient expanded to eps^0.  The pieces are
    1/eps^2, 1/eps  poles: UV, and IR (soft/collinear) from divergent C0, D0 written out,
    A0, B0      written out: logs, LogM(z) = log(z - i0), DiscB(s, m1, m2) (see B0?),
    C0(...), D0(...)  IR-finite ones stay symbols; .n() gives > 30 digits, explicit() the dilogarithms,
    mu          the renormalisation scale.
  .masters() gives the exact coefficients in d before expanding; .coefficients() a dict."""


def loop(numerator, *props, kin=None, euclidean=False, explain=False):
    r"""
    One-loop integral  Int numerator / prod D_i, reduced to A0, B0, C0, D0.

        loop("1", ["l", "m"])                                            # A0(m)
        loop("1", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})          # B0(s; m, m)
        loop("l.p", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})
        loop("l^mu l^nu", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})
        loop("l^2", ["l", "m"], ["l + p", "m"], kin={"p^2": "P2"}, euclidean=True)

    A propagator is [momentum, mass]; its momentum contains the loop momentum l once.
    The numerator may contain l^2, l.p, p.q, l^mu, p^mu (an index is a name after ^),
    products and numbers.  kin gives the scalar products of the external momenta,
    for example {"p^2": "s", "p1.p2": "t/2"}; missing ones become symbols (p1p2, p2 for p^2).
    Everything is exact: symbols stay symbols.  explain=True prints what the output means.
    """
    if explain:
        print(_INFO)
    text = numerator.replace('**', '^')
    # a metric factor g(mu,nu) in the numerator (for example from a Dirac trace)
    text = _TOKEN_G.sub(lambda mm: 'fsMET_%s_%s' % tuple(sorted((mm.group(1), mm.group(2)))), text)
    metric_idx = {i for mm in re.finditer(r'fsMET_(\w+?)_(\w+)', text) for i in mm.groups()}
    vec_names = {v for p in props for v in re.findall(r'[A-Za-z_]\w*', str(p[0]))}
    vec_names |= {a for mm in _TOKEN_DOT.finditer(text) for a in mm.groups()}
    vec_names |= {mm.group(1) for mm in _TOKEN_IDX.finditer(text)}
    for key in (kin or {}):
        vec_names |= set(re.findall(r'[A-Za-z_]\w*', key))
    vec_names.discard('l')
    vec_names = sorted(vec_names)
    indices = sorted({mm.group(2) for mm in _TOKEN_IDX.finditer(text)} | metric_idx)
    aux = {i: 'fsA' + i for i in indices}
    fixed, struct_names = {}, {}
    for a_i in indices:
        for a_j in indices:
            if a_i < a_j:
                nm = 'fsG%s_%s' % (a_i, a_j)
                fixed[(aux[a_i], aux[a_j])] = nm
                struct_names[nm] = ("delta^{%s %s}" if euclidean else "g^{%s %s}") % (a_i, a_j)
        for p in vec_names:
            nm = 'fsP%s_%s' % (p, a_i)
            fixed[(aux[a_i], p)] = nm
            struct_names[nm] = "%s^%s" % (p, a_i)
    masses = [SR(m) for _, m in props]
    msyms = {str(x) for m in masses for x in m.variables()}
    # scalar symbols written in the numerator (masses, couplings, d, ...)
    words = set(re.findall(r'\b[A-Za-z]\w*\b', re.sub(r'\^[A-Za-z]\w*', '', _TOKEN_DOT.sub('', text))))
    msyms |= {w for w in words if w not in vec_names and w != 'l' and not w.startswith('fs') and w not in indices}
    K = Kin(vec_names + list(aux.values()), kin, symbols=msyms, euclidean=euclidean, fixed=fixed)
    for v in vec_names:
        SR.var(v)
    fam = []
    for q, m in props:
        v, cl = _vec(str(q), vec_names)
        if cl == -1:
            v = {a: -c for a, c in v.items()}
        elif cl != 1:
            raise ValueError("each propagator must contain the loop momentum l once: %s" % q)
        fam.append((v, K.to_K(m)))
    q0 = fam[0][0]
    fam = tuple((_add(q, q0, -1), m) for q, m in fam)
    eng = _Engine(K)

    ph = {}

    def hold(a, b):
        name = '_dot_%s_%s' % (a, b)
        ph[name] = (a, b)
        return name
    s = re.sub(r'\s*([-+*/^(),])\s*', r'\1', text.strip())     # spaces next to operators
    s = re.sub(r'\s+', '*', s)                                   # the rest mean multiplication
    vecs_l = set(vec_names) | {'l'}
    s = _TOKEN_SQ.sub(lambda mm: ('%s.%s' % (mm.group(1), mm.group(1))) if mm.group(1) in vecs_l else mm.group(0), s)
    s = _TOKEN_DOT.sub(lambda mm: hold(*sorted(mm.groups())), s)
    s = _TOKEN_IDX.sub(lambda mm: hold(*sorted((mm.group(1), aux[mm.group(2)]))), s)
    s = s.replace('^', '**')
    names = {n: SR.var(n) for n in ph}
    names.update({v: SR.var(v) for v in vec_names})
    for w in msyms:
        names.setdefault(w, SR.var(w))
    for mm in re.finditer(r'fsMET_(\w+?)_(\w+)', s):
        a_i, a_j = mm.groups()
        names[mm.group(0)] = SR.var('fsG%s_%s' % tuple(sorted((a_i, a_j))))
    expr = SR(eval(s, {}, names)).expand()
    terms = expr.operands() if expr.operator() is not None and 'add' in str(expr.operator()) else [expr]
    total = {}
    shift = dict(q0)                      # l_old = l_new - q0 (poly_numerator shifts l -> l - qr)
    metric_syms = {names[n]: tuple(n[len('fsMET_'):].split('_', 1)) for n in names if n.startswith('fsMET_')}
    for mm in re.finditer(r'fsMET_(\w+?)_(\w+)', s):
        metric_syms[names[mm.group(0)]] = (mm.group(1), mm.group(2))
    aux_of = {v: k for k, v in aux.items()}           # fsAmu -> mu

    def contract(term):
        """Einstein summation over repeated indices: g(i,i) -> d, g(i,j) v^i -> v^j,
        g(i,j) g(j,k) -> g(i,k), v^i w^i -> v.w.  Returns the contracted term."""
        facs, coef = [], term
        for n, (a, b) in ph.items():
            k = int(term.degree(names[n]))
            if k:
                coef = coef.subs({names[n]: 1})
                facs += [('dot', a, b)] * k
        for sym, (i, j) in metric_syms.items():
            k = int(term.degree(sym))
            if k:
                coef = coef.subs({sym: 1})
                facs += [('g', i, j)] * k

        def idx_of(f):
            if f[0] == 'g':
                return [f[1], f[2]]
            return [aux_of[x] for x in f[1:] if x in aux_of]
        changed = True
        while changed:
            changed = False
            count = {}
            for f in facs:
                for i in idx_of(f):
                    count[i] = count.get(i, 0) + 1
            for i, c in count.items():
                if c < 2:
                    continue
                holders = [f for f in facs if i in idx_of(f)]
                f1 = holders[0]
                if f1[0] == 'g' and f1[1] == f1[2] == i:
                    facs.remove(f1); coef *= SR.var('d'); changed = True; break
                f2 = holders[1]
                facs.remove(f1); facs.remove(f2)

                def other(f):             # what an index-i factor attaches to on its other end
                    if f[0] == 'g':
                        return ('idx', f[2] if f[1] == i else f[1])
                    return ('vec', f[2] if f[1] == aux[i] else f[1])
                o1, o2 = other(f1), other(f2)
                if o1[0] == 'vec' and o2[0] == 'vec':
                    facs.append(('dot',) + tuple(sorted((o1[1], o2[1]))))
                elif o1[0] == 'idx' and o2[0] == 'idx':
                    if o1[1] == o2[1]:
                        coef *= SR.var('d')
                    else:
                        facs.append(('g',) + tuple(sorted((o1[1], o2[1]))))
                else:
                    v = o1[1] if o1[0] == 'vec' else o2[1]
                    j = o2[1] if o1[0] == 'vec' else o1[1]
                    facs.append(('dot',) + tuple(sorted((v, aux[j]))))
                changed = True
                break
        out = coef
        for f in facs:
            if f[0] == 'dot':
                nm = hold(f[1], f[2])
                if nm not in names:
                    names[nm] = SR.var(nm)
                out *= names[nm]
            else:
                a_i, a_j = sorted((f[1], f[2]))
                out *= SR.var('fsG%s_%s' % (a_i, a_j))
        return out

    new_terms = []
    for term in terms:
        t = contract(term).expand()
        new_terms += t.operands() if t.operator() is not None and 'add' in str(t.operator()) else [t]
    terms = new_terms
    for term in terms:
        coef, lv, l2 = term, [], 0
        for n, (a, b) in ph.items():
            k = int(term.degree(names[n]))
            if not k:
                continue
            coef = coef.subs({names[n]: 1})
            if a == 'l' and b == 'l':
                l2 += k
            elif 'l' in (a, b):
                lv += [{(b if a == 'l' else a): 1}] * k
            else:
                coef *= SR(str(K.dot({a: 1}, {b: 1}))) ** k
        c = K.to_K(coef)
        part = eng.poly_numerator(fam, lv, l2, shift) if q0 else eng.scalar(fam, lv, l2)
        _acc(total, part, c)
    # split by tensor structure: the structure symbols enter polynomially in the numerators
    R = K.R
    sgens = [R(n) for n in struct_names]
    split = {}
    sidx = [R.gens().index(g) for g in sgens]
    for key, c in total.items():
        num, den = c.numerator(), c.denominator()
        groups = {}
        uni = R.ngens() == 1                     # only d is a symbol: keys are plain integers
        for mon, cf in num.dict().items():
            mon = (mon,) if uni else mon
            exps = tuple(mon[i] for i in sidx)
            rest = list(mon)
            for i in sidx:
                rest[i] = 0
            groups.setdefault(exps, {})[rest[0] if uni else tuple(rest)] = cf
        for exps, mons in groups.items():
            split.setdefault(exps, {})[key] = K.K(R(mons)) / K.K(den)

    def label(exps):
        out = []
        for g, k in zip(struct_names.values(), exps):
            out += [g] * k
        return " ".join(out) if out else "1"
    parts, mparts = [], []
    for exps, dct in sorted(split.items(), key=lambda t: label(t[0])):
        val = _finalize(dct, K, euclidean)
        if not val.is_trivial_zero():
            parts.append((label(exps), val))
            mparts.append((label(exps), _masters_view(dct, euclidean)))
    return LoopResult(parts, mparts, euclidean, _INFO)


_TOKEN_G = re.compile(r'\bg\(\s*(\w+)\s*,\s*(\w+)\s*\)')
_TOKEN_IDX = re.compile(r'\b([A-Za-z]\w*)\^([A-Za-z]\w*)\b')
_TOKEN_SQ = re.compile(r'(?<![.\w])([A-Za-z]\w*)\^2\b')     # not after a dot: p.pp^2 is (p.pp)^2
_TOKEN_DOT = re.compile(r'\b([A-Za-z]\w*)\.([A-Za-z]\w*)\b')


# ---------------------------------------------------------------------------- Package-X style names
def _masses_syms(*xs):
    return {str(v) for x in xs for v in SR(x).variables()}


def PVB(r, n, s, m1, m2, explain=False):
    """PV coefficient B_{00..(2r) 11..(n)}(s; m1, m2) as Package-X's PVB[r, n, s, m1, m2]:
    PVB(0,0) = B0, PVB(0,1) = B1, PVB(1,0) = B00, PVB(0,2) = B11.  Exact, expanded to eps^0."""
    if explain:
        print(PVB.__doc__)
    K = Kin(['p'], {'p^2': s}, symbols=_masses_syms(s, m1, m2))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m1)), ({'p': 1}, K.to_K(m2))), 2 * r + n)
    return _finalize(T[(r, (1,) * n)], K, False)


def PVC(r, n1, n2, s1, s12, s2, m0, m1, m2, explain=False):
    """PV coefficient C_{00..1..2..}(s1, s12, s2; m0, m1, m2) as Package-X's PVC[r, n1, n2, ...].
    Propagators (l, m0), (l + p1, m1), (l + p2, m2); s1 = p1^2, s2 = p2^2, s12 = (p1 - p2)^2."""
    if explain:
        print(PVC.__doc__)
    K = Kin(['p1', 'p2'], {'p1^2': s1, 'p2^2': s2, 'p1.p2': (SR(s1) + SR(s2) - SR(s12)) / 2},
            symbols=_masses_syms(s1, s12, s2, m0, m1, m2))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m0)), ({'p1': 1}, K.to_K(m1)), ({'p2': 1}, K.to_K(m2))), 2 * r + n1 + n2)
    return _finalize(T[(r, (1,) * n1 + (2,) * n2)], K, False)


def PVD(r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, explain=False):
    """PV coefficient D_{00..1..2..3..} as Package-X's PVD[r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3].
    Propagators (l, m0), (l + p1, m1), (l + p2, m2), (l + p3, m3); s1 = p1^2, s2 = (p2 - p1)^2,
    s3 = (p3 - p2)^2, s4 = p3^2, s12 = p2^2, s23 = (p3 - p1)^2.  Exact, expanded to eps^0
    (IR poles explicit when the box is soft or collinear divergent)."""
    if explain:
        print(PVD.__doc__)
    s1, s2, s3, s4, s12, s23 = [SR(x) for x in (s1, s2, s3, s4, s12, s23)]
    K = Kin(['p1', 'p2', 'p3'], {'p1^2': s1, 'p2^2': s12, 'p3^2': s4, 'p1.p2': (s1 + s12 - s2) / 2,
                                 'p2.p3': (s12 + s4 - s3) / 2, 'p1.p3': (s1 + s4 - s23) / 2},
            symbols=_masses_syms(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m0)), ({'p1': 1}, K.to_K(m1)), ({'p2': 1}, K.to_K(m2)), ({'p3': 1}, K.to_K(m3))),
                   2 * r + n1 + n2 + n3)
    return _finalize(T[(r, (1,) * n1 + (2,) * n2 + (3,) * n3)], K, False)


def explicit(expr):
    r"""
    Replace every C0(...) and D0(...) with numerical arguments in expr by its closed form in
    dilogarithms and logarithms (scalar.c0_closed, scalar.d0_closed): exact algebraic
    arguments, branches fixed by the +i0.  Symbolic arguments, and the rare exceptional points
    without a regular closed form, are left as they are.
    """
    from .scalar import c0_closed, d0_closed

    def rep_c(*args):
        if all(SR(a).is_numeric() for a in args):
            try:
                return c0_closed(*_exact_args(args))
            except ZeroDivisionError:
                pass
        return _C0f(*args)

    def rep_d(*args):
        if all(SR(a).is_numeric() for a in args):
            try:
                return d0_closed(*_exact_args(args))
            except ZeroDivisionError:
                pass
        return _D0f(*args)
    return SR(expr).substitute_function(_C0f, rep_c).substitute_function(_D0f, rep_d)
