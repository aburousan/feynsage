r"""
One-loop tensor integrals and Passarino-Veltman functions, in the conventions of
Package-X and LoopTools.

Conventions
-----------
Minkowski metric (+,-,-,-).  A propagator (q, m) is 1/(q^2 - m^2 + i0).  The measure is

    mu^(2 eps) e^(eps gamma_E) * (1/(i pi^(d/2))) * Int d^d l ,   d = 4 - 2 eps
    (the same as (16 pi^2/i) mu^(2 eps) (e^gamma_E/(4 pi))^eps Int d^d l/(2 pi)^d, the MS-bar convention)

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
from ._parse import sr

eps = SR.var('eps', latex_name=r'\epsilon')
mu = SR.var('mu')
Dim = var('d')                      # the space-time dimension while reducing, d = 4 - 2 eps at the end

# ---------------------------------------------------------------------------- special functions
def _mp():
    import mpmath
    return mpmath


def _prec(parent):
    return getattr(parent, 'prec', lambda: 53)() if parent is not None else 53


def _mpf(x, prec):
    """An exact or numerical real argument as an mpmath number with prec bits, without passing
    through a double."""
    from sage.all import RealField
    mp = _mp()
    return mp.mpf(RealField(prec)(sr(x)).str(truncate=False)) if not isinstance(x, (int, float)) else mp.mpf(x)


def _cplx(z):
    return complex(z)


def _logm_evalf(self, z, parent=None, algorithm=None):
    """log(z - i0) for real z, to the precision of parent."""
    mp = _mp()
    prec = _prec(parent)
    with mp.workprec(prec + 20):
        try:
            x = _mpf(z, prec + 20)
            val = mp.log(-x) - 1j * mp.pi if x < 0 else mp.log(x)
        except (TypeError, ValueError):
            val = mp.log(mp.mpc(complex(z)))
        return _to_parent(val, parent)


LogM = function('LogM', nargs=1, evalf_func=_logm_evalf, derivative_func=lambda self, z, diff_param=None: 1 / z,
                print_latex_func=lambda self, z: r"\ln\left(%s - i0\right)" % z._latex_())


def _discb_num(s, m1, m2, prec=53):
    """DiscB(s; m1, m2) with s -> s + i0, to prec bits: the i0 is a positive imaginary part far
    below the working precision."""
    mp = _mp()
    with mp.workprec(prec + 40):
        s, m1, m2 = _mpf(s, prec + 40), _mpf(m1, prec + 40), _mpf(m2, prec + 40)
        if s == 0 and m1 ** 2 == m2 ** 2 and m1 != 0:
            return mp.mpc(-2, 0)                    # the limit s -> 0 of beta log((beta - 1)/(beta + 1))
        z = s + 1j * mp.mpf(2) ** (-prec - 30) * max(1, abs(s), m1 ** 2, m2 ** 2)
        lam = z ** 2 + m1 ** 4 + m2 ** 4 - 2 * z * m1 ** 2 - 2 * z * m2 ** 2 - 2 * m1 ** 2 * m2 ** 2
        r = mp.sqrt(lam)
        val = r / z * mp.log((m1 ** 2 + m2 ** 2 - z + r) / (2 * m1 * m2))
        # below the threshold (m1 + m2)^2 DiscB is real: drop the rounding left by the tiny i0
        return mp.mpc(val.real, 0) if s < (abs(m1) + abs(m2)) ** 2 else val


def _discb_evalf(self, s, m1, m2, parent=None, algorithm=None):
    return _to_parent(_discb_num(s, m1, m2, _prec(parent)), parent)


def _discb_deriv(self, s, m1, m2, diff_param=None):
    r"""Derivatives of DiscB = sqrt(lam)/s log((m1^2 + m2^2 - s + sqrt(lam))/(2 m1 m2)).  In the
    mass derivatives the terms odd in sqrt(lam) cancel, so they hold on every branch:
        d/dm1 DiscB = 2 m1 (m1^2 - s - m2^2)/lam DiscB + (m1^2 + s - m2^2)/(m1 s)."""
    lam = s ** 2 + m1 ** 4 + m2 ** 4 - 2 * s * m1 ** 2 - 2 * s * m2 ** 2 - 2 * m1 ** 2 * m2 ** 2
    if diff_param == 0:
        dlam = 2 * s - 2 * m1 ** 2 - 2 * m2 ** 2
        return (s * dlam - 2 * lam) / (2 * s * lam) * DiscB(s, m1, m2) - 1 / s
    a, b = (m1, m2) if diff_param == 1 else (m2, m1)
    return 2 * a * (a ** 2 - s - b ** 2) / lam * DiscB(s, m1, m2) + (a ** 2 + s - b ** 2) / (a * s)


DiscB = function('DiscB', nargs=3, evalf_func=_discb_evalf, derivative_func=_discb_deriv,
                 print_latex_func=lambda self, s, a, b: r"\Lambda(%s; %s, %s)" % (s._latex_(), a._latex_(), b._latex_()))


def _exact_args(args):
    """Numbers -> exact rationals (floats by their simplest rational), exact algebraic numbers kept."""
    out = []
    for a in args:
        a = sr(a)
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


def _master_deriv(kind, args, idx):
    r"""
    d/d(argument idx) of C0(s1, s12, s2, m0, m1, m2) or D0(s1, ..., s23, m0, ..., m3), written
    again in A0, B0, C0, D0.  A mass derivative is a raised propagator, d/dm_i = 2 m_i d/dm_i^2.
    An invariant derivative comes from q_k.d/dq_l I = -2 Int (l + q_l).q_k / (D_l^2 prod_(j != l) D_j),
    a linear system for the derivatives in the scalar products G_ab = q_a.q_b.  Both kinds of
    integral reduce to the masters by integration by parts.
    """
    N = 3 if kind == 'C' else 4
    args = [sr(a) for a in args]
    sub = {}
    safe = []
    for i, a in enumerate(args):                      # numbers outside Q(symbols) get a name
        try:
            if a.is_numeric():
                QQ(a)
            else:
                PolynomialRing(QQ, [str(v) for v in a.variables()]).fraction_field()(str(a))
            safe.append(a)
        except (TypeError, ValueError):
            z = SR.var('fsZ%d' % i)
            sub[z] = a
            safe.append(z)
    if kind == 'C':
        s1, s12, s2, *ms = safe
        vecs = ['p1', 'p2']
        kin = {'p1^2': s1, 'p2^2': s2, 'p1.p2': (s1 + s2 - s12) / 2}
        pair = {0: (0, 1), 1: (1, 2), 2: (0, 2)}
    else:
        s1, s2, s3, s4, s12, s23, *ms = safe
        vecs = ['p1', 'p2', 'p3']
        kin = {'p1^2': s1, 'p2^2': s12, 'p3^2': s4, 'p1.p2': (s1 + s12 - s2) / 2,
               'p2.p3': (s12 + s4 - s3) / 2, 'p1.p3': (s1 + s4 - s23) / 2}
        pair = {0: (0, 1), 1: (1, 2), 2: (2, 3), 3: (0, 3), 4: (0, 2), 5: (1, 3)}
    syms = {str(v) for a in safe for v in sr(a).variables()}
    kn = Kin(vecs, kin, symbols=syms)
    K = kn.K
    eng = _Engine(kn)
    props = tuple([({}, kn.to_K(ms[0]))] + [({v: 1}, kn.to_K(ms[i + 1])) for i, v in enumerate(vecs)])
    if idx >= N * (N - 1) // 2:                       # a mass
        i = idx - N * (N - 1) // 2
        dct = {k: 2 * kn.to_K(ms[i]) * v for k, v in eng.master(props + (props[i],)).items()}
    else:
        n = N - 1
        q = [props[a][0] for a in range(N)]
        G = lambda a, b: eng.dot(q[a], q[b])
        unk = [(a, b) for a in range(1, N) for b in range(a, N)]
        rows, rhs = [], []
        for k in range(1, N):
            for l in range(1, N):
                dup = props + (props[l],)
                o = {}
                _acc(o, eng.scalar(dup, [q[k]], 0), K(-2))
                _acc(o, eng.scalar(dup, [], 0), -2 * G(l, k))
                rows.append([(G(k, b) if a == l else 0) + (G(k, a) if b == l else 0) for a, b in unk])
                rhs.append(o)
        keys = sorted({kk for o in rhs for kk in o}, key=str)
        pick, sel = [], []                            # independent rows: each row is an exact identity,
        for i, row in enumerate(rows):                # but one integral can carry two master labels
            if matrix(K, sel + [row]).rank() == len(sel) + 1:
                pick.append(i)
                sel.append(row)
        M = matrix(K, sel)
        B = matrix(K, len(pick), len(keys), [rhs[i].get(kk, 0) for i in pick for kk in keys])
        X = M.solve_right(B)
        f = {u: {kk: X[j, c] for c, kk in enumerate(keys) if X[j, c] != 0} for j, u in enumerate(unk)}
        a0, b0 = pair[idx]
        dct = {}
        if a0 == 0:                                   # s_0b = G_bb, and G_bc = (s_0b + s_0c - s_bc)/2
            _acc(dct, f[(b0, b0)], K(1))
            for c in range(1, N):
                if c != b0:
                    _acc(dct, f[(min(b0, c), max(b0, c))], K(1) / 2)
        else:
            _acc(dct, f[(a0, b0)], K(-1) / 2)
    out = _merge_mu(_finalize(dct, kn, False))
    return out.subs(sub) if sub else out


_C0f = function('C0', nargs=6, evalf_func=_c0_evalf)
_D0f = function('D0', nargs=10, evalf_func=_d0_evalf)
_DERIV_CACHE = {}


def _occurrences(expr):
    """The distinct C0(...) and D0(...) inside expr, as (kind, function call)."""
    out = {}
    for kind, f, n in (('C', _C0f, 6), ('D', _D0f, 10)):
        for hit in sr(expr).find(f(*[SR.wild(i) for i in range(n)])):
            out[hit] = kind
    return out


def loop_diff(expr, x):
    r"""
    d expr / dx for an expression with A0, B0 (logs, LogM, DiscB), C0 and D0, written again in the
    same functions.  The derivatives of C0 and D0 are exact (integration by parts, _master_deriv):

        loop_diff(B0(s, 1, 2), s)       loop_diff(C0(s, 0, 5, 1, 1, 1), s)
    """
    expr = sr(expr)
    occ = _occurrences(expr)
    if not occ:
        return _tidy(expr.diff(x))
    syms = {f: SR.var('fsF%d' % i) for i, f in enumerate(occ)}
    plain = expr.subs(syms)
    out = plain.diff(x)
    for f, kind in occ.items():
        args = f.operands()
        df = 0
        for j, a in enumerate(args):
            da = sr(a).diff(x)
            if da.is_zero():
                continue
            key = (kind, tuple(args), j)
            if key not in _DERIV_CACHE:
                _DERIV_CACHE[key] = _master_deriv(kind, args, j)
            df += _DERIV_CACHE[key] * da
        out += plain.diff(syms[f]) * df
    return _tidy(_merge_mu(out.subs({v: f for f, v in syms.items()})))


def _tidy(expr, max_terms=400):
    """Collect expr by its special functions (DiscB, LogM, C0, D0) and simplify each rational coefficient,
    e.g. ((s - 5) s - s^2 + 10 s - 9) DiscB/... -> (5 s - 9) DiscB/....  Left as it is when large or if
    anything fails: this only changes how the result looks."""
    try:
        e = sr(expr).expand()
        terms = e.operands() if e.operator() is not None and 'add' in e.operator().__name__ else [e]
        if len(terms) > max_terms:
            return expr
        funcs = []
        for fn, n in ((DiscB, 3), (LogM, 1), (_C0f, 6), (_D0f, 10)):
            funcs += e.find(fn(*[SR.wild(i) for i in range(n)]))
        out, rest = sr(0), e
        for F in funcs:
            c = e.coefficient(F)
            if c.is_zero() or c.has(F):
                continue
            out += c.simplify_rational().factor() * F
            rest = rest - c * F
        return out + rest.expand().simplify_rational()
    except Exception:
        return expr


def _merge_mu(expr):
    """log(c mu^2) -> log(c) + log(mu^2), so that the mu dependence that cancels does cancel."""
    w = SR.wild(0)
    return sr(expr).subs(log(w * mu ** 2) == log(w) + log(mu ** 2)).expand()


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
    x = sr(x)
    return x.is_trivial_zero() or (x.is_numeric() and x == 0)


def A0(m):
    """Tadpole A0(m) = m^2 (1/eps + log(mu^2/m^2) + 1) in the Package-X normalisation; A0(0) = 0.
    Output: a symbolic expression in eps (the 1/eps UV pole) and mu (the renormalisation scale)."""
    m = sr(m)
    if _is_zero(m):
        return sr(0)
    return m ** 2 * (1 / eps + log(mu ** 2 / m ** 2) + 1)


def B0(s, m1, m2):
    """Two-point function B0(p^2 = s; m1, m2) in the Package-X / LoopTools normalisation.
    Output: 1/eps + finite part, the finite part in closed form (as Package-X's LoopRefine):
    logs, LogM(z) = log(z - i0) and DiscB(s, m1, m2) = sqrt(lam)/s * log((m1^2+m2^2-s+sqrt(lam))/(2 m1 m2))
    with s -> s + i0.  Numbers in, a number with .n(); symbols in, a formula (masses assumed nonzero)."""
    s, m1, m2 = sr(s), sr(m1), sr(m2)
    z1, z2, zs = _is_zero(m1), _is_zero(m2), _is_zero(s)
    if z1 and not z2:
        m1, m2, z1, z2 = m2, m1, z2, z1          # put the massless line second
    if zs:
        if z1 and z2:
            return sr(0)                          # scaleless
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
    return sr(expr).expand().coefficient(eps, -1)


def pole_parts(expr):
    """(coefficient of 1/eps^2, coefficient of 1/eps)."""
    e = sr(expr).expand()
    return e.coefficient(eps, -2), e.coefficient(eps, -1)


def finite_part(expr, mu_value=None):
    """The eps^0 part; with mu_value the scale mu is set to that number."""
    e = sr(expr).expand().coefficient(eps, 0)
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
            raw[tuple(sorted((a, b)))] = sr(val)
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
        return self.K(str(sr(x)))

    def dot(self, u, v):
        K = self.K
        return sum((K(cu) * K(cv) * self.sp[tuple(sorted((a, b)))] for a, cu in u.items() for b, cv in v.items()), K(0))


def _vec(s, vectors):
    """'l + 2 p1 - p2' -> ({'p1': 2, 'p2': -1}, coefficient of l)."""
    e = sr(s)
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
    return dd(tuple(sr(z) for z in nodes))


def _manifestly_nonneg(e):
    """Every term of the expanded e has a positive coefficient and even powers of all symbols."""
    e = sr(e).expand()
    if e.is_numeric():
        return bool(e >= 0)
    terms = e.operands() if 'add' in getattr(e.operator(), '__name__', '') else [e]
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
                return sr(0)                 # 0^(c - eps) = 0 in dimensional regularisation
            return z ** (s + Kd - j) / prod([s + i for i in range(1, Kd - j + 1)], sr(1))
        nodes = [Msq[i] for i in range(N) for _ in range(n[i] + 1)]
        return prod([factorial(x) for x in n], sr(1)) * _divided_difference(nodes, deriv)
    for G in groups:
        if any(not (Msq[i] - Msq[G[0]]).is_trivial_zero() for i in G):
            raise ZeroDivisionError("vanishing Gram determinant with unequal masses on kinematically identical "
                                    "lines: not covered; move the kinematics slightly")
    if len(groups) != 2:
        raise ZeroDivisionError("vanishing Gram determinant with more than two independent line groups: "
                                "not covered; move the kinematics slightly")
    weight, alpha = sr(1), []
    for G in groups:
        nG = sum(n[i] for i in G)
        weight *= prod([factorial(n[i]) for i in G], sr(1)) / factorial(nG + len(G) - 1)
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
    out = sr(0)
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
        self.use_pv = True                           # False: projection only (kept as a check)

    def dot(self, u, v):
        return self.kin.dot(u, v)

    def msq(self, m):
        return self.K(m) ** 2

    # masters ----------------------------------------------------------------
    def master(self, props):
        groups = {}
        for q, mm in props:
            groups.setdefault((_key(q), mm), []).append((q, mm))
        if len(groups) < len(props):                 # a propagator raised to a power
            return self.weighted(tuple(v[0] for v in groups.values()), tuple(len(v) for v in groups.values()))
        N = len(props)
        if N == 5:
            return self.pentagon(props)
        if N > 5:
            raise NotImplementedError("more than five propagators")
        q = [p[0] for p in props]
        m = [p[1] for p in props]
        sq = lambda a, b: self.dot(_add(q[a], q[b], -1), _add(q[a], q[b], -1))
        if N == 1:                                   # A0(0) = 0: scaleless
            return {('A', m[0]): self.K(1)} if m[0] != 0 else {}
        if N == 2:
            a, b = sorted([m[0], m[1]], key=str)
            if a == 0 and b == 0 and sq(0, 1) == 0:      # scaleless bubble
                return {}
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

    def pentagon(self, props):
        r"""
        The scalar pentagon in d = 4 - 2 eps as a sum of the five boxes with one line removed (Melrose;
        Bern, Dixon and Kosower):  E0 = -(1/2) sum_i b_i D0^(i) + O(eps),  b = S^-1 (1, ..., 1),
        S_ij = (m_i^2 + m_j^2 - s_ij)/2 the modified Cayley matrix.  The O(eps) term multiplies a
        finite six-dimensional pentagon, so it is dropped.
        """
        if self.e != 1:
            raise NotImplementedError("pentagons in Euclidean kinematics")
        key = ('P5', tuple((_key(q), m) for q, m in props))
        if key in self.cache:
            return self.cache[key]
        q = [p[0] for p in props]
        msq = [self.msq(p[1]) for p in props]
        sq = lambda a, b: self.dot(_add(q[a], q[b], -1), _add(q[a], q[b], -1))
        S = matrix(self.K, 5, 5, lambda i, j: (msq[i] + msq[j] - (sq(i, j) if i != j else 0)) / 2)
        if S.det() == 0:
            raise ZeroDivisionError("pentagon with a singular Cayley matrix")
        b = S.inverse() * matrix(self.K, 5, 1, [1] * 5)
        out = {}
        for i in range(5):
            sub = props[:i] + props[i + 1:]
            q0 = sub[0][0]
            sub = tuple((_add(qq, q0, -1), mm) for qq, mm in sub)
            _acc(out, self.master(sub), -b[i, 0] / 2)
        self.cache[key] = out
        return out

    def tensor_pent(self, props, r):
        r"""
        Pentagon tensor coefficients in four dimensions (Denner and Dittmaier, Nucl. Phys. B 658
        (2003) 175): the four momenta span space-time, so l^mu = sum_kl q_k^mu Zinv_kl 2 l.q_l and
        2 l.q_l = D_l - D_0 - f_l.  With a formal vector x and T(x) = T^{mu..} x_mu..:
            T_r(x) = sum_kl (q_k.x) Zinv_kl [T^(l)_(r-1)(x) - T^(0)_(r-1)(x) - f_l T_(r-1)(x)].
        Valid up to O(eps) while the pentagon is ultraviolet finite, rank <= 5.
        """
        if r == 0:
            return {(0, ()): self.master(props)}
        if r > 5:
            raise NotImplementedError("pentagon tensors of rank > 5 (ultraviolet divergent)")
        K, e = self.K, self.e
        Zi = self.zinv(props)
        if Zi is None:
            raise ZeroDivisionError("pentagon with a vanishing Gram determinant")
        m0sq = self.msq(props[0][1])
        f = {k: self.dot(props[k][0], props[k][0]) - e * self.msq(props[k][1]) + e * m0sq for k in range(1, 5)}

        def W(rr, n, ms):
            return QQ(factorial(rr)) / (2 ** n * factorial(n) * prod([factorial(ms.count(i)) for i in set(ms)]))

        def poly(T, rr):                     # {(n, ms)} -> {(n, m_1..m_4)} with the weights W
            out = {}
            for (n, ms), c in T.items():
                mon = (n,) + tuple(ms.count(i) for i in range(1, 5))
                _acc(out.setdefault(mon, {}), c, W(rr, n, ms))
            return out
        T1 = poly(self.tensor(props, r - 1), r - 1)
        R0 = poly(self.removed_0(props, r - 1), r - 1)
        acc = {}
        for l in range(1, 5):
            Pl = poly(self.removed_k(props, l, r - 1), r - 1)
            br = {}
            for mon, c in Pl.items():
                _acc(br.setdefault(mon, {}), c, K(1))
            for mon, c in R0.items():
                _acc(br.setdefault(mon, {}), c, K(-1))
            for mon, c in T1.items():
                _acc(br.setdefault(mon, {}), c, -f[l])
            for k in range(1, 5):
                z = Zi[k - 1, l - 1]
                if z == 0:
                    continue
                for mon, c in br.items():
                    new = list(mon)
                    new[k] += 1
                    _acc(acc.setdefault(tuple(new), {}), c, z)
        res = {}
        for mon, c in acc.items():
            if c:
                n = mon[0]
                ms = tuple(i for i in range(1, 5) for _ in range(mon[i]))
                res[(n, ms)] = {kk: v / W(r, n, ms) for kk, v in c.items()}
        return res

    def weighted(self, props, weights):
        r"""
        A scalar integral with propagator i raised to weights[i]: reduced by integration by parts
        (feynsage's own Laporta reducer, exact in d) to integrals with every power 1, which are
        then the ordinary masters A0, B0, C0, D0.
        """
        key = ('W', tuple((_key(q), m) for q, m in props), weights)
        if key in self.cache:
            return self.cache[key]
        from .easy import family, ibp_reduce
        names = sorted({a for q, _ in props for a in q})
        # write the momenta in a basis of their span, so the family is complete
        basis, rows = [], []
        for q, _ in props:
            row = [q.get(a, 0) for a in names]
            if any(row) and matrix(QQ, rows + [row]).rank() == len(rows) + 1:
                basis.append(q)
                rows.append(row)
        A = matrix(QQ, rows).transpose() if rows else None

        def mom(q):
            if not q:
                return "l"
            c = A.solve_right(matrix(QQ, [[q.get(a, 0)] for a in names]))
            return " + ".join(["l"] + ["(%s)*fsr%d" % (c[j, 0], j) for j in range(len(basis)) if c[j, 0] != 0])
        kin = {}
        for i in range(len(basis)):
            for j in range(i, len(basis)):
                kin[("fsr%d^2" % i) if i == j else ("fsr%d.fsr%d" % (i, j))] = str(self.dot(basis[i], basis[j]))
        fam = family([(mom(q), str(m)) for q, m in props], kin=kin or None, loops=["l"],
                     euclidean=(self.e == -1), name="W")
        red = ibp_reduce(fam, [tuple(weights)])
        out = {}
        for mst, c in red[tuple(weights)].items():
            if any(a not in (0, 1) for a in mst):
                raise NotImplementedError("IBP left a master with a raised power: %s" % (mst,))
            sub = [props[i] for i, a in enumerate(mst) if a == 1]
            q0 = sub[0][0]
            sub = tuple((_add(q, q0, -1), m) for q, m in sub)
            _acc(out, self.master(sub), self.K(str(c)))
        self.cache[key] = out
        return out

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

    # Passarino-Veltman recursion ----------------------------------------------------
    def zinv(self, props):
        """Inverse of the Gram matrix Z_kl = 2 q_k.q_l (k, l = 1..N-1), or None if it is singular."""
        key = ('Z', tuple((_key(q), m) for q, m in props))
        if key not in self.cache:
            n = len(props) - 1
            Z = matrix(self.K, n, n, lambda a, b: 2 * self.dot(props[a + 1][0], props[b + 1][0]))
            R = self.kin.R
            den = lcm([x.denominator() for x in Z.list()])
            ZR = matrix(R, n, n, [R(x * den) for x in Z.list()])
            det = ZR.det()
            if det == 0:
                self.cache[key] = None
            else:
                adj = ZR.adjugate()
                self.cache[key] = matrix(self.K, n, n, [self.K(a) * self.K(den) / self.K(det) for a in adj.list()])
        return self.cache[key]

    def removed_k(self, props, k, r):
        """Tensor coefficients of rank r of the integral without propagator k (k >= 1), labelled
        by the momenta of props."""
        key = ('Rk', tuple((_key(q), m) for q, m in props), k, r)
        if key in self.cache:
            return self.cache[key]
        lab = [i for i in range(1, len(props)) if i != k]
        T = self.tensor(props[:k] + props[k + 1:], r)
        res = {(n, tuple(sorted(lab[a - 1] for a in ms))): c for (n, ms), c in T.items()}
        self.cache[key] = res
        return res

    def removed_0(self, props, r):
        r"""Tensor coefficients of rank r of the integral without propagator 0, labelled by the
        momenta of props.  Its loop momentum is shifted, l = l' - q_1, so with a formal vector x
        T(x) = sum_j binom(r, j) (-q_1.x)^j T'_(r-j)(x), where T' has the momenta q_a - q_1.
        A coefficient c of {g^n q_i...} contributes c W(r, n, ms) (x^2)^n prod (q_i.x)^m_i,
        W = r!/(2^n n! prod m_i!), which turns coefficients into polynomials and back."""
        key = ('R0', tuple((_key(q), m) for q, m in props), r)
        if key in self.cache:
            return self.cache[key]
        N = len(props)
        q1 = props[1][0]
        new = tuple((_add(q, q1, -1), m) for q, m in props[1:])
        P = PolynomialRing(QQ, ['fsX'] + ['fsy%d' % i for i in range(1, N)])
        X, y = P.gen(0), P.gens()[1:]

        def W(rr, n, ms):
            return QQ(factorial(rr)) / (2 ** n * factorial(n) * prod([factorial(ms.count(i)) for i in set(ms)]))
        acc = {}
        for j in range(r + 1):
            for (n, ms), c in self.tensor(new, r - j).items():
                if not c:
                    continue
                poly = binomial(r, j) * (-y[0]) ** j * W(r - j, n, ms) * X ** n
                for a in ms:                                   # new label a is q_(a+1) - q_1
                    poly *= y[a] - y[0]
                for mon, cf in zip(poly.monomials(), poly.coefficients()):
                    _acc(acc.setdefault(mon, {}), c, QQ(cf))
        res = {}
        for mon, c in acc.items():
            ex = mon.exponents()[0]
            n = ex[0]
            ms = tuple(i for i in range(1, N) for _ in range(ex[i]))
            if c:
                res[(n, ms)] = {k: v / W(r, n, ms) for k, v in c.items()}
        self.cache[key] = res
        return res

    def tensor_pv(self, props, r):
        r"""
        Tensor coefficients by the Passarino-Veltman recursion (Denner and Dittmaier,
        Nucl. Phys. B 734 (2006) 62, eqs. (5.10) and (5.11)), with f_k = q_k^2 - e m_k^2 + e m_0^2:
            T_{00 i3..ir} = [2 e m0^2 T_{i3..ir} + T^(0)_{i3..ir} + sum_k f_k T_{k i3..ir}] / (2(d + r - N - 1))
            T_{k i2..ir}  = sum_l Zinv_kl [T^(l)_{i2..ir} - T^(0)_{i2..ir} - f_l T_{i2..ir}
                                           - 2 sum_s delta_{l i_s} T_{00 i2..(no i_s)..ir}]
        T^(l) is the integral without propagator l.  Needs a non-singular Gram matrix.
        """
        N = len(props)
        K, e, d = self.K, self.e, self.d
        if r == 0:
            return {(0, ()): self.master(props)}
        if N == 1 and r % 2:
            return {}
        Zi = self.zinv(props) if N > 1 else None
        m0sq = self.msq(props[0][1])
        f = {k: self.dot(props[k][0], props[k][0]) - e * self.msq(props[k][1]) + e * m0sq for k in range(1, N)}
        T1 = self.tensor(props, r - 1)
        T2 = self.tensor(props, r - 2) if r >= 2 else {}
        labels = list(range(1, N))
        res = {}
        structs = [(n, ms) for n in range(r // 2, -1, -1)
                   for ms in itertools.combinations_with_replacement(labels, r - 2 * n)]
        for n, ms in structs:                                # larger n first: n = 0 needs n = 1
            out = {}
            if n >= 1:
                base = (n - 1, ms)
                _acc(out, T2.get(base, {}), 2 * e * m0sq)
                if N > 1:
                    _acc(out, self.removed_0(props, r - 2).get(base, {}), K(1))
                for k in labels:
                    _acc(out, T1.get((n - 1, tuple(sorted(ms + (k,)))), {}), f[k])
                out = {key: v / (2 * (d + r - N - 1)) for key, v in out.items()}
            else:
                k, rest = ms[0], ms[1:]
                R0 = self.removed_0(props, r - 1)
                for l in labels:
                    z = Zi[k - 1, l - 1]
                    if z == 0:
                        continue
                    S = dict(self.removed_k(props, l, r - 1).get((0, rest), {}))
                    _acc(S, R0.get((0, rest), {}), K(-1))
                    _acc(S, T1.get((0, rest), {}), -f[l])
                    cnt = rest.count(l)
                    if cnt:
                        rr = list(rest)
                        rr.remove(l)
                        _acc(S, res.get((1, tuple(rr)), {}), K(-2 * cnt))
                    _acc(out, S, z)
            if out:
                res[(n, ms)] = out
        return res

    def tensor(self, props, r):
        key = ('T', tuple((_key(q), m) for q, m in props), r)
        if key in self.cache:
            return self.cache[key]
        if r < 0:
            return {}
        if r == 0:
            return {(0, ()): self.master(props)}
        N = len(props)
        if N == 5:
            res = self.tensor_pent(props, r)
            self.cache[key] = res
            return res
        if self.use_pv and (N == 1 or (all(props[i][0] for i in range(1, N)) and self.zinv(props) is not None)):
            res = self.tensor_pv(props, r)
            self.cache[key] = res
            return res
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
            F = sr(0)
            for coef, mon in zip(poly.coefficients(), poly.monomials()):
                n = list(mon.exponents()[0])
                F += coef * _simplex_integral(groups, n, Msq, lambda i, j: -e * _SR(sij(i, j)), N, k)
            pref = (-1) ** len(ms) / prod([(4 - 2 * eps) + 2 * i for i in range(k)], sr(1))
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
    """An element of the reduction field (or a number) as a symbolic expression.  Ring elements convert
    by their variable names directly, which is exact and much faster than going through a string."""
    if isinstance(x, str):
        return sr(x)
    try:
        return SR(x)
    except (TypeError, ValueError):
        return sr(str(x))


_MF_CACHE = {}


def _master_formula(key, euclidean):
    """(full value with its poles, (1/eps^2 coefficient, 1/eps coefficient)) in the Package-X
    normalisation.  IR-divergent C0 and D0 come from feynsage.ir.  Cached: one master appears in
    many tensor structures."""
    ck = (key, euclidean)
    if key[0] != 'X' and ck in _MF_CACHE:
        return _MF_CACHE[ck]
    out = _master_formula_raw(key, euclidean)
    if key[0] != 'X':
        _MF_CACHE[ck] = out
    return out


def _master_formula_raw(key, euclidean):
    from .ir import ir_coefficients, with_mu
    kind = key[0]
    sign = -1 if euclidean else 1                    # Minkowski invariant = sign * Euclidean one
    zero = sr(0)
    if kind == 'A':
        m = _SR(key[1])
        full, pole = A0(m), (zero, m ** 2)
        n = 1
    elif kind == 'B':
        s, m1, m2 = sign * _SR(key[1]), _SR(key[2]), _SR(key[3])
        full, pole = B0(s, m1, m2), (zero, sr(1))
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
    total = sr(0)
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
    total = sr(0)
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

    def __init__(self, parts, master_parts, euclidean, description, raw=None, K=None):
        self.parts = parts                      # [(structure, expression in eps)]
        self.master_parts = master_parts        # [(structure, sum c(d) * master)]
        self.euclidean = euclidean
        self.description = description
        self.raw = raw or []                    # [(structure, {master key: coefficient})]
        self.K = K

    def discontinuity(self, s, t=None):
        """The discontinuity across the normal threshold cut in the channel s, I(s + i0) - I(s - i0), for
        every tensor structure (Package-X's Part -> Discontinuity[s]).  s must be one of the invariants.
        With t: the Mandelstam double spectral function (Discontinuity[s, t]), from the boxes."""
        if self.euclidean:
            raise ValueError("discontinuities are defined for Minkowski kinematics")
        return {st: _disc_dict(dct, self.K, s, t) for st, dct in self.raw}

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

    def _repr_latex_(self):
        from .tex import structure, laurent
        rows = [r'%s &:\quad %s' % (structure(st), laurent(e)) for st, e in self.parts]
        return r'$$\begin{aligned} %s \end{aligned}$$' % (r' \\[6pt] '.join(rows) if rows else '0')

    def _latex_(self):
        # used by Sage's latex() and by %display latex in notebooks
        return self._repr_latex_()[2:-2]

    def __repr__(self):
        if not self.parts:
            return "0"
        return "\n".join("%s :  %s" % (st, c) for st, c in self.parts)


_INFO = r"""loop(): a one-loop integral reduced to scalar functions.
  Normalisation (Package-X / LoopTools): mu^(2 eps) e^(eps gamma_E) Int d^d l/(i pi^(d/2)),
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
    masses = [sr(m) for _, m in props]
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
    expr = sr(eval(s, {}, names)).expand()
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
                coef *= sr(str(K.dot({a: 1}, {b: 1}))) ** k
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
    parts, mparts, raw = [], [], []
    for exps, dct in sorted(split.items(), key=lambda t: label(t[0])):
        val = _finalize(dct, K, euclidean)
        if not val.is_trivial_zero():
            parts.append((label(exps), val))
            mparts.append((label(exps), _masters_view(dct, euclidean)))
            raw.append((label(exps), dct))
    return LoopResult(parts, mparts, euclidean, _INFO, raw=raw, K=K)


_TOKEN_G = re.compile(r'\bg\(\s*(\w+)\s*,\s*(\w+)\s*\)')
_TOKEN_IDX = re.compile(r'\b([A-Za-z]\w*)\^([A-Za-z]\w*)\b')
_TOKEN_SQ = re.compile(r'(?<![.\w])([A-Za-z]\w*)\^2\b')     # not after a dot: p.pp^2 is (p.pp)^2
_TOKEN_DOT = re.compile(r'\b([A-Za-z]\w*)\.([A-Za-z]\w*)\b')


# ---------------------------------------------------------------------------- Package-X style names
def _masses_syms(*xs):
    return {str(v) for x in xs for v in sr(x).variables()}


def PVB(r, n, s, m1, m2, explain=False, series=None, disc=None):
    """PV coefficient B_{00..(2r) 11..(n)}(s; m1, m2) as Package-X's PVB[r, n, s, m1, m2]:
    PVB(0,0) = B0, PVB(0,1) = B1, PVB(1,0) = B00, PVB(0,2) = B11.  Exact, expanded to eps^0."""
    if explain:
        print(PVB.__doc__)
    if series is not None:
        z = _series_or_none(2, r, (1,) * n, {(0, 1): s}, (m1, m2), series)
        return z if z is not None else loop_series(PVB(r, n, s, m1, m2), series)
    K = Kin(['p'], {'p^2': s}, symbols=_masses_syms(s, m1, m2))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m1)), ({'p': 1}, K.to_K(m2))), 2 * r + n)
    if disc is not None:
        return _disc_dict(T.get((r, (1,) * n), {}), K, *(disc if isinstance(disc, (tuple, list)) else (disc,)))
    return _finalize(T[(r, (1,) * n)], K, False)


def PVC(r, n1, n2, s1, s12, s2, m0, m1, m2, explain=False, series=None, disc=None):
    """PV coefficient C_{00..1..2..}(s1, s12, s2; m0, m1, m2) as Package-X's PVC[r, n1, n2, ...].
    Propagators (l, m0), (l + p1, m1), (l + p2, m2); s1 = p1^2, s2 = p2^2, s12 = (p1 - p2)^2."""
    if explain:
        print(PVC.__doc__)
    if series is not None:
        z = _series_or_none(3, r, (1,) * n1 + (2,) * n2, {(0, 1): s1, (1, 2): s12, (0, 2): s2}, (m0, m1, m2), series)
        return z if z is not None else loop_series(PVC(r, n1, n2, s1, s12, s2, m0, m1, m2), series)
    K = Kin(['p1', 'p2'], {'p1^2': s1, 'p2^2': s2, 'p1.p2': (sr(s1) + sr(s2) - sr(s12)) / 2},
            symbols=_masses_syms(s1, s12, s2, m0, m1, m2))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m0)), ({'p1': 1}, K.to_K(m1)), ({'p2': 1}, K.to_K(m2))), 2 * r + n1 + n2)
    if disc is not None:
        return _disc_dict(T.get((r, (1,) * n1 + (2,) * n2), {}), K, *(disc if isinstance(disc, (tuple, list)) else (disc,)))
    return _finalize(T[(r, (1,) * n1 + (2,) * n2)], K, False)


def PVD(r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, explain=False, series=None, disc=None):
    """PV coefficient D_{00..1..2..3..} as Package-X's PVD[r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3].
    Propagators (l, m0), (l + p1, m1), (l + p2, m2), (l + p3, m3); s1 = p1^2, s2 = (p2 - p1)^2,
    s3 = (p3 - p2)^2, s4 = p3^2, s12 = p2^2, s23 = (p3 - p1)^2.  Exact, expanded to eps^0
    (IR poles explicit when the box is soft or collinear divergent)."""
    if explain:
        print(PVD.__doc__)
    if series is not None:
        z = _series_or_none(4, r, (1,) * n1 + (2,) * n2 + (3,) * n3,
                            {(0, 1): s1, (1, 2): s2, (2, 3): s3, (0, 3): s4, (0, 2): s12, (1, 3): s23},
                            (m0, m1, m2, m3), series)
        return z if z is not None else loop_series(PVD(r, n1, n2, n3, s1, s2, s3, s4, s12, s23, m0, m1, m2, m3), series)
    s1, s2, s3, s4, s12, s23 = [sr(x) for x in (s1, s2, s3, s4, s12, s23)]
    K = Kin(['p1', 'p2', 'p3'], {'p1^2': s1, 'p2^2': s12, 'p3^2': s4, 'p1.p2': (s1 + s12 - s2) / 2,
                                 'p2.p3': (s12 + s4 - s3) / 2, 'p1.p3': (s1 + s4 - s23) / 2},
            symbols=_masses_syms(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m0)), ({'p1': 1}, K.to_K(m1)), ({'p2': 1}, K.to_K(m2)), ({'p3': 1}, K.to_K(m3))),
                   2 * r + n1 + n2 + n3)
    if disc is not None:
        return _disc_dict(T.get((r, (1,) * n1 + (2,) * n2 + (3,) * n3), {}), K, *(disc if isinstance(disc, (tuple, list)) else (disc,)))
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
        if all(sr(a).is_numeric() for a in args):
            try:
                return c0_closed(*_exact_args(args))
            except ZeroDivisionError:
                pass
        return _C0f(*args)

    def rep_d(*args):
        if all(sr(a).is_numeric() for a in args):
            try:
                return d0_closed(*_exact_args(args))
            except ZeroDivisionError:
                pass
        return _D0f(*args)
    return sr(expr).substitute_function(_C0f, rep_c).substitute_function(_D0f, rep_d)


# ---------------------------------------------------------------------------- series (LoopRefineSeries)
def _zero_momentum_series(N, k, ms, sij, Msq, v, v0, order):
    r"""
    Taylor series in v about v0 of the PV coefficient T_{00..(2k) ms} of an N-point function whose
    external invariants s_ij(v) all vanish at v0 (Package-X normalisation, Minkowski):
        T = pref * Int_simplex prod_(l in ms) x_l Delta^s,  Delta = sum_i x_i m_i^2 - sum_(i<j) x_i x_j s_ij,
    s = d/2 + k - N.  With M = sum_i x_i m_i^2 and S = sum x_i x_j s_ij = O(v - v0):
        Delta^s = sum_j binom(s, j) M^(s - j) (-S)^j,
    and every term is a simplex integral of a monomial times a power of M, an exact divided
    difference in the masses (_simplex_integral with one group of lines).
    """
    h = SR.var('fs_h')
    xs = [SR.var('fs_x%d' % i) for i in range(N)]
    ser = {key: sr(val).subs({v: v0 + h}).series(h, order + 1).truncate() for key, val in sij.items()}
    S = sum(xs[i] * xs[j] * ser[(i, j)] for (i, j) in ser)
    s = 2 - eps + k - N
    base = [0] * N
    for l in ms:
        base[l] += 1
    total = SR(0)
    for j in range(order + 1):
        poly = ((-S) ** j).expand().series(h, order + 1).truncate().expand() if j else SR(1)
        bj = prod([s - i for i in range(j)], SR(1)) / factorial(j)
        poly = sr(poly).expand()
        terms = poly.operands() if 'add' in getattr(poly.operator(), '__name__', '') else [poly]
        for t in terms:
            ex = [int(t.degree(x)) for x in xs] + [int(t.degree(h))]
            c = (t / prod([x ** e for x, e in zip(xs + [h], ex)], SR(1))).simplify_rational()
            nvec = [base[i] + ex[i] for i in range(N)]
            total += c * bj * h ** ex[N] * _simplex_integral([list(range(N))], nvec, Msq, None, N, k - j)
    dh = 2 - eps
    pref = (-1) ** len(ms) / prod([(4 - 2 * eps) + 2 * i for i in range(k)], SR(1))
    pref *= gamma(k + dh) * gamma(N - k - dh) / gamma(dh) * (-1) ** (N + k)
    pref *= exp(eps * euler_gamma) * mu ** (2 * eps)
    out = SR(0)
    full = (pref * total).expand()
    for t in range(order + 1):
        ct = full.coefficient(h, t) if t else full.subs({h: 0})
        ce = ct.series(eps, 1).truncate().expand()
        ce = sum(ce.coefficient(eps, a) * eps ** a for a in (-2, -1, 0))
        out += _merge_mu(ce).simplify_rational() * (v - v0) ** t
    return out


def _series_or_none(N, r_k, ms, sij, masses, series):
    """The zero-momentum route when every invariant vanishes at the expansion point, else None."""
    v, v0, n = series
    if all(sr(x).subs({v: v0}).is_zero() for x in sij.values()):
        return _zero_momentum_series(N, r_k, ms, sij, [sr(m) ** 2 for m in masses], v, v0, n)
    return None


def loop_series(expr, *specs):
    r"""
    Taylor series of a one-loop result about a point, as Package-X's LoopRefineSeries:
        loop_series(expr, (s, s0, n))                 to order (s - s0)^n
        loop_series(expr, (s, s0, ns), (t, t0, nt))   first in s, then in t
    expr may contain A0, B0 (logs, LogM, DiscB), C0 and D0; their derivatives are exact
    (integration by parts).  At a point where expr is analytic but its closed form is singular term
    by term (for example B0 at s = 0), use the series= option of PVB, PVC, PVD instead, which expands
    the Feynman-parameter integral about zero external momenta.
    """
    out = sr(expr)
    for v, v0, n in specs:
        h = SR.var('fs_h')
        terms, cur = [], out
        for t in range(n + 1):
            try:
                val = cur.subs({v: v0})
            except (ValueError, ZeroDivisionError, RuntimeError) as err:
                raise ValueError("no Taylor series from the closed form at %s = %s (%s); try the series= option "
                                 "of PVB/PVC/PVD for an expansion about zero momenta" % (v, v0, err))
            terms.append(_merge_mu(val) / factorial(t) * (v - v0) ** t)
            if t < n:
                cur = loop_diff(cur, v)
        out = sum(terms)
    return out


def PVA(r, m, explain=False, series=None):
    """PV coefficient A_{00..(2r)}(m) as Package-X's PVA[r, m]: PVA(0, m) = A0(m).  Exact, expanded to eps^0."""
    if explain:
        print(PVA.__doc__)
    K = Kin([], {}, symbols=_masses_syms(m))
    eng = _Engine(K)
    T = eng.tensor((({}, K.to_K(m)),), 2 * r)
    return _finalize(T.get((r, ()), {}), K, False)


# ---------------------------------------------------------------------------- kinematic helpers
def kallen(a, b, c):
    """Kallen function lambda(a, b, c) = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc (Package-X's Kallen-lambda)."""
    a, b, c = sr(a), sr(b), sr(c)
    return a ** 2 + b ** 2 + c ** 2 - 2 * a * b - 2 * a * c - 2 * b * c


def kibble(s1, s2, s3, s4, s12, s23):
    """Kibble polynomial phi(s1, s2, s3, s4; s12, s23) of a four-point function (Package-X's Kibble-phi);
    its zeros are the boundary of the physical region."""
    s1, s2, s3, s4, s12, s23 = [sr(x) for x in (s1, s2, s3, s4, s12, s23)]
    return (-s1 * s12 * s2 + s1 * s12 * s23 - s12 ** 2 * s23 + s12 * s2 * s23 - s12 * s23 ** 2 - s1 ** 2 * s3
            + s1 * s12 * s3 + s1 * s2 * s3 + s1 * s23 * s3 + s12 * s23 * s3 - s2 * s23 * s3 - s1 * s3 ** 2
            + s1 * s2 * s4 + s12 * s2 * s4 - s2 ** 2 * s4 - s1 * s23 * s4 + s12 * s23 * s4 + s2 * s23 * s4
            + s1 * s3 * s4 - s12 * s3 * s4 + s2 * s3 * s4 - s2 * s4 ** 2)


def mandelstam(momenta=("p1", "p2", "p3", "p4"), masses=("m1", "m2", "m3", "m4"), stu=("s", "t", "u")):
    r"""
    Scalar products of p1 + p2 -> p3 + p4 in Mandelstam variables, as a kin dictionary
    (Package-X's MandelstamRelations): s = (p1 + p2)^2, t = (p1 - p3)^2, u = (p1 - p4)^2.

        loop(..., kin=mandelstam())
    """
    p1, p2, p3, p4 = momenta
    m1, m2, m3, m4 = [sr(m) for m in masses]
    s, t, u = [sr(x) for x in stu]
    return {"%s^2" % p1: m1 ** 2, "%s^2" % p2: m2 ** 2, "%s^2" % p3: m3 ** 2, "%s^2" % p4: m4 ** 2,
            "%s.%s" % (p1, p2): (s - m1 ** 2 - m2 ** 2) / 2, "%s.%s" % (p3, p4): (s - m3 ** 2 - m4 ** 2) / 2,
            "%s.%s" % (p1, p3): (m1 ** 2 + m3 ** 2 - t) / 2, "%s.%s" % (p2, p4): (m2 ** 2 + m4 ** 2 - t) / 2,
            "%s.%s" % (p1, p4): (m1 ** 2 + m4 ** 2 - u) / 2, "%s.%s" % (p2, p3): (m2 ** 2 + m3 ** 2 - u) / 2}


def disc_expand(expr):
    """Replace DiscB(s, m1, m2) by its definition sqrt(lambda)/s log((m1^2 + m2^2 - s + sqrt(lambda))/(2 m1 m2))
    (Package-X's DiscExpand).  The +i0 of s is then implicit: use it where s is below threshold."""
    def rep(s, m1, m2):
        r = sqrt(kallen(s, m1 ** 2, m2 ** 2))
        return r / s * log((m1 ** 2 + m2 ** 2 - s + r) / (2 * m1 * m2))
    return sr(expr).substitute_function(DiscB, rep)


def _sgn(a):
    a = sr(a)
    if not a.is_numeric() or a == 0:
        raise ValueError("the side of the branch cut must be a nonzero number")
    return 1 if a > 0 else -1


def _ln_evalf(self, x, a, parent=None, algorithm=None):
    mp = _mp()
    prec = _prec(parent)
    with mp.workprec(prec + 20):
        v = _mpf(x, prec + 20)
        val = mp.log(-v) + _sgn(a) * 1j * mp.pi if v < 0 else mp.log(v)
        return _to_parent(val, parent)


def _dilog_evalf(self, x, a, parent=None, algorithm=None):
    mp = _mp()
    prec = _prec(parent)
    with mp.workprec(prec + 20):
        v = _mpf(x, prec + 20)
        val = mp.polylog(2, v)
        if v > 1:
            val = mp.re(val) + _sgn(a) * 1j * mp.pi * mp.log(v)
        return _to_parent(val, parent)


Ln = function('Ln', nargs=2, evalf_func=_ln_evalf,
              print_latex_func=lambda self, x, a: r"\ln\left(%s + i %s 0\right)" % (x._latex_(), a._latex_()))
Ln.__doc__ = "Ln(x, a): log(x + i a 0) for real x, the side of the cut set by the sign of a (Package-X's Ln)."
DiLog = function('DiLog', nargs=2, evalf_func=_dilog_evalf,
                 print_latex_func=lambda self, x, a: r"\mathrm{Li}_2\left(%s + i %s 0\right)" % (x._latex_(), a._latex_()))
DiLog.__doc__ = "DiLog(x, a): Li2(x + i a 0) for real x, the side of the cut set by the sign of a (Package-X's DiLog)."


def continued_dilog(x1, a1, x2, a2):
    r"""
    Beenakker-Denner continued dilogarithm Li2(x1 + i a1 0, x2 + i a2 0) = Li2(1 - x1 x2) + eta(x1, x2) log(1 - x1 x2)
    for real x1, x2 (Package-X's ContinuedDiLog[{x1, a1}, {x2, a2}]), as a number.
    eta = 2 pi i [theta(-Im x1) theta(-Im x2) theta(Im x1 x2) - theta(Im x1) theta(Im x2) theta(-Im x1 x2)].
    """
    from sage.all import CC
    mp = _mp()
    with mp.workprec(120):
        x1, x2 = _mpf(x1, 120), _mpf(x2, 120)
        s1, s2 = _sgn(a1), _sgn(a2)
        im12 = s1 * x2 + s2 * x1                      # sign of Im(x1 x2) to first order in the i0's
        s12 = 1 if im12 > 0 else -1
        eta = 0
        if s1 < 0 and s2 < 0 and s12 > 0:
            eta = 2j * mp.pi
        elif s1 > 0 and s2 > 0 and s12 < 0:
            eta = -2j * mp.pi
        w = 1 - x1 * x2                                # its i0 is -s12
        li = mp.polylog(2, w)
        lg = mp.log(-w) - s12 * 1j * mp.pi if w < 0 else mp.log(w)
        if w > 1:
            li = mp.re(li) - s12 * 1j * mp.pi * mp.log(w)
        return CC(mp.re(li + eta * lg)) + CC(0, 1) * CC(mp.im(li + eta * lg))


# ---------------------------------------------------------------------------- discontinuities (cuts)
def _master_lines(key):
    """(masses, invariant s_ij for every pair i < j) of a master key, in Minkowski space."""
    k = key[0]
    if k == 'B':
        return [_SR(key[2]), _SR(key[3])], {(0, 1): _SR(key[1])}
    if k == 'C':
        s1, s12, s2 = [_SR(x) for x in key[1:4]]
        return [_SR(x) for x in key[4:7]], {(0, 1): s1, (1, 2): s12, (0, 2): s2}
    if k == 'D':
        s01, s12, s23, s03, s02, s13 = [_SR(x) for x in key[1:7]]
        return [_SR(x) for x in key[7:11]], {(0, 1): s01, (1, 2): s12, (2, 3): s23, (0, 3): s03, (0, 2): s02, (1, 3): s13}
    return None, None


def _cut(masses, sij, i, j):
    r"""
    Discontinuity across the normal threshold of the channel s = s_ij (lines i and j cut), by
    Cutkosky's rule in the Package-X normalisation (Int d^4 l/(i pi^2)):
        Disc I = 4 i Int d^4 l delta+(D_i) delta+(D_j) prod_(k != i, j) 1/D_k.
    In the rest frame of P = q_j - q_i the remaining denominators are a_k + b_k.n on the sphere:
        bubble   2 pi i sqrt(lam_m)/s
        triangle (i pi sqrt(lam_m)/(s |b|)) log((a + |b|)/(a - |b|))
        box      (i pi sqrt(lam_m)/(s sqrt(X))) log((Y + sqrt(X))/(Y - sqrt(X))),
                 Y = a_k a_l - b_k.b_l,  X = Y^2 - (a_k^2 - b_k^2)(a_l^2 - b_l^2),
    each times theta(s - (m_i + m_j)^2).  Everything except the square roots is rational in the invariants.
    """
    from sage.all import heaviside
    N = len(masses)
    S = lambda a, b: sij[(min(a, b), max(a, b))] if a != b else sr(0)
    s = S(i, j)
    mi2, mj2 = masses[i] ** 2, masses[j] ** 2
    lam_m = kallen(s, mi2, mj2)
    theta = heaviside(s - (masses[i] + masses[j]) ** 2)
    others = [k for k in range(N) if k not in (i, j)]
    if not others:
        return 2 * pi * I * sqrt(lam_m) / s * theta
    QP = {k: (s + S(i, k) - S(j, k)) / 2 for k in others}          # Q_k.P with Q_k = q_k - q_i
    QQ_ = lambda k, l: (S(i, k) + S(i, l) - S(k, l)) / 2           # Q_k.Q_l
    l2 = lam_m / (4 * s)                                            # |l|^2
    a = {k: mi2 + S(i, k) - masses[k] ** 2 + (mj2 - mi2 - s) * QP[k] / s for k in others}
    bb = lambda k, l: 4 * l2 * (QP[k] * QP[l] / s - QQ_(k, l))      # b_k . b_l
    if len(others) == 1:
        k = others[0]
        b = sqrt(bb(k, k))
        return I * pi * sqrt(lam_m) / (s * b) * log((a[k] + b) / (a[k] - b)) * theta
    k, l = others
    Y = a[k] * a[l] - bb(k, l)
    X = Y ** 2 - (a[k] ** 2 - bb(k, k)) * (a[l] ** 2 - bb(l, l))
    return I * pi * sqrt(lam_m) / (s * sqrt(X)) * log((Y + sqrt(X)) / (Y - sqrt(X))) * theta


def _double_spectral(masses, sij):
    r"""Mandelstam double spectral function of a box with s and t on the two pairs of opposite lines
    (Package-X's Discontinuity[s, t]): 2 pi i theta(det Y)/sqrt(det Y), Y_ij = m_i^2 + m_j^2 - s_ij the
    modified Cayley matrix."""
    from sage.all import heaviside
    Y = matrix(SR, 4, 4, lambda i, j: masses[i] ** 2 + masses[j] ** 2 - (sij[(min(i, j), max(i, j))] if i != j else 0))
    dY = Y.det().expand()
    return 2 * pi * I * heaviside(dY) / sqrt(dY)


def _disc_dict(dct, K, s, t=None):
    """Discontinuity in the channel s of sum_M c_M(d) M: the coefficients at d = 4 times the cut masters.
    With t as well: the double spectral function (boxes only)."""
    K = K.K if isinstance(K, Kin) else K
    s = sr(s)
    total = sr(0)
    if t is not None:
        t = sr(t)
        for key, c in dct.items():
            masses, sij = _master_lines(key)
            if masses is None or len(masses) != 4:
                continue
            ok = [bool((sij[a] - s).expand() == 0) and bool((sij[b] - t).expand() == 0) or
                  bool((sij[a] - t).expand() == 0) and bool((sij[b] - s).expand() == 0)
                  for a, b in (((0, 2), (1, 3)), ((0, 1), (2, 3)), ((0, 3), (1, 2)))]
            if any(ok):
                total += _SR(_at4(c, K)[0]) * _double_spectral(masses, sij)
        return total
    for key, c in dct.items():
        masses, sij = _master_lines(key)
        if masses is None:
            if key[0] == 'X':
                raise NotImplementedError("discontinuity of a vanishing-Gram master: move the kinematics slightly")
            continue                                                   # A0: no cut
        for (i, j), v in sij.items():
            if bool((v - s).expand() == 0):
                total += _SR(_at4(c, K)[0]) * _cut(masses, sij, i, j)
    return total


def loop_many(jobs, nproc=None):
    r"""
    Many one-loop integrals at once, on every core: jobs is a list of (numerator, propagators,
    options), each as for loop(numerator, *propagators, **options).  Returns the LoopResults in order.

        loop_many([("1", (["l", "m"], ["l + p", "m"]), {"kin": {"p^2": "s"}}),
                   ("l^mu", (["l", "m"], ["l + p", "M"]), {"kin": {"p^2": "s"}})])
    """
    from .parallel import pmap
    def one(job):
        num, props, opts = (tuple(job) + ({},))[:3] if len(job) == 2 else job
        return loop(num, *props, **(opts or {}))
    return pmap(one, list(jobs), nproc=nproc)
