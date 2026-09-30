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

and C0, D0 are symbolic functions that evaluate numerically (see c0_numeric).

Quick start
-----------
    from feynsage.pv import *
    loop("1", ["l", "m"], ["l + p", "m"], kin={"p^2": "s"})             # B0(s, m, m)
    loop("l^mu l^nu", ["l", "m1"], ["l + p", "m2"], kin={"p^2": "s"})   # g^{mu nu} B00 + p^mu p^nu B11
    B0(s, m, m).explicit()   # closed form, like Package-X's LoopRefine
"""
import itertools
import re
from sage.all import (SR, var, function, log, sqrt, pi, I, gamma, matrix, factorial, PolynomialRing, lcm,
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


def _c0_evalf(self, *args, parent=None, algorithm=None):
    v = c0_numeric(*[float(a) for a in args])
    return parent(v) if parent is not None else v


def _d0_evalf(self, *args, parent=None, algorithm=None):
    v = d0_numeric(*[float(a) for a in args])
    return parent(v) if parent is not None else v


C0 = function('C0', nargs=6, evalf_func=_c0_evalf)
D0 = function('D0', nargs=10, evalf_func=_d0_evalf)


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
    """The coefficient of 1/eps (the UV pole, if the integral is IR finite)."""
    return SR(expr).expand().coefficient(eps, -1)


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
            raise ZeroDivisionError("vanishing Gram determinant (exceptional kinematics, e.g. collinear momenta); "
                                    "move the kinematics slightly or use symbols")
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
    """(full value with its pole, the 1/eps coefficient) in the Package-X normalisation."""
    kind = key[0]
    sign = -1 if euclidean else 1                    # Minkowski invariant = sign * Euclidean one
    if kind == 'A':
        m = _SR(key[1])
        full, pole = A0(m), m ** 2
        n = 1
    elif kind == 'B':
        s, m1, m2 = sign * _SR(key[1]), _SR(key[2]), _SR(key[3])
        full, pole = B0(s, m1, m2), SR(1)
        if _is_zero(s) and _is_zero(m1) and _is_zero(m2):
            pole = SR(0)
        n = 2
    elif kind == 'C':
        args = [sign * _SR(x) for x in key[1:4]] + [_SR(x) for x in key[4:]]
        full, pole = C0(*args), SR(0)
        n = 3
    else:
        args = [sign * _SR(x) for x in key[1:7]] + [_SR(x) for x in key[7:]]
        full, pole = D0(*args), SR(0)
        n = 4
    if euclidean and n % 2:
        full, pole = -full, -pole
    return full, pole


def _at4(c, K):
    """c(d = 4) and dc/dd(d = 4) for c in K."""
    R = K.ring()
    dgen = R('d')
    num, den = c.numerator(), c.denominator()
    n4, d4 = num.subs({dgen: 4}), den.subs({dgen: 4})
    if d4 == 0:
        raise ZeroDivisionError("a coefficient has a pole at d = 4; the eps expansion needs more terms")
    dn, dd = num.derivative(dgen).subs({dgen: 4}), den.derivative(dgen).subs({dgen: 4})
    return K(n4) / K(d4), (K(dn) * K(d4) - K(n4) * K(dd)) / K(d4) ** 2


def _finalize(dct, K, euclidean):
    """sum_M c_M(d) M  ->  expression in eps (up to eps^0), with the rational terms from d = 4 - 2 eps."""
    K = K.K if isinstance(K, Kin) else K
    total = SR(0)
    for key, c in dct.items():
        full, pole = _master_formula(key, euclidean)
        c4, dc4 = _at4(c, K)
        total += _SR(c4) * full
        if not pole.is_trivial_zero() and dc4 != 0:
            total += -2 * _SR(dc4) * pole
    return total


def _masters_view(dct, euclidean):
    """sum_M c_M(d) M with the masters as plain symbols (like Package-X's LoopIntegrate)."""
    total = SR(0)
    for key, c in dct.items():
        k = key[0]
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
    1/eps       UV pole (IR poles are not separated; C0 and D0 are assumed IR finite),
    A0, B0      written out: logs, LogM(z) = log(z - i0), DiscB(s, m1, m2) (see B0?),
    C0(...), D0(...)  exact symbols; .n() evaluates them numerically (two-grid check),
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
        for mon, cf in num.dict().items():
            exps = tuple(mon[i] for i in sidx)
            rest = list(mon)
            for i in sidx:
                rest[i] = 0
            groups.setdefault(exps, {})[tuple(rest)] = cf
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
_TOKEN_SQ = re.compile(r'\b([A-Za-z]\w*)\^2\b')
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
