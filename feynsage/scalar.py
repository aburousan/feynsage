r"""
The scalar three-point function C0 in closed form: dilogarithms and logarithms of exact
algebraic numbers, with every branch fixed by the +i0 of the propagators.

Conventions as in Package-X / LoopTools:  C0(s1, s12, s2; m0, m1, m2) with propagators
(l, m0), (l + p1, m1), (l + p2, m2), s1 = p1^2, s2 = p2^2, s12 = (p1 - p2)^2.

Method.  Write p_jk^2 for the invariant opposite to vertex i and
Q_i(y) = p_jk^2 y^2 - (p_jk^2 - m_j^2 + m_k^2) y + m_k^2 - i0 (the Feynman denominator on one
edge of the parameter triangle).  Then, with alpha = sqrt(lambda(s1, s12, s2)),

    C0 = -(1/alpha) sum_i Int_0^1 dy [log Q_i(y) - log Q_i(y0_i)] / (y - y0_i),

y0_i from A. Denner, Fortsch. Phys. 41 (1993) 307, eq. (4.27).  For lambda > 0 the integrals
are Denner's twelve dilogarithms (eq. 4.26, whose eta terms vanish then).  For lambda < 0,
y0_i is complex and each integral is done with explicit bookkeeping of the branch crossings
(at most one Li2 cut crossing per term, and constant 2 pi i shifts).  A zero invariant
p_jk^2 = 0 either removes its term (y0_i -> infinity) or makes Q_i linear with y0_i the
L'Hopital limit.  lambda = 0 with all invariants zero has its own elementary formula.

Everything is computed twice: as an mpmath number with a tiny explicit -i eps, and as an
exact Sage expression with eps -> 0.  Where an argument of log or Li2 lies on its cut, the
exact expression uses the continuation that the number selects:
    log(-x +- i0) = log x +- i pi,   Li2(x +- i0) = pi^2/3 - log^2(x)/2 - Li2(1/x) +- i pi log x  (x > 1).
c0_closed() checks at the end that the two agree.
"""
import mpmath as mp
from sage.all import SR, QQ, RR, sqrt, log, polylog, pi, I, Integer

DPS = 80
EPS = mp.mpf(10) ** -45


class V:
    """A number twice: n (mpmath, with -i eps) and e (exact Sage expression, eps -> 0).
    e = None means numbers only (the exact side is not built)."""
    __slots__ = ('n', 'e')

    def __init__(self, n, e):
        self.n, self.e = mp.mpc(n), (SR(e) if e is not None else None)

    @staticmethod
    def of(x):
        if isinstance(x, V):
            return x
        q = QQ(x) if not isinstance(x, (float,)) else RR(x).simplest_rational()
        return V(mp.mpf(int(q.numerator())) / int(q.denominator()), q)

    @staticmethod
    def ofx(x, exact=True):
        """Any real number (rational or algebraic) as a V pair; exact=False keeps the number only."""
        if isinstance(x, V):
            return x
        x = SR(x)
        try:
            q = QQ(x)
            n = mp.mpf(int(q.numerator())) / int(q.denominator())
            return V(n, q if exact else None)
        except (TypeError, ValueError):
            return V(mp.mpc(complex(x.n(prec=int(DPS * 3.4) + 20))), x if exact else None)

    def __add__(s, o): o = V.of(o); return V(s.n + o.n, _ex(s.e, o.e, lambda a, b: a + b))
    __radd__ = __add__
    def __sub__(s, o): o = V.of(o); return V(s.n - o.n, _ex(s.e, o.e, lambda a, b: a - b))
    def __rsub__(s, o): o = V.of(o); return V(o.n - s.n, _ex(o.e, s.e, lambda a, b: a - b))
    def __mul__(s, o): o = V.of(o); return V(s.n * o.n, _ex(s.e, o.e, lambda a, b: a * b))
    __rmul__ = __mul__
    def __truediv__(s, o): o = V.of(o); return V(s.n / o.n, _ex(s.e, o.e, lambda a, b: a / b))
    def __rtruediv__(s, o): o = V.of(o); return V(o.n / s.n, _ex(o.e, s.e, lambda a, b: a / b))
    def __neg__(s): return V(-s.n, -s.e if s.e is not None else None)

    def is_zero(s):
        return (s.e is not None and s.e.is_trivial_zero()) or abs(s.n) < mp.mpf(10) ** (-DPS // 2)


def _ex(a, b, f):
    return None if a is None or b is None else f(a, b)


def _tiny(z):
    return abs(mp.im(z)) <= mp.mpf(10) ** (-DPS // 3) * max(1, abs(z))


def _on_cut(z, side):
    """Is the exact value of z real and on the cut?  side = 'neg' (log, sqrt: (-oo, 0)) or 'gt1'
    (Li2: (1, oo)).  Decided from the exact expression (at high precision), not from the
    infinitesimal imaginary part of the number, which can be amplified by the algebra."""
    if not (_tiny(z.n) or abs(mp.im(z.n)) < mp.mpf(10) ** -8 * max(1, abs(z.n))):
        return False
    try:
        ev = complex(z.e.n(prec=int(DPS * 3.4)))
    except (TypeError, ValueError, ArithmeticError):
        return False
    if abs(ev.imag) > 1e-40 * max(1, abs(ev)):
        return False
    return ev.real < 0 if side == 'neg' else ev.real > 1


def _side(z):
    return 1 if mp.im(z.n) >= 0 else -1


def vsqrt(x):
    """sqrt with the same branch in both parts (principal, following the numerical value)."""
    n = mp.sqrt(x.n)
    ex = x.e
    if ex is None:
        return V(n, None)
    if _on_cut(x, 'neg'):
        return V(n, _side(x) * I * sqrt(-ex))
    return V(n, sqrt(ex))


def vlog(z):
    n = mp.log(z.n)
    if z.e is None:
        return V(n, None)
    if _on_cut(z, 'neg'):                      # on the cut: the side is set by the i eps
        return V(n, log(-z.e) + _side(z) * I * pi)
    return V(n, log(z.e))


def _li2(z):
    """Li2(z) for an mpmath complex at the working precision, computed by PARI (compiled, about 15 times
    faster than mpmath's polylog here) with mpmath as the fallback.  The numbers go through decimal
    strings at full precision: same value as mpmath to ~1e-84 at 80 digits and the same side of the
    cut [1, oo) for arguments with a tiny imaginary part."""
    try:
        from sage.all import ComplexField, RealField, pari
        bits = mp.mp.prec + 10
        RF, CF = RealField(bits), ComplexField(bits)
        dig = int(bits / 3.32) + 5
        zz = CF(RF(mp.nstr(mp.re(z), dig, min_fixed=-10**9, max_fixed=10**9)),
                RF(mp.nstr(mp.im(z), dig, min_fixed=-10**9, max_fixed=10**9)))
        w = CF(pari(zz).dilog(precision=bits).sage())
        return mp.mpc(mp.mpf(w.real().str(truncate=False)), mp.mpf(w.imag().str(truncate=False)))
    except Exception:
        return mp.polylog(2, z)


def vli2(z):
    n = _li2(mp.mpc(z.n))
    if z.e is None:
        return V(n, None)
    if _on_cut(z, 'gt1'):                      # on the cut [1, oo)
        x = z.e
        return V(n, pi ** 2 / 3 - log(x) ** 2 / 2 - polylog(2, 1 / x) + _side(z) * I * pi * log(x))
    if z.e.is_trivial_zero():
        return V(0, 0)
    return V(n, polylog(2, z.e))


def two_pi_i_times(k):
    return V(2j * mp.pi * k, 2 * pi * I * k)


def _int_of(z):
    """z = 2 pi i k numerically -> k."""
    k = mp.im(z) / (2 * mp.pi)
    kr = int(mp.nint(k))
    if abs(k - kr) > mp.mpf(10) ** -10 or abs(mp.re(z)) > mp.mpf(10) ** -10:
        raise ArithmeticError("branch bookkeeping failed (not a multiple of 2 pi i)")
    return kr


# ---------------------------------------------------------------------------- building blocks
def _roots(a, b, c):
    """Roots of a y^2 - b y + c - i eps (a may be 0: then one root), as V pairs."""
    eps = V(1j * EPS * max(1, abs(c.n), abs(b.n), abs(a.n)), 0)
    ce = c - eps
    if a.is_zero():
        if b.is_zero():
            return []
        return [ce / b]
    if a.e is None or b.e is None or c.e is None:
        dn = mp.sqrt(b.n * b.n - 4 * a.n * ce.n)
        return [V((b.n + dn) / (2 * a.n), None), V((b.n - dn) / (2 * a.n), None)]
    disc_exact = V(b.n * b.n - 4 * a.n * c.n, b.e ** 2 - 4 * a.e * c.e)
    r = vsqrt(disc_exact)
    exact = [(b.e + r.e) / (2 * a.e), (b.e - r.e) / (2 * a.e)]
    dn = mp.sqrt(b.n * b.n - 4 * a.n * ce.n)
    nums = [(b.n + dn) / (2 * a.n), (b.n - dn) / (2 * a.n)]
    # pair each numerical root with the exact root it approaches
    ev = [complex(SR(x).n(prec=200)) for x in exact]
    if abs(nums[0] - ev[0]) + abs(nums[1] - ev[1]) > abs(nums[0] - ev[1]) + abs(nums[1] - ev[0]):
        exact = exact[::-1]
    return [V(nums[0], exact[0]), V(nums[1], exact[1])]


def _li2_pair_term(y0, x):
    """Denner's pair  Li2((y0 - 1)/(y0 - x)) - Li2(y0/(y0 - x))  (= -R(y0, x))."""
    return vli2((y0 - 1) / (y0 - x)) - vli2(y0 / (y0 - x))


def _R_complex(y0, x):
    r"""R(y0, x) = Int_0^1 dy [log(y - x) - log(y0 - x)]/(y - y0) for complex y0 (not on [0,1]),
    with log(y - x) continuous along real y (Im x != 0)."""
    t = lambda yv: (y0 - yv) / (y0 - x)
    t0, t1 = t(V.of(0)), t(V.of(1))
    # crossing of 1 - t = (y - x)/(y0 - x) over the negative real axis, for y in (0, 1)
    w0 = (V.of(0) - x) / (y0 - x)
    w1 = (V.of(1) - x) / (y0 - x)
    dw = w1 - w0                                    # w(y) = w0 + y dw
    res = V(0, 0)
    yc = None
    if mp.im(dw.n) != 0:
        ycn = -mp.im(w0.n) / mp.im(dw.n)
        if 0 < ycn < 1 and mp.re(w0.n + ycn * dw.n) < 0:
            # y_c is the root of Im(w0 + y dw) = 0: exact form
            if w0.e is None or dw.e is None:
                yc = V(ycn, None)
            else:
                yce = -(w0.e.imag_part() if hasattr(w0.e, 'imag_part') else 0) / dw.e.imag_part()
                yc = V(ycn, yce.simplify_full() if hasattr(yce, 'simplify_full') else yce)
    if yc is None:
        res = vli2(t0) - vli2(t1)
        k = _int_of(mp.log(-x.n) - mp.log(y0.n - x.n) - mp.log(w0.n))      # branch count: numbers only
        if k:
            res = res + two_pi_i_times(k) * (vlog(V.of(1) - y0) - vlog(-y0))
        return res
    # split at y_c: approach the cut from each side
    h = mp.mpf(10) ** -(DPS - 15)                  # far below the result, still resolved at DPS digits
    tl = V(t(V(ycn - h, 0)).n, t(yc).e)
    tr = V(t(V(ycn + h, 0)).n, t(yc).e)
    res = vli2(t0) - vli2(tl) + vli2(tr) - vli2(t1)
    for (ya, yb, wa) in ((V.of(0), yc, w0), (yc, V.of(1), (V(1, 1) - x) / (y0 - x))):
        ym = V((ya.n + yb.n) / 2, _ex(ya.e, yb.e, lambda u, v: (u + v) / 2))
        wm = (ym - x) / (y0 - x)
        k = _int_of(mp.log(ym.n - x.n) - mp.log(y0.n - x.n) - mp.log(wm.n))
        if k:
            res = res + two_pi_i_times(k) * (vlog(yb - y0) - vlog(ya - y0))
    return res


# ---------------------------------------------------------------------------- C0
def _c0_terms(s1, s12, s2, m0, m1, m2, exact=True):
    P = {(1, 0): V.ofx(s1, exact), (2, 0): V.ofx(s2, exact), (2, 1): V.ofx(s12, exact)}
    p2 = lambda a, b: P[(max(a, b), min(a, b))]
    Ms = [V.ofx(m, exact) * V.ofx(m, exact) for m in (m0, m1, m2)]
    lam = p2(1, 0) * p2(1, 0) + p2(2, 1) * p2(2, 1) + p2(2, 0) * p2(2, 0) \
        - 2 * (p2(1, 0) * p2(2, 1) + p2(2, 1) * p2(2, 0) + p2(2, 0) * p2(1, 0))
    if lam.is_zero():
        return None, lam
    if mp.re(lam.n) > 0:
        alpha = vsqrt(lam)
    else:
        alpha = V(-1j * mp.sqrt(-lam.n), -I * sqrt(-lam.e) if lam.e is not None else None)
    total = V(0, 0)
    positive = mp.re(lam.n) > 0
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        pjk, pki, pij = p2(j, k), p2(k, i), p2(i, j)
        mi, mj, mk = Ms[i], Ms[j], Ms[k]
        if pjk.is_zero():
            N0 = (mk - mj) * (pki - pij + alpha)
            if not N0.is_zero():
                continue                                   # y0_i -> infinity: no contribution
            dalpha = -(pki + pij) / alpha
            y0 = ((-pki - pij + 2 * mi - mj - mk) + dalpha * (mk - mj) + alpha) / (2 * alpha)
        else:
            y0 = (pjk * (pjk - pki - pij + 2 * mi - mj - mk) - (pki - pij) * (mj - mk)
                  + alpha * (pjk - mj + mk)) / (2 * alpha * pjk)
        a, b, c = pjk, pjk - mj + mk, mk
        xs = _roots(a, b, c)
        if positive:
            for x in xs:
                total = total + _li2_pair_term(y0, x)
        else:
            S = V(0, 0)
            for x in xs:
                S = S + _R_complex(y0, x)
            lead = a if not a.is_zero() else (-b)
            # branch counts (numbers only); the i eps is the one used in _roots
            ieps = 1j * EPS * max(1, abs(c.n), abs(b.n), abs(a.n))
            K1 = mp.log(c.n - ieps) - mp.log(lead.n) - sum(mp.log(-x.n) for x in xs)
            Qy0 = lead.n
            for x in xs:
                Qy0 = Qy0 * (y0.n - x.n)
            K2 = mp.log(lead.n) + sum(mp.log(y0.n - x.n) for x in xs) - mp.log(Qy0)
            kk = _int_of(K1 + K2)
            if kk:
                S = S + two_pi_i_times(kk) * (vlog(V.of(1) - y0) - vlog(-y0))
            total = total - S
    return total / alpha, lam


def c0_closed(s1, s12, s2, m0, m1, m2, check=True):
    r"""
    C0(s1, s12, s2; m0, m1, m2) in closed form: an exact Sage expression in polylog(2, .)
    and log of algebraic numbers (inputs are made exact rationals).  Valid for all real
    invariants and real masses where C0 is finite, except lambda(s1, s12, s2) = 0 with some
    invariant nonzero (vanishing Gram determinant), which raises.  check=True compares the
    expression with an independent high-precision evaluation.
    """
    with mp.workdps(DPS):
        args = [V.of(x).e for x in (s1, s12, s2, m0, m1, m2)]
        if all(SR(x).is_trivial_zero() for x in args[:3]):
            return _c0_zero_momenta(*args[3:])
        val, lam = _c0_terms(*args)
        if val is None:
            raise ZeroDivisionError("lambda(s1, s12, s2) = 0 (vanishing Gram determinant): no closed form here")
        if check:
            num = complex(val.e.n(prec=200))
            if abs(num - complex(val.n)) > 1e-12 * max(1, abs(complex(val.n))):
                raise ArithmeticError("closed form and its numerical value disagree: %s vs %s" % (num, complex(val.n)))
        return val.e


def _chop(v):
    """Imaginary parts below 1e-20 of |v| are what is left of the infinitesimal i eps (it acts like
    sqrt(eps) exactly on a threshold): clear them."""
    v = mp.mpc(v)
    if abs(mp.im(v)) < mp.mpf(10) ** -20 * max(1, abs(v)):
        v = mp.mpc(mp.re(v), 0)
    return v


def c0_value(s1, s12, s2, m0, m1, m2, full=False):
    """High-precision value of C0 from the closed-form construction (complex)."""
    with mp.workdps(DPS):
        args = [V.of(x).e for x in (s1, s12, s2, m0, m1, m2)]
        if all(SR(x).is_trivial_zero() for x in args[:3]):
            e = _c0_zero_momenta(*args[3:]).n(prec=300)
            v = mp.mpc(mp.mpf(str(e.real())), mp.mpf(str(e.imag())))
            return _chop(v) if full else complex(_chop(v))
        val, lam = _c0_terms(*args, exact=False)
        if val is None:
            val = _c0_gram_zero(args)
            return _chop(val) if full else complex(_chop(val))
        return _chop(val.n) if full else complex(_chop(val.n))


def _c0_gram_zero(args, rel=None):
    r"""
    C0 where lambda(s1, s12, s2) = 0 (for example the vertex at q^2 = 0 with p1^2 = p2^2, the g - 2
    kinematics), where the closed form divides by sqrt(lambda).  C0 is analytic in the invariants
    there (unless one sits exactly on a threshold), so it is extrapolated from the closed form at one
    invariant shifted by +-delta and +-2 delta:  C0(0) = [4 C(delta) - C(2 delta)]/3 with C the mean
    of the two signs, an error of order delta^4.  delta is 10^GRAM_DELTA times the smallest non-zero
    scale (|s_i| or m_i^2), so it stays far below every threshold.
    """
    nums = [mp.mpf(str(SR(x).n(prec=300))) for x in args]
    scales = [abs(SR(x)) for x in args[:3] if not SR(x).is_trivial_zero()] + \
             [SR(x) ** 2 for x in args[3:] if not SR(x).is_trivial_zero()]
    small = min(scales, key=lambda z: z.n(prec=300)) if scales else SR(1)
    s1, s12, s2 = nums[:3]
    grads = [abs(2 * s1 - 2 * s12 - 2 * s2), abs(2 * s12 - 2 * s1 - 2 * s2), abs(2 * s2 - 2 * s1 - 2 * s12)]
    k = max(range(3), key=lambda i: grads[i])
    delta = small * QQ(10) ** (rel if rel is not None else GRAM_DELTA)
    def mean_at(h):
        vals = []
        for sign in (1, -1):
            a = list(args)
            a[k] = SR(a[k]) + sign * h
            v, lam = _c0_terms(*a, exact=False)
            if v is None:
                raise ZeroDivisionError("lambda(s1, s12, s2) = 0 and no shift helps here")
            vals.append(v.n)
        return (vals[0] + vals[1]) / 2
    return (4 * mean_at(delta) - mean_at(2 * delta)) / 3


GRAM_DELTA = -12                 # delta = 10^GRAM_DELTA times the smallest scale


def _c0_zero_momenta(m0, m1, m2):
    r"""C0(0,0,0; m0, m1, m2) = -Int_simplex 1/(x0 m0^2 + x1 m1^2 + x2 m2^2): elementary."""
    a, b, c = SR(m0) ** 2, SR(m1) ** 2, SR(m2) ** 2
    if (a - b).is_trivial_zero() and (b - c).is_trivial_zero():
        return -1 / (2 * a)
    def g(x, y, z):                               # distinct-mass formula, symmetric
        return x * log(x) / ((x - y) * (x - z))
    if (a - b).is_trivial_zero() or (b - c).is_trivial_zero() or (a - c).is_trivial_zero():
        # two equal masses: take the limit with a symbol
        t = SR.var('fs_t')
        if (a - b).is_trivial_zero():
            expr = -(g(a + t, b, c) + g(b, a + t, c) + g(c, a + t, b))
        elif (b - c).is_trivial_zero():
            expr = -(g(a, b + t, c) + g(b + t, a, c) + g(c, a, b + t))
        else:
            expr = -(g(a + t, b, c) + g(b, a + t, c) + g(c, a + t, b))
        return expr.limit(fs_t=0).simplify_full()
    return -(g(a, b, c) + g(b, a, c) + g(c, a, b))


# ---------------------------------------------------------------------------- D0
r"""
D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3), propagators (l, m0), (l + p1, m1), (l + p1 + p2, m2),
(l + p1 + p2 + p3, m3); p_ij^2 between lines i and j: p01 = s1, p12 = s2, p23 = s3, p03 = s4,
p02 = s12, p13 = s23.  The function is symmetric under all 24 relabellings of the lines.

All masses nonzero: A. Denner, Fortsch. Phys. 41 (1993) 307, eq. (4.43) (16 dilogarithms,
after relabelling so that r_02 is real) and its companion formula when all |r_ij| = 1.
One or more masses zero: A. Denner and S. Dittmaier, Nucl. Phys. B844 (2011) 199,
arXiv:1005.2076, eqs. (3.76), (3.78), (3.82), (3.84).
Infinitesimals: the i eps of those papers is ED; the -i0 of the internal masses is EDM << ED.
IR-divergent boxes are not handled here (see feynsage.ir).
"""
ED = mp.mpf(10) ** -50
EDM = mp.mpf(10) ** -64
_ZIM = mp.mpf(0)          # eta needs the sign of every infinitesimal imaginary part, however small


def _th(x):
    return 1 if x > _ZIM else 0


def _eta_k(a, b):
    """eta(a, b) = log(ab) - log a - log b as an integer k (eta = 2 pi i k), from the numbers."""
    an, bn = mp.mpc(a.n if isinstance(a, V) else a), mp.mpc(b.n if isinstance(b, V) else b)
    ab = an * bn
    return (_th(-an.imag) * _th(-bn.imag) * _th(ab.imag) - _th(an.imag) * _th(bn.imag) * _th(-ab.imag))


def _eta_t_k(a, bt, b):
    """Denner's eta-tilde(a, b~) with b = lim b~ (his eq. 4.42), as an integer."""
    an, btn, bn = mp.mpc(a.n), mp.mpc(bt.n), mp.mpc(b.n)
    if abs(mp.im(bn)) > mp.mpf(10) ** -20 * max(1, abs(bn)):
        return _eta_k(a, bt)
    if mp.re(bn) < 0:
        return _th(-an.imag) * _th(-btn.imag) - _th(an.imag) * _th(btn.imag)
    return 0


def _eta_term(k, v):
    return two_pi_i_times(k) * v if k else V(0, 0)


def _pair(n, cands, exact):
    """V(n, the candidate exact expression closest to n)."""
    if not exact:
        return V(n, None)
    ev = [complex(SR(c).n(prec=220)) for c in cands]
    j = min(range(len(cands)), key=lambda i: abs(ev[i] - complex(n)))
    return V(n, cands[j])


def _quad_roots(a, b, c, d, exact):
    """Roots of a x^2 + b x + c + i ED d: numbers with the i eps, exact roots of a x^2 + b x + c.
    A root at x = 0 (c = 0) makes log(-x_k) singular in the formulas: then another labelling is used."""
    scale = max(1, abs(a.n), abs(b.n))
    if (c.e is not None and c.e.is_trivial_zero()) or abs(c.n) < mp.mpf(10) ** -(DPS // 2) * scale:
        raise ZeroDivisionError("a root at x = 0 in this labelling")
    if abs(a.n) < mp.mpf(10) ** -(DPS // 2) * max(1, abs(b.n), abs(c.n)):
        raise ZeroDivisionError("a = 0 in this labelling")
    # a vanishing discriminant (x1 = x2) is a vanishing Cayley determinant, a leading Landau singularity of
    # the box: decided on the unperturbed numbers, since the i eps alone separates the roots by sqrt(eps)
    disc0 = b.n * b.n - 4 * a.n * c.n
    if abs(disc0) < mp.mpf(10) ** -(DPS // 2) * max(abs(b.n) ** 2, abs(a.n * c.n), mp.mpf(10) ** -300):
        raise ZeroDivisionError("x1 = x2 (vanishing Cayley determinant)")
    cc = c.n + 1j * ED * d.n
    D = mp.sqrt(b.n * b.n - 4 * a.n * cc)
    nums = [(-b.n + D) / (2 * a.n), (-b.n - D) / (2 * a.n)]
    nums0 = [(-b.n + mp.sqrt(b.n * b.n - 4 * a.n * c.n)) / (2 * a.n), (-b.n - mp.sqrt(b.n * b.n - 4 * a.n * c.n)) / (2 * a.n)]
    if abs(nums0[0] - nums[0]) > abs(nums0[1] - nums[0]):
        nums0 = nums0[::-1]
    if exact:
        # the root through vsqrt: a radicand that is exactly real and negative but written with
        # complex radicals becomes I*sqrt(-R), whose evaluation does not jump between branches
        De = vsqrt(b * b - 4 * a * c).e
        cands = [(-b.e + De) / (2 * a.e), (-b.e - De) / (2 * a.e)]
    else:
        cands = None
    xs = [_pair(n, cands, exact) for n in nums]
    x0 = [V(n0, x.e) for n0, x in zip(nums0, xs)]
    if abs(nums[0] - nums[1]) < mp.mpf(10) ** -(DPS // 3):
        raise ZeroDivisionError("x1 = x2 (vanishing Cayley determinant)")
    return xs, x0


def _sgn(x):
    x = mp.re(x)
    return 1 if x > 0 else (-1 if x < 0 else 0)


def _d0_denner(P, M, exact):
    """Denner (4.43), all masses nonzero; P[(i, j)] (i < j) and M[i] are V pairs, r_02 real
    (else all |r_ij| = 1 and the second formula is used)."""
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    K = {}
    for i in range(4):
        for j in range(4):
            if i != j:
                K[(i, j)] = (M[i] * M[i] + M[j] * M[j] - pp(i, j)) / (M[i] * M[j])
    k = lambda i, j: K[(i, j)]
    rt, r = {}, {}
    k02n = mp.re(k(0, 2).n)
    for ij in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        kk = k(*ij)
        kn = kk.n - 1j * ED
        rtn = (kn + mp.sqrt(kn * kn - 4)) / 2
        # With r_02 and r_13 both real, (4.43) holds only when |r_02| and |r_13| lie on the same side of
        # the unit circle (Im r~_02 and Im r~_13 of the same sign); this root formula has |r| > 1 for
        # k > 2, so for k_02 k_13 < 0 the other root is taken for r_13.  (With r_13 complex every
        # choice gives the same value.)  Found by comparing all labellings with Package-X.
        if ij == (1, 3) and abs(mp.re(kk.n)) > 2 and abs(k02n) >= 2 and mp.re(kk.n) * k02n < 0:
            rtn = (kn - mp.sqrt(kn * kn - 4)) / 2
        r0 = [(kk.n + mp.sqrt(kk.n * kk.n - 4)) / 2, (kk.n - mp.sqrt(kk.n * kk.n - 4)) / 2]
        r0n = min(r0, key=lambda z: abs(z - rtn))
        cands = [(kk.e + sqrt(kk.e ** 2 - 4)) / 2, (kk.e - sqrt(kk.e ** 2 - 4)) / 2] if exact else None
        rt[ij] = _pair(rtn, cands, exact)
        r[ij] = V(r0n, rt[ij].e)
    r02, r13, r02t, r13t = r[(0, 2)], r[(1, 3)], rt[(0, 2)], rt[(1, 3)]
    a = k(2, 3) / r13 + r02 * k(0, 1) - k(0, 3) * r02 / r13 - k(1, 2)
    b = (r13 - 1 / r13) * (r02 - 1 / r02) + k(0, 1) * k(2, 3) - k(0, 3) * k(1, 2)
    c = k(0, 1) / r02 + r13 * k(2, 3) - k(0, 3) * r13 / r02 - k(1, 2)
    d = k(1, 2) - r02 * k(0, 1) - r13 * k(2, 3) + r02 * r13 * k(0, 3)
    xs, x0 = _quad_roots(a, b, c, d, exact)
    gam = lambda p, q: _sgn(a.n * (xs[p].n - xs[q].n))
    s = [rt[(0, 3)], rt[(0, 1)], rt[(1, 2)], rt[(2, 3)]]
    one = V.of(1)
    total = V(0, 0)
    for kk in range(2):
        xk = xs[kk]
        xkj = [xk, xk / r13, xk * r02 / r13, xk * r02]
        for j in range(4):
            sg = (-1) ** (j + kk + 1)
            for sj in (s[j], 1 / s[j]):
                z = one + sj * xkj[j]
                total = total + sg * (vli2(z) + _eta_term(_eta_k(-xkj[j], sj), vlog(z)))
    ie = lambda v: V(v, 0)
    if abs(mp.im(r02.n)) > mp.mpf(10) ** -20:            # all |r_ij| = 1
        Pf = lambda y: (sum((k(i, j) * y[i] * y[j] for i in range(4) for j in range(i + 1, 4)), V(0, 0))
                        + sum((v * v for v in y), V(0, 0)))
        zero = V(0, 0)
        for kk in range(2):
            xk, xk0 = xs[kk], x0[kk]
            g = gam(kk, 1 - kk)
            sg = (-1) ** kk
            term = V(0, 0)
            u1 = xk0 / r13
            e1 = _eta_k(-xk, 1 / r13t)
            if e1:
                term = term + _eta_term(e1, vlog(Pf([one, u1, zero, zero]) / u1 - ie(u1.n * ED * b.n * g)) + vlog(u1))
            u2 = r02 * xk0
            e2 = _eta_k(-xk, r02t)
            if e2:
                term = term + _eta_term(e2, vlog(Pf([zero, zero, one, u2]) / u2 - ie(u2.n * ED * b.n * g)) + vlog(u2))
            u3 = r02 * xk0 / r13
            e3 = _eta_k(-xk, r02t / r13t) + _eta_k(r02t, 1 / r13t)
            if e3:
                term = term - _eta_term(e3, vlog(Pf([zero, one, u3, zero]) / u3 - ie(u3.n * ED * b.n * g)) + vlog(u3))
            e4 = (1 - g * _sgn(b.n)) * _eta_k(-xk, -r02t / r13t) * _eta_k(r02t, 1 / r13t)
            if e4:
                term = term + V(-4 * mp.pi ** 2 * e4, -4 * pi ** 2 * e4)
            total = total + sg * term
        return total / (M[0] * M[1] * M[2] * M[3] * a * (xs[0] - xs[1]))
    Q = lambda y0, y1, y3: (1 / r02 - r02) * y0 + (k(1, 2) - r02 * k(0, 1)) * y1 + (k(2, 3) - r02 * k(0, 3)) * y3
    Qb = lambda y0, y2, y3: (1 / r13 - r13) * y3 + (k(1, 2) - r13 * k(2, 3)) * y2 + (k(0, 1) - r13 * k(0, 3)) * y0
    zero = V(0, 0)
    for kk in range(2):
        xk, xk0 = xs[kk], x0[kk]
        g = gam(kk, 1 - kk)
        sg = (-1) ** kk
        e1 = _eta_t_k(-xk, r02t, r02)
        e2 = _eta_t_k(-xk, 1 / r13t, 1 / r13)
        e3 = _eta_t_k(-xk, r02t / r13t, r02 / r13) + _eta_k(r02t, 1 / r13t)
        e4 = _eta_k(r02t, 1 / r13t) * _eta_t_k(-xk, -r02t / r13t, -r02 / r13)
        sI = _sgn(mp.im(r13t.n))
        term = V(0, 0)
        if e1:
            L1 = (vlog(Q(1 / xk0, zero, one) - ie(1j * ED))
                  + vlog(Qb(zero, one, r02 * xk0) / d + ie(1j * ED * g * _sgn(r02.n) * sI)))
            term = term + _eta_term(e1, vlog(r02 * xk) + L1)
        if e2:
            L2 = (vlog(Q(r13 / xk0, one, zero) - ie(1j * ED))
                  + vlog(Qb(one, zero, xk0) / d + ie(1j * ED * g * sI)))
            term = term + _eta_term(e2, vlog(xk / r13) + L2)
        if e3:
            L3 = (vlog(Q(r13 / xk0, one, zero) - ie(1j * ED))
                  + vlog(Qb(zero, one, r02 * xk0) / d + ie(1j * ED * g * _sgn(r02.n) * sI)))
            term = term - _eta_term(e3, vlog(r02 * xk / r13) + L3)
        if e4:
            term = term + V(-4 * mp.pi ** 2 * e4, -4 * pi ** 2 * e4)
        total = total + sg * term
    return total / (M[0] * M[1] * M[2] * M[3] * a * (xs[0] - xs[1]))


def _L2(x1, x2):
    """Denner-Dittmaier's continued Li2(x1, x2) = Li2(1 - x1 x2) + eta(x1, x2) ln(1 - x1 x2)."""
    z = V.of(1) - x1 * x2
    return vli2(z) + _eta_term(_eta_k(x1, x2), vlog(z))


def _d0_dd(P, M, nz, exact):
    """Denner-Dittmaier sec. 3.3, lines arranged so that  nz = 1: m2 = 0;  2: m1 = m2 = 0;
    3: m1 = m2 = m3 = 0;  4: all zero.  P[(i, j)] and M[i] are V pairs (masses not squared)."""
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    M2 = [m * m for m in M]
    Ms = [V(m.n - (1j * EDM if not m.is_zero() else 0), m.e) for m in M2]      # m^2 - i0
    Y, Yr = {}, {}
    for i in range(4):
        for j in range(4):
            if i != j:
                re = M2[i] + M2[j] - pp(i, j)
                Yr[(i, j)] = re
                Y[(i, j)] = V(Ms[i].n + Ms[j].n - pp(i, j).n - 1j * EDM, re.e)
    y = lambda i, j: Y[(i, j)]
    yr = lambda i, j: Yr[(i, j)]

    def rroots(i, j):
        """r_ij,1/2 with m_i^2 + Y x + m_j^2 x^2 = m_i^2 (1 + r1 x)(1 + r2 x)."""
        Yn, Mi, Mj = y(i, j).n, Ms[i].n, Ms[j].n
        D = mp.sqrt(Yn * Yn - 4 * Mi * Mj)
        nums = [(Yn + D) / (2 * Mi), (Yn - D) / (2 * Mi)]
        if exact:
            De = sqrt(yr(i, j).e ** 2 - 4 * M2[i].e * M2[j].e)
            cands = [(yr(i, j).e + De) / (2 * M2[i].e), (yr(i, j).e - De) / (2 * M2[i].e)]
        else:
            cands = None
        return [_pair(n, cands, exact) for n in nums]

    one = V.of(1)
    total = V(0, 0)
    if nz == 1:
        r13 = rroots(1, 3)[0]
        Yn, m1r, m3r = yr(1, 3).n, M2[1].n, M2[3].n
        D = mp.sqrt(Yn * Yn - 4 * m1r * m3r)
        r13r = V(min([(Yn + D) / (2 * m1r), (Yn - D) / (2 * m1r)], key=lambda z: abs(z - r13.n)), r13.e)
        r03, r01 = rroots(0, 3), rroots(0, 1)
        m0, m1, m3 = M2[0], M2[1], M2[3]
        a = m1 * r13r * yr(2, 3) - m3 * yr(1, 2)
        b = yr(0, 2) * (m1 * r13r - m3 / r13r) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2)
        c = yr(0, 2) * (yr(0, 1) - yr(0, 3) / r13r) - m0 * (yr(1, 2) - yr(2, 3) / r13r)
        d = yr(1, 2) - yr(2, 3) / r13r
        xs, _ = _quad_roots(a, b, c, d, exact)
        for kk, xk in enumerate(xs):
            t = (_L2(-xk, y(2, 3) / y(0, 2)) - _L2(-xk, r03[0]) - _L2(-xk, r03[1])
                 - _L2(-xk * r13, y(1, 2) / y(0, 2)) + _L2(-xk * r13, r01[0]) + _L2(-xk * r13, r01[1]))
            e = _eta_k(-xk, r13)
            if e:
                Pn = Ms[0] + y(0, 3) * xk + Ms[3] * xk * xk - V(1j * ED, 0)
                Qn = y(0, 2) + y(2, 3) * xk                   # Q(1, 0, 0, x) at m2 = 0 (DD 3.72)
                t = t + _eta_term(e, vlog(Pn / Qn) + vlog(y(0, 2) / Ms[0]))
            total = total + (-1) ** kk * t
    elif nz in (2, 3):
        m0, m3 = M2[0], M2[3]
        a = yr(1, 3) * yr(2, 3) - m3 * yr(1, 2)
        b = yr(0, 2) * yr(1, 3) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2)
        c = yr(0, 1) * yr(0, 2) - m0 * yr(1, 2)
        d = yr(1, 2)
        xs, _ = _quad_roots(a, b, c, d, exact)
        r03 = rroots(0, 3) if nz == 2 else None
        if nz == 2 and yr(0, 1).is_zero() and yr(2, 3).is_zero():          # DD (3.80)
            for kk, xk in enumerate(xs):
                lx = vlog(-xk)
                t = (-_L2(-xk, r03[0]) - _L2(-xk, r03[1]) - lx * lx / 2
                     - lx * (vlog(y(1, 3) / y(1, 2)) + vlog(y(0, 2) / Ms[0])))
                total = total + (-1) ** kk * t
            return total / (a * (xs[0] - xs[1]))
        for kk, xk in enumerate(xs):
            t = _L2(-xk, y(2, 3) / y(0, 2))
            t = t - (_L2(-xk, r03[0]) + _L2(-xk, r03[1]) if nz == 2 else _L2(-xk, y(0, 3) / Ms[0]))
            t = t + _L2(-xk, y(1, 3) / y(0, 1)) - vlog(-xk) * (vlog(y(0, 1) / y(1, 2)) + vlog(y(0, 2) / Ms[0]))
            total = total + (-1) ** kk * t
    else:
        a = yr(1, 3) * yr(2, 3)
        b = yr(0, 2) * yr(1, 3) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2)
        c = yr(0, 1) * yr(0, 2)
        d = yr(1, 2)
        xs, _ = _quad_roots(a, b, c, d, exact)
        for kk, xk in enumerate(xs):
            t = _L2(-xk, y(2, 3) / y(0, 2)) - _L2(-1 / xk, y(0, 1) / y(1, 3))
            t = t - vlog(-xk) * (vlog(y(1, 3) / y(1, 2)) + vlog(y(0, 2) / y(0, 3)))
            total = total + (-1) ** kk * t
    return total / (a * (xs[0] - xs[1]))


_DD_PATTERN = {1: (1, 1, 0, 1), 2: (1, 0, 0, 1), 3: (1, 0, 0, 0), 4: (0, 0, 0, 0)}
_DD_NEEDS = {1: [(0, 2)], 2: [(0, 2), (0, 1), (1, 2)], 3: [(0, 2), (0, 1), (1, 2)],
             4: [(0, 2), (1, 3), (1, 2), (0, 3)]}


def _d0_terms(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, exact=True):
    import itertools
    inv = {(0, 1): s1, (1, 2): s2, (2, 3): s3, (0, 3): s4, (0, 2): s12, (1, 3): s23}
    inv = {k: V.ofx(v, exact) for k, v in inv.items()}
    M = [V.ofx(m, exact) for m in (m0, m1, m2, m3)]
    pp = lambda i, j: inv[(min(i, j), max(i, j))]
    zero = [SR(m).is_trivial_zero() if not isinstance(m, V) else m.is_zero() for m in (m0, m1, m2, m3)]
    nz = sum(zero)
    perms = list(itertools.permutations(range(4)))
    errors = []
    if nz == 0:
        def k02(pm):
            Mi, Mj = M[pm[0]].n, M[pm[2]].n
            return abs((Mi * Mi + Mj * Mj - pp(pm[0], pm[2]).n) / (Mi * Mj))
        real = [pm for pm in perms if k02(pm) >= 2]
        order = real + [pm for pm in perms if pm not in real] if real else perms
        for pm in order:
            P = {(i, j): pp(pm[i], pm[j]) for i in range(4) for j in range(i + 1, 4)}
            try:
                return _d0_denner(P, [M[q] for q in pm], exact)
            except (ZeroDivisionError, ArithmeticError) as e:
                errors.append(e)
    else:
        for pm in perms:
            if tuple(0 if zero[q] else 1 for q in pm) != _DD_PATTERN[nz]:
                continue
            Mp = [M[q] for q in pm]
            P = {(i, j): pp(pm[i], pm[j]) for i in range(4) for j in range(i + 1, 4)}
            Yr = lambda i, j: Mp[i] * Mp[i] + Mp[j] * Mp[j] - P[(min(i, j), max(i, j))]
            needs = _DD_NEEDS[nz]
            if nz == 2 and Yr(0, 1).is_zero() and Yr(2, 3).is_zero():
                needs = [(0, 2), (1, 2), (1, 3)]                          # DD (3.80)
            if any(Yr(i, j).is_zero() for i, j in needs):
                continue
            try:
                return _d0_dd(P, Mp, nz, exact)
            except (ZeroDivisionError, ArithmeticError) as e:
                errors.append(e)
    raise ZeroDivisionError("no labelling of the lines gives a regular closed form here "
                            "(vanishing Cayley/Gram determinant or a singular point)")


def d0_closed(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, check=True):
    r"""
    D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3) in closed form: an exact Sage expression in
    polylog(2, .) and log of algebraic numbers (16 dilogarithms for nonzero masses, fewer when
    masses vanish), every branch fixed by the +i0.  IR-finite boxes only (IR-divergent ones are
    in feynsage.ir).  check=True compares the expression with the high-precision number.
    """
    from .ir import classify_d0
    if classify_d0(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3) is not None:
        raise ValueError("this D0 is IR divergent: use feynsage.ir.ir_coefficients('D', args)")
    with mp.workdps(DPS):
        val = _d0_terms(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, exact=True)
        if check:
            for prec in (220, 400):
                num = complex(val.e.n(prec=prec))
                if abs(num - complex(val.n)) > 1e-12 * max(1, abs(complex(val.n))):
                    raise ArithmeticError("closed form and its numerical value disagree: %s vs %s" % (num, complex(val.n)))
        return val.e


def d0_value(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, full=False):
    """High-precision value of an IR-finite D0 from the same formulas (complex, ~30 digits)."""
    with mp.workdps(DPS):
        v = _chop(_d0_terms(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, exact=False).n)
        return v if full else complex(v)


def c0_values(arglist, nproc=None):
    """[c0_value(*args) for args in arglist], on every core when the list is long enough."""
    from .parallel import pmap
    return pmap(lambda a: c0_value(*a), list(arglist), nproc=nproc)


def d0_values(arglist, nproc=None):
    """[d0_value(*args) for args in arglist], on every core when the list is long enough."""
    from .parallel import pmap
    return pmap(lambda a: d0_value(*a), list(arglist), nproc=nproc)
