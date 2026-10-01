r"""
IR-divergent triangles and boxes in dimensional regularisation, d = 4 - 2 eps.

The formulas are those of R. K. Ellis and G. Zanderighi, "Scalar one-loop integrals for QCD",
JHEP 02 (2008) 002, arXiv:0712.1851, section 4.3 (six triangles) and 4.4 (sixteen boxes),
which together cover every IR-divergent C0 and D0 up to relabelling.  Their analytic
continuation (their sec. 2.2) is used throughout: every invariant gets +i eps, every mass
-i eps, a log of a product is the sum of the logs of its factors and
    Li2(1 - prod x_i)  ->  Li2(1 - P) + log(1 - P) [log P - sum log x_i]          (|P| <= 1)
                       ->  -Li2c(1/x_i) - (sum log x_i)^2 / 2                     (|P| > 1).
Normalisation: EZ divide by r_Gamma; Package-X (and feynsage) multiply by e^(eps gamma_E),
and e^(eps gamma_E) r_Gamma = 1 - pi^2 eps^2/12 + O(eps^3), so the eps^0 coefficient gets
-pi^2/12 times the 1/eps^2 coefficient.  Two typos of the paper are corrected (and checked
against Package-X): eq. (4.14) has ln(m2 m3/mu^2), eq. (4.27) has Li2.

Arguments follow Package-X:  C0(s1, s12, s2; m0, m1, m2) and
D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3).  The EZ arguments are the same with
masses squared (m_{i+1}^2 = m_i^2).

Each formula is written once against a backend:
  _Num  numbers in, an exact closed form (branches resolved) and its value out;
  _Sym  symbols in, a formula in LogM (log(z - i0)) and Li2c (the continued dilogarithm).
"""
import mpmath as mp
from sage.all import SR, QQ, RR, log, polylog, pi, I, sqrt, function, Integer

from . import scalar as S
from .scalar import V, vlog, vli2, two_pi_i_times, _int_of


# ---------------------------------------------------------------------------- symbolic special functions
def _li2c_num(args, delta=None):
    """Li2(1 - prod z_k^e_k) with every z_k -> z_k - i0 (mpmath)."""
    with mp.workdps(40):
        d = mp.mpf(10) ** -30 if delta is None else delta
        zs = [(mp.mpc(complex(args[2 * k])) - 1j * d * max(1, abs(complex(args[2 * k]))), int(args[2 * k + 1]))
              for k in range(len(args) // 2)]
        return complex(_li2c_mp(zs))


def _li2c_mp(zs):
    P = mp.mpc(1)
    for z, e in zs:
        P *= z ** e
    L = sum(e * mp.log(z) for z, e in zs)
    if abs(P) <= 1:
        return mp.polylog(2, 1 - P) + mp.log(1 - P) * (mp.log(P) - L)
    return -_li2c_mp([(z, -e) for z, e in zs]) - L * L / 2


def _li2c_evalf(self, *args, parent=None, algorithm=None):
    v = _li2c_num(args)
    return parent(v) if parent is not None else v


Li2c = function('Li2c', evalf_func=_li2c_evalf,
                print_latex_func=lambda self, *a: r"\mathrm{Li}_2^{c}\!\left(1 - %s\right)"
                % r"\,".join("(%s)^{%s}" % (a[2 * k]._latex_(), a[2 * k + 1]._latex_()) for k in range(len(a) // 2)))


# ---------------------------------------------------------------------------- backends
class _Num:
    """Numbers: each quantity is a V pair (mpmath value with explicit i eps, exact expression)."""

    def __init__(self):
        self.d = S.EPS

    def x(self, v):
        return V.of(v)

    def f(self, z):
        """The factor z - i0 for a real z."""
        z = V.of(z) if not isinstance(z, V) else z
        return V(z.n - 1j * self.d * max(1, abs(z.n)), z.e)

    def ln(self, *fac):
        """sum e * log(z) over (factor, e) pairs (EZ eq. 2.10)."""
        out = V(0, 0)
        for z, e in fac:
            out = out + e * vlog(z)
        return out

    def li2c(self, *fac):
        """Li2(1 - prod z^e) continued (EZ eqs. 2.11, 2.12)."""
        P = V(1, 1)
        for z, e in fac:
            P = P * (z if e == 1 else (1 / z if e == -1 else _pow(z, e)))
        L = self.ln(*fac)
        if abs(P.n) <= 1:
            res = vli2(V.of(1) - P)
            k = _int_of((vlog(P) - L).n)
            if k:
                res = res + two_pi_i_times(k) * vlog(V.of(1) - P)
            return res
        return -self.li2c(*[(z, -e) for z, e in fac]) - L * L / 2

    def pi2(self):
        return V(mp.pi ** 2, pi ** 2)

    def K(self, z, m, mp_):
        """x = -K(z + i eps, m - i eps, m' - i eps) (EZ eq. 4.15), as a factor with its own i eps."""
        zn, mn, mpn = V.of(z), V.of(m), V.of(mp_)
        if (zn - (mn - mpn) * (mn - mpn)).is_zero():
            return V(1, 1)                                 # -K = 1 at z = (m - m')^2
        d = self.d
        zc = mp.mpc(zn.n) + 1j * d
        mc, mpc = mp.mpc(mn.n) - 1j * d, mp.mpc(mpn.n) - 1j * d
        r = mp.sqrt(1 - 4 * mc * mpc / (zc - (mc - mpc) ** 2))
        num = -(1 - r) / (1 + r)
        re = sqrt(1 - 4 * mn.e * mpn.e / (zn.e - (mn.e - mpn.e) ** 2))
        ex = -(1 - re) / (1 + re)
        # the exact sqrt must be on the branch the number chose
        if abs(complex(ex.n(prec=120)) - complex(num)) > 1e-8 * max(1, abs(complex(num))):
            re = -re
            ex = -(1 - re) / (1 + re)
        return V(num, ex)

    def gamma_pm(self, i_m2, j_m2, p2, sign):
        """gamma_ij^+- of EZ eq. (4.32) with p^2 + i eps and masses - i eps."""
        a, b, p = V.of(i_m2), V.of(j_m2), V.of(p2)
        d = self.d
        an, bn, pn = mp.mpc(a.n) - 1j * d, mp.mpc(b.n) - 1j * d, mp.mpc(p.n) + 1j * d
        u = 1 - (an - bn) / pn
        rt = mp.sqrt(u * u - 4 * bn / pn)
        num = (u + sign * rt) / 2
        ue = 1 - (a.e - b.e) / p.e
        rte = sqrt(ue ** 2 - 4 * b.e / p.e)
        ex = (ue + sign * rte) / 2
        if abs(complex(ex.n(prec=120)) - complex(num)) > 1e-8 * max(1, abs(complex(num))):
            ex = (ue - sign * rte) / 2
        return V(num, ex)


def _pow(z, e):
    out = V(1, 1)
    for _ in range(abs(e)):
        out = out * z
    return out if e > 0 else 1 / out


class _Sym:
    """Symbols: LogM for logarithms of -i0 factors and Li2c for continued dilogarithms."""

    def x(self, v):
        return SR(v)

    def f(self, z):
        return SR(z)

    def ln(self, *fac):
        from .pv import LogM
        return sum((e * (log(z) if _positive(z) else LogM(z)) for z, e in fac), SR(0))

    def li2c(self, *fac):
        args = []
        for z, e in fac:
            args += [SR(z), SR(e)]
        return Li2c(*args)

    def pi2(self):
        return pi ** 2

    def K(self, *a):
        raise _NeedsNumbers()

    def gamma_pm(self, *a):
        raise _NeedsNumbers()


class _NeedsNumbers(Exception):
    pass


def _positive(z):
    """A mass or a mass squared (a symbol, its square, or a positive number): no i0 needed."""
    z = SR(z)
    if z.is_numeric():
        return bool(z > 0)
    if z.is_symbol():
        return True
    op = z.operator()
    return (op is not None and op.__name__ == 'pow' and z.operands()[0].is_symbol()
            and z.operands()[1] == 2)


# ---------------------------------------------------------------------------- helpers
def _z(x):
    """Exact zero test for numbers and symbols (structural, no Maxima)."""
    return SR(x).is_trivial_zero()


def _zero(x):
    return _z(x.e if isinstance(x, V) else x)


def _eq(a, b):
    return SR(a - b).is_trivial_zero() if not isinstance(a, V) else _z(a.e - b.e)


def _laurent_times_power(c, L):
    """(c[-2]/e^2 + c[-1]/e + c[0]) * exp(-e L), to order e^0."""
    m2, m1, z0 = c
    return (m2, m1 - L * m2, z0 - L * m1 + L * L * m2 / 2)


# ---------------------------------------------------------------------------- triangles (EZ 4.3)
def _tri(B, t, p1, p2, p3, M1, M2, M3, m1, m2, m3):
    """EZ triangle t with arguments (p_i^2; m_i^2), masses m_i (not squared) for the K cases."""
    f, ln, li2c = B.f, B.ln, B.li2c
    zero = B.x(0)
    if t == 1:
        L = ln((f(-p3), 1))
        return (1 / p3, -L / p3, L * L / (2 * p3))
    if t == 2:
        L2, L3 = ln((f(-p2), 1)), ln((f(-p3), 1))
        if _eq(p2, p3):                                    # EZ (4.7) at r = 0
            return (zero, -1 / p2, L2 / p2)
        D = p2 - p3
        return (zero, (L3 - L2) / D, (L2 * L2 - L3 * L3) / (2 * D))
    if t == 3:
        Lm = ln((f(M3), 1))
        A2, A3 = ln((f(M3 - p2), 1)), ln((f(M3 - p3), 1))
        if _eq(p2, p3) and _zero(p2):                      # C0(0,0,0; 0,0,m) = B0(0; 0, m)/m^2 - 0
            return (zero, 1 / M3, (1 - Lm) / M3)
        if _eq(p2, p3):                                    # EZ (4.9) at r = 0
            c = (zero, 1 / (M3 - p2), (-Lm - (M3 + p2) / p2 * (A2 - Lm)) / (M3 - p2))
            return c
        D = p2 - p3
        A = A3 - A2
        Bt = (li2c((f(M3 - p2), 1), (f(M3), -1)) - li2c((f(M3 - p3), 1), (f(M3), -1))
              + (A2 - Lm) * (A2 - Lm) - (A3 - Lm) * (A3 - Lm))
        return (zero, A / D, (Bt - Lm * A) / D)
    if t == 4:
        Lm = ln((f(M3), 1))
        X = Lm - ln((f(M3 - p2), 1))
        inner = (B.x(QQ(1) / 2), X, B.pi2() / 12 + X * X / 2 - li2c((f(M3), 1), (f(M3 - p2), -1)))
        c = _laurent_times_power(inner, Lm)
        D = p2 - M3
        return tuple(v / D for v in c)
    if t == 5:
        Lm = ln((f(M3), 1))
        c = _laurent_times_power((zero, B.x(-QQ(1) / 2), B.x(1)), Lm)
        return tuple(v / M3 for v in c)
    if t == 6:                                             # (m2^2, s, m3^2; 0, m2^2, m3^2), s = p2
        xs = B.K(p2, m2, m3)
        lm2, lm3 = ln((f(m2), 1)), ln((f(m3), 1))
        if _is_one(xs):                                    # EZ (4.16)
            if _eq(m2, m3):
                last = B.x(-2)                        # limit of (m3+m2)/(m3-m2) ln(m2/m3)
            else:
                last = (m3 + m2) / (m3 - m2) * (lm2 - lm3)
            return (zero, 1 / (2 * m2 * m3), (-(lm2 + lm3) - 2 - last) / (2 * m2 * m3))
        one = B.x(1)
        pref = xs / (m2 * m3 * (one - xs * xs))
        lx = ln((xs, 1))
        l1mx = ln((one - xs, 1), (one + xs, 1))
        c0 = (lx * (-lx / 2 + 2 * l1mx + lm2 + lm3) - B.pi2() / 6
              + li2c((one - xs, 1), (one + xs, 1)) + (lm2 - lm3) * (lm2 - lm3) / 2
              + li2c((xs, 1), (f(m2), 1), (f(m3), -1)) + li2c((xs, 1), (f(m3), 1), (f(m2), -1)))
        return (zero, -lx * pref, c0 * pref)
    raise ValueError(t)


def _is_one(x):
    if isinstance(x, V):
        return _z(x.e - 1)
    return _z(SR(x) - 1)


# ---------------------------------------------------------------------------- boxes (EZ 4.4)
def _box(B, t, p1, p2, p3, p4, s12, s23, M1, M2, M3, M4, m1, m2, m3, m4):
    f, ln, li2c = B.f, B.ln, B.li2c
    zero, one = B.x(0), B.x(1)
    L = lambda v: ln((f(v), 1))
    if t in (1, 2, 3, 4, 5):
        L12, L23 = L(-s12), L(-s23)
        sq = (L12 - L23) * (L12 - L23)
        if t == 1:
            c = (B.x(4), -2 * (L12 + L23), L12 * L12 + L23 * L23 - sq - B.pi2())
            return tuple(v / (s12 * s23) for v in c)
        if t == 2:
            L4 = L(-p4)
            c = (B.x(2), -2 * (L12 + L23 - L4),
                 L12 * L12 + L23 * L23 - L4 * L4 - 2 * li2c((f(-p4), 1), (f(-s12), -1))
                 - 2 * li2c((f(-p4), 1), (f(-s23), -1)) - sq - B.pi2() / 3)
            return tuple(v / (s12 * s23) for v in c)
        if t in (3, 5) and _zero(s23 * s12 - p2 * p4):
            # EZ (4.21) and (4.25) at r = 1 - p2 p4/(s12 s23) = 0: no residue when s12, s23 (and so
            # p2, p4) have opposite signs; otherwise this is the leading Landau singularity
            sg = lambda v: (v.n.real if isinstance(v, V) else SR(v).n())
            try:
                regular = sg(s12) * sg(s23) < 0
            except TypeError:
                regular = False
            if not regular:
                raise ZeroDivisionError("leading Landau singularity (s12 s23 = p2 p4 with s12 s23 > 0): the box diverges")
            L4 = L(-p4)
            L0 = lambda num, den, zv: (B.x(-1) if _eq(num, den) else (L(-num) - L(-den)) / (one - zv))
            if t == 3:
                c0 = 2 + 2 * (L12 + L23 - L4) + 2 * (L0(p4, s23, p4 / s23) + L0(p4, s12, p4 / s12))
                return (zero, -2 / (s12 * s23), c0 / (s12 * s23))
            L3 = L(-p3)
            c0 = -(-L12 + L3 - L23 - 2 - (one + p4 / s23) * L0(p4, s23, p4 / s23))
            return (zero, -1 / (s12 * s23), c0 / (s12 * s23))
        if t == 3:
            L2, L4 = L(-p2), L(-p4)
            c = (zero, -2 * (L12 + L23 - L2 - L4),
                 L12 * L12 + L23 * L23 - L2 * L2 - L4 * L4
                 - 2 * li2c((f(-p2), 1), (f(-s12), -1)) - 2 * li2c((f(-p2), 1), (f(-s23), -1))
                 - 2 * li2c((f(-p4), 1), (f(-s12), -1)) - 2 * li2c((f(-p4), 1), (f(-s23), -1))
                 + 2 * li2c((f(-p2), 1), (f(-p4), 1), (f(-s12), -1), (f(-s23), -1)) - sq)
            return tuple(v / (s23 * s12 - p2 * p4) for v in c)
        if t == 4:
            L3, L4 = L(-p3), L(-p4)
            X = L3 + L4 - L12
            c = (one, -2 * (L12 + L23 - L3 - L4) - X,
                 L12 * L12 + L23 * L23 - L3 * L3 - L4 * L4 + X * X / 2
                 - 2 * li2c((f(-p3), 1), (f(-s23), -1)) - 2 * li2c((f(-p4), 1), (f(-s23), -1)) - sq)
            return tuple(v / (s12 * s23) for v in c)
        if t == 5:
            L2, L3, L4 = L(-p2), L(-p3), L(-p4)
            X, Y = L2 + L3 - L23, L3 + L4 - L12
            c = (zero, -2 * (L12 + L23 - L2 - L3 - L4) - X - Y,
                 L12 * L12 + L23 * L23 - L2 * L2 - L3 * L3 - L4 * L4 + X * X / 2 + Y * Y / 2
                 - 2 * li2c((f(-p2), 1), (f(-s12), -1)) - 2 * li2c((f(-p4), 1), (f(-s23), -1))
                 + 2 * li2c((f(-p2), 1), (f(-p4), 1), (f(-s12), -1), (f(-s23), -1)) - sq)
            return tuple(v / (s23 * s12 - p2 * p4) for v in c)
    if t in (6, 7, 8, 9, 10):
        Msq = M4
        Lm = L(Msq)
        if t == 6:
            A, Bq = L(Msq - s23) - Lm, L(-s12) - Lm
            c = _laurent_times_power((B.x(2), -(2 * A + Bq), 2 * A * Bq - B.pi2() / 2), Lm)
            return tuple(-v / (s12 * (Msq - s23)) for v in c)
        if t == 7:
            A, Bq, C = L(Msq - s23) - Lm, L(-s12) - Lm, L(Msq - p4) - Lm
            inner = (B.x(QQ(3) / 2), -(2 * A + Bq - C),
                     -2 * li2c((f(Msq - p4), 1), (f(Msq - s23), -1)) + 2 * Bq * A - C * C - 5 * B.pi2() / 12)
            c = _laurent_times_power(inner, Lm)
            return tuple(v / (s12 * (s23 - Msq)) for v in c)
        if t == 8:
            Ls, A3, A4, Ad = L(-s12), L(Msq - p3), L(Msq - p4), L(Msq - s23)
            c = (one, -(Ls + 2 * Ad - A3 - A4),
                 -2 * li2c((f(Msq - p3), 1), (f(Msq - s23), -1)) - 2 * li2c((f(Msq - p4), 1), (f(Msq - s23), -1))
                 - li2c((f(Msq - p3), 1), (f(Msq - p4), 1), (f(-s12), -1), (f(Msq), -1))
                 - B.pi2() / 6 + Ls * Ls / 2 - (Ls - Lm) * (Ls - Lm) / 2 + 2 * Ls * (Ad - Lm)
                 - A3 * (A3 - Lm) - A4 * (A4 - Lm))
            return tuple(v / (s12 * (s23 - Msq)) for v in c)
        if t == 9:
            X = L(-s12) + L(Msq - s23) - L(-p2) - ln((f(m4), 1))
            c = (B.x(QQ(1) / 2), -X,
                 li2c((f(Msq - p3), 1), (f(Msq - s23), 1), (f(Msq), -1), (f(-p2), -1))
                 + 2 * li2c((f(-s12), 1), (f(-p2), -1)) + B.pi2() / 12 + X * X)
            return tuple(v / (s12 * (s23 - Msq)) for v in c)
        if t == 10:
            Y = L(Msq - p4) + L(-p2) - L(Msq - s23) - L(-s12)
            c = (zero, Y,
                 li2c((f(Msq - p3), 1), (f(Msq - s23), 1), (f(-p2), -1), (f(Msq), -1))
                 - li2c((f(Msq - p3), 1), (f(Msq - p4), 1), (f(-s12), -1), (f(Msq), -1))
                 + 2 * li2c((f(Msq - s23), 1), (f(Msq - p4), -1))
                 - 2 * li2c((f(-p2), 1), (f(-s12), -1))
                 + 2 * li2c((f(-p2), 1), (f(Msq - p4), 1), (f(-s12), -1), (f(Msq - s23), -1))
                 + 2 * (ln((f(m4), 1)) - L(Msq - s23)) * Y)
            return tuple(v / (s12 * s23 - Msq * s12 - p2 * p4 + Msq * p2) for v in c)
    if t in (11, 12, 13):
        l3, l4 = ln((f(m3), 1)), ln((f(m4), 1))
        if t == 11:
            A, Bq = L(M3 - s12), L(M4 - s23)
            c0 = 2 * (A - l3) * (Bq - l4) - B.pi2() / 2
            if _zero(p3):                                              # EZ (4.34)
                c0 = c0 - (l3 - l4) * (l3 - l4)
            else:
                c0 = c0 + (l3 - l4) * (l3 - l4)
                for sg in (1, -1):
                    g = B.gamma_pm(M3, M4, p3, sg)
                    r = ln((g, 1), (g - one, -1))
                    c0 = c0 - r * r / 2
            c = (one, -(A + Bq - l3 - l4), c0)
            return tuple(v / ((M3 - s12) * (M4 - s23)) for v in c)
        if t == 12:
            A, Bq, P4 = L(M4 - s23), L(M3 - s12), L(M4 - p4)
            Lm34 = 2 * (l4 - l3)
            c0 = (2 * (A - l3) * (Bq - l3) - (P4 - l3) * (P4 - l3) - B.pi2() / 12
                  + (P4 - Bq) * Lm34 - 2 * li2c((f(M4 - p4), 1), (f(M4 - s23), -1)))
            if _zero(p3):              # EZ (4.36)
                c0 = c0 - Lm34 * Lm34 / 2 - li2c((f(M4 - p4), 1), (f(M3 - s12), -1)) \
                    - li2c((f(m3), 2), (f(m4), -2), (f(M4 - p4), 1), (f(M3 - s12), -1))
            else:
                for sg in (1, -1):
                    g34 = B.gamma_pm(M3, M4, p3, sg)
                    r = ln((g34, 1), (g34 - one, -1))
                    c0 = c0 - r * r / 2
                    g43 = B.gamma_pm(M4, M3, p3, sg)
                    c0 = c0 - li2c((f(M4 - p4), 1), (f(M3 - s12), -1), (g43, 1), (g43 - one, -1))
            c = (B.x(QQ(1) / 2), -(A + Bq - P4 - l3), c0)
            return tuple(v / ((s12 - M3) * (s23 - M4)) for v in c)
        if t == 13:
            P2, P4, A3, A4 = L(M3 - p2), L(M4 - p4), L(M3 - s12), L(M4 - s23)
            Lm3, Lm4 = 2 * l3, 2 * l4
            c0 = (-2 * li2c((f(M3 - p2), 1), (f(M3 - s12), -1)) - 2 * li2c((f(M4 - p4), 1), (f(M4 - s23), -1))
                  + 2 * li2c((f(M3 - p2), 1), (f(M4 - p4), 1), (f(M3 - s12), -1), (f(M4 - s23), -1))
                  + 2 * A3 * A4 - P2 * P2 - P4 * P4 + (P2 - A4) * Lm3 + (P4 - A3) * Lm4)
            if _zero(p3):              # EZ (4.39)
                c0 = (c0 - li2c((f(M3 - p2), 1), (f(M4 - s23), -1))
                      - li2c((f(m4), 2), (f(m3), -2), (f(M3 - p2), 1), (f(M4 - s23), -1))
                      - li2c((f(M4 - p4), 1), (f(M3 - s12), -1))
                      - li2c((f(m3), 2), (f(m4), -2), (f(M4 - p4), 1), (f(M3 - s12), -1))
                      - (Lm4 - Lm3) * (Lm4 - Lm3) / 2)
            else:
                for sg in (1, -1):
                    g34 = B.gamma_pm(M3, M4, p3, sg)
                    g43 = B.gamma_pm(M4, M3, p3, sg)
                    c0 = (c0 - li2c((f(M3 - p2), 1), (f(M4 - s23), -1), (g34, 1), (g34 - one, -1))
                          - li2c((f(M4 - p4), 1), (f(M3 - s12), -1), (g43, 1), (g43 - one, -1)))
                    r = ln((g34, 1), (g34 - one, -1))
                    c0 = c0 - r * r / 2
            D = (M3 - s12) * (M4 - s23) - (M3 - p2) * (M4 - p4)
            c = (zero, P2 + P4 - A3 - A4, c0)
            return tuple(v / D for v in c)
    if t == 14:
        x = B.K(s23, m2, m4)
        Ls = L(-s12)
        if _is_one(x):
            return (zero, 1 / (m2 * m4 * s12), -Ls / (m2 * m4 * s12))
        pref = -2 * x * ln((x, 1)) / (m2 * m4 * s12 * (one - x * x))
        return (zero, pref, -Ls * pref)
    if t == 15:
        x = B.K(s23, m2, m4)
        yf = [(f(m2), 1), (f(M4 - p3), 1), (f(m4), -1), (f(M2 - p2), -1)]
        ly = ln(*yf)
        l2, l4 = ln((f(m2), 1)), ln((f(m4), 1))
        a = L(M2 - p2) - L(-s12)
        b = L(M4 - p3) - L(-s12)
        if _is_one(x):                                      # EZ (4.46)
            y = B.x(1)
            for z, e in yf:
                y = y * (z if e == 1 else 1 / z)
            last = B.x(-2) if _is_one(y) else (one + y) / (one - y) * ly     # limit y -> 1: -2
            return (zero, 1 / (2 * m2 * m4 * s12), (-(l2 + l4) + a + b - 2 - last) / (2 * m2 * m4 * s12))
        lx = ln((x, 1))
        pref = x / (m2 * m4 * s12 * (one - x * x))
        c0 = (lx * (-lx / 2 + l2 + l4 - a - b) - li2c((x, 1), (x, 1)) + ly * ly / 2
              + li2c((x, 1), *yf) + li2c((x, 1), *[(z, -e) for z, e in yf]))
        return (zero, -lx * pref, c0 * pref)
    if t == 16:
        x23, x2, x3 = B.K(s23, m2, m4), B.K(p2, m2, m3), B.K(p3, m3, m4)
        l3 = ln((f(m3), 1))
        if _is_one(x23):                                    # EZ (4.48)
            lx2, lx3 = ln((x2, 1)), ln((x3, 1))
            t1 = (one + x2 * x3) / (one - x2 * x3) * (lx2 + lx3)
            t2 = (x3 + x2) / (x3 - x2) * (lx2 - lx3) if not _eq(x2, x3) else B.x(-2)
            c0 = 2 * (l3 - L(M3 - s12)) - t1 - t2 - 2
            return (zero, 1 / (2 * m2 * m4 * (s12 - M3)), c0 / (2 * m2 * m4 * (s12 - M3)))
        lx = ln((x23, 1))
        pref = x23 / (m2 * m4 * (s12 - M3) * (one - x23 * x23))
        c0 = (-2 * lx * (l3 - L(M3 - s12)) + ln((x2, 1)) * ln((x2, 1)) + ln((x3, 1)) * ln((x3, 1))
              - li2c((x23, 1), (x23, 1))                  # Li2(1 - x23^2)
              + li2c((x23, 1), (x2, 1), (x3, 1)) + li2c((x23, 1), (x2, -1), (x3, -1))
              + li2c((x23, 1), (x2, 1), (x3, -1)) + li2c((x23, 1), (x3, 1), (x2, -1)))
        return (zero, -lx * pref, c0 * pref)
    raise ValueError(t)


# ---------------------------------------------------------------------------- classification
_TRI_ORDER = (1, 2, 5, 4, 3, 6)
_BOX_ORDER = (1, 2, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)


def _tri_labellings(p, m):
    """The 6 relabellings of EZ (sec. 4.3): p = (p1, p2, p3), m = (m1, m2, m3)."""
    out = []
    for pp, mm in ((p, m), ((p[0], p[2], p[1]), (m[1], m[0], m[2]))):
        for k in range(3):
            out.append((pp[k:] + pp[:k], mm[k:] + mm[:k]))
    return out


def _box_labellings(p, s, m):
    """The 8 relabellings of EZ eq. (4.17): p = (p1..p4), s = (s12, s23), m = (m1..m4)."""
    out = []
    refl = ((p[3], p[2], p[1], p[0]), s, (m[0], m[3], m[2], m[1]))
    for pp, ss, mm in ((p, s, m), refl):
        for k in range(4):
            out.append((pp[k:] + pp[:k], ss if k % 2 == 0 else (ss[1], ss[0]), mm[k:] + mm[:k]))
    return out


def _tri_matches(t, p, M):
    z, e = _z, (lambda a, b: _z(a - b))
    p1, p2, p3 = p
    M1, M2, M3 = M
    if t in (1, 2):
        if not (z(M1) and z(M2) and z(M3) and z(p1)):
            return False
        return z(p2) and not z(p3) if t == 1 else (not z(p2) and not z(p3))
    if t in (3, 4, 5):
        if not (z(M1) and z(M2) and not z(M3) and z(p1)):
            return False
        if t == 5:
            return e(p2, M3) and e(p3, M3)
        if t == 4:
            return e(p3, M3) and not e(p2, M3)
        return not e(p2, M3) and not e(p3, M3)
    if t == 6:
        return z(M1) and not z(M2) and not z(M3) and e(p1, M2) and e(p3, M3)


def _box_matches(t, p, s, M):
    z, e = _z, (lambda a, b: _z(a - b))
    p1, p2, p3, p4 = p
    M1, M2, M3, M4 = M
    if t <= 5:
        if not all(z(x) for x in M) or not z(p1):
            return False
        zp = (z(p2), z(p3), z(p4))
        return {1: (True, True, True), 2: (True, True, False), 3: (False, True, False),
                4: (True, False, False), 5: (False, False, False)}[t] == zp
    if t <= 10:
        if not (z(M1) and z(M2) and z(M3) and not z(M4) and z(p1)):
            return False
        if t == 6:
            return z(p2) and e(p3, M4) and e(p4, M4)
        if t == 7:
            return z(p2) and e(p3, M4) and not e(p4, M4)
        if t == 8:
            return z(p2) and not e(p3, M4) and not e(p4, M4)
        if t == 9:
            return not z(p2) and e(p4, M4)
        return not z(p2) and not e(p4, M4)
    if t <= 13:
        if not (z(M1) and z(M2) and not z(M3) and not z(M4) and z(p1)):
            return False
        if t == 11:
            return e(p2, M3) and e(p4, M4)
        if t == 12:
            return e(p2, M3) and not e(p4, M4)
        return not e(p2, M3) and not e(p4, M4)
    if t <= 15:
        if not (z(M1) and not z(M2) and z(M3) and not z(M4) and e(p1, M2) and e(p4, M4)):
            return False
        both = e(p2, M2) and e(p3, M4)
        return both if t == 14 else not both
    return z(M1) and not z(M2) and not z(M3) and not z(M4) and e(p1, M2) and e(p4, M4)


def _ir_divergent(P, M):
    """Soft or collinear divergence (Kinoshita's rules, as in Denner-Dittmaier sec. 4), with
    P[(i, j)] = (q_i - q_j)^2 for every pair of lines (adjacent and, for boxes, diagonal):
    collinear if P_ij = 0 with m_i = m_j = 0; soft at a massless line i if P_ij = m_j^2 and
    P_ik = m_k^2 for two other lines j, k."""
    N = len(M)
    z, e = _z, (lambda a, b: _z(a - b))
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    for i in range(N):
        for j in range(i + 1, N):
            if z(pp(i, j)) and z(M[i]) and z(M[j]):
                return True
    for i in range(N):
        if not z(M[i]):
            continue
        others = [j for j in range(N) if j != i and e(pp(i, j), M[j])]
        if len(others) >= 2:
            return True
    return False


def classify_c0(s1, s12, s2, m0, m1, m2):
    """None for an IR-finite C0, else (EZ triangle number, (p1, p2, p3), (m1, m2, m3)) in EZ labels."""
    p, m = tuple(SR(x) for x in (s1, s12, s2)), tuple(SR(x) for x in (m0, m1, m2))
    M = tuple(x ** 2 for x in m)
    if not _ir_divergent({(0, 1): p[0], (1, 2): p[1], (0, 2): p[2]}, M):
        return None
    if all(_z(x) for x in p + M):
        return ('0', p, m)
    for t in _TRI_ORDER:
        for pp, mm in _tri_labellings(p, m):
            if _tri_matches(t, pp, tuple(x ** 2 for x in mm)):
                return (t, pp, mm)
    raise NotImplementedError("IR-divergent C0%s not in the Ellis-Zanderighi list" % ((s1, s12, s2, m0, m1, m2),))


def classify_d0(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3):
    """None for an IR-finite D0, else (EZ box number, p, (s12, s23), masses) in EZ labels."""
    p = tuple(SR(x) for x in (s1, s2, s3, s4))
    s = (SR(s12), SR(s23))
    m = tuple(SR(x) for x in (m0, m1, m2, m3))
    M = tuple(x ** 2 for x in m)
    if not _ir_divergent({(0, 1): p[0], (1, 2): p[1], (2, 3): p[2], (0, 3): p[3], (0, 2): s[0], (1, 3): s[1]}, M):
        return None
    if all(_z(x) for x in p + s + M):
        return ('0', p, s, m)
    for t in _BOX_ORDER:
        for pp, ss, mm in _box_labellings(p, s, m):
            if _box_matches(t, pp, ss, tuple(x ** 2 for x in mm)):
                return (t, pp, ss, mm)
    raise NotImplementedError("IR-divergent D0%s at an exceptional point outside the Ellis-Zanderighi list "
                              "(for example a vanishing s12 or s23 between massless lines)"
                              % ((s1, s2, s3, s4, s12, s23, m0, m1, m2, m3),))


# ---------------------------------------------------------------------------- evaluation
def _vnum(x):
    x = SR(x)
    try:
        return V.of(QQ(x))
    except (TypeError, ValueError):
        return V(mp.mpc(complex(x.n(prec=220))), x)


def _laurent(kind, cls, B):
    """(c_-2, c_-1, c_0) at mu = 1 in the Package-X normalisation, in the backend B."""
    if kind == 'C':
        t, p, m = cls
        conv = (lambda x: _vnum(x)) if isinstance(B, _Num) else SR
        p, m = [conv(x) for x in p], [conv(x) for x in m]
        M = [x * x for x in m]
        c = _tri(B, t, *p, *M, *m)
    else:
        t, p, s, m = cls
        conv = (lambda x: _vnum(x)) if isinstance(B, _Num) else SR
        p, s, m = [conv(x) for x in p], [conv(x) for x in s], [conv(x) for x in m]
        M = [x * x for x in m]
        c = _box(B, t, *p, *s, *M, *m)
    c = [B.x(v) if not isinstance(v, (V,)) and isinstance(B, _Num) else v for v in c]
    return (c[0], c[1], c[2] - B.pi2() / 12 * c[0])         # EZ (r_Gamma) -> Package-X (e^(eps gamma))


def _numeric_args(args):
    return all(SR(a).is_numeric() or SR(a).is_constant() for a in args)


def ir_coefficients(kind, args, exact=True):
    r"""
    For an IR-divergent C0 (kind 'C', six Package-X arguments) or D0 (kind 'D', ten):
    the Laurent coefficients (c_-2, c_-1, c_0) at mu = 1, as Sage expressions.
    Numbers in: closed forms in log and polylog(2, .) with exact arguments.
    Symbols in: logs and Li2c (when the formula needs the roots K or gamma, the coefficients
    are the named functions C0IR(k, ...) / D0IR(k, ...) instead).  None if the integral is finite.
    """
    cls = (classify_c0 if kind == 'C' else classify_d0)(*args)
    if cls is None:
        return None
    if cls[0] == '0':
        return (SR(0), SR(0), SR(0))
    if _numeric_args(args):
        with mp.workdps(S.DPS):
            c = _laurent(kind, cls, _Num())
            return tuple(v.e if exact else SR(complex(v.n)) for v in c)
    try:
        return _laurent(kind, cls, _Sym())
    except _NeedsNumbers:
        F = C0IR if kind == 'C' else D0IR
        return tuple(F(k, *args) for k in (-2, -1, 0))


def ir_value(kind, args):
    """The three Laurent coefficients at mu = 1 as complex numbers (high precision underneath)."""
    cls = (classify_c0 if kind == 'C' else classify_d0)(*args)
    if cls is None:
        return None
    if cls[0] == '0':
        return (0j, 0j, 0j)
    with mp.workdps(S.DPS):
        return tuple(complex(v.n) for v in _laurent(kind, cls, _Num()))


def _irf_evalf(kind):
    def ev(self, k, *args, parent=None, algorithm=None):
        v = ir_value(kind, [RR(a).simplest_rational() if not SR(a).is_integer() else Integer(a) for a in args])
        out = v[{-2: 0, -1: 1, 0: 2}[int(k)]]
        return parent(out) if parent is not None else out
    return ev


C0IR = function('C0IR', nargs=7, evalf_func=_irf_evalf('C'))
D0IR = function('D0IR', nargs=11, evalf_func=_irf_evalf('D'))


def with_mu(c):
    """c_-2/eps^2 + c_-1/eps + c_0 with the scale mu restored (the integral carries mu^(2 eps))."""
    from .pv import eps, mu
    L = log(mu ** 2)
    return c[0] / eps ** 2 + (c[1] + L * c[0]) / eps + c[2] + L * c[1] + L ** 2 / 2 * c[0]
