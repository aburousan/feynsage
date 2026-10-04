r"""
C0 and D0 in closed form with symbols (Package-X's C0Expand and D0Expand).

The numerical closed forms of feynsage.scalar decide every branch of a log, a square root or a
dilogarithm from the numbers.  Here the invariants and masses may be symbols, so the side of each
cut is kept in the formula itself, the way the papers write it: the infinitesimals are written in
as two symbols,

    fs_eps   the i0 of the Feynman denominators (and of Denner's k_ij - i eps),
    fs_del   the -i0 of the internal masses (needed only with massless lines), fs_del << fs_eps,

and the formula holds in the limit fs_eps, fs_del -> 0+.  Conditional.n() evaluates it there,
with fs_eps = 1e-50 and fs_del = 1e-64 at 600 bits, as feynsage.scalar does with numbers.
Square roots, logarithms and dilogarithms are the ordinary principal ones.  Four small functions
complete the formulas:

    Eta(a, b)        eta(a, b) = log(ab) - log a - log b  (0 or +-2 pi i)
    EtaT(a, b~, b)   Denner's eta-tilde, his eq. (4.42), b = lim b~
    SqrtL(z0, z)     lim sqrt(z) as z -> z0, the root of z0 closest to sqrt(z) (for r = lim r~ and
                     the other quantities that the formulas use without their i0)
    ReSign(z)        the sign of Re z
    RootSign(k02, k13)  which root r_13 is (Denner's (4.43) needs |r_02| and |r_13| on the same side
                     of the unit circle when both are real)

As in the papers, the formula is written with named intermediate quantities (the roots x_k, the
r_ij, ...), whose definitions come with it (Conditional.defs); Conditional.inline() substitutes
them and gives one expression.

The formulas are those of the numerical code, with the same infinitesimals:
    C0: A. Denner, Fortsch. Phys. 41 (1993) 307, eqs. (4.26), (4.27) (lambda(s1, s12, s2) > 0).
    D0, all masses nonzero: Denner eq. (4.43), lines labelled so that r_02 is real.
    D0, one or more masses zero: A. Denner and S. Dittmaier, Nucl. Phys. B 844 (2011) 199,
        eqs. (3.76), (3.78), (3.80), (3.82), (3.84).

    c0_expand(s1, s12, s2, m0, m1, m2)          -> Conditional
    d0_expand(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, pair=(1, 3))
    expand_c0d0(expr)                            every C0(...) and D0(...) inside expr

A Conditional holds the expression, its definitions and the conditions under which it is valid,
like Package-X's ConditionalExpression.  Its value at a point is .n(subs).
"""
import itertools
import mpmath as mp
from sage.all import SR, I, sqrt, log, polylog, function, QQ, ComplexField

from ._parse import sr

__all__ = ['c0_expand', 'd0_expand', 'expand_c0d0', 'Conditional', 'SqrtL', 'Eta', 'EtaT', 'ReSign', 'RootSign']


# ---------------------------------------------------------------------------- the infinitesimals
EPS = SR.var('fs_eps', latex_name=r'\epsilon')       # the i0 of the Feynman denominators
DEL = SR.var('fs_del', latex_name=r'\delta')         # the -i0 of the internal masses
EPS_VALUE, DEL_VALUE = QQ(10) ** -50, QQ(10) ** -64   # as in feynsage.scalar: EPS^2 << DEL << EPS
WORK_PREC = 600


def _prec(parent):
    return getattr(parent, 'prec', lambda: 53)() if parent is not None else 53


def _m(x, prec):
    """A Sage number (any parent) as an mpmath mpc at prec bits, through decimal strings: every digit is
    kept and no integer type has to be shared with mpmath (some installations mix gmpy versions)."""
    try:
        z = ComplexField(prec)(x)
        return mp.mpc(mp.mpf(z.real().str(truncate=False)), mp.mpf(z.imag().str(truncate=False)))
    except (TypeError, ValueError):
        return mp.mpc(complex(x))


def _out(v, parent):
    from .pv import _to_parent
    return _to_parent(v, parent)


def _tiny(prec):
    """Below this relative size an imaginary part is rounding, far below the infinitesimals."""
    return mp.mpf(2) ** (-int(prec * 0.75))


def _imsign(z, prec):
    if abs(mp.im(z)) <= _tiny(prec) * max(1, abs(z)):
        return 0
    return 1 if mp.im(z) > 0 else -1


def _ev(fn):
    def evalf(self, *args, parent=None, algorithm=None):
        prec = _prec(parent)
        with mp.workprec(prec + 30):
            return _out(fn(*[_m(a, prec + 30) for a in args], prec=prec), parent)
    return evalf


def _sqrtl(z0, zf, prec):
    """lim sqrt(z) for z = zf -> z0: the root of z0 closest to the principal sqrt(zf)."""
    s, t = mp.sqrt(z0), mp.sqrt(zf)
    return s if abs(s - t) <= abs(s + t) else -s


def _eta(a, b, prec):
    """eta(a, b)/(2 pi i) = theta(-Im a) theta(-Im b) theta(Im ab) - theta(Im a) theta(Im b) theta(-Im ab)."""
    sa, sb, sab = _imsign(a, prec), _imsign(b, prec), _imsign(a * b, prec)
    th = lambda s: 1 if s > 0 else 0
    return th(-sa) * th(-sb) * th(sab) - th(sa) * th(sb) * th(-sab)


def _eta_t(a, bt, b, prec):
    """Denner's eta-tilde(a, b~)/(2 pi i), his eq. (4.42); b = lim b~."""
    if abs(mp.im(b)) > mp.mpf(10) ** -20 * max(1, abs(b)):
        return _eta(a, bt, prec)
    if mp.re(b) < 0:
        th = lambda s: 1 if s > 0 else 0
        sa, sb = _imsign(a, prec), _imsign(bt, prec)
        return th(-sa) * th(-sb) - th(sa) * th(sb)
    return 0


def _resign(z, prec):
    if abs(mp.re(z)) <= _tiny(prec) * max(1, abs(z)):
        return 0
    return 1 if mp.re(z) > 0 else -1


def _lat(x):
    return x._latex_() if hasattr(x, '_latex_') else str(x)


SqrtL = function('SqrtL', nargs=2, evalf_func=_ev(lambda a, b, prec: _sqrtl(a, b, prec)),
                 derivative_func=lambda self, z0, zf, diff_param=None: 1 / (2 * SqrtL(z0, zf)) if diff_param == 0 else 0,
                 print_latex_func=lambda self, z0, zf: r'\sqrt{%s}' % _lat(z0))
Eta = function('Eta', nargs=2, evalf_func=_ev(lambda a, b, prec: 2j * mp.pi * _eta(a, b, prec)),
               derivative_func=lambda self, *a, diff_param=None: 0,
               print_latex_func=lambda self, a, b: r'\eta\left(%s,\ %s\right)' % (_lat(a), _lat(b)))
EtaT = function('EtaT', nargs=3, evalf_func=_ev(lambda a, bt, b, prec: 2j * mp.pi * _eta_t(a, bt, b, prec)),
                derivative_func=lambda self, *a, diff_param=None: 0,
                print_latex_func=lambda self, a, bt, b: r'\tilde\eta\left(%s,\ %s\right)' % (_lat(a), _lat(bt)))
def _rootsign(k02, k13, prec):
    a, b = mp.re(k02), mp.re(k13)
    return -1 if abs(b) > 2 and abs(a) >= 2 and a * b < 0 else 1


RootSign = function('RootSign', nargs=2, evalf_func=_ev(lambda a, b, prec: _rootsign(a, b, prec)),
                    derivative_func=lambda self, a, b, diff_param=None: 0,
                    print_latex_func=lambda self, a, b: r'\sigma\left(%s,\ %s\right)' % (_lat(a), _lat(b)))
ReSign = function('ReSign', nargs=1, evalf_func=_ev(lambda z, prec: _resign(z, prec)),
                  derivative_func=lambda self, z, diff_param=None: 0,
                  print_latex_func=lambda self, z: r'\mathrm{sgn}\,\mathrm{Re}\left(%s\right)' % _lat(z))


# ---------------------------------------------------------------------------- values with their i0
class J:
    """A quantity twice: e0 its limit (EPS, DEL -> 0, the side of every cut kept by SqrtL) and ef
    the full expression with the infinitesimals EPS and DEL written in."""
    __slots__ = ('e0', 'ef')

    def __init__(self, e0, ef=None):
        self.e0 = SR(e0)
        self.ef = SR(e0) if ef is None else SR(ef)

    @staticmethod
    def of(x):
        return x if isinstance(x, J) else J(x)

    @staticmethod
    def inf(e0, eps=0, dl=0):
        """e0 plus i0 terms: eps and dl are the coefficients of EPS and DEL."""
        return J(e0, SR(e0) + eps * EPS + dl * DEL)

    def __add__(s, o): o = J.of(o); return J(s.e0 + o.e0, s.ef + o.ef)
    __radd__ = __add__
    def __sub__(s, o): o = J.of(o); return J(s.e0 - o.e0, s.ef - o.ef)
    def __rsub__(s, o): return J.of(o) - s
    def __neg__(s): return J(-s.e0, -s.ef)
    def __mul__(s, o): o = J.of(o); return J(s.e0 * o.e0, s.ef * o.ef)
    __rmul__ = __mul__
    def __truediv__(s, o): o = J.of(o); return J(s.e0 / o.e0, s.ef / o.ef)
    def __rtruediv__(s, o): return J.of(o) / s

    def limit(s):
        return J(s.e0)


def jsqrt(z):
    return J(SqrtL(z.e0, z.ef), sqrt(z.ef))


def jlog(z):
    return log(z.ef)


def jli2(z):
    return polylog(2, z.ef)


def jL2(x1, x2):
    """Li2(1 - x1 x2) + eta(x1, x2) log(1 - x1 x2)."""
    w = 1 - x1.ef * x2.ef
    return polylog(2, w) + Eta(x1.ef, x2.ef) * log(w)


def _zero(x):
    x = SR(x)
    return x.is_trivial_zero() or bool(x.expand().is_trivial_zero())


_COUNTER = [0]


class _Ctx:
    """Names the intermediate quantities of one formula; the definitions are kept in order."""

    def __init__(self, tag):
        _COUNTER[0] += 1
        self.tag = "%s%d" % (tag, _COUNTER[0])
        self.defs = []

    def name(self, label, latex, q):
        """A name for q (J or expression): returns J(symbol for the limit, symbol for the full value)."""
        q = J.of(q)
        s0 = SR.var('fs%s_%s' % (self.tag, label), latex_name=latex)
        self.defs.append((s0, q.e0))
        if bool((q.ef - q.e0).is_trivial_zero()):
            return J(s0)
        sf = SR.var('fs%s_%s_t' % (self.tag, label), latex_name=r'\tilde{%s}' % latex)
        self.defs.append((sf, q.ef))
        return J(s0, sf)


# ---------------------------------------------------------------------------- Conditional
class Conditional:
    """An expression with its definitions and the conditions under which it holds (Package-X's
    ConditionalExpression).  The expression contains the infinitesimals fs_eps and fs_del (the
    i0's, -> 0+) and named intermediates (.defs, in order).  .n() gives values; .inline() one
    expression with every definition put in."""

    def __init__(self, value, conditions=(), defs=()):
        self.value = SR(value)
        self.defs = list(defs)
        self.conditions = []
        for c in conditions:
            if not any(str(c) == str(d) for d in self.conditions):
                self.conditions.append(c)

    def __repr__(self):
        cond = " and ".join(str(c) for c in self.conditions) if self.conditions else "always"
        return "Conditional(<expression, %d leaves, %d definitions>, valid if %s)" % (
            self.leaves(), len(self.defs), cond)

    def _latex_(self):
        cond = r",\ ".join(c._latex_() if hasattr(c, '_latex_') else str(c) for c in self.conditions)
        return self.value._latex_() + (r"\qquad \left(%s\right)" % cond if cond else "")

    @staticmethod
    def _count(e):
        n = [0]
        def walk(x):
            ops = x.operands()
            if not ops:
                n[0] += 1
            for o in ops:
                walk(o)
        walk(SR(e))
        return n[0]

    def leaves(self):
        return self._count(self.value) + sum(self._count(d) for _, d in self.defs)

    def inline(self):
        """The value (and conditions) with every definition substituted: one expression."""
        v, conds = self.value, list(self.conditions)
        for s, d in reversed(self.defs):
            v = v.subs({s: d})
            conds = [c.subs({s: d}) if hasattr(c, 'subs') else c for c in conds]
        return Conditional(v, conds)

    def subs(self, *a, **k):
        return Conditional(self.value.subs(*a, **k), [c.subs(*a, **k) if hasattr(c, 'subs') else c for c in self.conditions],
                           [(s, d.subs(*a, **k)) for s, d in self.defs])

    def _point(self, a, k, prec):
        """The definitions evaluated at the point, in order (numbers at prec bits)."""
        pt = {}
        for x, v in list(dict(a[0]).items() if a else []) + list(k.items()):
            pt[SR.var(x) if isinstance(x, str) else x] = v
        pt[EPS], pt[DEL] = EPS_VALUE, DEL_VALUE
        CF = ComplexField(prec)
        vals = dict(pt)
        bad = []
        for s, d in self.defs:
            try:
                z = CF(SR(d).subs(vals).n(prec=prec))
            except (ValueError, ZeroDivisionError, TypeError, ArithmeticError) as e:
                bad.append((s, e))                    # singular here; the conditions say why
                continue
            vals[s] = z if z.imag() != 0 else z.real()
        self._bad = bad
        return vals

    def holds(self, *a, **k):
        """The conditions at a point: True, False, or None if one cannot be decided."""
        vals = self._point(a, k, WORK_PREC)
        CF = ComplexField(WORK_PREC)
        res = True
        for c in self.conditions:
            if not hasattr(c, 'lhs'):
                res = None
                continue
            try:
                lhs = CF(c.lhs().subs(vals).n(prec=WORK_PREC))
                rhs = CF(c.rhs().subs(vals).n(prec=WORK_PREC))
                d = lhs - rhs
                tol = 1e-40 * max(1, abs(lhs), abs(rhs))
                op = c.operator().__name__
                if op == 'ne':
                    ok = abs(d) > tol
                elif op in ('gt', 'ge', 'lt', 'le'):
                    x = d.real()
                    ok = {'gt': x > tol, 'ge': x > -tol, 'lt': x < -tol, 'le': x < tol}[op]
                else:
                    ok = abs(d) <= tol
                if not ok:
                    return False
            except (TypeError, ValueError, AttributeError, ZeroDivisionError):
                res = None
        return res

    def n(self, *a, prec=53, digits=None, check=True, **k):
        """Value at a point: C.n({s: 3, m: 1}) or C.n(s=3, m=1)."""
        if digits is not None:
            prec = int(digits * 3.33) + 4
        if check and (a or k) and self.holds(*a, **k) is False:
            raise ValueError("the point is outside the region where this closed form holds: %s" % self.conditions)
        w = max(WORK_PREC, 4 * prec)
        vals = self._point(a, k, w)
        if self._bad:
            raise ValueError("singular point of this closed form (%s): %s" % (self._bad[0][0], self._bad[0][1]))
        z = ComplexField(w)(self.value.subs(vals).n(prec=w))
        if abs(z.imag()) < 1e-20 * max(1, abs(z)):          # what is left of the i0 (sqrt(eps) on a threshold)
            z = ComplexField(w)(z.real())
        z = ComplexField(prec)(z)
        return z if z.imag() != 0 else z.real()

    def __call__(self, *a, **k):
        return self.n(*a, **k)

    def n_many(self, points, prec=53, check=True, nproc=None):
        """Values at many points (a list of dicts), on every core when that pays."""
        from .parallel import pmap
        return pmap(lambda p: self.n(p, prec=prec, check=check), list(points), nproc=nproc)


# ---------------------------------------------------------------------------- C0
def _roots_j(a, b, c):
    """Roots of a y^2 - b y + c - i eps (a, b, c exact; a may be literally 0) as J."""
    ce = J.inf(c, -I)
    if _zero(a):
        if _zero(b):
            return []
        return [ce / b]
    d = jsqrt(J(b) * J(b) - 4 * a * ce)
    return [(J(b) + d) / (2 * a), (J(b) - d) / (2 * a)]


def c0_expand(s1, s12, s2, m0, m1, m2, _alpha_sign=1):
    r"""
    C0(s1, s12, s2; m0, m1, m2) as an explicit formula in the arguments (Package-X's C0Expand):
    the twelve dilogarithms of Denner (1993), eq. (4.26), with y0_i and the roots y_i+- of his
    eq. (4.27) as named intermediates.  The i0 that picks the side of a cut is written in as
    fs_eps.  Valid where lambda(s1, s12, s2) > 0, the region of Package-X's C0Expand.  Arguments
    in Package-X order, masses not squared.  Returns a Conditional.

        C = c0_expand(s1, s12, s2, m0, m1, m2)
        C.n(s1=-1, s12=2, s2=-3, m0=1, m1=2, m2=3)
    """
    from .scalar import _c0_zero_momenta
    s1, s12, s2, m0, m1, m2 = [sr(x) for x in (s1, s12, s2, m0, m1, m2)]
    if all(_zero(x) for x in (s1, s12, s2)):
        return Conditional(_c0_zero_momenta(m0, m1, m2), [])
    ctx = _Ctx('C')
    P = {(1, 0): s1, (2, 0): s2, (2, 1): s12}
    p2 = lambda a, b: P[(max(a, b), min(a, b))]
    Ms = [m0 ** 2, m1 ** 2, m2 ** 2]
    lam = s1 ** 2 + s12 ** 2 + s2 ** 2 - 2 * (s1 * s12 + s12 * s2 + s2 * s1)
    if _zero(lam):
        raise ZeroDivisionError("lambda(s1, s12, s2) = 0 (vanishing Gram determinant): no closed form of this kind")
    conds = [lam > 0]
    # The formula holds for either root alpha = +-sqrt(lambda).  With an invariant p_jk literally zero,
    # lambda = (p_ki - p_ij)^2 and alpha = p_ij - p_ki makes that term's y0 regular everywhere (no
    # case split on the sign of p_ki - p_ij, which sqrt(lambda) = |p_ki - p_ij| would need).
    alpha = _alpha_sign * sqrt(lam)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        if _zero(p2(j, k)):
            alpha = p2(i, j) - p2(k, i)
            break
    alpha_raw = alpha
    alpha = ctx.name('alpha', r'\alpha', alpha).e0
    total = SR(0)
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        pjk, pki, pij = p2(j, k), p2(k, i), p2(i, j)
        mi, mj, mk = Ms[i], Ms[j], Ms[k]
        if _zero(pjk):
            N0 = (mk - mj) * (pki - pij + alpha_raw)
            if not _zero(N0):
                conds.append(N0 != 0)
                continue                                   # y0_i -> infinity: no contribution
            dalpha = -(pki + pij) / alpha
            y0 = ((-pki - pij + 2 * mi - mj - mk) + dalpha * (mk - mj) + alpha) / (2 * alpha)
        else:
            y0 = (pjk * (pjk - pki - pij + 2 * mi - mj - mk) - (pki - pij) * (mj - mk)
                  + alpha * (pjk - mj + mk)) / (2 * alpha * pjk)
        y0 = ctx.name('y%d' % i, r'y_{0%d}' % i, y0).e0
        for n, x in enumerate(_roots_j(pjk, pjk - mj + mk, mk)):
            x = ctx.name('x%d%s' % (i, 'pm'[n]), r'x_{%d%s}' % (i, '+-'[n]), x)
            y = J(y0) - x
            total += jli2((y0 - 1) / y) - jli2(y0 / y)
    return Conditional(total / alpha, conds, ctx.defs)


# ---------------------------------------------------------------------------- D0
def _quad_j(ctx, a, b, c, d):
    """Roots of a x^2 + b x + c + i eps d (a, b, c, d exact) as named J, and the conditions."""
    for v, what in ((c, "c"), (a, "a")):
        if _zero(v):
            raise ZeroDivisionError("%s = 0 in this labelling" % what)
    disc = b * b - 4 * a * c
    if _zero(disc):
        raise ZeroDivisionError("x1 = x2 (vanishing Cayley determinant)")
    D = jsqrt(J.inf(disc, -4 * I * a * d))
    xs = [ctx.name('x%d' % (n + 1), r'x_{%d}' % (n + 1), (-J(b) + sg * D) / (2 * a)) for n, sg in enumerate((1, -1))]
    return xs, [a != 0, c != 0, disc != 0]


def _d0_denner_sym(ctx, P, M):
    """Denner (4.43), all masses nonzero, r_02 real (that condition is returned first)."""
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    Kv = {}
    for i, j in itertools.combinations(range(4), 2):
        Kv[(i, j)] = Kv[(j, i)] = ctx.name('k%d%d' % (i, j), r'k_{%d%d}' % (i, j),
                                           (M[i] ** 2 + M[j] ** 2 - pp(i, j)) / (M[i] * M[j])).e0
    k = lambda i, j: Kv[(i, j)]
    rt, r = {}, {}
    for ij in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        kt = J.inf(k(*ij), -I)                                # k - i eps
        # r_13 on the same side of the unit circle as r_02 when both are real (see scalar._d0_denner)
        sg = RootSign(k(0, 2), k(1, 3)) if ij == (1, 3) else 1
        rt[ij] = ctx.name('r%d%d' % ij, r'r_{%d%d}' % ij, (kt + sg * jsqrt(kt * kt - 4)) / 2)
        r[ij] = rt[ij].limit()
    r02, r13, r02t, r13t = r[(0, 2)], r[(1, 3)], rt[(0, 2)], rt[(1, 3)]
    K = lambda i, j: J(k(i, j))
    a = ctx.name('a', 'a', K(2, 3) / r13 + r02 * K(0, 1) - K(0, 3) * r02 / r13 - K(1, 2))
    b = ctx.name('b', 'b', (r13 - 1 / r13) * (r02 - 1 / r02) + K(0, 1) * K(2, 3) - K(0, 3) * K(1, 2))
    c = ctx.name('c', 'c', K(0, 1) / r02 + r13 * K(2, 3) - K(0, 3) * r13 / r02 - K(1, 2))
    d = ctx.name('d', 'd', K(1, 2) - r02 * K(0, 1) - r13 * K(2, 3) + r02 * r13 * K(0, 3))
    xs, conds = _quad_j(ctx, a.e0, b.e0, c.e0, d.e0)
    conds = [(pp(0, 2) - (M[0] - M[2]) ** 2) * (pp(0, 2) - (M[0] + M[2]) ** 2) >= 0] + conds
    s = [rt[(0, 3)], rt[(0, 1)], rt[(1, 2)], rt[(2, 3)]]
    total = SR(0)
    for kk in range(2):
        xk = xs[kk]
        xkj = [xk, xk / r13, xk * r02 / r13, xk * r02]
        for j in range(4):
            sg = (-1) ** (j + kk + 1)
            for sj in (s[j], 1 / s[j]):
                total += sg * jL2(-xkj[j], sj)
    Q = lambda y0, y1, y3: (1 / r02 - r02) * y0 + (K(1, 2) - r02 * K(0, 1)) * y1 + (K(2, 3) - r02 * K(0, 3)) * y3
    Qb = lambda y0, y2, y3: (1 / r13 - r13) * y3 + (K(1, 2) - r13 * K(2, 3)) * y2 + (K(0, 1) - r13 * K(0, 3)) * y0
    one, zero = J(1), J(0)
    sI = -ReSign(k(1, 3))                                  # sign of Im r~_13
    sr02 = ReSign(r02.e0)
    ie = lambda cf: J.inf(0, I * cf)
    for kk in range(2):
        xk, xk0 = xs[kk], xs[kk].limit()
        g = ReSign(a.e0 * (xs[kk].e0 - xs[1 - kk].e0))
        sg = (-1) ** kk
        L1 = jlog(Q(1 / xk0, zero, one) + ie(-1)) + jlog(Qb(zero, one, r02 * xk0) / d + ie(g * sr02 * sI))
        L2 = jlog(Q(r13 / xk0, one, zero) + ie(-1)) + jlog(Qb(one, zero, xk0) / d + ie(g * sI))
        L3 = jlog(Q(r13 / xk0, one, zero) + ie(-1)) + jlog(Qb(zero, one, r02 * xk0) / d + ie(g * sr02 * sI))
        e1 = EtaT((-xk).ef, r02t.ef, r02t.e0)
        e2 = EtaT((-xk).ef, (1 / r13t).ef, (1 / r13t).e0)
        e3 = EtaT((-xk).ef, (r02t / r13t).ef, (r02t / r13t).e0) + Eta(r02t.ef, (1 / r13t).ef)
        e4 = Eta(r02t.ef, (1 / r13t).ef) * EtaT((-xk).ef, (-r02t / r13t).ef, (-r02t / r13t).e0)
        total += sg * (e1 * (jlog(r02 * xk) + L1) + e2 * (jlog(xk / r13) + L2)
                       - e3 * (jlog(r02 * xk / r13) + L3) + e4)
    return total / (M[0] * M[1] * M[2] * M[3] * a.e0 * (xs[0].e0 - xs[1].e0)), conds


def _d0_dd_sym(ctx, P, M, nz):
    """Denner-Dittmaier sec. 3.3 with symbols; lines arranged as in scalar._d0_dd."""
    pp = lambda i, j: P[(min(i, j), max(i, j))]
    M2 = [m ** 2 for m in M]
    Ms = [J.inf(M2[i], 0, 0 if _zero(M[i]) else -I) for i in range(4)]     # m^2 - i del
    yr = lambda i, j: M2[i] + M2[j] - pp(i, j)
    y = lambda i, j: Ms[i] + Ms[j] - pp(i, j) + J.inf(0, 0, -I)

    def rroots(i, j):
        Yj = y(i, j)
        D = jsqrt(Yj * Yj - 4 * Ms[i] * Ms[j])
        return [ctx.name('r%d%d%s' % (i, j, sfx), r'r_{%d%d,%d}' % (i, j, n + 1), (Yj + sg * D) / (2 * Ms[i]))
                for n, (sg, sfx) in enumerate(((1, 'a'), (-1, 'b')))]

    def named(a, b, c, d):
        return [ctx.name(nm, nm, v).e0 for nm, v in (('a', a), ('b', b), ('c', c), ('d', d))]

    total = SR(0)
    if nz == 1:
        r13 = rroots(1, 3)[0]
        rr = r13.e0
        r03, r01 = rroots(0, 3), rroots(0, 1)
        m0, m1, m3 = M2[0], M2[1], M2[3]
        a, b, c, d = named(m1 * rr * yr(2, 3) - m3 * yr(1, 2),
                           yr(0, 2) * (m1 * rr - m3 / rr) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2),
                           yr(0, 2) * (yr(0, 1) - yr(0, 3) / rr) - m0 * (yr(1, 2) - yr(2, 3) / rr),
                           yr(1, 2) - yr(2, 3) / rr)
        xs, conds = _quad_j(ctx, a, b, c, d)
        for kk, xk in enumerate(xs):
            t = (jL2(-xk, y(2, 3) / y(0, 2)) - jL2(-xk, r03[0]) - jL2(-xk, r03[1])
                 - jL2(-xk * r13, y(1, 2) / y(0, 2)) + jL2(-xk * r13, r01[0]) + jL2(-xk * r13, r01[1]))
            Pn = Ms[0] + y(0, 3) * xk + Ms[3] * xk * xk + J.inf(0, -I)
            Qn = y(0, 2) + y(2, 3) * xk
            t += Eta((-xk).ef, r13.ef) * (jlog(Pn / Qn) + jlog(y(0, 2) / Ms[0]))
            total += (-1) ** kk * t
        return total / (a * (xs[0].e0 - xs[1].e0)), conds
    if nz in (2, 3):
        m0, m3 = M2[0], M2[3]
        a, b, c, d = named(yr(1, 3) * yr(2, 3) - m3 * yr(1, 2),
                           yr(0, 2) * yr(1, 3) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2),
                           yr(0, 1) * yr(0, 2) - m0 * yr(1, 2), yr(1, 2))
        xs, conds = _quad_j(ctx, a, b, c, d)
        r03 = rroots(0, 3) if nz == 2 else None
        if nz == 2 and _zero(yr(0, 1)) and _zero(yr(2, 3)):                 # DD (3.80)
            for kk, xk in enumerate(xs):
                lx = jlog(-xk)
                t = (-jL2(-xk, r03[0]) - jL2(-xk, r03[1]) - lx * lx / 2
                     - lx * (jlog(y(1, 3) / y(1, 2)) + jlog(y(0, 2) / Ms[0])))
                total += (-1) ** kk * t
            return total / (a * (xs[0].e0 - xs[1].e0)), conds
        for kk, xk in enumerate(xs):
            t = jL2(-xk, y(2, 3) / y(0, 2))
            t -= (jL2(-xk, r03[0]) + jL2(-xk, r03[1])) if nz == 2 else jL2(-xk, y(0, 3) / Ms[0])
            t += jL2(-xk, y(1, 3) / y(0, 1)) - jlog(-xk) * (jlog(y(0, 1) / y(1, 2)) + jlog(y(0, 2) / Ms[0]))
            total += (-1) ** kk * t
        return total / (a * (xs[0].e0 - xs[1].e0)), conds
    a, b, c, d = named(yr(1, 3) * yr(2, 3), yr(0, 2) * yr(1, 3) + yr(0, 1) * yr(2, 3) - yr(0, 3) * yr(1, 2),
                       yr(0, 1) * yr(0, 2), yr(1, 2))
    xs, conds = _quad_j(ctx, a, b, c, d)
    for kk, xk in enumerate(xs):
        t = jL2(-xk, y(2, 3) / y(0, 2)) - jL2(-1 / xk, y(0, 1) / y(1, 3))
        t -= jlog(-xk) * (jlog(y(1, 3) / y(1, 2)) + jlog(y(0, 2) / y(0, 3)))
        total += (-1) ** kk * t
    return total / (a * (xs[0].e0 - xs[1].e0)), conds


def d0_expand(s1, s2, s3, s4, s12, s23, m0, m1, m2, m3, pair=(1, 3)):
    r"""
    D0(s1, s2, s3, s4; s12, s23; m0, m1, m2, m3) as an explicit formula in the arguments
    (Package-X's D0Expand).  Arguments in Package-X order, masses not squared.  Returns a
    Conditional.

    All masses nonzero: Denner's sixteen dilogarithms, eq. (4.43), which need one pair of lines
    (i, j) with r_ij real, p_ij^2 <= (m_i - m_j)^2 or p_ij^2 >= (m_i + m_j)^2.  pair chooses it:
    the default (1, 3) is the diagonal s23, the condition of Package-X's D0Expand; (0, 2) is s12.
    Some masses zero: Denner and Dittmaier's formulas, valid for all real kinematics (the
    conditions only exclude singular points).  IR-divergent boxes are in feynsage.ir.
    """
    from .scalar import _DD_PATTERN, _DD_NEEDS
    args = [sr(x) for x in (s1, s2, s3, s4, s12, s23, m0, m1, m2, m3)]
    inv = dict(zip(((0, 1), (1, 2), (2, 3), (0, 3), (0, 2), (1, 3)), args[:6]))
    M = args[6:]
    pp = lambda i, j: inv[(min(i, j), max(i, j))]
    zero = [_zero(m) for m in M]
    nz = sum(zero)
    errors = []
    if nz == 0:
        i0, j0 = pair
        for pm in [pm for pm in itertools.permutations(range(4)) if {pm[0], pm[2]} == {i0, j0}]:
            Pq = {(i, j): pp(pm[i], pm[j]) for i in range(4) for j in range(i + 1, 4)}
            ctx = _Ctx('D')
            try:
                v, conds = _d0_denner_sym(ctx, Pq, [M[q] for q in pm])
                return Conditional(v, conds, ctx.defs)
            except (ZeroDivisionError, ArithmeticError) as e:
                errors.append(e)
    else:
        for pm in itertools.permutations(range(4)):
            if tuple(0 if zero[q] else 1 for q in pm) != _DD_PATTERN[nz]:
                continue
            Mp = [M[q] for q in pm]
            Pq = {(i, j): pp(pm[i], pm[j]) for i in range(4) for j in range(i + 1, 4)}
            Yr = lambda i, j: Mp[i] ** 2 + Mp[j] ** 2 - Pq[(min(i, j), max(i, j))]
            needs = _DD_NEEDS[nz]
            if nz == 2 and _zero(Yr(0, 1)) and _zero(Yr(2, 3)):
                needs = [(0, 2), (1, 2), (1, 3)]
            if any(_zero(Yr(i, j)) for i, j in needs):
                continue
            ctx = _Ctx('D')
            try:
                v, conds = _d0_dd_sym(ctx, Pq, Mp, nz)
                return Conditional(v, [Yr(i, j) != 0 for i, j in needs] + conds, ctx.defs)
            except (ZeroDivisionError, ArithmeticError) as e:
                errors.append(e)
    raise ZeroDivisionError("no labelling of the lines gives a regular closed form here: %s" % errors)


def expand_c0d0(expr, pair=(1, 3)):
    r"""
    Replace every C0(...) and D0(...) inside expr (for example a loop() coefficient) by its closed
    form (Package-X's C0Expand and D0Expand together).  Returns a Conditional that collects the
    definitions and conditions of every replaced function.
    """
    from .pv import _C0f, _D0f
    expr = SR(expr)
    conds, defs = [], []

    def c0(*a):
        r = c0_expand(*a)
        conds.extend(r.conditions)
        defs.extend(r.defs)
        return r.value

    def d0(*a):
        r = d0_expand(*a, pair=pair)
        conds.extend(r.conditions)
        defs.extend(r.defs)
        return r.value

    out = expr.substitute_function(_C0f, c0).substitute_function(_D0f, d0)
    return Conditional(out, conds, defs)
