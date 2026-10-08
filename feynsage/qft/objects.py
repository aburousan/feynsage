r"""
Dirac matrices and spinors as ordered (non-commuting) Python objects.

    gamma(mu), slash(p), gamma5(), PL(), PR(), sigma(mu, nu), one()
    u(p, m), v(p, m), ubar(p, m), vbar(p, m)

They multiply with * in the order written and mix with ordinary Sage scalars and tensors:
    (slash(p) + m) * gamma(mu) * (slash(k) + m)       one product of three factors
    ubar(p2, m) * gamma(mu) * u(p1, m)                 a spinor bilinear (a number)
    J1 * metric(mu, nu) * J2                           two bilinears times a tensor

Internally an expression is a sum of terms.  A term is
    Sage coefficient  x  [one open chain]  x  [closed chains]
A chain is (left spinor or None, factors, right spinor or None) and every factor is a linear
combination of single matrices 1, gamma5, gamma^mu, slash(p) with Sage coefficients, so slash(p) + m
stays one factor and a product of propagators is not expanded in Python.  A chain with spinors at
both ends is a number (a bilinear): it commutes with everything and is kept in the closed list.
"""
from sage.all import SR, Integer, latex
from sage.structure.element import Element
from . import tensors as T

# ---------------------------------------------------------------------------- one factor
# atoms: ('1',) unit matrix, ('5',) gamma5, ('i', mu) gamma^mu, ('p', p) slash(p); mu, p Sage symbols


def _akey(a):
    return (a[0], str(a[1]) if len(a) > 1 else '')


class Factor:
    """A linear combination sum_a c_a A_a of single Dirac matrices A_a (immutable)."""
    __slots__ = ('items', '_key')

    def __init__(self, items):
        acc = {}
        for a, c in items:
            k = _akey(a)
            if k in acc:
                acc[k] = (a, acc[k][1] + c)
            else:
                acc[k] = (a, SR(c))
        self.items = tuple(sorted(((a, c) for a, c in acc.values() if not c.is_trivial_zero()),
                                  key=lambda t: _akey(t[0])))
        self._key = None

    def key(self):
        if self._key is None:
            self._key = tuple((_akey(a), str(c)) for a, c in self.items)
        return self._key

    def is_zero(self):
        return not self.items

    def scale(self, c):
        return Factor([(a, c * x) for a, x in self.items])

    def odd_even(self):
        """(even part, odd part): even = 1 and gamma5, odd = gamma^mu and slash(p)."""
        ev = Factor([(a, c) for a, c in self.items if a[0] in ('1', '5')])
        od = Factor([(a, c) for a, c in self.items if a[0] in ('i', 'p')])
        return ev, od

    def flip(self):
        """F' with gamma5 F = F' gamma5 (anticommuting gamma5): even part minus odd part."""
        return Factor([(a, -c if a[0] in ('i', 'p') else c) for a, c in self.items])

    def bar(self, conj):
        """gamma^0 F^dagger gamma^0: gamma^mu, slash p, 1 unchanged, gamma5 -> -gamma5, c -> c*."""
        return Factor([(a, -conj(c) if a[0] == '5' else conj(c)) for a, c in self.items])

    def gamma5_coefficient(self):
        for a, c in self.items:
            if a[0] == '5':
                return c
        return SR(0)

    def without_gamma5(self):
        return Factor([(a, c) for a, c in self.items if a[0] != '5'])

    def _shown(self):
        """Display order of textbooks: (p/ + m), (1 - gamma5)/2."""
        return sorted(self.items, key=lambda t: {'p': 0, 'i': 1, '1': 2, '5': 3}[t[0][0]])

    def __repr__(self):
        return _join([(c, _atom_repr(a)) for a, c in self._shown()])

    def _latex_(self):
        return _join([(c, _atom_latex(a)) for a, c in self._shown()], latex_mode=True)


def _atom_repr(a):
    return {'1': '1', '5': 'gamma5()'}.get(a[0]) or ('gamma(%s)' % a[1] if a[0] == 'i' else 'slash(%s)' % a[1])


def _atom_latex(a):
    if a[0] == '1':
        return '1'
    if a[0] == '5':
        return r'\gamma^{5}'
    if a[0] == 'i':
        return (r'\gamma_{%s}' if T._CONV['euclidean'] else r'\gamma^{%s}') % latex(a[1])
    return r'{\not{%s}}' % latex(a[1])


def _coef_str(c, latex_mode):
    if latex_mode:
        s = latex(c)
        return s if c.operator() is None or c.operator().__name__ not in ('add_vararg', 'add') else r'\left(%s\right)' % s
    s = str(c)
    return s if c.operator() is None or c.operator().__name__ not in ('add_vararg', 'add') else '(%s)' % s


def _join(pairs, latex_mode=False):
    if not pairs:
        return '0'
    parts = []
    for c, a in pairs:
        if a == '1':
            parts.append(_coef_str(c, latex_mode))
        elif (c - 1).is_trivial_zero():
            parts.append(a)
        elif (c + 1).is_trivial_zero():
            parts.append('-' + a)
        else:
            parts.append(_coef_str(c, latex_mode) + (' ' if latex_mode else '*') + a)
    s = ' + '.join(parts).replace('+ -', '- ')
    return s


# ---------------------------------------------------------------------------- spinors
class Spinor:
    """u(p, m), v(p, m) (columns, at the right end) or ubar(p, m), vbar(p, m) (rows, at the left end)."""
    __slots__ = ('kind', 'p', 'm')

    def __init__(self, kind, p, m=0):
        self.kind, self.p, self.m = kind, SR(p), SR(m)
        if not T.is_momentum(self.p):
            raise TypeError("%s(p, m): p must be one momentum made with momenta(...), got %s" % (kind, p))

    @property
    def row(self):
        return self.kind in ('ubar', 'vbar')

    @property
    def particle(self):
        return self.kind[0]          # 'u' or 'v'

    def key(self):
        return (self.kind, str(self.p), str(self.m))

    def partner_key(self):
        """The key of the spinor this one is summed with: u(p) with ubar(p), v(p) with vbar(p)."""
        k = {'u': 'ubar', 'ubar': 'u', 'v': 'vbar', 'vbar': 'v'}[self.kind]
        return (k, str(self.p), str(self.m))

    def dagger_bar(self):
        """The spinor at the same end after complex conjugation of a bilinear: u <-> ubar, v <-> vbar."""
        return Spinor({'u': 'ubar', 'ubar': 'u', 'v': 'vbar', 'vbar': 'v'}[self.kind], self.p, self.m)

    def __repr__(self):
        return '%s(%s%s)' % (self.kind, self.p, '' if self.m.is_trivial_zero() else ', %s' % self.m)

    def _latex_(self):
        base = r'\bar{%s}' % self.kind[0] if self.row else self.kind
        return r'%s(%s)' % (base, latex(self.p))


class Chain:
    """(left spinor | None) factor_1 ... factor_n (right spinor | None)."""
    __slots__ = ('left', 'factors', 'right')

    def __init__(self, left, factors, right):
        self.left, self.factors, self.right = left, tuple(factors), right

    @property
    def closed(self):
        return self.left is not None and self.right is not None

    def key(self):
        return (self.left.key() if self.left else None, tuple(f.key() for f in self.factors),
                self.right.key() if self.right else None)

    def __repr__(self):
        parts = ([repr(self.left)] if self.left else []) + ['(%s)' % f if len(f.items) > 1 else repr(f)
                                                             for f in self.factors] + ([repr(self.right)] if self.right else [])
        return '*'.join(parts) if parts else '1'

    def _latex_(self):
        parts = ([self.left._latex_()] if self.left else []) + [r'\left(%s\right)' % f._latex_() if len(f.items) > 1 else f._latex_()
                                                                 for f in self.factors] + ([self.right._latex_()] if self.right else [])
        s = ' '.join(parts) if parts else '1'
        return r'\left[%s\right]' % s if self.closed else s


# ---------------------------------------------------------------------------- the expression
class DiracError(TypeError):
    pass


def _is_scalar(x):
    if isinstance(x, (int, float, complex)):
        return True
    if isinstance(x, Element):
        return True
    return False


class DiracExpr:
    """A sum of terms (coefficient, open chain or None, closed chains).  Build it with gamma(), slash(),
    u(), ... and *, +, -; take its trace with dirac_trace(expr) or expr.trace()."""
    __array_priority__ = 1000

    def __init__(self, terms):
        self.terms = [t for t in terms if not SR(t[0]).is_trivial_zero()]

    # -- kind: 'scalar', 'matrix', 'row', 'column'
    @staticmethod
    def _tkind(t):
        o = t[1]
        if o is None:
            return 'scalar'
        if o.left is None and o.right is None:
            return 'matrix'
        return 'row' if o.left is not None else 'column'

    def kind(self):
        ks = {self._tkind(t) for t in self.terms}
        if len(ks) > 1:
            raise DiracError("mixed expression: %s" % sorted(ks))
        return ks.pop() if ks else None

    # -- arithmetic
    def _scalar_as(self, c, kind):
        c = SR(c)
        if kind == 'matrix':
            return DiracExpr([(c, Chain(None, (), None), ())])
        if kind in ('scalar', None):
            return DiracExpr([(c, None, ())])
        raise DiracError("cannot add the number %s to a %s spinor expression" % (c, kind))

    def _simple_factor(self):
        """If every term is c * (one factor or the unit matrix), the merged Factor; else None."""
        items = []
        for c, o, cl in self.terms:
            if cl or o is None or o.left is not None or o.right is not None or len(o.factors) > 1:
                return None
            if o.factors:
                items.extend((a, c * x) for a, x in o.factors[0].items)
            else:
                items.append((('1',), c))
        return Factor(items)

    def __add__(self, other):
        if isinstance(other, DiracExpr):
            ka, kb = self.kind(), other.kind()
            if ka and kb and ka != kb:
                raise DiracError("cannot add a %s and a %s (for example a Dirac matrix and a spinor bilinear)" % (ka, kb))
            fa, fb = self._simple_factor(), other._simple_factor()
            if fa is not None and fb is not None:
                f = Factor(list(fa.items) + list(fb.items))
                return DiracExpr([(SR(1), Chain(None, (f,), None), ())]) if not f.is_zero() else DiracExpr([])
            return DiracExpr(_combine(self.terms + other.terms))
        if _is_scalar(other):
            return self + self._scalar_as(other, self.kind())
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return DiracExpr([(-c, o, cl) for c, o, cl in self.terms])

    def __sub__(self, other):
        if isinstance(other, DiracExpr) or _is_scalar(other):
            return self + (-other)
        return NotImplemented

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if isinstance(other, DiracExpr):
            out = []
            for ta in self.terms:
                for tb in other.terms:
                    out.append(_term_mul(ta, tb))
            return DiracExpr(_combine(out))
        if _is_scalar(other):
            c = SR(other)
            return DiracExpr([(c * x, o, cl) for x, o, cl in self.terms])
        return NotImplemented

    def __rmul__(self, other):
        if _is_scalar(other):
            c = SR(other)
            return DiracExpr([(c * x, o, cl) for x, o, cl in self.terms])
        return NotImplemented

    def __truediv__(self, other):
        if _is_scalar(other):
            return self * (1 / SR(other))
        return NotImplemented

    def __pow__(self, n):
        n = int(n)
        if n < 0:
            raise DiracError("negative powers of Dirac matrices are not supported; write the propagator numerator and denominator separately")
        out = one() if self.kind() == 'matrix' else DiracExpr([(SR(1), None, ())])
        for _ in range(n):
            out = out * self
        return out

    def __matmul__(self, other):
        return self * other

    # -- operations
    def trace(self, **kw):
        """Dirac trace (see dirac_trace)."""
        from .operations import dirac_trace
        return dirac_trace(self, **kw)

    def conjugate(self):
        from .operations import conjugate
        return conjugate(self)

    def bar(self):
        from .operations import dirac_bar
        return dirac_bar(self)

    def is_zero(self):
        return not self.terms

    def __len__(self):
        return len(self.terms)

    def __iter__(self):
        return iter(self.terms)

    # -- display
    def _term_str(self, t, latex_mode=False):
        c, o, cl = t
        pieces = [x._latex_() if latex_mode else repr(x) for x in cl]
        if o is not None:
            s = o._latex_() if latex_mode else repr(o)
            if s != '1' or not pieces:
                pieces.append(s)
        body = (' ' if latex_mode else '*').join(p for p in pieces if p != '1') or '1'
        if body == '1':
            return _coef_str(c, latex_mode)
        if (c - 1).is_trivial_zero():
            return body
        if (c + 1).is_trivial_zero():
            return '-' + body
        return _coef_str(c, latex_mode) + (' ' if latex_mode else '*') + body

    def __repr__(self):
        if not self.terms:
            return '0'
        return ' + '.join(self._term_str(t) for t in self.terms).replace('+ -', '- ')

    def _latex_(self):
        if not self.terms:
            return '0'
        return ' + '.join(self._term_str(t, True) for t in self.terms).replace('+ -', '- ')

    def _repr_latex_(self):
        return '$' + self._latex_() + '$'


def _combine(terms):
    """Add the coefficients of identical (open chain, closed chains)."""
    acc, order = {}, []
    for c, o, cl in terms:
        k = (o.key() if o is not None else None, tuple(sorted(x.key() for x in cl)))
        if k in acc:
            acc[k][0] += c
        else:
            acc[k] = [SR(c), o, tuple(cl)]                 # kept in the order written; the key is sorted
            order.append(k)
    return [tuple(acc[k]) for k in order]


def _concat(a, b):
    if a.right is not None:
        raise DiracError("%r already ends with the spinor %r; nothing can follow it in the same line" % (a, a.right))
    if b.left is not None:
        raise DiracError("%r is a row spinor; it must stand at the left end of a line" % b.left)
    return Chain(a.left, a.factors + b.factors, b.right)


def _term_mul(ta, tb):
    ca, oa, cla = ta
    cb, ob, clb = tb
    closed = list(cla) + list(clb)
    if oa is None:
        o = ob
    elif ob is None:
        o = oa
    else:
        if oa.right is not None and oa.left is None and ob.left is not None:
            raise DiracError("%r * %r is an outer product of spinors; write the spin sum with spin_sum()" % (oa.right, ob.left))
        o = _concat(oa, ob)
    if o is not None and o.closed:
        closed.append(o)
        o = None
    return (ca * cb, o, tuple(closed))


# ---------------------------------------------------------------------------- constructors
def _single(atom, c=1):
    return DiracExpr([(SR(1), Chain(None, (Factor([(atom, SR(c))]),), None), ())])


def one():
    """The 4x4 unit matrix in Dirac space."""
    return DiracExpr([(SR(1), Chain(None, (), None), ())])


def gamma(*args, **kw):
    """gamma(mu): the Dirac matrix gamma^mu for a Lorentz index mu (from lorentz_indices).
    Called with anything else it is Sage's Euler Gamma function, so gamma(5) = 24 still works."""
    if len(args) == 1 and not kw and not isinstance(args[0], DiracExpr):
        x = args[0]
        if isinstance(x, str) and x in T._KIND:
            x = SR.var(x)
        try:
            xs = SR(x)
        except (TypeError, ValueError):
            xs = None
        if xs is not None and T.is_index(xs):
            return _single(('i', xs))
        if xs is not None and T.is_momentum(xs):
            raise DiracError("gamma(%s): %s is a momentum; write slash(%s) for gamma^mu %s_mu" % (xs, xs, xs, xs))
    from sage.all import gamma as _euler_gamma
    return _euler_gamma(*args, **kw)


def slash(p):
    """slash(p) = gamma^mu p_mu for a momentum or a combination such as p - k."""
    parts = T.vector_parts(p, "slash(p) needs a momentum")
    if not parts:
        return DiracExpr([])
    return DiracExpr([(SR(1), Chain(None, (Factor([(('p', v), c) for v, c in parts]),), None), ())])


def gamma5():
    """gamma5 (= i g^0 g^1 g^2 g^3 in Minkowski space, g_1 g_2 g_3 g_4 in Euclidean space)."""
    return _single(('5',))


def PL():
    """The chiral projector (1 - gamma5)/2."""
    return DiracExpr([(SR(1), Chain(None, (Factor([(('1',), Integer(1) / 2), (('5',), -Integer(1) / 2)]),), None), ())])


def PR():
    """The chiral projector (1 + gamma5)/2."""
    return DiracExpr([(SR(1), Chain(None, (Factor([(('1',), Integer(1) / 2), (('5',), Integer(1) / 2)]),), None), ())])


def chiral(expr, p, hand):
    r"""
    Keep one chirality of the external spinor with momentum p: u(p) -> P_hand u(p) (and v(p) the same),
    ubar(p) -> ubar(p) P_other.  For a massless fermion this picks one helicity, so that
        spin_sum(chiral(M, p1, "L") * conjugate(chiral(M, p1, "L")))
    is |M|^2 for a left-handed electron with momentum p1 (summed over the other spins).  For a massless
    antifermion v(p) with P_L v = v the helicity is +1/2 (right-handed positron).
    """
    if hand not in ("L", "R"):
        raise ValueError('hand must be "L" or "R"')
    pr = (PL() if hand == "L" else PR()).terms[0][1].factors
    pl = (PR() if hand == "L" else PL()).terms[0][1].factors
    key = str(SR(p))

    def fix(ch):
        if ch is None:
            return ch
        f = ch.factors
        if ch.right is not None and str(ch.right.p) == key:
            f = f + pr
        if ch.left is not None and str(ch.left.p) == key:
            f = pl + f
        return Chain(ch.left, f, ch.right)
    return DiracExpr([(c, fix(o), tuple(fix(x) for x in cl)) for c, o, cl in expr.terms])


def sigma(*args, **kw):
    """sigma(mu, nu) = (i/2) [gamma^mu, gamma^nu] for Lorentz indices or momenta (sigma(mu, p) =
    sigma^{mu nu} p_nu).  With other arguments it is Sage's divisor function sigma(n, k)."""
    from sage.all import I
    if len(args) != 2 or kw or not all(_is_slot(x) for x in args):
        from sage.all import sigma as _sage_sigma
        return _sage_sigma(*args, **kw)
    mu, nu = args
    return (I / 2) * (_g_or_slash(mu) * _g_or_slash(nu) - _g_or_slash(nu) * _g_or_slash(mu))


def _is_slot(x):
    try:
        x = SR(x)
    except (TypeError, ValueError):
        return False
    if T.is_index(x) or T.is_momentum(x):
        return True
    try:
        return bool(T.vector_parts(x))
    except TypeError:
        return False


def _g_or_slash(x):
    return gamma(x) if T.is_index(SR(x)) else slash(x)


def _sp(kind, p, m):
    s = Spinor(kind, p, m)
    if s.row:
        return DiracExpr([(SR(1), Chain(s, (), None), ())])
    return DiracExpr([(SR(1), Chain(None, (), s), ())])


def u(p, m=0):
    """The particle spinor u(p) with mass m (a column: it stands at the right end of a line)."""
    return _sp('u', p, m)


def v(p, m=0):
    """The antiparticle spinor v(p) with mass m (a column)."""
    return _sp('v', p, m)


def ubar(p, m=0):
    """u-bar(p) = u(p)^dagger gamma^0 (a row: it stands at the left end of a line)."""
    return _sp('ubar', p, m)


def vbar(p, m=0):
    """v-bar(p) = v(p)^dagger gamma^0 (a row)."""
    return _sp('vbar', p, m)
