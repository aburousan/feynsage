r"""
Feynman rules from a Lagrangian (FeynCalc's FeynRule, FieldStrength, CovariantD).

    phi = field("phi")                                  # real scalar
    psi = field("psi", "dirac");  A = field("A", "vector")
    L = -lam/24 * phi()**4 - e * psi.bar() * gamma(mu) * psi() * A(mu)
    feynman_rule(L, [(A, p1, mu1), (psi, p2), (psi.bar, p3)])   # -> I e^... gamma^mu1 ... (all momenta incoming)

The vertex is i times the derivative of the action with respect to the fields of the legs, in momentum space
with every momentum incoming: a derivative d_mu acting on a field whose leg carries momentum p becomes -i p_mu.
Identical fields are summed over all their assignments to the legs.  Indices that the Lagrangian sums over
are renamed to the indices of the legs (a metric or colour delta is written out when a summed index sits on two
fields).  For fermions the vertex is the Dirac matrix between psibar and psi (bilinears only).

Fields: field(name, kind="scalar" | "complex" | "dirac" | "vector" | "ghost", colour=None | "adjoint" |
"fundamental").  A field is called with its indices in the order Lorentz then colour: A(mu, a), psi(i),
phi().  A complex, Dirac or ghost field has a partner .bar (psibar, phi^dagger, the antighost).  d(mu, X)
is the derivative of a field or of a product (Leibniz rule).  field_strength and covariant_derivative build
the non-abelian ones with FeynCalc's conventions.
"""
from itertools import permutations, product

from sage.all import SR, I

from .qft import tensors as T
from .qft.objects import DiracExpr


_counter = [0]


def _fresh(kind):
    _counter[0] += 1
    name = "fr%s%d" % ("c" if kind == "colour" else "l", _counter[0])
    if kind == "colour":
        from .qft.color import color_indices
        return SR(color_indices(name))
    return SR(T.lorentz_indices(name))


def _fresh_quark():
    from .qft.color import quark_colors
    _counter[0] += 1
    return SR(quark_colors("frq%d" % _counter[0]))


class Field:
    """A field of a Lagrangian (made by field()); call it with its indices to get a term, A(mu, a)."""

    def __init__(self, name, kind="scalar", colour=None, partner=None, is_bar=False):
        self.name, self.kind, self.colour, self.is_bar = name, kind, colour, is_bar
        self._partner = partner

    @property
    def bar(self):
        """psibar for psi, phi^dagger for a complex phi, the antighost for a ghost (and back)."""
        if self.kind in ("scalar", "vector"):
            return self
        if self._partner is None:
            self._partner = Field(self.name, self.kind, self.colour, self, not self.is_bar)
        return self._partner

    @property
    def key(self):
        return (self.name, self.is_bar)

    def __call__(self, *idx):
        lor = None
        col = None
        idx = list(idx)
        if self.kind == "vector":
            if not idx:
                raise ValueError("%s is a vector field: give its Lorentz index, %s(mu)" % (self.name, self.name))
            lor = SR(idx.pop(0))
        if self.colour is not None:
            if not idx:
                raise ValueError("%s has a colour index: %s(..., a)" % (self.name, self.name))
            col = SR(idx.pop(0))
        if idx:
            raise ValueError("too many indices for %s" % self.name)
        return Lag([(SR(1), None, (_Factor(self, lor, col, ()),))])

    def __repr__(self):
        return self.name + ("bar" if self.is_bar and self.kind == "dirac" else "^dagger" if self.is_bar else "")


def field(name, kind="scalar", colour=None):
    """A field: kind "scalar" (real), "complex", "dirac", "vector" or "ghost"; colour None, "adjoint" or
    "fundamental".  Indices when called: Lorentz (vectors) then colour."""
    if kind not in ("scalar", "complex", "dirac", "vector", "ghost"):
        raise ValueError("unknown kind of field %r" % kind)
    return Field(name, kind, colour)


class _Factor:
    def __init__(self, f, lor, col, derivs):
        self.f, self.lor, self.col, self.derivs = f, lor, col, tuple(derivs)

    def renamed(self, sub):
        s = lambda x: None if x is None else SR(x).subs(sub)
        return _Factor(self.f, s(self.lor), s(self.col), tuple(SR(m).subs(sub) for m in self.derivs))


class Lag:
    """A sum of terms coefficient * (Dirac matrix) * fields; built from fields with +, -, * and d()."""

    def __init__(self, terms):
        self.terms = [t for t in terms if not (isinstance(t[0], (int,)) and t[0] == 0)]

    def __add__(self, other):
        return Lag(self.terms + _lag(other).terms)

    __radd__ = __add__

    def __neg__(self):
        return Lag([(-c, D, F) for c, D, F in self.terms])

    def __sub__(self, other):
        return self + (-_lag(other))

    def __rsub__(self, other):
        return _lag(other) - self

    def __mul__(self, other):
        if isinstance(other, DiracExpr):
            return Lag([(c, other if D is None else D * other, F) for c, D, F in self.terms])
        o = _lag(other)
        out = []
        for c1, D1, F1 in self.terms:
            for c2, D2, F2 in o.terms:
                D = D1 if D2 is None else (D2 if D1 is None else D1 * D2)
                out.append((c1 * c2, D, F1 + F2))
        return Lag(out)

    def __rmul__(self, other):
        if isinstance(other, DiracExpr):
            return Lag([(c, other if D is None else other * D, F) for c, D, F in self.terms])
        return _lag(other) * self

    def __truediv__(self, x):
        return Lag([(c / x, D, F) for c, D, F in self.terms])

    def __pow__(self, n):
        out = Lag([(SR(1), None, ())])
        for _ in range(int(n)):
            out = out * self
        return out

    def __repr__(self):
        parts = []
        for c, D, F in self.terms:
            fs = " ".join("%s%s%s" % ("".join("d%s " % m for m in f.derivs), f.f, "" if f.lor is None and f.col is None
                                          else "(%s)" % ", ".join(str(x) for x in (f.lor, f.col) if x is not None)) for f in F)
            parts.append("(%s)%s %s" % (c, "" if D is None else " " + repr(D), fs))
        return " + ".join(parts) or "0"


def _lag(x):
    if isinstance(x, Lag):
        return x
    return Lag([(SR(x), None, ())])


def d(mu, X):
    """The derivative d_mu of a field or of a product of fields (Leibniz rule)."""
    X = _lag(X)
    out = []
    for c, D, F in X.terms:
        for k, f in enumerate(F):
            G = list(F)
            G[k] = _Factor(f.f, f.lor, f.col, f.derivs + (SR(mu),))
            out.append((c, D, tuple(G)))
    return Lag(out)


def field_strength(A, mu, nu, a, g):
    """F^a_{mu nu} = d_mu A^a_nu - d_nu A^a_mu + g f^{abc} A^b_mu A^c_nu (FeynCalc's FieldStrength); for an
    abelian field (colour None) without the last term."""
    from .qft.color import f_color
    F = d(mu, A(nu, a)) - d(nu, A(mu, a)) if A.colour else d(mu, A(nu)) - d(nu, A(mu))
    if A.colour:
        b, c = _fresh("colour"), _fresh("colour")
        F = F + g * f_color(a, b, c) * A(mu, b) * A(nu, c)
    return F


def covariant_derivative(mu, psi, A, g, i=None):
    """D_mu psi_i = d_mu psi_i - i g (T^a)_{ij} A^a_mu psi_j for a fundamental psi and a gluon A (FeynCalc's
    CovariantD); for a field without colour and an abelian A: d_mu psi + i g A_mu psi (the QED sign of the
    electron with charge g)."""
    from .qft.color import T_color
    if psi.colour == "fundamental":
        j, a = _fresh_quark(), _fresh("colour")
        return d(mu, psi(i)) - I * g * T_color(a, i, j) * A(mu, a) * psi(j)
    return d(mu, psi()) + I * g * A(mu) * psi()


def _fresh_term(c, D, F):
    """Every index the term sums over gets a fresh name, so that it never clashes with an index of a leg."""
    from .qft.color import is_adjoint, is_fundamental
    syms = set()
    for f in F:
        for x in (f.lor, f.col) + f.derivs:
            if x is not None:
                syms.add(SR(x))
    for v in SR(c).variables():
        if T.is_index(v) or is_adjoint(v) or is_fundamental(v):
            syms.add(v)
    sub = {}
    for x in syms:
        sub[x] = _fresh_quark() if is_fundamental(x) else _fresh("colour") if is_adjoint(x) else _fresh("lorentz")
    D2 = None if D is None else _rename_dirac(D, sub)
    return SR(c).subs(sub), D2, tuple(f.renamed(sub) for f in F)


def _split_shared(c, D, F):
    """An index that sits on two fields gets a fresh name on the second one, with a metric or delta."""
    from .qft.color import color_delta
    seen = {}
    G = []
    for f in F:
        sub = {}
        for x, kind in [(f.lor, "lorentz"), (f.col, "colour")] + [(m, "lorentz") for m in f.derivs]:
            if x is None:
                continue
            if str(x) in seen and seen[str(x)] is not f:
                if kind == "colour" and f.f.colour == "fundamental":
                    y = _fresh_quark()
                else:
                    y = _fresh(kind)
                c = c * (T.metric(x, y) if kind == "lorentz" else color_delta(x, y))
                sub[x] = y
            else:
                seen[str(x)] = f
        G.append(f.renamed(sub) if sub else f)
    return c, D, tuple(G)


def _rename_dirac(D, sub):
    from .qft.operations import _rename_factor
    from .qft.objects import Chain
    terms = []
    for c, o, cl in D.terms:
        o2 = None if o is None else Chain(o.left, tuple(_rename_factor(F, sub) for F in o.factors), o.right)
        terms.append((SR(c).subs(sub), o2, cl))
    return DiracExpr(terms)


def feynman_rule(L, legs, factor_i=True):
    """
    The vertex of the Lagrangian L for the legs [(field, momentum, Lorentz index, colour index), ...] (indices
    only as the field has them; field.bar for psibar, phi^dagger or an antighost; all momenta incoming).
    Returns i times the sum over the terms of L with exactly these fields; a Dirac matrix (DiracExpr) when
    the legs are a fermion and an antifermion, else a Sage expression.  factor_i=False leaves out the i.
    FeynCalc's FeynRule.
    """
    want = sorted([(lg[0].name, lg[0].is_bar) for lg in legs])
    total = SR(0)
    dirac_total = None
    for term in _lag(L).terms:
        c, D, F = _split_shared(*_fresh_term(*term))
        if sorted(f.f.key for f in F) != want:
            continue
        slots = {}
        for k, f in enumerate(F):
            slots.setdefault(f.f.key, []).append(k)
        legs_by = {}
        for n, lg in enumerate(legs):
            legs_by.setdefault((lg[0].name, lg[0].is_bar), []).append(n)
        keys = sorted(slots)
        for choice in product(*[permutations(legs_by[k]) for k in keys]):
            assign = {}
            for k, perm in zip(keys, choice):
                for fpos, legn in zip(slots[k], perm):
                    assign[fpos] = legn
            sub, w = {}, c
            for fpos, f in enumerate(F):
                lg = legs[assign[fpos]]
                p = SR(lg[1])
                rest = list(lg[2:])
                if f.lor is not None:
                    sub[f.lor] = SR(rest.pop(0))
                if f.col is not None:
                    sub[f.col] = SR(rest.pop(0))
                for m in f.derivs:
                    w = w * (-I) * T.comp(p, m)
            w = SR(w).subs(sub)
            if D is None:
                total += w
            else:
                Dr = _rename_dirac(D, sub) * w
                dirac_total = Dr if dirac_total is None else dirac_total + Dr
    from .qft.operations import contract, dummy_indices
    pref = I if factor_i else 1
    if dirac_total is not None:
        return contract(dirac_total * pref)
    total = (pref * total).expand()
    if dummy_indices(total):
        total = _unhide(contract(_hide(total)))      # Lorentz contraction only: f^{abc} stays as it is
    return _canon_sym(_canon_f(total)).expand()


_HID = None


def _hide(e):
    """f_color -> an opaque function during the Lorentz contraction (contract would expand it in traces)."""
    global _HID
    from sage.all import function
    from .qft.color import COLF
    if _HID is None:
        _HID = function("fsHiddenF", nargs=3)
    w = [SR.wild(k) for k in range(3)]
    return SR(e).substitute_function(COLF, _HID)


def _unhide(e):
    from .qft.color import COLF
    return SR(e).substitute_function(_HID, COLF) if _HID is not None else SR(e)


def _canon_sym(e):
    """Metric(nu, mu) -> Metric(mu, nu) and Dot(q, p) -> Dot(p, q): symmetric, so equal terms combine."""
    e = SR(e)
    sub = {}
    for op in {x.operator() for x in [SR(T.metric(*T.lorentz_indices("fsa fsb"))),
                                      SR(T.dot(*T.momenta("fsp fsq")))]}:
        for F in e.find(op(SR.wild(0), SR.wild(1))):
            a_, b_ = F.operands()
            if str(b_) < str(a_):
                sub[F] = op(b_, a_)
    return e.subs(sub) if sub else e


def _canon_f(e):
    """Each f^{abc} with its indices in a fixed order (f is totally antisymmetric), so equal terms combine."""
    from .qft.color import COLF
    e = SR(e)
    w = [SR.wild(k) for k in range(3)]
    sub = {}
    for F in e.find(COLF(*w)):
        args = list(F.operands())
        order = sorted(range(3), key=lambda k: str(args[k]))
        perm_sign = 1
        p = list(order)
        for x in range(3):                             # the sign of the sorting permutation
            for y in range(x + 1, 3):
                if p[x] > p[y]:
                    perm_sign = -perm_sign
        sub[F] = perm_sign * COLF(*[args[k] for k in order])
    return e.subs(sub) if sub else e
