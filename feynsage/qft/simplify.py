r"""
Simplify open Dirac lines and spinor bilinears (no trace).

    simplify_dirac(gamma(mu) * gamma(nu) * gamma(mu))          -> -2 gamma(nu)          (dim=4)
    simplify_dirac(metric(mu, nu) * gamma(mu) * gamma(nu))     -> 4
    simplify_dirac(ubar(p2, m) * slash(p2) * gamma(mu) * u(p1, m))  -> m ubar(p2) gamma(mu) u(p1)

Two steps.  FORM contracts the tensors into the gamma matrices (g^{mu nu} gamma_nu -> gamma^mu,
p_mu gamma^mu -> slash(p)).  Then the normal ordering of feynsage.dirac (Package-X's
FermionLineExpand) puts each line in a fixed order with {g^a, g^b} = 2 g^{ab}: repeated indices give
the dimension, repeated momenta their square, gamma5 is anticommuted to the right end (NDR, as in
feynsage.dirac), and at a spinor the Dirac equation removes its own momentum,
    slash(p) u(p) = m u(p),  slash(p) v(p) = -m v(p),  ubar(p) slash(p) = m ubar(p),  vbar(p) slash(p) = -m vbar(p).
"""
from sage.all import SR, Integer
from .. import dirac as _d
from . import tensors as T
from .objects import DiracExpr, Factor, Chain, _combine
from .compiler import Program
from .parser import Reader
from .operations import _execute, _check


class _Order(_d._Order):
    """The ordering of feynsage.dirac with the tensors of feynsage.qft."""

    def __init__(self, left, right, dimv, euclid, names):
        super().__init__(left=left, right=right)
        self.dimv, self.euclid, self.names = dimv, euclid, names

    def half_anti(self, a, b):
        N = self.names
        if a[0] == 'i' and b[0] == 'i':
            if a[1] == b[1]:
                return SR(self.dimv)
            return T.metric(N[a[1]], N[b[1]], euclidean=self.euclid)
        if a[0] == 'p' and b[0] == 'p':
            return T.dot(N[a[1]], N[b[1]])
        idx, mom = (a[1], b[1]) if a[0] == 'i' else (b[1], a[1])
        return T.comp(N[mom], N[idx])


_CH_FACTOR = {'1': None,
              '5': lambda: Factor([(('5',), Integer(1))]),
              'L': lambda: Factor([(('1',), Integer(1) / 2), (('5',), -Integer(1) / 2)]),
              'R': lambda: Factor([(('1',), Integer(1) / 2), (('5',), Integer(1) / 2)])}


def _expand_line(c, atoms):
    """Coefficient and a list of atoms (FORM output of one line) -> [(coefficient, dirac.py items)]."""
    out = [(SR(c), [])]
    for a in atoms:
        if a[0] == '5':
            out = [(x, it + [('c', '5')]) for x, it in out]
        elif a[0] == '16':                            # g6_ = 1 + g5 = 2 PR
            out = [(2 * x, it + [('c', 'R')]) for x, it in out]
        elif a[0] == '17':                            # g7_ = 1 - g5 = 2 PL
            out = [(2 * x, it + [('c', 'L')]) for x, it in out]
        else:
            out = [(x, it + [(a[0], str(a[1]))]) for x, it in out]
    return out


def _simplify_line(c, atoms, left, right, dimv, euclid, names):
    """[(coefficient, factors)] for c * (left) atoms (right), normal ordered with the Dirac equation."""
    pl = str(left.p) if left is not None else None
    pr = str(right.p) if right is not None else None
    order = _Order(pl, pr, dimv, euclid, names)
    res = {}
    for c0, items in _expand_line(c, atoms):
        c1, slots, ch = _d._chiral_right(SR(c0), items)
        if c1.is_trivial_zero():
            continue
        for s, cc in _d._normal_order(c1, slots, order).items():
            s, chs = list(s), ch
            if right is not None and s and s[-1] == ('p', pr):
                s.pop()
                cc *= right.m if right.particle == 'u' else -right.m
                if chs == '5':
                    cc = -cc
                elif chs in ('L', 'R'):
                    chs = 'R' if chs == 'L' else 'L'
            if left is not None and s and s[0] == ('p', pl):
                s.pop(0)
                cc *= left.m if left.particle == 'u' else -left.m
            k = (tuple(s), chs)
            res[k] = res.get(k, SR(0)) + cc
    out = []
    for (s, chs), cc in res.items():
        cc = cc.expand()
        if cc.is_trivial_zero():
            continue
        fs = [Factor([((kind, names[n]), Integer(1))]) for kind, n in s]
        if _CH_FACTOR[chs] is not None:
            fs.append(_CH_FACTOR[chs]())
        out.append((cc, fs))
    return out


def simplify_dirac(expr, dim=4, rules=None, euclidean=None, eps_in_d=False, debug=False):
    r"""
    Simplify every Dirac line of expr: contract the tensors into the gamma matrices (FORM), order the
    matrices, sum repeated indices, use p/ p/ = p^2 and the Dirac equation at spinors.  Returns a
    DiracExpr (or a Sage number when nothing Dirac is left).  rules (dot(p, p): m^2, ...) are put in
    at the end.
    """
    dim = T.normalize_dim(dim)
    if not isinstance(expr, DiracExpr):
        from .operations import contract
        return contract(expr, rules=rules, dim=dim, euclidean=euclidean, eps_in_d=eps_in_d, debug=debug)
    from .operations import _guard_indices
    _guard_indices(expr)
    eu = T._euclid(euclidean)
    dimv = 4 if dim == 4 else SR(dim)
    prog = Program(dim, eu)
    texts, shapes = [], []
    for c, o, cl in expr.terms:
        chains = ([o] if o is not None else []) + list(cl)
        parts = [prog.scalar(c)]
        for n, ch in enumerate(chains, 1):
            parts.append(prog.chain(ch.factors, n))
        texts.append('*'.join(parts))
        shapes.append((o is not None, chains))
    _check(prog, dim, eps_in_d)
    bodies = _execute(prog, texts, 0, contract=True, debug=debug)
    rd = Reader(prog)
    names = {}
    for kind, x in prog.sage_of.values():
        if kind in ('index', 'vector'):
            names[str(x)] = x
    out = []
    for body, (has_open, chains) in zip(bodies, shapes):
        for c, lines in rd.terms(body):
            for at in lines.values():
                for a in at:
                    if a[0] in ('i', 'p'):
                        names[str(a[1])] = a[1]
            pieces = [[(c, None)]]
            for n, ch in enumerate(chains, 1):
                alts = _simplify_line(1, lines.get(n, []), ch.left, ch.right, dimv, eu, names)
                pieces.append([(x, Chain(ch.left, fs, ch.right)) for x, fs in alts])
            combos = [(SR(1), [])]
            for alts in pieces:
                combos = [(cx * x, chs + ([chn] if chn is not None else [])) for cx, chs in combos for x, chn in alts]
            for cx, chs in combos:
                op = chs[0] if has_open else None
                closed = tuple(chs[1:] if has_open else chs)
                out.append((cx, op, closed))
    res = DiracExpr(_combine(out))
    if rules:
        sub = {SR(k): SR(v) for k, v in rules.items() if not T.is_momentum(SR(k))}
        res = DiracExpr(_combine([(SR(c).subs(sub).expand(), o, cl) for c, o, cl in res.terms]))
    res = DiracExpr([(SR(c).expand(), o, cl) for c, o, cl in res.terms])
    if all(o is not None and not o.factors and not cl for _, o, cl in res.terms) and res.terms:
        return sum((c for c, _, _ in res.terms), SR(0))              # only multiples of the unit matrix
    if all(o is None and not cl for _, o, cl in res.terms):
        return sum((c for c, _, _ in res.terms), SR(0))
    return res
