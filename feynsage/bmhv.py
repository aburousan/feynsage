r"""
Dirac traces in the BMHV (Breitenlohner-Maison-'t Hooft-Veltman) scheme of gamma5, the scheme in which gamma5
needs no care with two or more fermion lines, at the price of 4- and (d-4)-dimensional pieces.

    bmhv_trace(G, four_dim=())        Tr[G] for a product G of gamma(mu), slash(p), gamma5() (Tr 1 = 4)

Space-time is split, d = 4 + (d - 4): every index, metric and momentum has a 4-dimensional (bar) and a
(d-4)-dimensional (hat) part, g = gbar + ghat, gbar.ghat = 0.  gamma5 = i g^0 g^1 g^2 g^3 anticommutes with the
barred gammas and commutes with the hatted ones, so  g5 g^mu = -g^mu g5 + 2 ghat^mu g5.  The trace is done in
three steps:
  1. every gamma5 is moved to the right end with that rule (g5 g5 = 1);
  2. every gamma is split into its bar and hat parts; an index summed inside the trace must be bar on both
     gammas or hat on both (a bar index contracted with a hat one is zero); the hats are moved to the right of
     the bars (they anticommute);
  3. Tr[bars hats (g5)] = Tr_4[bars (g5)] Tr_(d-4)[hats] / 4, each piece by FORM in its own dimension.
The result uses MetricBar, MetricHat, CompBar, CompHat, DotBar, DotHat and Epsilon (four-dimensional);
to_d() writes the bar objects as full minus hat (gbar = g - ghat), the form FeynCalc prints with $BreitMaison.
Momenta in four_dim have no hat part (external momenta are usually taken four-dimensional).
"""
from itertools import product

from sage.all import SR, function, Integer

from .qft import tensors as T

MetricBar = function("MetricBar", nargs=2, print_latex_func=lambda s, a, b: r"\bar g^{%s %s}" % (a, b))
MetricHat = function("MetricHat", nargs=2, print_latex_func=lambda s, a, b: r"\hat g^{%s %s}" % (a, b))
CompBar = function("CompBar", nargs=2, print_latex_func=lambda s, p, a: r"\bar{%s}^{%s}" % (p, a))
CompHat = function("CompHat", nargs=2, print_latex_func=lambda s, p, a: r"\hat{%s}^{%s}" % (p, a))
DotBar = function("DotBar", nargs=2, print_latex_func=lambda s, p, q: r"\bar{%s}\cdot\bar{%s}" % (p, q))
DotHat = function("DotHat", nargs=2, print_latex_func=lambda s, p, q: r"\hat{%s}\cdot\hat{%s}" % (p, q))


def _atom_sequences(G):
    """G (a DiracExpr without spinors) as [(coefficient, [atoms])], every factor expanded into single matrices."""
    out = []
    for c, o, cl in G.terms:
        if cl:
            raise ValueError("bmhv_trace takes a product of Dirac matrices (no spinors)")
        seqs = [(SR(c), [])]
        for F in (o.factors if o is not None else ()):
            new = []
            for cc, atoms in seqs:
                for atom, w in F.items:
                    if atom[0] == '1':
                        new.append((cc * w, atoms))
                    elif atom[0] in ('i', 'p', '5'):
                        new.append((cc * w, atoms + [atom]))
                    else:
                        raise ValueError("bmhv_trace: unknown Dirac matrix %r" % (atom,))
            seqs = new
        out += seqs
    return out


def _g5_right(seq):
    """Move every gamma5 to the right end: [(sign, slots, n5 mod 2)]; slots are (atom, 'd' or 'h')."""
    work = [(Integer(1), [(a, 'd') if a[0] != '5' else (a, '5') for a in seq])]
    done = []
    while work:
        c, s = work.pop()
        # the leftmost gamma5 that has something other than gamma5s to its right
        k = next((j for j in range(len(s) - 1) if s[j][1] == '5' and any(x[1] != '5' for x in s[j + 1:])), None)
        if k is None:
            n5 = sum(1 for x in s if x[1] == '5')
            done.append((c, [x for x in s if x[1] != '5'], n5 % 2))
            continue
        nxt = s[k + 1]
        if nxt[1] == '5':                                  # g5 g5 = 1
            work.append((c, s[:k] + s[k + 2:]))
        elif nxt[1] == 'h':                                # g5 ghat = ghat g5
            work.append((c, s[:k] + [nxt, s[k]] + s[k + 2:]))
        else:                                              # g5 g = -g g5 + 2 ghat g5
            work.append((-c, s[:k] + [nxt, s[k]] + s[k + 2:]))
            work.append((2 * c, s[:k] + [(nxt[0], 'h'), s[k]] + s[k + 2:]))
    return done


def bmhv_trace(G, four_dim=()):
    """Tr[G] in the BMHV scheme (Tr 1 = 4) for a DiracExpr G made of gamma(mu), slash(p), gamma5() and numbers.
    four_dim: momenta without a (d-4)-dimensional part."""
    from .qft.operations import dirac_trace
    from .qft.objects import DiracExpr, Chain, Factor
    dh = SR.var('fsdh')
    fourd = {str(SR(p)) for p in four_dim}
    bar_jobs, hat_jobs, combos = [], [], []

    def matrix(atoms):
        return DiracExpr([(SR(1), Chain(None, tuple(Factor([(a, SR(1))]) for a in atoms), None), ())])
    for c, seq in _atom_sequences(G):
        for c5, slots, odd in _g5_right(seq):
            free = [k for k, (a, t) in enumerate(slots) if t == 'd']
            idx_pos = {}
            for k, (a, t) in enumerate(slots):
                if a[0] == 'i':
                    idx_pos.setdefault(str(a[1]), []).append(k)
            for choice in product('bh', repeat=len(free)):
                kind = [t for a, t in slots]
                for k, ch in zip(free, choice):
                    kind[k] = ch
                if any(kind[k] == 'h' and slots[k][0][0] == 'p' and str(slots[k][0][1]) in fourd for k in range(len(slots))):
                    continue                               # a four-dimensional momentum has no hat part
                if any(len(ps) == 2 and kind[ps[0]] != kind[ps[1]] for ps in idx_pos.values()):
                    continue                               # a summed index: bar.hat = 0
                bars = [slots[k][0] for k in range(len(slots)) if kind[k] == 'b']
                hats = [slots[k][0] for k in range(len(slots)) if kind[k] == 'h']
                if len(hats) % 2 or (len(bars) % 2) or (odd and len(bars) < 4):
                    continue
                # move the hats to the right of the bars: one sign for every (hat, bar) pair out of order
                sign, nh = 1, 0
                for k in range(len(slots)):
                    if kind[k] == 'h':
                        nh += 1
                    else:
                        sign *= (-1) ** nh
                bar_jobs.append(matrix(bars + ([('5',)] if odd else [])) if bars or odd else None)
                hat_jobs.append(matrix(hats) if hats else None)
                combos.append(c * c5 * sign)
    if not combos:
        return SR(0)
    bt = dirac_trace([b for b in bar_jobs if b is not None], dim=4) if any(b is not None for b in bar_jobs) else []
    ht = dirac_trace([h for h in hat_jobs if h is not None], dim=dh) if any(h is not None for h in hat_jobs) else []
    bt, ht = list(bt), list(ht)
    rename_bar = lambda e: SR(e).substitute_function(_op(T.metric), MetricBar).substitute_function(
        _op(T.comp), CompBar).substitute_function(_op(T.dot), DotBar)
    rename_hat = lambda e: SR(e).substitute_function(_op(T.metric), MetricHat).substitute_function(
        _op(T.comp), CompHat).substitute_function(_op(T.dot), DotHat).subs({dh: SR.var('d') - 4})
    total = SR(0)
    for c, b, h in zip(combos, bar_jobs, hat_jobs):
        tb = rename_bar(bt.pop(0)) if b is not None else SR(4)
        th = rename_hat(ht.pop(0)) if h is not None else SR(4)
        total += c * tb * th / 4
    return total.expand()


def _op(f):
    """The Sage function behind metric, comp, dot (whose Python wrappers check their arguments)."""
    a, b = T.lorentz_indices("fsbma fsbmb")
    p, q = T.momenta("fsbmp fsbmq")
    return {T.metric: SR(T.metric(a, b)), T.comp: SR(T.comp(p, a)), T.dot: SR(T.dot(p, q))}[f].operator()


def to_d(e):
    """The bar objects written as d-dimensional minus hat: gbar = g - ghat, pbar = p - phat, pbar.qbar = p.q - phat.qhat."""
    e = SR(e)
    w = [SR.wild(k) for k in range(2)]
    sub = {}
    for F in e.find(MetricBar(w[0], w[1])):
        sub[F] = SR(T.metric(*F.operands())) - MetricHat(*F.operands())
    for F in e.find(CompBar(w[0], w[1])):
        sub[F] = SR(T.comp(*F.operands())) - CompHat(*F.operands())
    for F in e.find(DotBar(w[0], w[1])):
        sub[F] = SR(T.dot(*F.operands())) - DotHat(*F.operands())
    return e.subs(sub).expand() if sub else e


# ---------------------------------------------------------------------------- contraction and one-loop integrals
def _ops():
    a, b, c, e = T.lorentz_indices("fsbca fsbcb fsbcc fsbce")
    p, q = T.momenta("fsbcp fsbcq")
    return (SR(T.metric(a, b)).operator(), SR(T.comp(p, a)).operator(), SR(T.dot(p, q)).operator(),
            SR(T.epsilon(a, b, c, e)).operator())


def _flat(m):
    """A monomial as a list of factors (positive integer powers written out)."""
    op = m.operator()
    fs = m.operands() if (op is not None and "mul" in getattr(op, "__name__", str(op))) else [m]
    out = []
    for f in fs:
        fo = f.operator()
        if fo is not None and "pow" in getattr(fo, "__name__", str(fo)) and f.operands()[1].is_integer() and f.operands()[1] > 0:
            out += [f.operands()[0]] * int(f.operands()[1])
        else:
            out.append(f)
    return out


def bmhv_contract(e, four_dim=(), four_indices=(), dim=None):
    """Contract the repeated Lorentz indices of an expression in Metric, Comp, Dot (d dimensions), MetricHat,
    CompHat, DotHat ((d-4) dimensions) and Epsilon (4 dimensions), after writing the bar objects as full minus hat.
    four_dim: momenta without a hat part, four_indices: free indices that are four-dimensional (external)."""
    MET, COMP, DOT, EPS = _ops()
    d = SR.var('d') if dim is None else SR(dim)
    fourd = {str(SR(p)) for p in four_dim}
    four_i = {str(SR(i)) for i in four_indices}
    e = to_d(SR(e)).expand()
    terms = e.operands() if e.operator() is not None and "add" in getattr(e.operator(), "__name__", "") else [e]
    work = [(SR(1), _flat(t)) for t in terms if not t.is_zero()]
    out = SR(0)
    while work:
        c, fs = work.pop()
        r = _rule(fs, MET, COMP, DOT, EPS, d, fourd, four_i)
        if r is None:
            t = c
            for F in fs:
                t *= F
            out += t
        else:
            for cc, gs in r:
                if not SR(cc).is_zero():
                    work.append((c * cc, gs))
    return out.expand()


def _is_idx(x):
    return T.is_index(SR(x))


def _rule(fs, MET, COMP, DOT, EPS, d, fourd, four_i):
    """One contraction step on a monomial, or None when nothing is left to do."""
    hat = lambda F: F.operator() in (MetricHat, CompHat, DotHat)
    for k, F in enumerate(fs):
        op = F.operator()
        rest = fs[:k] + fs[k + 1:]
        if op in (CompHat,) and (str(F.operands()[0]) in fourd or str(F.operands()[1]) in four_i):
            return [(SR(0), [])]
        if op is DotHat and any(str(x) in fourd for x in F.operands()):
            return [(SR(0), [])]
        if op is MetricHat and any(str(x) in four_i for x in F.operands()):
            return [(SR(0), [])]
        if op is EPS:
            args = [str(x) for x in F.operands()]
            if len(set(args)) < 4:
                return [(SR(0), [])]
        if op in (MET, MetricHat):
            a, b = F.operands()
            if str(a) == str(b):
                return [(d if op is MET else d - 4, rest)]
            for x, y in ((a, b), (b, a)):
                if not _is_idx(x):
                    continue
                for j, G in enumerate(rest):
                    gop = G.operator()
                    if str(x) not in [str(z) for z in G.operands()]:
                        continue
                    others = rest[:j] + rest[j + 1:]
                    if op is MET:                          # g^{x y} X(..x..) -> X(..y..)
                        return [(SR(1), others + [G.subs({x: y})])]
                    # ghat^{x y} with:
                    if gop is MET:
                        z = [u for u in G.operands() if str(u) != str(x)][0]
                        return [(SR(1), others + [MetricHat(y, z)])]
                    if gop is MetricHat:
                        z = [u for u in G.operands() if str(u) != str(x)][0]
                        return [(SR(1), others + [MetricHat(y, z)])]
                    if gop is COMP:
                        p = G.operands()[0]
                        return [(SR(0), [])] if str(p) in fourd else [(SR(1), others + [CompHat(p, y)])]
                    if gop is CompHat:
                        return [(SR(1), others + [CompHat(G.operands()[0], y)])]
                    if gop is EPS:
                        return [(SR(0), [])]
        if op in (COMP, CompHat):
            p, a = F.operands()
            if not _is_idx(a):
                continue
            for j, G in enumerate(rest):
                gop = G.operator()
                if gop in (COMP, CompHat) and str(G.operands()[1]) == str(a):
                    q = G.operands()[0]
                    others = rest[:j] + rest[j + 1:]
                    h = (op is CompHat) or (gop is CompHat)
                    if h:
                        if str(p) in fourd or str(q) in fourd:
                            return [(SR(0), [])]
                        return [(SR(1), others + [DotHat(p, q)])]
                    return [(SR(1), others + [SR(T.dot(p, q))])]
                if gop is EPS and str(a) in [str(z) for z in G.operands()]:
                    others = rest[:j] + rest[j + 1:]
                    if op is CompHat:
                        return [(SR(0), [])]
                    return [(SR(1), others + [G.subs({a: p})])]
    return None


def bmhv_loop(N, loop_momentum, props, kin, four_dim=(), four_indices=()):
    """Int d^d l N / prod(props) in loop()'s normalisation, for a numerator N from bmhv_trace (bar and hat objects,
    the loop momentum l in d dimensions, the external momenta four_dim without hat part).  Every l in an epsilon
    tensor or a hat object becomes l^a with a new index; loop() does the tensor integral; the result is contracted
    in d dimensions (ghat.g = ghat, ghat^a_a = d - 4, ghat.p = 0, eps.ghat = 0) and expanded in eps (d = 4 - 2 eps).
    Returns the coefficient expression (poles and finite part)."""
    from .pv import loop, eps
    MET, COMP, DOT, EPS = _ops()
    l = SR(loop_momentum)
    N = bmhv_contract(N, four_dim, four_indices)
    terms = N.operands() if N.operator() is not None and "add" in getattr(N.operator(), "__name__", "") else [N]
    counter = [0]

    def new():
        counter[0] += 1
        return SR(T.lorentz_indices("fsbl%d" % counter[0]))
    total = SR(0)
    for t in terms:
        outside, lpart = SR(1), []
        for F in _flat(t):
            op = F.operator()
            args = F.operands() if op is not None else []
            if op is DOT and str(l) in [str(x) for x in args]:
                o = [x for x in args if str(x) != str(l)]
                lpart.append("l^2" if not o else "l.%s" % o[0])
            elif op is COMP and str(args[0]) == str(l):
                lpart.append("l^%s" % args[1])
            elif op is EPS and str(l) in [str(x) for x in args]:
                a = new()
                outside *= F.subs({l: a})
                lpart.append("l^%s" % a)
            elif op is DotHat and all(str(x) == str(l) for x in args):
                a, b = new(), new()
                outside *= MetricHat(a, b)
                lpart += ["l^%s" % a, "l^%s" % b]
            elif op is CompHat and str(args[0]) == str(l):
                a = new()
                outside *= MetricHat(args[1], a)
                lpart.append("l^%s" % a)
            else:
                if l in SR(F).variables():
                    raise ValueError("bmhv_loop: cannot integrate the factor %s" % F)
                outside *= F
        r = loop(" ".join(lpart) if lpart else "1", *props, kin=kin)
        for st, coef in r.parts:
            import re
            tensor = SR(1)
            for g_a, g_b, pnm, a in re.findall(r"g\^\{(\w+) (\w+)\}|(\w+)\^(\w+)", str(st)):
                if g_a:
                    tensor *= SR(T.metric(SR(T.lorentz_indices(g_a)), SR(T.lorentz_indices(g_b))))
                else:
                    tensor *= SR(T.comp(SR(T.momenta(pnm)), SR(T.lorentz_indices(a))))
            total += bmhv_contract(outside * tensor, four_dim, four_indices) * SR(coef)
    return canon_epsilon(total.subs({SR.var('d'): 4 - 2 * eps})).expand()


def canon_epsilon(e):
    """Every Epsilon with its arguments in a fixed order, times the sign of the permutation (so that the same
    tensor written in two orders adds up); an Epsilon with a repeated argument is zero."""
    from sage.all import Permutation
    EPS = _ops()[3]
    e = SR(e)
    sub = {}
    for F in e.find(EPS(*[SR.wild(k) for k in range(4)])):
        args = list(F.operands())
        if len({str(a) for a in args}) < 4:
            sub[F] = 0
            continue
        order = sorted(range(4), key=lambda k: str(args[k]))
        sub[F] = Permutation([k + 1 for k in order]).sign() * EPS(*[args[k] for k in order])
    return e.subs(sub) if sub else e
